#ifndef EFG_HPP_CHAINX
#define EFG_HPP_CHAINX

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <filesystem>
#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <iterator>
#include <unordered_map>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp> // init_support for rank queries
#include "concurrentqueue.h" // https://github.com/cameron314/concurrentqueue

//#define EFG_HPP_DEBUG

using std::vector, std::map, std::unordered_map, std::set, std::pair, std::get, std::string, std::istringstream, std::cerr, sdsl::bit_vector;

namespace chainx_block_graph {

class Elasticfoundergraph;
class GAFHit;

class Elasticfoundergraph {
	private:
		int m = 0, n = 0; // rows, cols	
		vector<int> cuts, heights, cumulative_height;
		vector<string> ordered_node_ids;
		unordered_map<string,int> node_indexes;
		vector<string> ordered_node_labels;
		map<int,vector<int>> edges; // adjacency lists // TODO decide this map
		vector<string> walk_ids;
		vector<vector<int>> walks;
		vector<vector<bool>> orientations; // true for +, false for -
		bit_vector leaders;
		//sdsl::rank_support_v5<> block_rank_support;
		vector<int> block;
		// following data structures can be empty

		// shortest_paths[i] is the length of a shortest path from the beginning of first block to the end of the (i-1)-th block, extremes included
		sdsl::int_vector<> shortest_paths; 

	friend bool are_colinear_eds(const GAFHit &a1, const GAFHit &a2, const Elasticfoundergraph &graph);
	friend int max_gap_eds(const GAFHit &a1, const GAFHit &a2, const Elasticfoundergraph &graph);
	friend int overlap_eds(const GAFHit &a1, const GAFHit &a2, const Elasticfoundergraph &graph);
	friend GAFHit dummy_start(const GAFHit &a, const Elasticfoundergraph &graph);
	friend GAFHit dummy_end(const GAFHit &a, const Elasticfoundergraph &graph);

	public:
		Elasticfoundergraph(std::ifstream &graphstream)
		{
			int nodes = 0;

			for (std::string line; std::getline(graphstream, line); ) {
				if (line[0] == 'M') {
					assert(m == 0 && n == 0);

					char c;
					istringstream(line) >> c >> m >> n;
				} else if (line[0] == 'X') {
					assert(cuts.size() == 0);

					istringstream buffer(line.substr(2)); // ignore starting "X\t"
					vector<int> ccuts({
						std::istream_iterator<int>(buffer), std::istream_iterator<int>()
					});
					istringstream().swap(buffer); // clear buffer
					std::swap(ccuts, cuts);
				} else if (line[0] == 'B') {
					assert(heights.size() == 0);

					istringstream buffer(line.substr(2)); // ignore starting "X\t"
					vector<int> hheights({
						std::istream_iterator<int>(buffer), std::istream_iterator<int>()
					});
					istringstream().swap(buffer); // clear buffer
					std::swap(hheights, heights);

					vector<int> cumulative_h(heights.size() + 1, 0);
					for (int i = 1; i <= heights.size(); i++) {
						cumulative_h[i] = cumulative_h[i-1] + heights[i-1];
					}
					std::swap(cumulative_h, cumulative_height);
				} else if (line[0] == 'S') {
					char c;
					string id, label;

					istringstream(line) >> c >> id >> label;
					if (node_indexes.contains(id)) {
						ordered_node_ids[node_indexes[id]] = id;
						ordered_node_labels[node_indexes[id]] = label;
					} else {
						node_indexes[id] = nodes++;
						ordered_node_ids.push_back(id);
						ordered_node_labels.push_back(label);
					}
				} else if (line[0] == 'L') {
					char c;
					std::string id1, id2;
					istringstream(line) >> c >> id1 >> c >> id2;

					if (!node_indexes.contains(id1)) {
						node_indexes[id1] = nodes++;
						ordered_node_ids.push_back("");
						ordered_node_labels.push_back("");
					}
					if (!node_indexes.contains(id2)) {
						node_indexes[id2] = nodes++;
						ordered_node_ids.push_back("");
						ordered_node_labels.push_back("");
					}
					edges[node_indexes[id1]].push_back(node_indexes[id2]);
				} else if (line[0] == 'P') {
					// do nothing
				} else {
					std::cerr << "Unrecognized line " << line[0] << ": skipping..." << std::endl;
				}
			}
			if (!check()) {
				exit(1);
			}
			for (const auto &id    : ordered_node_ids)    assert(id.size() > 0);
			for (const auto &label : ordered_node_labels) assert(label.size() > 0);

			block = vector<int>(ordered_node_ids.size());
			for (int b = 0, i = 0; b < heights.size(); b++) {
				for (int j = 0; j < heights[b]; j++) {
					block[i] = b;
					i += 1;
				}
			}
		}

		bool check()
		{
			if (heights.size() != cuts.size()) {
				std::cerr << "Error: number of cuts and blocks mismatch!" << std::endl;
				return false;
			}

			int block_sum = 0;
			for (auto h : heights) block_sum += h;
			if (block_sum != ordered_node_ids.size()) {
				std::cerr << "Error: sum of block heights does not correspond to node number!" << std::endl;
				return false;
			}

			return true;
		}

		int get_node(string node_id) const
		{
			//TODO: warn the user
			assert(node_indexes.contains(node_id));
			return node_indexes.at(node_id);
		}

		string get_id(int node) const
		{
			//TODO: warn the user
			assert(node >= -1 and node <= (int)ordered_node_ids.size());
			if (node == -1) {
				return "dummystart";
			} else if (node == ordered_node_ids.size()) {
				return "dummyend";
			} else {
				return ordered_node_ids[node];
			}
		}

		int get_label_length(int node) const
		{
			assert(node >= 0 && node < (int)ordered_node_ids.size());
			return ordered_node_labels[node].size();
		}

		int get_block(int node) const
		{
			if (node == -1)
				return -1; // dummy start block
			if (node == ordered_node_ids.size())
				return heights.size(); // dummy end block

			return block.at(node);
		}

		void init_eds_support()
		{

			int shortest_path = 0;
			for (int node = 0, j = 0; j < heights.size(); j++) {
				long unsigned int minlength = ordered_node_labels[node].length();
				for (int k = 0; k < heights[j]; k++) {
					minlength = std::min(minlength, ordered_node_labels[node].length());
					node += 1;
				}
				shortest_path += minlength;
			}

			shortest_paths = sdsl::int_vector(heights.size(), 0, shortest_path);
			shortest_path = 0;
			for (int node = 0, j = 0; j < heights.size(); j++) {
				long unsigned int minlength = ordered_node_labels[node].length();
				for (int k = 0; k < heights[j]; k++) {
					minlength = std::min(minlength, ordered_node_labels[node].length());
					node += 1;
				}
				shortest_paths[j] = shortest_path + minlength;
				shortest_path += minlength;
			}
			sdsl::util::bit_compress(shortest_paths);
#ifdef EFG_HPP_DEBUG
			std::cerr << "shortest path lengths from the start of EDS to the i-th block are ";
			for (auto l : shortest_paths)
				std::cerr << l << " ";
			std::cerr << std::endl;
#endif
		}

		int shortest_path_eds(int block1, int block2) const
		{
			if (block1 > block2)
				return 0;
			else if (block1 == 0)
				return shortest_paths[block2];
			else
				return shortest_paths[block2] - shortest_paths[block1 - 1];
		}

		void to_stream(std::ostream *out)
		{
			*out << "M\t" << m << "\t" << n << std::endl;

			*out << "X\t";
			for (auto x : cuts) {
				if (x == -1)
					continue;

				*out << x << "\t";
			}
			*out << std::endl;

			*out << "B\t";
			for (auto h : heights) {
				if (h == 0)
					continue;

				*out << h << "\t";
			}
			*out << std::endl;

			for (int i = 0; i < ordered_node_ids.size(); i++) {
				if (ordered_node_labels[i] == "")
					continue;

				*out << "S\t" << ordered_node_ids[i] << "\t" << ordered_node_labels[i] << std::endl;
			}

			for (int i = 0; i < ordered_node_ids.size(); i++) {
				if (ordered_node_labels[i] == "")
					continue;

				for (auto j : edges[i]) {
					*out << "L\t" << ordered_node_ids[i] << "\t+\t" << ordered_node_ids[j] << "\t+\t0M" << std::endl;
				}
			}

			// TODO: this part does unnecessary things, remove
			for (int i = 0; i < walks.size(); i++) {
				*out << "P\t" << walk_ids[i] << "\t";

				bool modified = false;
				int startingnode = walks[i][0];
				if (ordered_node_labels[startingnode] != "") {
					*out << ordered_node_ids[startingnode] << ((orientations[i][0]) ? "+" : "-");
				} else {
					modified = true;
					// go backwards until first non-simplified node that reaches this one
					for (int j = startingnode - 1; j >= 0; j -= 1) {
						// we assume reachability
						if (edges[j][0] = startingnode) {
							// update?
							if (ordered_node_labels[j] != "") {
								*out << ordered_node_ids[j] << ((orientations[i][0]) ? "+" : "-");
								break;
							} else {
								startingnode = j;
							}
						}
					}
				}

				for (int j = 1; j < walks[i].size(); j++) {
					if (j == walks[i].size() - 1 && ordered_node_labels[walks[i][j]] == "")
						modified = true;
					if (ordered_node_labels[walks[i][j]] == "")
						continue;

					*out << "," << ordered_node_ids[walks[i][j]] << ((orientations[i][j]) ? "+" : "-");
				}
				*out << "\t*" << std::endl;

				if (modified) {
					std::cerr << "Warning: string spelled by path " << walk_ids[i] << " gets longer." << std::endl;
				}
			}
		}
};

class GAFHit {
	private:
		// query info
		int qlength, qstart, qend;
		// path info
		bool pstrand; // true for +, false for -
		vector<int> path; // ids are the 0-based node index in the graph
		vector<bool> orientations; // true for +, false for -
		int plength, pstart, pend;
		// we ignore residue matches, alignment block length, mapping quality


	friend bool are_colinear_eds(const GAFHit &a1, const GAFHit &a2, const Elasticfoundergraph &graph);
	friend int max_gap_eds(const GAFHit &a1, const GAFHit &a2, const Elasticfoundergraph &graph);
	friend int overlap_eds(const GAFHit &a1, const GAFHit &a2, const Elasticfoundergraph &graph);
	friend GAFHit dummy_start(const GAFHit &a, const Elasticfoundergraph &graph);
	friend GAFHit dummy_end(const GAFHit &a, const Elasticfoundergraph &graph);

	public:
		GAFHit(istringstream &descr, const Elasticfoundergraph &efg, string &qname)
		{
			// query info
			descr >> qname >> qlength >> qstart >> qend;

			// path info
			char strandc;
			descr >> strandc;
			pstrand = (strandc == '+');
			//TODO: warn user that we assume only forward matches
			assert(pstrand);

			string pathdescr;
			descr >> std::ws;
			std::getline(descr, pathdescr, '\t');
			istringstream pathstream(pathdescr);

			char orientation;
			pathstream >> orientation;
			//TODO: generalize orientation or warn user
			assert(orientation == '>');
			string s;
			while (std::getline(pathstream, s, '>')) {
				path.push_back(efg.get_node(s));
				orientations.push_back(true);
			}

			descr >> plength >> pstart >> pend;
		}

		GAFHit(int qlength, int qstart, int qend, vector<int> path, int plength, int pstart, int pend, bool reverse = false)
		{
			this->qlength = qlength;
			this->qstart = qstart;
			this->qend = qend;
			this->pstrand = true;
			this->path = path;
			if (reverse)
				this->orientations = vector<bool>(path.size(), false);
			else
				this->orientations = vector<bool>(path.size(), true);
			this->plength = plength;
			this->pstart = pstart;
			this->pend = pend;
		}

		bool operator==(const GAFHit &a) const {
			return (qstart == a.qstart and 
				qend == a.qend and
				path == a.path and
				pstart == a.pstart and
				pend == a.pend and
				pstrand == a.pstrand and
				orientations == a.orientations);
		}

		int get_query_start() const
		{
			return qstart;
		}

		int get_query_end() const
		{
			return qend;
		}

		int get_query_length() const
		{
			return qlength;
		}

		int start_distance_query(GAFHit &a)
		{
			if (qstart <= a.qstart)
				return a.qstart - qstart - 1;
			else
				return qstart - a.qstart - 1;
		}

		int gap_query(GAFHit &a)
		{
			//assert(qend <= a.qstart);
			return (a.qstart - qend);
		}

		static int gap_query(GAFHit &a1, GAFHit &a2)
		{
			assert(a1.qend <= a2.qstart);
			return (a2.qstart - a1.qend);
		}

		// for reverse complement versions of queries, reverse the match and remove the "_rev" prefix
		void reverse()
		{
			//qname = qname.substr(4);
			int newqstart = qlength - qend;
			int newqend = qlength - qstart + 1;

			qstart = newqstart;
			qend = newqend;

			int newpstart = plength - pend;
			int newpend = plength - pstart + 1;

			pstart = newpstart;
			pend = newpend;

			std::reverse(path.begin(), path.end());
			for (int i = 0; i < orientations.size(); i++)
				orientations[i] = !orientations[i];
		}

		// split long matches into a perfect chain of matches spanning one node
		vector<GAFHit> split(const Elasticfoundergraph &efg)
		{
			vector<GAFHit> sol;
			int qlength = this->qlength;
			int qstart = this->qstart;
			//int qend = this->qend;
			int plength = this->plength;
			int pstart = this->pstart;
			int pend = this->pend;

			for (int i = 0; i < path.size(); i++) {
				int node_length = efg.get_label_length(path[i]);
				if (pstart >= node_length || pend <= 0)
					continue;
				// the match must involve node efg.path[i]

				int nodeend = std::min(node_length, pend);
				int match_length = nodeend - pstart;
				if (this->orientations[i])
					sol.push_back(GAFHit(
								qlength,
								qstart,
								qstart + match_length,
								vector<int>( {path[i]} ),
								node_length,
								pstart,
								nodeend
							       ));
				else
					sol.push_back(GAFHit(
								qlength,
								qstart,
								qstart + match_length,
								vector<int>( {path[i]} ),
								node_length,
								pstart,
								nodeend,
								true
							       ));

				qstart += match_length;
				plength -= node_length;
				pstart = 0;
				pend -= node_length;
			}

			return sol;
		}

		// split long matches into an perfect chain of matches spanning one node, but remove the resulting matches of length 1 
		vector<GAFHit> split_single_graphaligner(const Elasticfoundergraph &efg)
		{
			vector<GAFHit> sol;
			int qlength = this->qlength;
			int qstart = this->qstart;
			//int qend = this->qend;
			int plength = this->plength;
			int pstart = this->pstart;
			int pend = this->pend;

			for (int i = 0; i < path.size(); i++) {
				int node_length = efg.get_label_length(path[i]);
				if (pstart >= node_length || pend <= 0)
					continue;
				// the match must involve node efg.path[i]
				assert(pstart >= 0);

				int nodeend = std::min(node_length, pend);
				int match_length = nodeend - pstart;
				if (match_length > 1) {
					sol.push_back(GAFHit(
								qlength,
								qstart,
								qstart + match_length,
								vector<int>( {path[i]} ),
								node_length,
								pstart,
								nodeend,
								!this->orientations[i]
							       ));
				}

				qstart += match_length;
				plength -= node_length;
				pstart = 0;
				pend -= node_length;
			}

			return sol;
		}

		string to_string(const Elasticfoundergraph &efg, const string &qname)
		{
			string out;
			// GAF format
			out += qname + "\t";                  // query name
			out += std::to_string(qlength) + "\t";                // query length
			out += std::to_string(qstart) + "\t";                 // query start (0-based, closed)
			out += std::to_string(qend) + "\t";                   // query end (open)
			out += std::string() + ((pstrand) ? "+" : "-") + "\t"; // strand

			for (int i = 0; i < path.size(); i++) { // path
				out += std::string() + ((orientations[i]) ? ">" : "<"); 
				out += efg.get_id(path[i]);
			}
			out += "\t";

			out += std::to_string(plength) + "\t"; // path length
			out += std::to_string(pstart) + "\t";  // start position on the path
			out += std::to_string(pend) + "\t";    // end position on the path
			out += "0\t0\t255";    // FIXME

			return out;
		}

		string to_string_node(Elasticfoundergraph &efg, int node, const string &qname)
		{
			string out;
			int qlength_corrected;
			int qstart_corrected;
			int qend_corrected;
			// GAF format
			out += qname + "\t";                  // query name
			out += std::to_string(qlength) + "\t";                // query length
			out += std::to_string(qstart) + "\t";                 // query start (0-based, closed)
			out += std::to_string(qend) + "\t";                   // query end (open)
			out += std::string() + ((pstrand) ? "+" : "-") + "\t"; // strand

			for (int i = 0; i < path.size(); i++) { // path
				out += std::string() + ((orientations[i]) ? ">" : "<"); 
				out += efg.get_id(path[i]);
			}
			out += "\t";

			out += std::to_string(plength) + "\t"; // path length
			out += std::to_string(pstart) + "\t";  // start position on the path
			out += std::to_string(pend) + "\t";    // end position on the path
			out += "0\t0\t255";    // FIXME

			return out;
		}
};

string read_gaf_query_id(string const &gafline) {
	string id;
	istringstream(gafline) >> id;
	return id;
}

std::unordered_map<string, vector<GAFHit>> read_gaf_chaining(std::ifstream &anchorsstream, Elasticfoundergraph &efg, moodycamel::ConcurrentQueue<string> &anchorsqueue)
{
	unordered_map<string,vector<GAFHit>> buckets;
	// read file
	// sort anchors in a bucket for each corresponding query
	{
		string line;
		while (std::getline(anchorsstream, line)) {
			istringstream linestream(line);
			string qname;
			const GAFHit a(linestream, efg, qname);
			if (buckets.contains(qname)) {
				buckets[qname].push_back(a);
			} else {
				anchorsqueue.enqueue(qname);
				buckets[qname] = vector<GAFHit>({ dummy_start(a, efg) });
				buckets[qname].push_back(a);
			}
		}
	}

	return buckets;
}

std::unordered_map<string, vector<GAFHit>> read_gaf_chaining(std::stringstream &anchorsstream, const Elasticfoundergraph &efg)
{
	unordered_map<string,vector<GAFHit>> buckets;
	// read file
	// sort anchors in a bucket for each corresponding query
	{
		string line;
		while (std::getline(anchorsstream, line)) {
			istringstream linestream(line);
			string qname;
			const GAFHit a(linestream, efg, qname);
			if (buckets.contains(qname)) {
				buckets[qname].push_back(a);
			} else {
				buckets[qname] = vector<GAFHit>({ dummy_start(a, efg) });
				buckets[qname].push_back(a);
			}
		}
	}

	return buckets;
}

// colinear here means that startpoints (endpoints) strictly precede one another
// (i_a < j_a && i_b < j_b && i_c < j_c && i_d < j_d in the original algorithm)
bool are_colinear_eds(const GAFHit &a1, const GAFHit &a2, const Elasticfoundergraph &graph)
{
	assert(a1.path.size() <= 2 && a2.path.size() <= 2);

	// TODO: implement short circuit differently?
	if (a1.qstart >= a2.qstart || a1.qend >= a2.qend)
		return false;

	// no path overlap
	if (graph.get_block(a1.path.back()) < graph.get_block(a2.path[0]))
		return true;

	// the paths match
	if (a1.path == a2.path)
		return (a1.pstart < a2.pstart && a1.pend < a2.pend);

	// one node overlap, paths do not match
	if (a1.path.back() == a2.path[0]) {
		if (a1.path.size() == 1)
			return (a1.pstart < a2.pstart && a1.pend < a2.pend);
		else // a1 starts from previous block
			return (a1.pstart < graph.get_label_length(a1.path[0]) + a2.pstart && a1.pend < graph.get_label_length(a1.path[0]) + a2.pend);
	}

	return false;
}

int max_gap_eds(const GAFHit &a1, const GAFHit &a2, const Elasticfoundergraph &graph)
{
	// we assume that the GAFHit utilizes all nodes in the path 
	assert(a1.path.size() <= 2 && a2.path.size() <= 2);
	assert(graph.shortest_paths.size() > 0);
	int graphdistance;

	const int b1 = graph.get_block(a1.path.back());
	const int b2 = graph.get_block(a2.path.at(0));

	if (b1 < b2) { // paths do not overlap
		const int a1suf	= a1.plength - a1.pend;
		const int a2pre = a2.pstart;

		graphdistance = a1suf + graph.shortest_path_eds(b1 + 1, b2 - 1) + a2pre;
	} else if (a1.path.back() == a2.path.at(0)) { // paths overlap in one node
		const int a1suf = a1.plength - a1.pend;
		const int a1end = graph.get_label_length(a1.path.back()) - a1suf;

		graphdistance = std::max(0, a2.pstart - a1end);
	} else if (a1.path == a2.path) {
		// paths overlap in two nodes
		graphdistance = std::max(0, a2.pstart - a1.pend);
	} else {
		// they are not co-linear?
		assert(false);
	}

	return std::max(graphdistance, std::max(0, a2.qstart - a1.qend));
}

int overlap_eds(const GAFHit &a1, const GAFHit &a2, const Elasticfoundergraph &graph)
{
	const int b1 = graph.get_block(a1.path.back());
	const int b2 = graph.get_block(a2.path[0]);
	int graphoverlap = 0;

	if (b1 < b2) { // paths do not overlap
		graphoverlap = 0;
	} else if (a1.path.back() == a2.path[0]) { // paths overlap in one node
		const int a1suf = a1.plength - a1.pend;
		const int a1end = graph.get_label_length(a1.path.back()) - a1suf;

		graphoverlap = std::max(0, a1end - a2.pstart);
	} else if (a1.path == a2.path) {
		// paths overlap in two nodes
		graphoverlap = std::max(0, a1.pend - a2.pstart);
	} else {
		// they are not co-linear?
		assert(false);
	}

	return std::abs(graphoverlap - std::max(0, a1.qend - a2.qstart));
}

GAFHit dummy_start(const GAFHit &a, const Elasticfoundergraph &graph)
{
	return GAFHit(
			a.qlength,
			-1, // qstart
			0, // qend
			vector<int>( {-1} ), // dummy node -1
			0, // plength
			0, // pstart,
			0 // pend
		);
}

GAFHit dummy_end(const GAFHit &a, const Elasticfoundergraph &graph)
{
	return GAFHit(
			a.qlength,
			a.qlength, // qstart
			a.qlength+1, // qend
			vector<int>( {graph.ordered_node_ids.size()} ), // dummy node n, where 0..n-1 are the proper nodes
			0, // plength
			0, // pstart,
			0 // pend
		);
}


bool is_sorted(const vector<GAFHit> &hits)
{
	for (long int i = 0; i < hits.size() - 1; i++) {
		if (hits[i].get_query_start() > hits[i+1].get_query_start())
			return false;
	}
	return true;
}

int compute_coverage_greedy(const vector<GAFHit> &anchorlist)
{
	// we assume anchorlist is sorted by starting position in the query
	// we assume the anchors were computed greedily and are NOT nested
	int coverage = 0;
	int processed = 0;
	for (const GAFHit &a : anchorlist) {
		if (a.get_query_end() <= processed)
			continue;

		coverage += a.get_query_end() - std::max(a.get_query_start(), processed + 1);
		processed = a.get_query_end();
	}

	assert(coverage <= anchorlist.at(0).get_query_length());
	return coverage;
}

} // Namespace chainx_block_graph

#endif
