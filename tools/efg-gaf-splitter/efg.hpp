#ifndef EFG_HPP
#define EFG_HPP

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
#include <utility>
#include <algorithm>

//#define EFG_HPP_DEBUG

using std::vector, std::map, std::unordered_map, std::set, std::pair, std::get, std::string, std::istringstream, std::cerr;

class Elasticfoundergraph;
class GAFAnchor;

class Elasticfoundergraph {

	private:
		int m = 0, n = 0; // rows, cols	
		vector<int> cuts, heights, cumulative_height;
		vector<string> ordered_node_ids;
		unordered_map<string,int> node_indexes;
		vector<string> ordered_node_labels;
		map<int,vector<int>> edges; // adjacency lists // TODO decide this map

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
					node_indexes[id] = nodes;
					ordered_node_ids.push_back(id);
					ordered_node_labels.push_back(label);

					nodes += 1;
				} else if (line[0] == 'L') {
					char c;
					std::string id1, id2;
					istringstream(line) >> c >> id1 >> c >> id2;

					edges[node_indexes[id1]].push_back(node_indexes[id2]);
				} else if (line[0] == 'P') {
				} else {
					std::cerr << "Unrecognized line " << line[0] << ": skipping..." << std::endl;
				}
			}
			if (!check()) {
				exit(1);
			}
		}

		bool check() const
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

		int get_node(const string node_id) const
		{
			assert(node_indexes.contains(node_id));
			return node_indexes.at(node_id);
		}

		string get_label(int node) const
		{
			assert(node >= 0 && node < (int)ordered_node_ids.size());
			return ordered_node_labels[node];
		}

		string get_id(int node) const
		{
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

		void to_stream(std::ostream *out) const
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
				*out << "S\t" << ordered_node_ids[i] << "\t" << ordered_node_labels[i];
				*out << std::endl;
			}

			for (int i = 0; i < ordered_node_ids.size(); i++) {
				for (auto j : edges.at(i)) {
					*out << "L\t" << ordered_node_ids[i] << "\t+\t" << ordered_node_ids[j] << "\t+\t0M" << std::endl;
				}
			}
		}
};

class GAFAnchor {
	private:
		// query info
		string qname;
		int qlength, qstart, qend;
		// path info
		bool pstrand; // true for +, false for -
		vector<int> path; // ids are the 0-based node index in the graph
		vector<bool> orientations; // true for +, false for -
		int plength, pstart, pend;
		// we ignore residue matches, alignment block length, mapping quality


	public:
		GAFAnchor(istringstream &descr, Elasticfoundergraph const &efg)
		{
			// query info
			descr >> qname >> qlength >> qstart >> qend;

			// path info
			char strandc;
			descr >> strandc;
			pstrand = ((strandc == '+') ? true : false);
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

		GAFAnchor(string qname, int qlength, int qstart, int qend, vector<int> path, int plength, int pstart, int pend, bool reverse = false)
		{
			this->qname = qname;
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

			assert(pend <= plength);
		}

		// constructor that FIXES the input by checking the graph
		GAFAnchor(const Elasticfoundergraph &efg, string qname, int qlength, int qstart, int qend, vector<int> path, int plength, int pstart, int pend, bool reverse = false)
		{
			this->qname = qname;
			this->qlength = qlength;
			this->qstart = qstart;
			this->qend = qend;
			this->pstrand = true;
			// 
			int imin = 0;
			while (pstart >= efg.get_label_length(path[imin])) {
				pstart  -= efg.get_label_length(path[imin]);
				pend    -= efg.get_label_length(path[imin]);
				plength -= efg.get_label_length(path[imin]);
				imin += 1;
			}
			int imax = imin, pplength = 0; //plength gets ignored?
			while (imax < path.size() and pplength < pend) {
				pplength += efg.get_label_length(path[imax]);
				imax += 1;
			}

			this->path.clear();
			for (int i = imin; i < imax; i++)
				this->path.push_back(path[i]);

			assert(this->path.size() <= 2);
			if (reverse)
				this->orientations = vector<bool>(this->path.size(), false);
			else
				this->orientations = vector<bool>(this->path.size(), true);
			this->plength = pplength;
			this->pstart = pstart;
			this->pend = pend;
		}

		// empty constructor
		GAFAnchor()
		{
		}

		int get_query_start()
		{
			return qstart;
		}

		string get_query_id()
		{
			return qname;
		}

		int get_path_length()
		{
			return path.size();
		}

		int get_length()
		{
			return qend - qstart;
		}

		int start_distance_query(GAFAnchor &a)
		{
			if (qstart <= a.qstart)
				return a.qstart - qstart - 1;
			else
				return qstart - a.qstart - 1;
		}

		static int gap_query(GAFAnchor &a1, GAFAnchor &a2)
		{
			assert(a1.qend <= a2.qstart);
			return (a2.qstart - a1.qend);
		}

		void reverse()
		{
			qname = qname.substr(4);
			int newqstart = qlength - qend;
			int newqend = qlength - qstart;

			qstart = newqstart;
			qend = newqend;

			int newpstart = plength - pend;
			int newpend = plength - pstart;

			assert(newpstart >= 0);
			pstart = newpstart;
			pend = newpend;

			std::reverse(path.begin(), path.end());
			for (int i = 0; i < orientations.size(); i++)
				orientations[i] = !orientations[i];
		}

		// split long matches into a perfect chain of matches spanning one node
		vector<GAFAnchor> split_single(const Elasticfoundergraph &efg) const
		{
			vector<GAFAnchor> sol;
			int qlength = this->qlength;
			int qstart = this->qstart;
			//int qend = this->qend;
			int plength = this->plength;
			int pstart = this->pstart;
			int pend = this->pend;

			for (int i = 0; i < path.size(); i++) {
				int node_length = efg.get_label_length(path[i]);
				if (pstart >= node_length || pend <= 0) {
					pstart -= node_length;
					plength -= node_length;
					pend -= node_length;
					continue;
				}
				// the match must involve node efg.path[i]

				int nodeend = std::min(node_length, pend);
				int match_length = nodeend - pstart;
				if (this->orientations[i])
					sol.push_back(GAFAnchor(
								qname,
								qlength,
								qstart,
								qstart + match_length,
								vector<int>( {path[i]} ),
								node_length,
								pstart,
								nodeend
							       ));
				else
					sol.push_back(GAFAnchor(
								qname,
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
				//plength -= node_length;
				pstart = 0;
				pend -= node_length;
			}

			return sol;
		}

		// split long matches into a perfect chain of matches spanning one node, but remove the resulting matches of length 1 
		vector<GAFAnchor> split_single_graphaligner(const Elasticfoundergraph &efg) const
		{
			vector<GAFAnchor> sol;
			int qlength = this->qlength;
			int qstart = this->qstart;
			//int qend = this->qend;
			//int plength = this->plength;
			int pstart = this->pstart;
			int pend = this->pend;

			for (int i = 0; i < path.size(); i++) {
				int node_length = efg.get_label_length(path[i]);
				if (pstart >= node_length || pend <= 0) {
					pstart -= node_length;
					pend -= node_length;
					continue;
				}
				// the match must involve node efg.path[i]
				assert(pstart >= 0);

				int nodeend = std::min(node_length, pend);
				int match_length = nodeend - pstart;
				if (match_length > 1) {
					sol.push_back(GAFAnchor(
								qname,
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
				//plength -= node_length;
				pstart = 0;
				pend -= node_length;
			}

			return sol;
		}

		string to_string(const Elasticfoundergraph &efg)
		{
			string out;
			// GAF format
			out += qname + "\t";                                  // query name
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

		string to_string_node(Elasticfoundergraph &efg, int node)
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

void anchors_to_stream(std::ostream *out, const Elasticfoundergraph &graph, vector<GAFAnchor> &matches)
{
	for (GAFAnchor m : matches) {
		*out << m.to_string(graph) << std::endl;
	}
}

void anchors_to_stream_split_single(std::ostream *out, const Elasticfoundergraph &graph, const vector<GAFAnchor> &matches)
{
	for (GAFAnchor m : matches) {
		for (GAFAnchor n : m.split_single(graph)) {
			*out << n.to_string(graph) << std::endl;
		}
	}
}
void anchors_to_stream_split_single_graphaligner(std::ostream *out, const Elasticfoundergraph &graph, vector<GAFAnchor> &matches)
{
	for (GAFAnchor m : matches) {
		for (GAFAnchor n : m.split_single_graphaligner(graph)) {
			*out << n.to_string(graph) << std::endl;
		}
	}
}

vector<vector<GAFAnchor>> read_gaf(std::ifstream &anchorsstream, Elasticfoundergraph const &efg)
{
	vector<GAFAnchor> unsorted_anchors;
	// read file
	{
		string line;
		while (std::getline(anchorsstream, line)) {
			istringstream linestream(line);
			unsorted_anchors.push_back(GAFAnchor(linestream, efg));
		}
	}

	// sort anchors in a bucket for each corresponding query
	unordered_map<string,vector<GAFAnchor>> buckets;
	while (!unsorted_anchors.empty()) {
		GAFAnchor a = unsorted_anchors.back();
		unsorted_anchors.pop_back();
		const string query = a.get_query_id();
		if (buckets.contains(query))
			buckets[query].push_back(a);
		else
			buckets[query] = vector<GAFAnchor>({ a });
	}

	// unload the buckets into a vector<vector<GAFAnchor>>
	vector<vector<GAFAnchor>> anchors;
	for (pair<string,vector<GAFAnchor>> bucket : buckets) {
		anchors.push_back(vector<GAFAnchor>(std::move(bucket.second)));
	}

	return anchors;
}

bool read_gaf_single(std::ifstream &anchorsstream, Elasticfoundergraph const &efg, GAFAnchor &seed)
{
	string line;
	while (std::getline(anchorsstream, line)) {
		if (line == "")
			continue;
		istringstream linestream(line);
		seed = GAFAnchor(linestream, efg);
		return 1;
	}
	return 0; // file is consumed
}

#endif
