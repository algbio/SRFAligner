#ifndef EFG_HPP_LOCATE
#define EFG_HPP_LOCATE

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
#include <utility> // std::tie, <
#include <sdsl/bit_vectors.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/util.hpp> // init_support for rank queries
#include <sdsl/config.hpp> // util things

#include "efg-locate.hpp" // input parameters (Param)

using namespace sdsl;
typedef sdsl::csa_wt<wt_hutu<>, 16, 64, text_order_sa_sampling<>, isa_sampling<>, succinct_byte_alphabet<>> csa_type;

//#define EFG_HPP_DEBUG

using std::vector, std::map, std::unordered_map, std::set, std::pair, std::get, std::string, std::istringstream, std::cerr, sdsl::bit_vector, std::tie;
typedef sdsl::csa_wt<>::size_type size_type;

namespace efg_locate {
class Elasticfoundergraph;
class GAFAnchor;

class Elasticfoundergraph {
	friend int efg_backward_search(Elasticfoundergraph &efg, const string &pattern, vector<int> &path);
	friend int efg_backward_search_old(Elasticfoundergraph &efg, const string &pattern, vector<int> &path);
	friend int efg_backward_search_ignorechars(Elasticfoundergraph &efg, const string &ignorechars, const string &pattern, vector<int> &path);
	friend int efg_backward_search_greedy(const Elasticfoundergraph &efg, const string &pattern_id, const string &pattern, const int edgemincount, int &q, vector<GAFAnchor> &match, const string &ignorechars);
	friend int efg_backward_search_greedy_old(const Elasticfoundergraph &efg, const string &pattern_id, const string &pattern, const int edgemincount, int &q, vector<GAFAnchor> &match, const string &ignorechars);
	friend void first_search(const Elasticfoundergraph &, const string &, int &, size_type &, size_type &, int &, size_type &, size_type &, const string &);
	friend void first_search(const Elasticfoundergraph &, const string &, int, int &, size_type &, size_type &, int &, size_type &, size_type &, int &, size_type &, size_type &, const string &);
	friend void simple_search(const Elasticfoundergraph &, const string &, int &, vector<int> &, int &, int &, int &, const string &);
	friend int find_connecting_vertex(const Elasticfoundergraph &, const string &, const int, const int, const int, const int, const int, vector<int> &, int &);

	private:
		int m = 0, n = 0; // rows, cols	
		vector<int> cuts, heights, cumulative_height;
		vector<string> ordered_node_ids;
		unordered_map<string,int> node_indexes;
		vector<string> ordered_node_labels;
		map<int,vector<int>> edges; // adjacency lists
		//vector<string> walk_ids;
		//vector<vector<int>> walks;
		//vector<vector<bool>> orientations; // true for +, false for -

		// following data structures can be empty
		bit_vector is_source;
		csa_type edge_index;
		bit_vector node_leaders, edge_leaders;
		sdsl::rank_support_v5<> node_leaders_rank_support, edge_leaders_rank_support;
		sdsl::select_support_mcl<> node_leaders_select_support, edge_leaders_select_support;


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
				// continue anyway
			}
			for (const auto &id    : ordered_node_ids)    assert(id.size() > 0);
			for (const auto &label : ordered_node_labels) assert(label.size() > 0);
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
			//TODO: warn the user
			assert(node_indexes.contains(node_id));
			return node_indexes.at(node_id);
		}

		string get_label(int node) const
		{
			//TODO: warn the user
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

		void init_pattern_matching_support()
		{
			is_source = bit_vector(ordered_node_ids.size() + 1, true);
			bit_vector is_sink(ordered_node_ids.size() + 1, true);

			for (int i = 0; i < ordered_node_ids.size(); i++) {
				for (int j : edges[i]) {
					is_sink[i] = false;
					is_source[j] = false;
				}
			}

#ifdef EFG_HPP_DEBUG
			for (int i = 0; i < ordered_node_ids.size(); i++) {
				if (is_source[i])
					cerr << "DEBUG: node " << ordered_node_ids[i] << " is a source" << std::endl;
				if (is_sink[i])
					cerr << "DEBUG: node " << ordered_node_ids[i] << " is a sink" << std::endl;
			}
#endif

			// add one supersource with dummy node label that is NOT going to be queried
			int supersource = ordered_node_ids.size();
			ordered_node_ids.push_back("supersource");
			ordered_node_labels.push_back("0"); // TODO ignorechar here?
			for (int i = 0; i < ordered_node_ids.size() - 1; i++) {
				if (is_source[i]) {
					edges[supersource].push_back(i);
				}
			}
			is_sink[supersource] = false;

			// TODO: parameterize/check special characters!
			string edge_concat = "#";
			for (int i = 0; i < ordered_node_ids.size(); i++) {
				if (is_sink[i]) {
					edge_concat += "#";
				} // else
				for (int j : edges[i]) {
					if (is_source[i])
						edge_concat += "$";

					edge_concat += ordered_node_labels[i];
					edge_concat += ordered_node_labels[j];

					if (is_sink[j])
						edge_concat += "$";
					edge_concat += "#";
				}
			}
#ifdef EFG_HPP_DEBUG
			cerr << "DEBUG: edge_concat is " << std::endl;
			cerr << edge_concat << std::endl;
#endif

			node_leaders = bit_vector(edge_concat.size(), 0);
			edge_leaders = bit_vector(edge_concat.size(), 0);
			for (int i = 0, k = 1; i < ordered_node_ids.size(); i++) {
				if (is_source[i])
					node_leaders[k] = 1;
				else
					node_leaders[k] = 1;

				if (is_sink[i]) {
					k += 1;
				} // else
				for (int j : edges[i]) {
					edge_leaders[k] = 1;
					k += ordered_node_labels[i].size() + ordered_node_labels[j].size();
					if (is_source[i])
						k += 1;

					if (is_sink[j])
						k += 1;
					k += 1;
				}
			}
			// a string occurring in position x (0-indexed) 
			node_leaders_rank_support = sdsl::rank_support_v5<>(&node_leaders);
			node_leaders_select_support = sdsl::select_support_mcl<>(&node_leaders);
			// and from the edge_leaders_rank_support()
			edge_leaders_rank_support = sdsl::rank_support_v5<>(&edge_leaders);
			edge_leaders_select_support = sdsl::select_support_mcl<>(&edge_leaders);
#ifdef EFG_HPP_DEBUG
			cerr << "DEBUG: node_leaders is " << std::endl;
			for (auto b : node_leaders)
				cerr << ((b == 0) ? ' ' : '*');
			cerr << std::endl;
			cerr << "DEBUG: edge_leaders is " << std::endl;
			for (auto b : edge_leaders)
				cerr << ((b == 0) ? ' ' : '*');
			cerr << std::endl;
#endif

			sdsl::construct_im(edge_index, edge_concat, 1);

#ifdef EFG_HPP_DEBUG
			cerr << "DEBUG: compressed suffix array is " << std::endl;
			cerr << sdsl::extract(edge_index, 0, edge_index.size()-1) << std::endl;
			cerr << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << std::endl;
			csXprintf(cerr, "%2I %2S %3s %3P %2p %3B   %:1T", edge_index);
#endif
		}

		// stub for storing and loading index, TODO check
		void store_pattern_matching_support(const string &basename)
		{
			cerr << "Storing indexes...";
			sdsl::store_to_file(is_source, basename + ".is_source");
			sdsl::store_to_file(edge_index, basename + ".edge_index");
			sdsl::store_to_file(node_leaders, basename + ".node_leaders");
			sdsl::store_to_file(edge_leaders, basename + ".edge_leaders");
			sdsl::store_to_file(node_leaders_rank_support,   basename + ".node_leaders_rank_support");
			sdsl::store_to_file(node_leaders_select_support, basename + ".node_leaders_select_support");
			sdsl::store_to_file(edge_leaders_rank_support, basename + ".edge_leaders_rank_support");
			sdsl::store_to_file(edge_leaders_select_support, basename + ".edge_leaders_select_support");
			cerr << " done.\n";
		}

		void load_pattern_matching_support(const string &basename)
		{
			if (std::filesystem::exists(std::filesystem::path({ basename + ".edge_index" }))) {
				cerr << "Loading stored indexes...";
				sdsl::load_from_file(is_source, basename + ".is_source");
				sdsl::load_from_file(edge_index, basename + ".edge_index");
				sdsl::load_from_file(node_leaders, basename + ".node_leaders");
				sdsl::load_from_file(edge_leaders, basename + ".edge_leaders");
				sdsl::load_from_file(node_leaders_rank_support, basename + ".node_leaders_rank_support");
				sdsl::util::init_support(node_leaders_rank_support, &node_leaders);
				sdsl::load_from_file(node_leaders_select_support, basename + ".node_leaders_select_support");
				sdsl::util::init_support(node_leaders_select_support, &node_leaders);
				sdsl::load_from_file(edge_leaders_rank_support, basename + ".edge_leaders_rank_support");
				sdsl::util::init_support(edge_leaders_rank_support, &edge_leaders);
				sdsl::load_from_file(edge_leaders_select_support, basename + ".edge_leaders_select_support");
				sdsl::util::init_support(edge_leaders_select_support, &edge_leaders);
				cerr << " done.\n";
			} else {
				cerr << "Indexes not found, computing...";
				init_pattern_matching_support();
				cerr << " done.\n";
				store_pattern_matching_support(basename);
			}
		}

		std::pair<int,int> locate_edge(size_type lex_rank) const
		{
			int startnode = node_leaders_rank_support(edge_index[lex_rank]+1)-1;
			int startnodeindex = node_leaders_select_support(startnode+1);
			int endnode = edges.at(startnode).at(edge_leaders_rank_support(edge_index[lex_rank]+1) - edge_leaders_rank_support(startnodeindex+1));
			return std::pair(startnode, endnode);
		}

		std::tuple<int,int,int> locate_edge_and_position(size_type lex_rank) const
		{
			int pos = edge_index[lex_rank]; // position in the text via SA
			int startnode = node_leaders_rank_support(pos+1)-1;
			int startnodeindex = node_leaders_select_support(startnode+1);
			int endnode = edges.at(startnode).at(edge_leaders_rank_support(pos+1) - edge_leaders_rank_support(startnodeindex+1));
			int edgepos = edge_leaders_select_support(edge_leaders_rank_support(pos+1)); // position in the text of the edge
			return std::tuple(startnode, endnode, pos - edgepos + ((is_source[startnode]) ? -1 : 0));
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

pair<vector<string>,vector<string>> read_patterns(std::ifstream &patternsfs)
{
	vector<string> pattern_ids, patterns;

	//TODO proper FASTA parsing
	string pattern;
	for (string line; getline(patternsfs, line);) {
		if (line.size() >= 1 and line[0] == '>') {
			pattern_ids.push_back(line.substr(1));
			if (pattern.size() > 0) {
				patterns.emplace_back(pattern);
				pattern.clear();
			}
		} else if (line.size() >= 1 and line[0] != '>') {
			pattern += line;
		}
	}
	if (pattern.size() > 0) {
		patterns.emplace_back(pattern);
		pattern.clear();
	}

	return pair<vector<string>, vector<string>>(pattern_ids, patterns);
}

void path_to_stream(std::ostream *out, Elasticfoundergraph &graph, string &path_id, vector<int> &path)
{
	*out << "P\t" << path_id << "\t";
	*out << graph.get_id(path[0]) << "+";
	for (int n = 1; n < path.size(); n++) {
		*out << "," << graph.get_id(path[n]) << "+";
	}
	*out << "\t*" << std::endl;
}

void path_to_stream(std::ostream *out, Elasticfoundergraph &graph, string &path_id, vector<vector<int>> &paths)
{
	for (int p = 0; p < paths.size(); p++) {
		*out << "P\t" << path_id << "-" << p+1 << "\t";
		*out << graph.get_id(paths[p][0]) << "+";
		for (int n = 1; n < paths[p].size(); n++) {
			*out << "," << graph.get_id(paths[p][n]) << "+";
		}
		*out << "\t*" << std::endl;
	}
}

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
		GAFAnchor(istringstream &descr, Elasticfoundergraph &efg)
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

		bool check(const Elasticfoundergraph &efg, const string &pattern) const
		{
			assert(path.size() > 0);
			assert(orientations.at(0)); // only + nodes, for now

			if (qend - qstart != pend - pstart)
				return false;

			string pathsequence;
			for (auto &n : path) {
				pathsequence += efg.get_label(n);
			}

			return (pattern.substr(qstart, qend - qstart) == pathsequence.substr(pstart, pend - pstart));
		}

		bool operator<(const GAFAnchor &a) const
		{
			// TODO: does this order matter?
			return tie(qstart, pstart, qend, pend, path) <
				tie(a.qstart, a.pstart, a.qend, a.pend, a.path);
		}

		bool operator==(const GAFAnchor &a) const
		{
			// TODO: does this order matter?
			return tie(qstart, pstart, qend, pend, path) ==
				tie(a.qstart, a.pstart, a.qend, a.pend, a.path);
		}

		int get_query_start() const
		{
			return qstart;
		}

		string get_query_id() const
		{
			return qname;
		}

		int get_path_length() const
		{
			return path.size();
		}

		int get_length() const
		{
			return qend - qstart;
		}

		int start_distance_query(GAFAnchor &a) const
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

void anchors_to_stream_split_single(std::ostream *out, const Elasticfoundergraph &graph, const vector<GAFAnchor> &matches, const bool keepedgematches)
{
	// matches can be their rev_ version
	vector<GAFAnchor> buffer;
	long currentqstart = -1;
	bool currentreverse = false;
	for (long i = 0; i < matches.size(); i++) {
		if (matches[i].get_query_start() == currentqstart and
				currentreverse == (matches[i].get_query_id().substr(0,4) == "rev_")) {
			if (keepedgematches and matches[i].get_path_length() <= 2) {
				buffer.push_back(matches[i]);
			} else {
				for (GAFAnchor n : matches[i].split_single(graph))
					buffer.push_back(std::move(n));
			}
		} else {
			std::sort(buffer.begin(), buffer.end());
			for (int n = 0; n < (int)buffer.size() - 1; n++) {
				if (buffer[n] != buffer[n+1])
					*out << buffer[n].to_string(graph) << "\n";
			}
			if (buffer.size() > 0)
				*out << buffer.back().to_string(graph) << "\n";
			buffer.clear();

			currentqstart = matches[i].get_query_start();
			currentreverse = (matches[i].get_query_id().substr(0,4) == "rev_");
			if (keepedgematches and matches[i].get_path_length() <= 2) {
				buffer.push_back(matches[i]);
			} else {
				for (GAFAnchor n : matches[i].split_single(graph))
					buffer.push_back(std::move(n));
			}
		}
	}
	std::sort(buffer.begin(), buffer.end());
	for (int n = 0; n < (int)buffer.size() - 1; n++) {
		if (buffer[n] != buffer[n+1])
			*out << buffer[n].to_string(graph) << "\n";
	}
	if (buffer.size() > 0)
		*out << buffer.back().to_string(graph) << "\n";
}

void anchors_to_stream_split_single_graphaligner(std::ostream *out, const Elasticfoundergraph &graph, vector<GAFAnchor> const &matches, const bool keepedgematches)
{
	// TODO: remove duplicate anchors?
	for (GAFAnchor m : matches) {
		if (keepedgematches and m.get_path_length() <= 2) {
			*out << m.to_string(graph) << "\n";
		} else {
			for (GAFAnchor n : m.split_single_graphaligner(graph)) {
				*out << n.to_string(graph) << "\n";
			}
		}
	}
}

} // Namespace efg_locate

#endif
