#ifndef ALGO_HPP
#define ALGO_HPP

#include <vector>
#include <iostream>
#include <string>
#include <sdsl/csa_wt.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <syncstream>
#include <thread>
#include <unordered_set>
#include <unordered_map>
#include "concurrentqueue.h" // https://github.com/cameron314/concurrentqueue
#include <algorithm>

#include "efg-locate.hpp" // input parameters (Param)
#include "efg.hpp"

//#define ALGO_DEBUG

using std::vector, std::cerr, std::endl, sdsl::backward_search, std::tie;
typedef sdsl::csa_wt<>::size_type size_type;

std::mutex mutp, mutoutput, mutcerr;
moodycamel::ConcurrentQueue<std::string*> outputqueue { 50, 64, 64 }; // TODO: parameterize this!
moodycamel::ConcurrentQueue<std::pair<std::string,std::string>> readqueue { 50, 64, 64 }; // TODO: parameterize this!

namespace efg_locate {

// TODO: find a better place for this
char complement(char n)
{
	switch(n)
	{   
		case 'A':
			return 'T';
		case 'T':
			return 'A';
		case 'G':
			return 'C';
		case 'C':
			return 'G';
		case 'N':
			return 'N';
	}   
	assert(false);
	return ' ';
}   

/*
 * find longest suffix pattern[f..q] of pattern[0..q] that is prefix of some
 * edge label l(u)l(v). If f exists, then f_l and f_r are the lex range of
 * pattern[f..q] in the edge index, otherwise f = -1. lastq_l lastq_r are the lex range of
 * pattern[q'+1..q] where q' is the final state of q': pattern[q'+1..q] is
 * longest suffix of pattern[0..q] that does not contain chars in ignorechars
 * and occurs in some edge label
 */
void inline first_search(
	const Elasticfoundergraph &efg,
	const string &pattern,
	int &q,
	size_type &lastq_l,
	size_type &lastq_r,
	int &f,
	size_type &f_l,
	size_type &f_r,
	const string &ignorechars = "")
{
	const sdsl::csa_wt<> &edge_index = efg.edge_index;
	lastq_l = 0; // lex bounds in edge_index
	lastq_r = edge_index.size() - 1;
	size_type l_res, r_res; // temporary results

	f = -1;
	while (q >= 0 and ignorechars.find(pattern[q]) == std::string::npos) {
		const int res = backward_search(edge_index, lastq_l, lastq_r, pattern[q], l_res, r_res);

		if (res == 0) { // no match
			break;
		} else {
			lastq_l = l_res;
			lastq_r = r_res;
			q -= 1;

			// test if we can read the edge separator char '#' from [lastq_l..lastq_r]
			if (backward_search(edge_index, lastq_l, lastq_r, '#', l_res, r_res) != 0) {
				f = q + 1;
				f_l = lastq_l;
				f_r = lastq_r;
			}
		}
	}
}
/*
 * version of the above first_search that additionally finds the shortest suffix
 * pattern[qq..q] of pattern[0..q] occurring at most count times in the
 * edge index. If such suffix exists, qq_l and qq_r correspond to the lex
 * range of pattern[qq..q], otherwise qq = -1
 */
void inline first_search(
	const Elasticfoundergraph &efg,
	const string &pattern,
	const int count,
	int &q,
	size_type &lastq_l,
	size_type &lastq_r,
	int &f,
	size_type &f_l,
	size_type &f_r,
	int &qq,
	size_type &qq_l,
	size_type &qq_r,
	const string &ignorechars = "")
{
	const sdsl::csa_wt<> &edge_index = efg.edge_index;
	lastq_l = 0; // lex bounds in edge_index
	lastq_r = edge_index.size() - 1;
	size_type l_res, r_res; // temporary results

	f = -1;
	qq = -1;
	while (q >= 0 and ignorechars.find(pattern[q]) == std::string::npos) {
		const int res = backward_search(edge_index, lastq_l, lastq_r, pattern[q], l_res, r_res);

		if (res == 0) { // no match
			break;
		} else {
			lastq_l = l_res;
			lastq_r = r_res;
			q -= 1;

			if (qq == -1 and lastq_l <= lastq_r and count >= 1 and lastq_r - lastq_l + 1 <= count) {
				qq = q + 1;
				qq_l = lastq_l;
				qq_r = lastq_r;
			}

			// test if we can read the edge separator char '#' from [lastq_l..lastq_r]
			if (backward_search(edge_index, lastq_l, lastq_r, '#', l_res, r_res) != 0) {
				f = q + 1;
				f_l = lastq_l;
				f_r = lastq_r;
			}
		}
	}
}

/*
 * find ONE match in iEFG efg of the longest suffix of pattern[0..q] that does
 * not contain any char of ignorechars and that ends at a node boundary, i.e. it
 * is a suffix of l(u_1)...l(u_k) with u_1, ... , u_k a path in efg. After the
 * execution q holds value q' such that pattern[q'+1..q] is such suffix, path
 * contains nodes u_1, ..., u_k-1, first_node_pos is such that
 * u_1[first_node_pos..] is the matched part of u_1, y is such that
 * pattern[y..q] corresponds to the matched part of u_k.
 * If no such suffix exists, q' = q, y = -1, and first_node_pos = -1.
 *
 * Prerequisites:
 * efg.init_pattern_matching_support() must have been called
 *
 * Completeness:
 * to simplify edge cases, we assume that pattern[0..q] cannot occur in efg
 * starting from a source of the graph by adding a supersource in efg with a
 * dummy node label in init_pattern_matching_support() TODO
 */
void inline simple_search(
	const Elasticfoundergraph &efg,
	const string &pattern,
	int &q,
	vector<int> &path,
	int &first_node_pos,
	int &y,
	int &u_k,
	const string &ignorechars = "")
{
	const sdsl::csa_wt<> &edge_index = efg.edge_index;
	const int startq = q;

	path.clear(); // path will be filled in backwards (thus reverse) order
	size_type lastq_l = 0; // lex bounds in edge_index
	size_type lastq_r = edge_index.size() - 1;
	size_type l_res, r_res; // temporary results
	first_node_pos = -1;
	y = -1;
	int last_edge_boundary = -1;

	backward_search(edge_index, 0, edge_index.size() - 1, '#', lastq_l, lastq_r);
	while (q >= 0 and ignorechars.find(pattern[q]) == std::string::npos) {
		const int res = backward_search(edge_index, lastq_l, lastq_r, pattern[q], l_res, r_res);
		if (res == 0) {
#ifdef ALGO_DEBUG
			cerr << "DEBUG: simple search failed at index " << q << " and after collecting " << path.size() << "(+1) nodes\n";
#endif
			break;
		}

		lastq_l = l_res;
		lastq_r = r_res;
		q -= 1;

		// test if we can read the edge separator char '#' from [lastq_l..lastq_r]
		if (backward_search(edge_index, lastq_l, lastq_r, '#', l_res, r_res) != 0) {
			// full edge matched, first node is thus unique
			int startnode, endnode;
			tie(startnode, endnode) = efg.locate_edge(lastq_l);
#ifdef ALGO_DEBUG
			cerr << "DEBUG: read full edge " << efg.get_label(startnode) << " -> " << efg.get_label(endnode) << "in pattern[" << q+1 << ".." << startq << "]\n";
#endif

			path.push_back(startnode);
			last_edge_boundary = q + 1;
			q += efg.get_label_length(startnode);
			if (y == -1) {
				u_k = endnode;
				y = q + 1;
			}
			backward_search(edge_index, 0, edge_index.size() - 1, '#', lastq_l, lastq_r);
		}
	}

	if (q == startq) {
		// no suffix of pattern[0..q] ends at some node label l(u)
		assert(y == -1);
		assert(first_node_pos == -1);
		return;
	}

	if (path.size() == 0) {
		// no full edge was identified
		int startnode, endnode, position;
		tie(startnode, endnode, position) = efg.locate_edge_and_position(lastq_l);

		if (position < efg.get_label_length(startnode)) {
			// full node endnode was matched
			y = q + 1 + efg.get_label_length(startnode) - position;
			path.push_back(startnode);
			u_k = endnode;
			first_node_pos = position;
		} else {
			// no full node was matched
			y = q + 1;
		}
	} else if (q + 1 < last_edge_boundary) { // path.size() > 0
		// some full edges were identified, and pattern[q+1..qstart] starts with a suffix of l(u_1)
		int startnode, endnode, position;
		tie(startnode, endnode, position) = efg.locate_edge_and_position(lastq_l);

		assert(path.back() == endnode);
		path.push_back(startnode);
		first_node_pos = position;
	} else { // path.size() > 0 and q + 1 >= last_edge_boundary
		// pattern[q+1..qstart] starts with a full edge
		first_node_pos = 0;
		assert(q + 1 == last_edge_boundary);
	}

	std::reverse(path.begin(), path.end());
}

/*
 * find ONE connecting vertex v such that l(v) = pattern[f..x] is a proper prefix of pattern[f..qstart]
 * and suffix of pattern[y..x]; if one such vertex exists, append the resulting nodes u v w to path and store in first_node_pos the matched part of l(u), and return value 1; otherwise return value 0 and first_node_pos contains an undefined value
 *
 * Prerequisites:
 * efg.init_pattern_matching_support() must have been called
 *
 * Completeness:
 * to simplify edge cases, we assume that pattern[0..q] cannot occur in efg
 * starting from a source of the graph by adding a supersource in efg with a
 * dummy node label in init_pattern_matching_support()
 */
int inline find_connecting_vertex(
	const Elasticfoundergraph &efg,
	const string &pattern,
	const int y,
	const int f,
	const int f_l,
	const int f_r,
	const int qstart,
	vector<int> &path,
	int &first_node_pos)
{
	const sdsl::csa_wt<> &edge_index = efg.edge_index;

	std::unordered_set<int> middle_block;
	std::unordered_map<int,int> middle_edges; // just one candidate last edge

	int q = qstart;
	size_type l = 0; // lex bounds in edge_index
	size_type r = edge_index.size() - 1;
	size_type l_res, r_res; // temporary results

	// if it is a semi-repeat-free match, all occurrences of pattern[f..qstart] must be prefix of some l(u)l(v)
	// TODO is this correct for sources?
	backward_search(edge_index, f_l, f_r, '#', l_res, r_res);
	if (f_r - f_l != r_res - l_res)
		return 0;

	for (int x = f; x < qstart; x++) {
		// 1. check if pattern[y..x] is suffix of some l(u)l(v)
		backward_search(edge_index, 0, edge_index.size() - 1, '#', l, r);
		q = x;
		while (q >= y and backward_search(edge_index, l, r, pattern[q], l_res, r_res) != 0) {
			l = l_res;
			r = r_res;
			q -= 1;
		}
		if (q != y - 1) // l(v) is not suffix of l(u)l(v) = pattern[y..x]
			continue;

		// 2. check if pattern[f..x] contains a full node label
		/*backward_search(edge_index, 0, edge_index.size() - 1, '#', l, r);
		q = x;
		while (q >= f and backward_search(edge_index, l, r, pattern[q], l_res, r_res) != 0) {
			l = l_res;
			r = r_res;
			q -= 1;
		}
		if (q != f - 1) // cannot be node full label
			continue;*/

		int startnode, endnode, position;
		tie(startnode, endnode, position) = efg.locate_edge_and_position(l);
		if (f - y != efg.get_label_length(startnode) - position) // pattern[f..x] is not full node label
			continue;
		const int v = endnode;
#ifdef ALGO_DEBUG
		cerr << "DEBUG: found candidate node " << efg.get_id(v) << "(pattern[" << y << ".." << x << "] as connecting vertex, testing\n";
#endif

		if (middle_block.size() == 0) {
			// find all nodes l(vv) such that l(vv) is prefix of pattern[f..qstart] (<= H^2)
			l = f_l;
			r = f_r;
			for (size_type k = l; k <= r; k++) {
				int snode, enode;
				tie(snode, enode) = efg.locate_edge(k);

				middle_block.insert(snode);
				if (snode == v and f + efg.get_label_length(v) - 1 < qstart) { // TODO investigate
					// connecting vertex is v
					path.push_back(startnode); // u
					first_node_pos = position;
					path.push_back(v);
					assert(f + efg.get_label_length(v) - 1 < qstart); // this edge case should not be possible
					path.push_back(enode); // w
					return 1;
				}

				if (!middle_edges.contains(snode)) {
					middle_edges[snode] = enode;
				}
			}
		} else if (middle_block.contains(v)) {
			path.push_back(startnode);
			first_node_pos = position;
			path.push_back(v);
			assert(f + efg.get_label_length(v) - 1 < qstart); // this edge case should not be possible
			path.push_back(middle_edges[v]);
			return 1;
		}
	}

	// no connecting vertex
	return 0;
}

/*
 * Find ONE occurrence of pattern in iEFG efg. If successful, return value is
 * >0 and path contains the node indices of such path, otherwise return value is
 * 0 and path is empty.
 *
 * Prerequisites:
 * efg.init_pattern_matching_support() must have been called
 *
 * Note: 
 * The implementation follow a "simplified" version of the algorithms in
 * https://doi.org/10.1016/j.tcs.2023.114269 working in
 * O(|Q| + min(|Q|,L)^2 + H^2) time, where L is the maximum node label
 * length and H is the maximum block height of the graph. 
 */
int efg_backward_search(Elasticfoundergraph &efg, const string &pattern, vector<int> &path)
{
	path.clear();
	int q = pattern.size() - 1, f;
	size_type lastq_l, lastq_r, f_l, f_r;

	// find largest suffix pattern[f..q] of pattern[..q] that is prefix of some edge label l(u)l(v)
	first_search(efg, pattern, q, lastq_l, lastq_r, f, f_l, f_r);

	// case 1a: pattern occurs in some edge label l(u)l(v)
	if (q == -1) {
		// pick arbitrary occurrence
		int startnode, endnode;
		int position;
		tie(startnode, endnode, position) = efg.locate_edge_and_position(lastq_l); 
#ifdef ALGO_DEBUG
		cerr << "DEBUG: occurrence found in edge " << efg.get_id(startnode) << "->" << efg.get_id(endnode) << endl;
#endif

		// trim first or last node
		if (position < efg.get_label_length(startnode))
			path.push_back(startnode);
		if (position + pattern.length() >= efg.get_label_length(endnode))
			path.push_back(endnode);

		return 1;
	}

	// case 1b: pattern does not occur in some edge and f is not well-defined
	if (f == -1)
		return 0;

#ifdef ALGO_DEBUG
	cerr << "DEBUG: pattern[f..] is pattern[" << f << "..]" << endl;
#endif

	int first_node_pos, y, u_k;
	q = f - 1;
	simple_search(efg, pattern, q, path, first_node_pos, y, u_k);

	// case 2a: no match of pattern[0..f-1] in efg ending at node boundary
	if (q != -1)
		return 0;

#ifdef ALGO_DEBUG
	cerr << "DEBUG: simple_search matched pattern[" << q+1 << ".." << f - 1 << "] and found nodes ";
	for (auto n : path) {
		cerr << efg.get_id(n) << ",";
	}
	cerr << ", first_node_pos is " << first_node_pos << ", y is " << y << "\n";
#endif

	if (find_connecting_vertex(efg, pattern, y, f, f_l, f_r, pattern.size()-1, path, first_node_pos) > 0) {
		// case 3a: pattern occurs in efg
		return 1;
	} else {
		// case 3b: no connecting vertex exists
		return 0;
	}
}

int efg_backward_search_ignorechars(Elasticfoundergraph &efg, const string &ignorechars, const string &pattern, vector<vector<int>> &paths)
{
	paths.clear();
	int startq = 0;
	for (int q = 0; q < pattern.size(); q++) {
		if (ignorechars.find(pattern[q]) != std::string::npos) {
			// search for current pattern without special characters, if any
			if (startq < q - 1) {
#ifdef ALGO_DEBUG
				cerr << "searching for pattern substring " << pattern.substr(startq, q - startq) << endl;
#endif
				vector<int> partial_path;
				if (efg_backward_search(efg, pattern.substr(startq, q - startq), partial_path) == 0)
					return 0; // not found
				paths.push_back(std::move(partial_path));
			}
			startq = q + 1;
		}
	}
	if (startq < pattern.size()) {
#ifdef ALGO_DEBUG
		cerr << "searching for pattern substring " << pattern.substr(startq) << endl;
#endif
		vector<int> partial_path;
		if (efg_backward_search(efg, pattern.substr(startq), partial_path) == 0)
			return 0; // not found
		paths.push_back(std::move(partial_path));
	}

	return 1;
}

void inline output_edge_count_matches(const Elasticfoundergraph &efg, const string &pattern, const string &pattern_id, const int qq, const int startq, const int qq_l, const int qq_r, vector<GAFAnchor> &match)
{
		// output all edge-count matches
		for (int i = qq_l; i <= qq_r; i++) {
			int startnode, endnode, position;
			tie(startnode, endnode, position) = efg.locate_edge_and_position(i);

			match.push_back(GAFAnchor(
						efg,
						pattern_id,
						pattern.size(),
						qq,
						startq+1,
						std::vector({startnode,endnode}),
						efg.get_label_length(startnode) + efg.get_label_length(endnode),
						position,
						position + (startq+1-qq)));
			//assert(match.back().check(efg, pattern)); // TODO implement faster sanity check
		}
}

/*
 * copy of efg_backward_search that matches pattern[0..q] until failure, saves result in a GAFAnchor, and updates q such that pattern[0..q] is the rest of the pattern that was not matched
 */
int efg_backward_search_greedy(const Elasticfoundergraph &efg, const string &pattern_id, const string &pattern, const int edgemincount, int &q, vector<GAFAnchor> &match, const string &ignorechars = "")
{
	const int startq = q;
	int f, qq;
	size_type lastq_l, lastq_r, f_l, f_r, qq_l, qq_r;

	// find largest suffix pattern[f..q] of pattern[..q] that is prefix of some edge label l(u)l(v)
	first_search(efg, pattern, edgemincount, q, lastq_l, lastq_r, f, f_l, f_r, qq, qq_l, qq_r, ignorechars);

	// case 1a: the whole pattern occurs in some edge label l(u)l(v)
	if (q == -1) {
		if (qq == -1) { // not a semi-repeat-free (probably) or edge-count match
			// policy 1: do nothing
			// policy 2: f
			//q = std::max(q, f - 1);
			return 0;
		}

		output_edge_count_matches(efg, pattern, pattern_id, qq, startq, qq_l, qq_r, match);
		// policy 1: do nothing
		// policy 2: f
		//q = std::max(q, f - 1);
		return 1;
	}

	// case 1b: pattern does not occur in some edge and f is not well-defined
	if (f == -1) {
		if (q == startq)
			q -= 1;
		if (qq == -1) { // not a semi-repeat-free (probably) or edge-count match
			// policy 1: do nothing
			// policy 2: f, do nothing in this case
			return 0;
		}
		//std::cerr << "case 1b and after is " << q;
		output_edge_count_matches(efg, pattern, pattern_id, qq, startq, qq_l, qq_r, match);
		// policy 1: do nothing
		// policy 2: f, do nothing in this case
		return 1;
	}

	int first_node_pos, y, u_k;
	vector<int> path;
	q = f - 1;
	simple_search(efg, pattern, q, path, first_node_pos, y, u_k, ignorechars);

	// case 2a: no suffix of pattern[0..f-1] is a suffix of some edge label l(u)l(v)
	if (q == f - 1) {
		if (qq == -1) { // not a semi-repeat-free (probably) or edge-count match
			// policy 1: do nothing
			// policy 2: f
			//q = std::max(q, f - 1);
			return 0;
		}

		output_edge_count_matches(efg, pattern, pattern_id, qq, startq, qq_l, qq_r, match);
		// policy 1: we consumed pattern[qq..]
		q = std::min(q, qq);
		// policy 2: f
		//q = std::max(q, f - 1);
		return 1;
	}

#ifdef ALGO_DEBUG
	cerr << "DEBUG: simple_search matched pattern[" << q+1 << ".." << f - 1 << "] and found nodes ";
	for (auto n : path) {
		cerr << efg.get_id(n) << ",";
	}
	cerr << ", first_node_pos is " << first_node_pos << ", y is " << y << "\n";
#endif

	if (path.size() == 0 and find_connecting_vertex(efg, pattern, y, f, f_l, f_r, startq, path, first_node_pos) > 0) {
		assert(first_node_pos >= 0);
		int pathlength = 0;
		for (int i : path)
			pathlength += efg.get_label_length(i);
		match.push_back(GAFAnchor(
					pattern_id,
					pattern.size(),
					q+1,
					startq+1,
					path,
					pathlength,
					first_node_pos,
					first_node_pos + (startq+1-(q+1))));
		//assert(match.back().check(efg, pattern));
		return 1;
	} else {
		int u_node_pos;
		if (find_connecting_vertex(efg, pattern, y, f, f_l, f_r, startq, path, u_node_pos) > 0) {
			assert(u_node_pos == 0);
			assert(first_node_pos >= 0);
			int pathlength = 0;
			for (int i : path)
				pathlength += efg.get_label_length(i);

			match.push_back(GAFAnchor(
						pattern_id,
						pattern.size(),
						q+1,
						startq+1,
						path,
						pathlength,
						first_node_pos,
						first_node_pos + (startq+1-(q+1))));
			//assert(match.back().check(efg, pattern));
			return 1;
		}
	}

	if (path.size() > 0) {
		// lucky full node match in simple_search
		path.push_back(u_k);
		int pathlength = 0;
		for (int i : path)
			pathlength += efg.get_label_length(i);
		match.push_back(GAFAnchor(
					pattern_id,
					pattern.size(),
					q+1,
					f,
					path,
					pathlength,
					first_node_pos,
					first_node_pos + (f-(q+1))));
		//assert(match.back().check(efg, pattern));
		return 1;
	}

	// no semi-repeat-free match was found
	if (qq == -1) { // not a semi-repeat-free (probably) or edge-count match
		// policy 1: do nothing
		// policy 2: f
		//q = std::max(q, f - 1);
		return 0;
	}

	output_edge_count_matches(efg, pattern, pattern_id, qq, startq, qq_l, qq_r, match);
	// policy 1: we consumed pattern[qq..]
	q = std::min(q, qq);
	// policy 2: f
	//q = std::max(q, f - 1);
	return 1;
}

int approx_efg_backward_search(const Elasticfoundergraph &efg, const string &pattern_id, const string &pattern, Params &params, vector<GAFAnchor> &matches)
{
	//TODO: skip coverage computation if not needed?
	long coverage = 0, rev_coverage = 0, full_node_matches = 0, rev_full_node_matches = 0;

#ifdef ALGO_DEBUG
	cerr << "Searching pattern " << pattern << endl;
#endif
	matches.clear();
	int q = pattern.size() - 1;
	while (q >= 0) {
		while (params.ignorechars.find(pattern[q]) != std::string::npos)
			q--;
		if (q < 0)
			break;

		int startq = q;
		vector<GAFAnchor> match;
#ifdef ALGO_DEBUG
		cerr << "q before greedy search is " << q << endl;
#endif
		int res = efg_backward_search_greedy(efg, pattern_id, pattern, ((params.edgemincountheuristic) ? 1 : params.edgemincount), q, match, params.ignorechars);
#ifdef ALGO_DEBUG
		cerr << "q after greedy search is " << q << endl;
#endif
		if (res > 0) {
			coverage += startq - q;
			for (auto &m : match)
				matches.push_back(m);
		}
	}

	vector<GAFAnchor> reversecompl_matches;
	if (params.reversecompl) {
		string reverse_pattern(pattern.size(),0);
		const string pattern_id_rev = ((params.renamereversecomplement) ? "rev_" + pattern_id : pattern_id);
#ifdef ALGO_DEBUG
		cerr << "Searching reverse complement pattern " << reverse_pattern << endl;
#endif
		for (int i = 0; i < pattern.size(); i++)
			reverse_pattern[i] = complement(pattern[pattern.size()-1-i]);
		q = pattern.size() - 1;
		while (q >= 0) {
			while (params.ignorechars.find(pattern[q]) != std::string::npos)
				q--;
			if (q < 0)
				break;
			int startq = q;
			vector<GAFAnchor> match;
			int res = efg_backward_search_greedy(efg, pattern_id_rev, reverse_pattern, ((params.edgemincountheuristic) ? 1 : params.edgemincount), q, match, params.ignorechars);
			if (res > 0) {
				rev_coverage += startq - q;
				for (auto &m : match) {
					if (!params.renamereversecomplement)
						m.reverse();
					reversecompl_matches.push_back(m);
				}
			}
		}
	}

	for (auto m : matches) {
		full_node_matches += std::max(0, m.get_path_length() - 2); // probably correct
	}
	for (auto m : reversecompl_matches) {
		rev_full_node_matches += std::max(0, m.get_path_length() - 2); // probably correct
	}
	if (params.reportstats) {
		std::osyncstream oss(std::cout);
		oss << pattern_id << "\t"; // pattern name
		oss << pattern.size() << "\t"; // length
		oss << ((long)coverage*100)/pattern.size() << "\t"; // coverage %
		oss << ((long)rev_coverage*100)/pattern.size() << "\t"; // reverse coverage %
		oss << full_node_matches << "\t"; // lower bound on full node matches
		oss << rev_full_node_matches << "\t"; // lower bound on full reverse node matches
		oss << "\n";
	}

	if ((matches.size() == 0 and reversecompl_matches.size() == 0) or 
			(((long)coverage*100)/pattern.size() < params.mincoverage) && (((long)rev_coverage*100)/pattern.size() < params.mincoverage))
		return 0;
	// DO NOT TRUST EXACT MATCHES THAT MUCH
	/*if (((long)coverage*100)/pattern.size() < mincoverage)
		matches.clear();
	if (((long)rev_coverage*100)/pattern.size() >= mincoverage)*/
	std::reverse(matches.begin(), matches.end()); // to have forward (and reverse) matches sorted
	if (params.renamereversecomplement)
		std::reverse(reversecompl_matches.begin(), reversecompl_matches.end());
	matches.insert(matches.end(), reversecompl_matches.begin(), reversecompl_matches.end());

	reversecompl_matches.clear();
	if (params.edgemincountheuristic and full_node_matches == 0 and rev_full_node_matches == 0) {
#ifdef ALGO_DEBUG
		cerr << "Did not find long matches, actually computing edge matches with min count equal to " << params.edgemincount << endl;
#endif
		q = pattern.size() - 1;
		while (q >= 0) {
			int startq = q;
			vector<GAFAnchor> match;
			int res = efg_backward_search_greedy(efg, pattern_id, pattern, params.edgemincount, q, match);
			if (res > 0) {
				for (auto &m : match)
					matches.push_back(m);
			}
		}
		if (params.reversecompl) {
			string reverse_pattern(pattern.size(),0);
			const string pattern_id_rev = ((params.renamereversecomplement) ? "rev_" + pattern_id : pattern_id);
			for (int i = 0; i < pattern.size(); i++)
				reverse_pattern[i] = complement(pattern[pattern.size()-1-i]);
			q = pattern.size() - 1;
			while (q >= 0) {
				int startq = q;
				vector<GAFAnchor> match;
				int res = efg_backward_search_greedy(efg, pattern_id_rev, reverse_pattern, params.edgemincount, q, match);
				if (res > 0) {
					rev_coverage += startq - q;
					for (auto &m : match) {
						if (!params.renamereversecomplement)
							m.reverse();
						reversecompl_matches.push_back(m);
					}
				}
			}
		}
		matches.insert(matches.end(), reversecompl_matches.begin(), reversecompl_matches.end());
	}

	return 1;
}

int approx_efg_backward_search_ignorechars(const Elasticfoundergraph &efg, const string &ignorechars, const string &pattern_id, Params &params, const string &pattern, vector<GAFAnchor> &matches)
{
	matches.clear();
	vector<string> maximal_substrings;
	int startq = 0;
	for (int q = 0; q < pattern.size(); q++) {
		if (ignorechars.find(pattern[q]) != std::string::npos) {
			// search for current pattern without special characters, if any
			if (startq < q - 1) {
				maximal_substrings.push_back(pattern.substr(startq, q - startq));
			}
			startq = q + 1;
		}
	}
	if (startq < pattern.size()) {
		maximal_substrings.push_back(pattern.substr(startq));
	}

	for (int s = 0; s < maximal_substrings.size(); s++) {
		vector<GAFAnchor> substring_matches;
		int res = approx_efg_backward_search(efg, pattern_id + "-" + std::to_string(s), maximal_substrings[s], params, substring_matches);
		//TODO: more efficient move?
		if (res > 0)
			matches.insert(matches.end(), substring_matches.begin(), substring_matches.end());
	}

	if (matches.size() > 0)
		return 1;
	else
		return 0;
}

void reader_worker(std::ifstream &patternsfs, std::atomic<bool> &input_done)
{
	//TODO proper FASTA parsing
	string patternid, pattern;
	for (string line; getline(patternsfs, line);) {
		if (line.size() >= 1 and line[0] == '>') {
			if (pattern.size() > 0) {
				std::pair<std::string, std::string> p(patternid, std::move(pattern));
				while (!readqueue.try_enqueue(p))
					std::this_thread::sleep_for(std::chrono::milliseconds(10));
				pattern.clear();
			}
			patternid = line.substr(1);
		} else if (line.size() >= 1 and line[0] != '>') {
			pattern += line;
		}
	}
	if (pattern.size() > 0) {
		std::pair<std::string, std::string> p(patternid, std::move(pattern));
		while (!readqueue.try_enqueue(p))
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
	}

	input_done = true;
}

void writer_worker(std::atomic<bool>& done, Params &params)
{
	string *ptr;
	while (true) {
		while (outputqueue.try_dequeue(ptr)) {
			params.outputfs << *ptr;
			delete(ptr);
		}
		// queue is empty?
		if (!done)
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
		else
			break;
	}
}

void approx_worker(const Elasticfoundergraph &graph, const vector<string> &pattern_ids, const vector<string> &patterns, Params &params, std::atomic<bool> &input_done)
{
	std::osyncstream oss(cerr);
	int currentp;
	vector<GAFAnchor> matches;
	while (true) {
		std::pair<std::string, std::string> p;
		while (readqueue.try_dequeue(p)) {
			/*{
			  std::scoped_lock lck {mutp};
			  if (p == patterns.size())
			  return;

			  currentp = p++;
			  }*/

			if (approx_efg_backward_search(graph, p.first, p.second, params, matches) != 0) {
				//std::scoped_lock lck {mutoutput};
				std::stringstream localoutput;
				if (params.splitoutputmatches)
					//anchors_to_stream_split_single(&params.outputfs, graph, matches);
					anchors_to_stream_split_single(&localoutput, graph, matches);
				else if (params.splitoutputmatchesgraphaligner)
					//anchors_to_stream_split_single_graphaligner(&params.outputfs, graph, matches);
					anchors_to_stream_split_single_graphaligner(&localoutput, graph, matches);
				else
					//anchors_to_stream(&params.outputfs, graph, matches);
					anchors_to_stream(&localoutput, graph, matches);
				std::string *s = new std::string(localoutput.str());
				while (!outputqueue.try_enqueue(s))
					std::this_thread::sleep_for(std::chrono::milliseconds(10));
			} else {
				//std::scoped_lock lck {mutcerr};
				oss << "Cannot find any semi-repeat-free match of " << p.first << "\n";
			}
		}
		if (input_done)
			break;
		else
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
	}
}

} // Namespace efg_locate

#endif
