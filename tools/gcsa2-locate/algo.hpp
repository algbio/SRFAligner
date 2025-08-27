#ifndef ALGO_HPP
#define ALGO_HPP

#include <fstream>
#include <string>
#include <chrono>
#include <tuple>
#include <atomic>
#include <gcsa/algorithms.h>

#include "concurrentqueue.h" // https://github.com/cameron314/concurrentqueue
#include "gafanchor.hpp"
#include "gcsa2-locate.hpp" // Params

//#define ALGO_HPP_DEBUG

using std::pair, std::string, graph::SimpleGraph, gcsa2_locate::Params, gafanchor::GAFAnchor, std::cerr, std::endl;

moodycamel::ConcurrentQueue<pair<string,string>> readqueue { 50, 64, 64 }; // TODO: parameterize this!
moodycamel::ConcurrentQueue<string*> outputqueue { 50, 64, 64 }; // TODO: parameterize this!

// TODO: find a better place for these
char complement(const char n)
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
	return 'N';
}
string reverse_complement(const string &s)
{
	string reverse_pattern(s.size(),0);
	for (long unsigned int i = 0; i < s.size(); i++)
		reverse_pattern[i] = complement(s[s.size()-1-i]);
	return reverse_pattern;
}

void reader_worker(std::ifstream &patternsfs, std::atomic<bool> &input_done)
{
	//TODO proper FASTA parsing
	string patternid, pattern;
	for (string line; getline(patternsfs, line);) {
		if (line.size() >= 1 and line[0] == '>') {
			if (pattern.size() > 0) {
				pair<string, string> p(patternid, std::move(pattern));
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
		pair<string, string> p(patternid, std::move(pattern));
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

void inline output_gcsa_matches(const SimpleGraph &graph, const string &pattern, const string &pattern_id, const long unsigned int qstart, const long unsigned int qend, const vector<gcsa::node_type> &positions, vector<GAFAnchor> &match, bool splitoutputmatchesgraphaligner = false)
{
	for (auto const &n : positions) {
		const auto node = gcsa::Node::id(n);
		const auto offset = gcsa::Node::offset(n);
		const auto graph_node_length = graph.get_label_length(node);
		const auto match_length = std::min(qend - qstart + 1, graph_node_length - offset);
		if (splitoutputmatchesgraphaligner and match_length == 1) continue;
		match.push_back(GAFAnchor(
					pattern_id,
					pattern.length(),
					qstart,
					qstart + match_length,
					vector<long unsigned int>({ node }),
					graph_node_length,
					offset,
					offset + match_length,
					gcsa::Node::rc(n)
				));
	}
}

int gcsa_backward_search(
		const SimpleGraph &graph,
		const gcsa::GCSA &index,
		const string &query_id,
		const string &query,
		const Params &params,
		vector<GAFAnchor> &matches)
{
	matches.clear();

	gcsa::range_type r = index.find(query);
	if (gcsa::Range::empty(r)) {
		return 0;
	} else {
		vector<gcsa::node_type> positions;
		index.locate(r, positions);
		output_gcsa_matches(graph, query, query_id, 0, query.size() - 1, positions, matches, params.splitoutputmatchesgraphaligner);
		return matches.size();
	}
}

void exact_worker(const SimpleGraph &graph, const gcsa::GCSA &index, Params &params, std::atomic<bool> &input_done)
{
	vector<GAFAnchor> matches;
	while (true) {
		pair<string, string> p;
		while (readqueue.try_dequeue(p)) {
			const string pattern_id_rev = ((params.reversecompl) ? p.first : "");
			if (gcsa_backward_search(graph, index, p.first, p.second, params, matches) != 0) {
				std::stringstream localoutput;
				anchors_to_stream(&localoutput, matches);
				string *s = new string(localoutput.str());
				while (!outputqueue.try_enqueue(s))
					std::this_thread::sleep_for(std::chrono::milliseconds(10));
			}
			if (params.reversecompl and gcsa_backward_search(graph, index, pattern_id_rev, reverse_complement(p.second), params, matches) != 0) {
				std::stringstream localoutput;
				for (auto &m : matches)
					m.reverse();
				anchors_to_stream(&localoutput, matches);
				string *s = new string(localoutput.str());
				while (!outputqueue.try_enqueue(s))
					std::this_thread::sleep_for(std::chrono::milliseconds(10));
			}
		}
		if (input_done)
			break;
		else
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
	}
}

void gcsa_backward_search_greedy(
		const gcsa::GCSA &index, 
		const string &query,
		int &q,
		gcsa::range_type &range)
{
	assert(q >= 0);

	gcsa::range_type new_range = index.charRange(index.alpha.char2comp[query[q]]);
	while (q >= 0) {
		const unsigned long int lex_length = gcsa::Range::length(new_range);
		if (lex_length == 1) {
			range = new_range;
			q -= 1;
			break;
		} else if (lex_length >= 1) {
			range = new_range;
			q -= 1;
			if (q >= 0)
				new_range = index.LF(range, index.alpha.char2comp[query[q]]);
		} else {
			break;
		}
	}
}

struct exactedgematch {
	// compact description of substring Q[qstart..qend] (0-indexed) in graph edges
	int qstart;
	int qend;
	gcsa::range_type range; // GCSA lex range
};

void update_longest_match(vector<exactedgematch> &v, int qstart, int qend, gcsa::range_type range)
{
	assert(v.size() > 0);
	int minpos = 0;
	int minqstart = v[0].qstart;
	int minqend = v[0].qend;

	for (long unsigned int i = 1; i < v.size(); i++) {
		if (v[i].qend - v[i].qstart < minqend - minqstart) {
			minpos = i;
			minqstart = v[i].qstart;
			minqend = v[i].qend;
		}
	}

	if (qend - qstart > minqend - minqstart) {
#ifdef ALGO_HPP_DEBUG
		std::cerr << "DEBUG: substituting longest match Q[" << v[minpos].qstart << ".." << v[minpos].qend << "] ([" << v[minpos].range.first << ".." << v[minpos].range.second << " with longer match Q[" << qstart << ".." << qend << "] ([" << range.first << ".." << range.second << "\n";
#endif
		v[minpos].qstart = qstart;
		v[minpos].qend = qend;
		v[minpos].range = range;
	}
}

int approx_gcsa_backward_search(
		const SimpleGraph &graph,
		const gcsa::GCSA &index,
		const string &query_id,
		const string &query,
		const Params &params,
		vector<GAFAnchor> &matches)
{
	matches.clear();
	vector<exactedgematch> longest_matches;
	longest_matches.reserve(params.edgelongestcount);
	for (int i = 0; i < params.edgelongestcount; i++)
		longest_matches.push_back(exactedgematch({0, 0, gcsa::Range::empty_range()}));

	// find semi-repeat-free-equivalent matches
	int q = query.length() - 1;
	while (q >= 0) {
		while (params.ignorechars.find(query[q]) != std::string::npos)
			q--;
		if (q < 0)
			break;

		int startq = q;
		gcsa::range_type res;
		gcsa_backward_search_greedy(index, query, q, res);
#ifdef ALGO_HPP_DEBUG
				cerr << "DEBUG: greedy match processed Q[" << q+1 << ".." << startq << "]";
				cerr << " and found a string with lex range of size " << gcsa::Range::length(res) << "." << endl;
#endif
		if (gcsa::Range::length(res) == 1) {
			// unique occurrence in GCSA index
			vector<gcsa::node_type> pos;
			index.locate(res, pos);
#ifdef ALGO_HPP_DEBUG
				cerr << "lex range corresponds to " << pos.size() << " positions." << endl;
#endif
			if (pos.size() == 1) {
				output_gcsa_matches(graph, query, query_id, q+1, startq, pos, matches, params.splitoutputmatchesgraphaligner);
			}
		} else if (gcsa::Range::length(res) > 1 and params.edgelongestcount > 0 and gcsa::Range::length(res) < params.edgelongestcountmax) {
			update_longest_match(longest_matches, q + 1, startq, res);
		}
	}
	// output longest non-repeat-free matches
	if (params.edgelongestcount > 0) {
		for (exactedgematch &m : longest_matches) {
			if ((!gcsa::Range::empty(m.range)) and (gcsa::Range::length(m.range) <= params.edgelongestcountmax)) {
				vector<gcsa::node_type> pos;
				index.locate(m.range, pos);
				if (pos.size() <= params.edgelongestcountmax) {
					output_gcsa_matches(graph, query, query_id, m.qstart, m.qend, pos, matches, params.splitoutputmatchesgraphaligner);
				}
			}
		}
	}

	if (params.reversecompl) {
		vector<GAFAnchor> reversecompl_matches;
		assert(longest_matches.size() == (unsigned long int)params.edgelongestcount);
		for (int i = 0; i < params.edgelongestcount; i++)
			longest_matches.push_back(exactedgematch({0, 0, gcsa::Range::empty_range()}));

		string reverse_query(query.size(),0);
		const string query_id_rev = query_id;
		for (long unsigned int i = 0; i < query.size(); i++)
			reverse_query[i] = complement(query[query.size()-1-i]);
		int q = reverse_query.length() - 1;
		while (q >= 0) {
			while (params.ignorechars.find(query[q]) != std::string::npos)
				q--;
			if (q < 0)
				break;

			int startq = q;
			gcsa::range_type res;
			gcsa_backward_search_greedy(index, reverse_query, q, res);
			if (gcsa::Range::length(res) == 1) {
				// unique occurrence in GCSA index
				vector<gcsa::node_type> pos;
				index.locate(res, pos);
				assert(pos.size() == 1);
				output_gcsa_matches(graph, reverse_query, query_id_rev, q+1, startq, pos, reversecompl_matches, params.splitoutputmatchesgraphaligner);
			} else if (gcsa::Range::length(res) > 1 and params.edgelongestcount > 0) {
				update_longest_match(longest_matches, q + 1, startq, res);
			}
		}
		// output longest non-repeat-free matches
		if (params.edgelongestcount > 0) {
			for (exactedgematch &m : longest_matches) {
				if ((!gcsa::Range::empty(m.range)) and (gcsa::Range::length(m.range) <= params.edgelongestcountmax)) {
					vector<gcsa::node_type> pos;
					index.locate(m.range, pos);
					output_gcsa_matches(graph, reverse_query, query_id_rev, m.qstart, m.qend, pos, reversecompl_matches, params.splitoutputmatchesgraphaligner);
				}
			}
		}
		// merge matches
		for (auto &m : reversecompl_matches)
			m.reverse();

		matches.insert(matches.end(), reversecompl_matches.begin(), reversecompl_matches.end());
	}
	return matches.size();
}

void approx_worker(const SimpleGraph &graph, const gcsa::GCSA &index, Params &params, std::atomic<bool> &input_done)
{
	vector<GAFAnchor> matches;
	while (true) {
		std::pair<std::string, std::string> p;
		while (readqueue.try_dequeue(p)) {
			if (approx_gcsa_backward_search(graph, index, p.first, p.second, params, matches) != 0) {
				std::stringstream localoutput;
				anchors_to_stream(&localoutput, matches);
				std::string *s = new std::string(localoutput.str());
				while (!outputqueue.try_enqueue(s))
					std::this_thread::sleep_for(std::chrono::milliseconds(10));
			}
		}
		if (input_done)
			break;
		else
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
	}
}


#endif // ALGO_CPP
