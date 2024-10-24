//code adapted from https://github.com/at-cg/ChainX
//Chirag Jain, Daniel Gibney and Sharma Thankachan. Algorithms for Colinear Chaining with Overlaps and Gap Costs. Journal of Computational Biology, 2022
//license: https://www.apache.org/licenses/LICENSE-2.0.txt

#ifndef CHAIN_HPP 
#define CHAIN_HPP

#include "efg.hpp"

//#define CHAIN_HPP_DEBUG

using std::vector, std::swap;

namespace chainx_block_graph {
/**
 * Compute optimal chain based on anchor-restricted edit distance using
 * strong precedence criteria optimized to run faster using engineering
 * trick(s), comparison mode: global.
 * We assume the anchors are sorted by the starting positions in the
 * linear text and that there are two dummy anchors marking the
 * beginning and end of the references.
 * graph.init_eds_support() must have been called before this function
 **/
vector<GAFHit> chain_global_eds(vector<GAFHit> &anchors, const Elasticfoundergraph &graph, const int initial_guess, const double ramp_up_factor, Stats &stats, const bool removesol = false)
{
	//graph.init_eds_support();

	int n = anchors.size();
	vector<int> costs(n, 0);
	vector<int> backtrack(n, 0);

	int bound_redit = initial_guess; //distance assumed to be <= initial_guess
	int revisions = 0;
	//with this assumption on upper bound of distance, a gap of >bound_redit will not be allowed between adjacent anchors

	while (true) {
		int inner_loop_start = 0;

		for(int j=1; j<n; j++) {
			//compute cost[j] here
			int find_min_cost = std::numeric_limits<int>::max();
			int backtrack_min_cost = std::numeric_limits<int>::max();

			// anchor i < anchor j 

			while (anchors[inner_loop_start].gap_query(anchors[j]) > bound_redit)
				inner_loop_start++;

			for(int i=j-1; i>=inner_loop_start; i--) {
				if (costs[i] < std::numeric_limits<int>::max() and are_colinear_eds(anchors[i], anchors[j], graph)) {
					int g = max_gap_eds(anchors[i], anchors[j], graph);
					int o = overlap_eds(anchors[i], anchors[j], graph);
#ifdef CHAIN_HPP_DEBUG
					std::cerr << "anchors[" << i << "] -> anchors[" << j << "]: g = " << g << ", o = " << o << std::endl;
#endif
					if (costs[i] + g + o < find_min_cost) {
						find_min_cost = costs[i] + g + o;
						backtrack_min_cost = i;
					}
				}
			}
			//save optimal cost at offset j
			costs[j] = find_min_cost;
			backtrack[j] = backtrack_min_cost;
		}

		if (costs[n-1] > bound_redit) {
			bound_redit = bound_redit * ramp_up_factor;
			revisions++;
		} else {
			break;
		}
	}

#ifdef CHAIN_HPP_DEBUG
		std::cerr << "Cost DP array = ";
		for (auto c : costs)
			std::cerr << c << " ";
		std::cerr << std::endl;
		std::cerr << "Backtrack array = ";
		for (auto b : backtrack)
			std::cerr << b << " ";
		std::cerr << std::endl;
		std::cerr << "Chaining cost computed " << revisions + 1 << " times" << "\n";
#endif

	//TODO: consider freeing here the space of costs array
	// backtrack optimal solution
	vector<GAFHit> solution;
	for (int j = backtrack[n - 1]; j > 0; j = backtrack[j])
	{
		solution.push_back(anchors[j]);
	}
	std::reverse(solution.begin(), solution.end());
	if (removesol)
	{
		vector<GAFHit> newanchors;
		newanchors.reserve(anchors.size() - solution.size());
		for (int i = 0, j = 0; i < anchors.size(); i += 1)
		{
			if (j < solution.size() and anchors[i] == solution[j])
			{
				j += 1;
			} else {
				newanchors.push_back(anchors[i]);
			}
		}
		swap(anchors, newanchors);
	}

	//std::cout << "distance = " << costs[n-1] << std::endl;
	if (solution.size() > 0) {
		stats.maxiterations = std::max(stats.maxiterations, revisions);
		stats.miniterations = std::min(stats.miniterations, revisions);
		stats.totaliterations += revisions;
		stats.maxcost = std::max(stats.maxcost, costs[n-1]);
		stats.mincost = std::min(stats.mincost, costs[n-1]);
		stats.totalcost += costs[n-1];
		const double relativecost = (double)costs[n-1] / solution.at(0).get_query_length();
		stats.maxrelativecost = std::max(stats.maxrelativecost, relativecost);
		stats.minrelativecost = std::min(stats.minrelativecost, relativecost);
		stats.totalrelativecost += relativecost;
	}

	return solution;
}

/**
 * See chain_global_eds, comparison mode: semiglobal.
 **/
vector<GAFHit> chain_semiglobal_eds(vector<GAFHit> &anchors, const Elasticfoundergraph &graph, const int initial_guess, const double ramp_up_factor, Stats &stats, const bool removesol = false)
{
	//graph.init_eds_support();

	int n = anchors.size();
	vector<int> costs(n, 0);
	vector<int> backtrack(n, 0);

	int bound_redit = initial_guess; //distance assumed to be <= initial_guess
	int revisions = 0;
	//with this assumption on upper bound of distance, a gap of >bound_redit will not be allowed between adjacent anchors

	while (true) {
		int inner_loop_start = 0;

		for(int j=1; j<n; j++) {
			//compute cost[j] here
			int find_min_cost = std::numeric_limits<int>::max();
			int backtrack_min_cost = std::numeric_limits<int>::max();

			// anchor i < anchor j 

			while (anchors[inner_loop_start].gap_query(anchors[j]) > bound_redit)
				inner_loop_start++;

			{
				//always consider the first dummy anchor 
				//connection to first dummy anchor is done with modified cost to allow free gaps
				//int i_d = std::get<1>(anchors[0]) + std::get<2>(anchors[0]) - 1;
				//int qry_gap = j_c - i_d - 1;
				int queryg = GAFHit::gap_query(anchors[0], anchors[j]);
				find_min_cost = std::min(find_min_cost, costs[0] + queryg);
				backtrack_min_cost = 0;
#ifdef CHAIN_HPP_DEBUG
					std::cerr << "anchors[" << "0" << "] -> anchors[" << j << "]: queryg = " << queryg << std::endl;
#endif
			}

			//process all anchors in array for the final last dummy anchor
			if (j == n-1)
				inner_loop_start=0;

			for(int i=j-1; i>=inner_loop_start; i--) {
				if (costs[i] < std::numeric_limits<int>::max() and are_colinear_eds(anchors[i], anchors[j], graph)) {
					int g;
					if (j == n-1) //modified cost for the last dummy anchor to allow free gaps
						g = GAFHit::gap_query(anchors[i], anchors[j]);
					else
						g = max_gap_eds(anchors[i], anchors[j], graph);

					int o = overlap_eds(anchors[i], anchors[j], graph);
#ifdef CHAIN_HPP_DEBUG
					std::cerr << "anchors[" << i << "] -> anchors[" << j << "]: g = " << g << ", o = " << o << std::endl;
#endif
					if (costs[i] + g + o < find_min_cost) {
						find_min_cost = costs[i] + g + o;
						backtrack_min_cost = i;
					}
				}
			}
			//save optimal cost at offset j
			costs[j] = find_min_cost;
			backtrack[j] = backtrack_min_cost;
		}

		if (costs[n-1] > bound_redit) {
			bound_redit = bound_redit * ramp_up_factor;
			revisions++;
		} else {
			break;
		}
	}

#ifdef CHAIN_HPP_DEBUG
		std::cerr << "Cost DP array = ";
		for (auto c : costs)
			std::cerr << c << " ";
		std::cerr << std::endl;
		std::cerr << "Backtrack array = ";
		for (auto b : backtrack)
			std::cerr << b << " ";
		std::cerr << std::endl;
		std::cerr << "Chaining cost computed " << revisions + 1 << " times" << "\n";
#endif

	//TODO: consider freeing here the space of costs array
	// backtrack optimal solution
	vector<GAFHit> solution;
	for (int j = backtrack[n - 1]; j > 0; j = backtrack[j])
	{
		solution.push_back(anchors[j]);
	}
	std::reverse(solution.begin(), solution.end());
	if (removesol)
	{
		vector<GAFHit> newanchors;
		newanchors.reserve(anchors.size() - solution.size());
		for (int i = 0, j = 0; i < anchors.size(); i += 1)
		{
			if (j < solution.size() and anchors[i] == solution[j])
			{
				j += 1;
			} else {
				newanchors.push_back(anchors[i]);
			}
		}
		swap(anchors, newanchors);
	}

	//std::cerr << "distance = " << costs[n-1] << std::endl;
	//std::cerr << "length - distance = " << anchors[0].get_query_length() - costs[n-1]  << std::endl;
	if (solution.size() > 0) {
		stats.maxiterations = std::max(stats.maxiterations, revisions);
		stats.miniterations = std::min(stats.miniterations, revisions);
		stats.totaliterations += revisions;
		stats.maxcost = std::max(stats.maxcost, costs[n-1]);
		stats.mincost = std::min(stats.mincost, costs[n-1]);
		stats.totalcost += costs[n-1];
		const double relativecost = (double)costs[n-1] / solution.at(0).get_query_length();
		stats.maxrelativecost = std::max(stats.maxrelativecost, relativecost);
		stats.minrelativecost = std::min(stats.minrelativecost, relativecost);
		stats.totalrelativecost += relativecost;
	}

	return solution;
}

} // Namespace chainx_block_graph

#endif //CHAIN_HPP
