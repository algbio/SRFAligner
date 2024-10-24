#ifndef CHAINX_BLOCK_GRAPH_HPP
#define CHAINX_BLOCK_GRAPH_HPP

#include <iostream>
#include <fstream>
#include <limits>

using std::string, std::ifstream, std::ofstream;

namespace chainx_block_graph {
struct Params {
	ifstream graphfs;
	ifstream anchorsfs;
	ofstream outputfs;
	string ignorechars;
	bool unsorted_anchors;
	bool global;
	bool semiglobal;
	bool nosplit;
	bool splitgraphaligner;
	int threads;
	int alternativealignments;
	int initialguess;
	double initialguesscov;
	double rampupfactor;
};

struct Stats {
	unsigned long long seeds = 0;
	unsigned long long reads = 0;
	int maxiterations = 0;
	int miniterations = std::numeric_limits<int>::max();
	unsigned long long totaliterations = 0;

	int maxcost = 0;
	int mincost = std::numeric_limits<int>::max();
	unsigned long long totalcost = 0;

	double maxrelativecost = 0;
	double minrelativecost = std::numeric_limits<double>::max();
	double totalrelativecost = 0;
};

struct Stats mergestats(const struct Stats &s1, const struct Stats &s2)
{
	struct Stats s;
	s.seeds = s1.seeds + s2.seeds;
	s.reads = s1.reads + s2.reads;
	s.maxiterations = std::max(s1.maxiterations, s2.maxiterations);
	s.miniterations = std::min(s1.miniterations, s2.miniterations);
	s.totaliterations = s1.totaliterations + s2.totaliterations;
	s.maxcost = std::max(s1.maxcost, s2.maxcost);
	s.mincost = std::min(s1.mincost, s2.mincost);
	s.totalcost = s1.totalcost + s2.totalcost;
	s.maxrelativecost = std::max(s1.maxrelativecost, s2.maxrelativecost);
	s.minrelativecost = std::min(s1.minrelativecost, s2.minrelativecost);
	s.totalrelativecost = s1.totalrelativecost + s2.totalrelativecost;
	return s;
}

}

#endif
