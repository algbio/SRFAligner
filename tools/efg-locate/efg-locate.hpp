#ifndef EFG_LOCATE_HPP
#define EFG_LOCATE_HPP

#include <iostream>
#include <fstream>

using std::string, std::ifstream, std::ofstream;

namespace efg_locate {
struct Params {
	ifstream graphfs;
	ifstream patternsfs;
	ofstream outputfs;
	string ignorechars;
	bool reversecompl;
	int threads;
	int mincoverage;
	bool reportstats;
	bool renamereversecomplement;
	bool splitoutputmatches;
	bool splitoutputmatchesgraphaligner;
	bool splitkeepedgematches;
	int edgemincount;
	int edgelongestcount;
	int edgelongestcountmax;
};
}

#endif
