#ifndef EFG_LOCATE_HPP
#define EFG_LOCATE_HPP

#include <iostream>
#include <fstream>

using std::string, std::ifstream, std::ofstream;

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
	int edgemincount;
	bool edgemincountheuristic;
};
#endif
