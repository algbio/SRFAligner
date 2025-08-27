#ifndef EFG_LOCATE_HPP
#define EFG_LOCATE_HPP

#include <iostream>
#include <fstream>
#include <filesystem>

using std::string, std::ifstream, std::ofstream, std::filesystem::path;

namespace gcsa2_locate {
struct Params {
	ifstream graphfs;
	path gcsapath;
	ifstream patternsfs;
	ofstream outputfs;
	string ignorechars;
	bool reversecompl;
	int threads;
	bool splitoutputmatchesgraphaligner;
	int edgelongestcount;
	long unsigned int edgelongestcountmax;
};
}

#endif
