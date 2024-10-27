#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <filesystem>
#include <vector>
#include <cassert>
#include <unordered_map>

#include "command-line-parsing/cmdline.h" // gengetopt-generated parser
#include "efg.hpp"

using std::string;

int main(int argc, char* argv[])
{
	gengetopt_args_info argsinfo;
	if (cmdline_parser(argc, argv, &argsinfo) != 0) exit(1);

	if (argsinfo.inputs_num == 0)
		{cmdline_parser_print_help(); exit(1);};
	if (argsinfo.inputs_num == 1)
		{std::cerr << argv[0] << ": missing GAF file" << std::endl; exit(1);};
	if (argsinfo.inputs_num > 2)
		{std::cerr << argv[0] << ": too many arguments" << std::endl; exit(1);};

	// open files
	std::filesystem::path graphpath {argsinfo.inputs[0]};
	std::ifstream graphfs {graphpath};
	if (!graphfs) {std::cerr << "Error opening graph file " << graphpath << "." << std::endl; exit(1);};

	std::filesystem::path gafpath {argsinfo.inputs[1]};
	std::ifstream gaffs {gafpath};
	if (!gaffs) {std::cerr << "Error opening GAF file " << gafpath << "." << std::endl; exit(1);};

	std::cerr << "Reading the graph..." << std::flush;
	Elasticfoundergraph graph(graphfs);
	std::cerr << " done." << std::endl;

	std::cerr << "Reading the seeds..." << std::flush;
	vector<vector<GAFAnchor>> seeds = read_gaf(gaffs, graph);
	std::cerr << " done." << std::endl;

	std::cerr << "Splitting the seeds..." << std::flush;
	for (auto &patternseeds : seeds) {
		for (auto &a : patternseeds) {
			for (auto &b : a.split_single_graphaligner(graph)) {
				if (b.get_query_id().find("rev_") != std::string::npos) {
					b.reverse();
				}
				std::cout << b.to_string(graph) << std::endl;
			}
		}
	}
	std::cerr << " done." << std::endl;
}
