#include <gcsa/algorithms.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <filesystem>
#include <atomic>

#include "gcsa2-locate.hpp"
#include "graph.hpp"
#include "command-line-parsing/cmdline.h" // gengetopt-generated parser
#include "algo.hpp"

using namespace gcsa2_locate;
using std::string, std::max;

int main(int argc, char* argv[])
{
	gengetopt_args_info argsinfo;
	if (cmdline_parser(argc, argv, &argsinfo) != 0) exit(1);

	if (argsinfo.inputs_num == 0)
		{cmdline_parser_print_help(); exit(1);};
	if (argsinfo.inputs_num == 1)
		{std::cerr << argv[0] << ": missing GCSA2 index file" << std::endl; exit(1);};
	if (argsinfo.inputs_num == 2)
		{std::cerr << argv[0] << ": missing patterns file" << std::endl; exit(1);};
	if (argsinfo.inputs_num == 3)
		{std::cerr << argv[0] << ": missing output file" << std::endl; exit(1);};
	if (argsinfo.inputs_num > 4)
		{std::cerr << argv[0] << ": too many arguments" << std::endl; exit(1);};

	Params params;
	params.gcsapath = path(argsinfo.inputs[1]);
	params.ignorechars = ((argsinfo.ignore_chars_arg != NULL) ? string(argsinfo.ignore_chars_arg): "");
	params.reversecompl = argsinfo.reverse_complement_flag;
	params.threads = argsinfo.threads_arg;
	params.splitoutputmatchesgraphaligner = argsinfo.split_output_matches_graphaligner_flag;
	params.edgelongestcount = argsinfo.approximate_edge_match_longest_arg;
	params.edgelongestcountmax = (long unsigned int)argsinfo.approximate_edge_match_longest_max_count_arg;

	// open and load graph
	std::filesystem::path graphpath {argsinfo.inputs[0]};
	params.graphfs = std::ifstream {graphpath};
	if (!params.graphfs) {std::cerr << "Error opening graph file " << graphpath << "." << std::endl; exit(1);};
	graph::SimpleGraph graph(params.graphfs);

	// load GCSA2 index file
	gcsa::GCSA index;
	if(!sdsl::load_from_file(index, params.gcsapath.c_str())) { std::cerr << "ERROR: cannot load the index from " << params.gcsapath << std::endl; exit(1); }

	// check and open output file
	std::filesystem::path outputpath {argsinfo.inputs[3]};
	if (std::filesystem::exists(outputpath)) {
		if (argsinfo.overwrite_flag) {
			params.outputfs = std::ofstream(outputpath, std::ios::out | std::ios::trunc);
		} else {
			std::cerr << "Error: output file already exists." << std::endl;
			exit(1);
		}
	} else {
		params.outputfs = std::ofstream(outputpath);
	}
	if (!params.outputfs) {std::cerr << "Error opening output file " << outputpath << "." << std::endl; exit(1);};

	// check and open patterns file
	std::filesystem::path patternspath {argsinfo.inputs[2]};
	params.patternsfs = std::ifstream {patternspath};
	if (!params.patternsfs) {std::cerr << "Error opening patterns file " << patternspath << "." << std::endl; exit(1);};

	std::atomic<bool> input_done = false;
	std::thread inputworker;
	vector<string> pattern_ids, patterns;
	inputworker = std::thread(reader_worker, std::ref(params.patternsfs), std::ref(input_done));

	// exact pattern matching
	if (!argsinfo.approximate_flag) {
		std::atomic<bool> workers_done = false;
		std::thread outputworker(writer_worker, std::ref(workers_done), std::ref(params));
		vector<std::thread> workers;
		for (int i = 0; i < max(1,params.threads); i++)
			workers.push_back(std::thread(exact_worker, std::ref(graph), std::ref(index), std::ref(params), std::ref(input_done)));
		for (long unsigned int i = 0; i < workers.size(); i++)
			workers[i].join();
		workers_done = true;
		inputworker.join();
		outputworker.join();
		// sanity check?
		outputworker = std::thread(writer_worker, std::ref(workers_done), std::ref(params));
		outputworker.join();
		return 0;
	} else {
		std::atomic<bool> workers_done = false;
		std::thread outputworker(writer_worker, std::ref(workers_done), std::ref(params));
		vector<std::thread> workers;
		for (int i = 0; i < std::max(1,params.threads); i++)
			workers.push_back(std::thread(approx_worker, std::ref(graph), std::ref(index), std::ref(params), std::ref(input_done)));
		for (long unsigned int i = 0; i < workers.size(); i++)
			workers[i].join();
		workers_done = true;
		inputworker.join();
		outputworker.join();
		// sanity check?
		outputworker = std::thread(writer_worker, std::ref(workers_done), std::ref(params));
		outputworker.join();
	}
}
