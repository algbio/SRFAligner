#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <filesystem>
#include <vector>
#include <cassert>
#include <unordered_map>

#include "efg-locate.hpp"
#include "command-line-parsing/cmdline.h" // gengetopt-generated parser
#include "efg.hpp"
#include "algo.cpp"

//#define LOCATE_DEBUG

using namespace efg_locate;
using std::string;

int main(int argc, char* argv[])
{
	gengetopt_args_info argsinfo;
	if (cmdline_parser(argc, argv, &argsinfo) != 0) exit(1);

	if (argsinfo.inputs_num == 0)
		{cmdline_parser_print_help(); exit(1);};
	if (argsinfo.inputs_num == 1)
		{std::cerr << argv[0] << ": missing patterns file" << std::endl; exit(1);};
	if (argsinfo.inputs_num == 2)
		{std::cerr << argv[0] << ": missing output file" << std::endl; exit(1);};
	if (argsinfo.inputs_num > 3)
		{std::cerr << argv[0] << ": too many arguments" << std::endl; exit(1);};

	Params params;
	params.ignorechars = ((argsinfo.ignore_chars_arg != NULL) ? string(argsinfo.ignore_chars_arg): "");
	params.reversecompl = argsinfo.reverse_complement_flag;
	params.threads = argsinfo.threads_arg;
	params.mincoverage = argsinfo.approximate_min_coverage_arg;
	params.reportstats = argsinfo.approximate_stats_flag;
	params.renamereversecomplement = argsinfo.rename_reverse_complement_flag;
	params.splitoutputmatches = argsinfo.split_output_matches_flag;
	params.splitoutputmatchesgraphaligner = argsinfo.split_output_matches_graphaligner_flag;
	params.edgemincount = argsinfo.approximate_edge_match_min_count_arg;
	params.edgemincountheuristic = argsinfo.approximate_edge_match_min_count_heuristic_flag; // default false

	// open graph file
	std::filesystem::path graphpath {argsinfo.inputs[0]};
	params.graphfs = std::ifstream {graphpath};
	if (!params.graphfs) {std::cerr << "Error opening graph file " << graphpath << "." << std::endl; exit(1);};

	// check and open output file
	std::filesystem::path outputpath {argsinfo.inputs[2]};
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

	std::cerr << "Reading the graph..." << std::flush;
	Elasticfoundergraph graph(params.graphfs);
	std::cerr << " done." << std::endl;

	std::cerr << "Indexing the graph..." << std::flush;
	graph.init_pattern_matching_support();
	std::cerr << " done." << std::endl;

#ifdef LOCATE_DEBUG 
	std::cerr << "DEBUG graph is " << std::endl;
	graph.to_stream(&std::cerr);
#endif

	// check and open patterns file
	std::filesystem::path patternspath {argsinfo.inputs[1]};
	params.patternsfs = std::ifstream {patternspath};
	if (!params.patternsfs) {std::cerr << "Error opening patterns file " << patternspath << "." << std::endl; exit(1);};

	std::atomic<bool> input_done = false;
	std::thread inputworker;
	vector<string> pattern_ids, patterns;
	if (argsinfo.approximate_flag && params.threads > 1) {
		std::cerr << "Locate" << std::endl;
		inputworker = std::thread(reader_worker, std::ref(params.patternsfs), std::ref(input_done));
	} else {
		std::cerr << "Reading the patterns..." << std::flush;
		std::tie(pattern_ids, patterns) = read_patterns(params.patternsfs);
		std::cerr << " done." << std::endl;
	}

#ifdef LOCATE_DEBUG 
	std::cerr << std::endl;
	for (int i = 0; i < pattern_ids.size(); i++) {
		cerr << "DEBUG pattern:" << pattern_ids[i] << std::endl << patterns[i] << std::endl;
	}
#endif

	int returnvalue = 0;
	if (!argsinfo.approximate_flag and argsinfo.ignore_chars_arg == NULL) {
		if (params.reversecompl) {std::cerr << "Reverse complement not implemented in this mode!" << std::endl; exit(1);};
		if (params.threads != -1) {std::cerr << "Multithreading not implemented in this mode!" << std::endl; exit(1);};
		if (params.splitoutputmatches or params.splitoutputmatchesgraphaligner) {std::cerr << "Splitting of output matches not implemented in this mode!" << std::endl; exit(1);};

		for (int p = 0; p < patterns.size(); p++) {
			vector<int> path;

			if (efg_backward_search(graph, patterns[p], path) != 0) {
				path_to_stream(&params.outputfs, graph, pattern_ids[p], path);
			} else {
				cerr << "Pattern " << pattern_ids[p] << " does not occur!" << std::endl;
				returnvalue = 1;
			}
		}
		return returnvalue;
	}

	if (!argsinfo.approximate_flag and argsinfo.ignore_chars_arg != NULL) {
		if (params.reversecompl) {std::cerr << "Reverse complement not implemented in this mode!"<<  std::endl; exit(1);};
		if (params.threads != -1) {std::cerr << "Multithreading not implemented in this mode!" << std::endl; exit(1);};
		if (params.splitoutputmatches or params.splitoutputmatchesgraphaligner) {std::cerr << "Splitting of output matches not implemented in this mode!" << std::endl; exit(1);};

		for (int p = 0; p < patterns.size(); p++) {
			vector<vector<int>> paths;

			if (efg_backward_search_ignorechars(graph, params.ignorechars, patterns[p], paths) != 0) {
				path_to_stream(&params.outputfs, graph, pattern_ids[p], paths);
			} else {
				cerr << "Pattern " << pattern_ids[p] << " does not occur!" << std::endl;
				returnvalue = 1;
			}
		}
		return returnvalue;
	}

	if (argsinfo.approximate_flag) {
		if (params.threads > 0) {
			std::atomic<bool> workers_done = false;
			std::thread outputworker(writer_worker, std::ref(workers_done), std::ref(params));
			vector<std::thread> workers;
			for (int i = 0; i < params.threads; i++)
				workers.push_back(std::thread(approx_worker, std::ref(graph), std::ref(pattern_ids), std::ref(patterns), std::ref(params), std::ref(input_done)));
			for (int i = 0; i < workers.size(); i++)
				workers[i].join();
			workers_done = true;
			inputworker.join();
			outputworker.join();
			// sanity check?
			outputworker = std::thread(writer_worker, std::ref(workers_done), std::ref(params));
			outputworker.join();
		} else {
			for (int p = 0; p < patterns.size(); p++) {
				vector<GAFAnchor> matches;

				if (approx_efg_backward_search(graph, pattern_ids[p], patterns[p], params, matches) != 0) {
					if (params.splitoutputmatches)
						anchors_to_stream_split_single(&params.outputfs, graph, matches);
					else if (params.splitoutputmatchesgraphaligner)
						anchors_to_stream_split_single_graphaligner(&params.outputfs, graph, matches);
					else
						anchors_to_stream(&params.outputfs, graph, matches);
				} else {
					cerr << "Cannot find any semi-repeat-free match of " << pattern_ids[p] << std::endl;
				}
			}
		}

		return 0;
	}

	cerr << "Mode not implemented!" << std::endl;
	return 1;
}
