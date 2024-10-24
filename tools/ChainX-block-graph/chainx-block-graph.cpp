#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <syncstream>
#include <filesystem>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <thread>
#include "concurrentqueue.h" // https://github.com/cameron314/concurrentqueue

#include "chainx-block-graph.hpp"
#include "command-line-parsing/cmdline.h" // gengetopt-generated parser
#include "efg.hpp"
#include "chaining.hpp"

//#define CHAINX_DEBUG

using namespace chainx_block_graph;
using std::string;
using std::move, std::back_inserter;
std::mutex mutp;
moodycamel::ConcurrentQueue<std::string*> outputqueue { 50, 64, 64 }; // TODO: parameterize this!
moodycamel::ConcurrentQueue<string> anchorsqueue { 50, 64, 64 }; // TODO: parameterize this!
moodycamel::ConcurrentQueue<std::pair<string,vector<string>>*> taskqueue { 50, 64, 64 }; // TODO: parameterize this!

void writer_worker(std::atomic<bool> &input_done, std::atomic<bool>& done, Params &params)
{
	string *ptr;
	while (true) {
		while (outputqueue.try_dequeue(ptr)) {
			params.outputfs << *ptr << "\n";
			delete(ptr);
		}
		// queue is empty?
		if (!done or !input_done)
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
		else
			break;
	}
}

void reader_worker(const Elasticfoundergraph &efg, std::atomic<bool> &input_done, Params &params)
{
	std::pair<string,vector<string>> *ptr;
	std::unordered_set<string> anchor_ids;
	string current_id = "", line;
	vector<string> gafhits; // string description of anchors
	string qname;
	while (std::getline(params.anchorsfs, line)) {
		if (line == "")
			continue;
	
		//istringstream linestream(line);
		//const GAFHit a(linestream, efg, qname);
		qname = read_gaf_query_id(line);

		if (qname == current_id) {
			//gafhits.push_back(a);
			gafhits.push_back(line);
		} else {
			if (current_id != "" and anchor_ids.contains(qname)) {
				std::cerr << "Fatal error for anchors of query " << qname << ": are the anchors sorted? Sort them, or try flag --unsorted_input" << std::endl;
				exit(1);
			}
			if (current_id != "") {
				//std::pair<string,vector<GAFHit>> *ptr = new std::pair(current_id,vector<GAFHit>(std::move(gafhits)));
				std::pair<string,vector<string>> *ptr = new std::pair(current_id,vector<string>(std::move(gafhits)));
				while (!taskqueue.try_enqueue(ptr))
					std::this_thread::sleep_for(std::chrono::milliseconds(10));
				gafhits.clear(); // TODO assert?
				////gafhits.push_back(dummy_start(a, efg));
				//istringstream liness(line);
				//GAFHit a(liness, efg, qname);
				//gafhits.push_back(dummy_start(a, efg).to_string(efg, qname));
			}
			if (current_id == "") {
				//istringstream liness(line);
				//GAFHit a(liness, efg, qname);
				//gafhits.push_back(dummy_start(a, efg).to_string(efg, qname));
			}

			anchor_ids.insert(qname);
			current_id = qname;
			//gafhits.push_back(a);
			gafhits.push_back(line);
		}
	}
	if (gafhits.size() > 0) {
		//taskqueue.enqueue(new std::pair(current_id,vector<GAFHit>(std::move(gafhits))));
		taskqueue.enqueue(new std::pair(current_id,vector<string>(std::move(gafhits))));
	}
	input_done = true;
}

void chain_worker(const Elasticfoundergraph &graph, unordered_map<string,vector<GAFHit>> &anchors, std::atomic<bool> &input_done, Params &params, vector<Stats> &stats, int statsindex)
{
	vector<GAFHit> matches;
	string query;
	bool warning = false;
	if (params.unsorted_anchors) {
		while (anchorsqueue.try_dequeue(query)) {
			stats[statsindex].reads += 1;
			/*{
			  std::scoped_lock lck {mutp};
			  if (p == anchors.end())
			  return;

			  query = p.first;
			  anchorlist = p.second;
			  ++p;
			  }*/
			vector<GAFHit> &anchorlist = anchors.at(query);
			stats[statsindex].seeds += anchorlist.size() - 1;

			//std::cerr << "Chaining for query " << anchorlist[0].get_query_id() << "..." << std::endl;
			//anchorlist.push_back(dummy_start(anchorlist.at(0), graph)); // dummy start is already in the list
			anchorlist.push_back(dummy_end(anchorlist.at(0), graph));
			if (!is_sorted(anchorlist)) {
				if (!warning) {
					std::cerr << "Anchors do not seem to be sorted by starting position in the query! Sorting...\n";
					warning = true;
				}
				std::sort(anchorlist.begin(), anchorlist.end(), [](GAFHit a1, GAFHit a2)
						{
						return (a1.get_query_start() < a2.get_query_start());
						});
			}

#ifdef CHAINX_DEBUG
			{
				std::osyncstream osscerr(cerr);
				std::cerr << "DEBUG Sorted anchors for query " << query << " are" << std::endl;
				for (auto &a : anchorlist)
					cerr << "DEBUG " << a.to_string(graph, query) << std::endl;
				std::cerr << std::endl;
			}
#endif

			vector<GAFHit> solution;
			int initial_guess;
			if (params.initialguesscov == 0) {
				initial_guess = params.initialguess;
			} else {
				initial_guess = anchorlist.at(0).get_query_length() - (compute_coverage_greedy(anchorlist) * params.initialguesscov);
			}
			if (params.global) {
				if (params.alternativealignments == 0) {
					solution = chain_global_eds(anchorlist, graph, initial_guess, params.rampupfactor, stats[statsindex]);
				} else {
					vector<GAFHit> chain;
					for (int i = params.alternativealignments; i > 0; i--) {
						chain = chain_global_eds(anchorlist, graph, initial_guess, params.rampupfactor,  stats[statsindex], true);
						solution.reserve(solution.size() + chain.size());
						move(chain.begin(), chain.end(), back_inserter(solution));
						chain.clear();
					}
					chain = chain_global_eds(anchorlist, graph, initial_guess, params.rampupfactor,  stats[statsindex]);
					solution.reserve(solution.size() + chain.size());
					move(chain.begin(), chain.end(), back_inserter(solution));
					//chain.clear();
				}
			} else if (params.semiglobal) {
				if (params.alternativealignments == 0) {
					solution = chain_semiglobal_eds(anchorlist, graph, initial_guess, params.rampupfactor,  stats[statsindex]);
				} else {
					vector<GAFHit> chain;
					for (int i = params.alternativealignments; i > 0; i--) {
						chain = chain_semiglobal_eds(anchorlist, graph, initial_guess, params.rampupfactor,  stats[statsindex], true);
						solution.reserve(solution.size() + chain.size());
						move(chain.begin(), chain.end(), back_inserter(solution));
						chain.clear();
					}
					chain = chain_semiglobal_eds(anchorlist, graph, initial_guess, params.rampupfactor,  stats[statsindex]);
					solution.reserve(solution.size() + chain.size());
					move(chain.begin(), chain.end(), back_inserter(solution));
					//chain.clear();
				}
			}
			//std::cerr << "...done." << std::endl;

			bool reverse = false;
			if (query.find("rev_") != std::string::npos) {
				for (auto &a : solution)
					a.reverse();
				reverse = true;
			}
			if (params.nosplit) {
				std::osyncstream ossoutput(params.outputfs);
				for (auto &a : solution) {
					string *ptr = new std::string(a.to_string(graph, (reverse) ? query.substr(4) : query));
					while (!outputqueue.try_enqueue(ptr))
						std::this_thread::sleep_for(std::chrono::milliseconds(10));
				}
			} else {
				std::osyncstream ossoutput(params.outputfs);
				for (auto &a : solution) {
					if (params.splitgraphaligner) {
						for (auto &b : a.split_single_graphaligner(graph)) {
							string *ptr = new std::string(b.to_string(graph, (reverse) ? query.substr(4) : query));
							while(!outputqueue.try_enqueue(ptr))
								std::this_thread::sleep_for(std::chrono::milliseconds(10));
						}
					} else {
						for (auto &b : a.split(graph)) {
							string *ptr = new std::string(b.to_string(graph, (reverse) ? query.substr(4) : query));
							while (!outputqueue.try_enqueue(ptr))
								std::this_thread::sleep_for(std::chrono::milliseconds(10));
						}
					}
				}
			}
		}
	} else {
		//std::pair<string,vector<GAFHit>>* task;
		std::pair<string,vector<string>>* task;
		while (true) {
			while (taskqueue.try_dequeue(task)) {
				stats[statsindex].reads += 1;
				vector<GAFHit> anchors;
				//task->second.push_back(dummy_start(task->second.at(0), graph)); // BUG?
				//task->second.push_back(dummy_end(task->second.at(0), graph));
				string query = string(task->first);
				istringstream dummyss(task->second.at(0));
				anchors.push_back(dummy_start(GAFHit(dummyss, graph, query), graph));

				for (string &line : task->second) {
					istringstream ss(line);
					anchors.push_back(GAFHit(ss, graph, query));
				}

				dummyss = istringstream(task->second.at(0));
				anchors.push_back(dummy_end(GAFHit(dummyss, graph, query), graph));
				stats[statsindex].seeds += anchors.size() - 2;
				if (!is_sorted(anchors)) {
					std::cerr << "Anchors do not seem to be sorted for read " << query << "...\nSorting...\n";
					std::sort(anchors.begin(), anchors.end(), [](GAFHit a1, GAFHit a2)
							{
							return (a1.get_query_start() < a2.get_query_start());
							});
				}

				/*//std::sort(task->second.begin(), task->second.end(), [](GAFHit a1, GAFHit a2)
				std::sort(anchors.begin(), anchors.end(), [](GAFHit a1, GAFHit a2)
						{
						return (a1.get_query_start() < a2.get_query_start());
						});*/

				vector<GAFHit> solution;
				int initial_guess;
				if (params.initialguesscov == 0) {
					initial_guess = params.initialguess;
				} else {
					initial_guess = anchors.at(0).get_query_length() - (compute_coverage_greedy(anchors) * params.initialguesscov);
				}
				if (params.global) {
					if (params.alternativealignments == 0) {
						solution = chain_global_eds(anchors, graph, initial_guess, params.rampupfactor,  stats[statsindex]);
					} else {
						vector<GAFHit> chain;
						for (int i = params.alternativealignments; i > 0; i--) {
							chain = chain_global_eds(anchors, graph, initial_guess, params.rampupfactor,  stats[statsindex], true);
							solution.reserve(solution.size() + chain.size());
							move(chain.begin(), chain.end(), back_inserter(solution));
							chain.clear();
						}
						chain = chain_global_eds(anchors, graph, initial_guess, params.rampupfactor,  stats[statsindex]);
						solution.reserve(solution.size() + chain.size());
						move(chain.begin(), chain.end(), back_inserter(solution));
						//chain.clear();
					}
				} else if (params.semiglobal) {
					if (params.alternativealignments == 0) {
						solution = chain_semiglobal_eds(anchors, graph, initial_guess, params.rampupfactor,  stats[statsindex]);
					} else {
						vector<GAFHit> chain;
						for (int i = params.alternativealignments; i > 0; i--) {
							chain = chain_semiglobal_eds(anchors, graph, initial_guess, params.rampupfactor,  stats[statsindex], true);
							solution.reserve(solution.size() + chain.size());
							move(chain.begin(), chain.end(), back_inserter(solution));
							chain.clear();
						}
						chain = chain_semiglobal_eds(anchors, graph, initial_guess, params.rampupfactor,  stats[statsindex]);
						solution.reserve(solution.size() + chain.size());
						move(chain.begin(), chain.end(), back_inserter(solution));
						//chain.clear();
					}
				}
				//std::cerr << "...done." << std::endl;

				bool reverse = false;
				if (query.find("rev_") != std::string::npos) {
					for (auto &a : solution)
						a.reverse();
					reverse = true;
				}
				if (params.nosplit) {
					std::osyncstream ossoutput(params.outputfs);
					for (auto &a : solution) {
						string *ptr = new std::string(a.to_string(graph, (reverse) ? query.substr(4) : query));
						while (!outputqueue.try_enqueue(ptr))
							std::this_thread::sleep_for(std::chrono::milliseconds(10));
					}
				} else {
					std::osyncstream ossoutput(params.outputfs);
					for (auto &a : solution) {
						if (params.splitgraphaligner) {
							for (auto &b : a.split_single_graphaligner(graph)) {
								string *ptr = new std::string(b.to_string(graph, (reverse) ? query.substr(4) : query));
								while(!outputqueue.try_enqueue(ptr))
									std::this_thread::sleep_for(std::chrono::milliseconds(10));
							}
						} else {
							for (auto &b : a.split(graph)) {
								string *ptr = new std::string(b.to_string(graph, (reverse) ? query.substr(4) : query));
								while (!outputqueue.try_enqueue(ptr))
									std::this_thread::sleep_for(std::chrono::milliseconds(10));
							}
						}
					}
				}
				delete(task);
			}
			if (input_done)
				break;
			else
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
		}
	}
}

int main(int argc, char* argv[])
{
	gengetopt_args_info argsinfo;
	if (cmdline_parser(argc, argv, &argsinfo) != 0) exit(1);

	if (argsinfo.inputs_num == 0)
		{cmdline_parser_print_help(); exit(1);};
	if (argsinfo.inputs_num == 1)
		{std::cerr << argv[0] << ": missing anchors in input and output file" << std::endl; exit(1);};
	if (argsinfo.inputs_num == 2)
		{std::cerr << argv[0] << ": missing output file" << std::endl; exit(1);};
	if (argsinfo.inputs_num > 3)
		{std::cerr << argv[0] << ": too many arguments" << std::endl; exit(1);};

	Params params;
	params.threads = argsinfo.threads_arg;
	params.global = argsinfo.global_flag;
	params.semiglobal = argsinfo.semi_global_flag;
	params.nosplit = argsinfo.no_split_output_matches_flag;
	params.unsorted_anchors = argsinfo.unsorted_input_flag;
	params.splitgraphaligner = argsinfo.split_output_matches_graphaligner_flag;
	params.alternativealignments = argsinfo.alternative_chains_arg;
	params.initialguess = argsinfo.initial_guess_arg;
	params.initialguesscov = argsinfo.initial_guess_coverage_arg;
	params.rampupfactor = argsinfo.ramp_up_factor_arg;

	Stats stats;

	if (!params.global and !params.semiglobal)
		{std::cerr << argv[0] << ": select mode (global or semiglobal)" << std::endl; exit(1);};
	if (params.global and params.semiglobal)
		{std::cerr << argv[0] << ": select only one mode (global or semiglobal)" << std::endl; exit(1);};
	// TODO check ramp-up factor

	// open files
	std::filesystem::path graphpath {argsinfo.inputs[0]};
	params.graphfs = std::ifstream {graphpath};
	if (!params.graphfs) {std::cerr << "Error opening graph file " << graphpath << "." << std::endl; exit(1);};

	std::filesystem::path anchorspath {argsinfo.inputs[1]};
	params.anchorsfs = std::ifstream {anchorspath};
	if (!params.anchorsfs) {std::cerr << "Error opening anchors file " << anchorspath << "." << std::endl; exit(1);};

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
	graph.init_eds_support();
	std::cerr << " done." << std::endl;

	std::unordered_map<string,vector<GAFHit>> anchors; 
	if (params.unsorted_anchors) {
		std::cerr << "Reading the anchors..." << std::flush;
		anchors = read_gaf_chaining(params.anchorsfs, graph, anchorsqueue);
		{
			long unsigned int n = 0;
			for (auto &v : anchors) {
				n += v.second.size() - 1;
			}
			std::cerr << " done (read " << n << " anchors)." << std::endl;
		}

#ifdef CHAINX_DEBUG
		std::cerr << std::endl;
		for (auto &bucket : anchors) {
			cerr << "DEBUG Set of anchors:" << std::endl;
			for (auto &a : bucket.second)
				cerr << "DEBUG " << a.to_string(graph, bucket.first) << std::endl;
		}
		std::cerr << std::endl;
#endif
	}


	if (params.threads == -1) {
		if (!params.unsorted_anchors) {
			std::cerr << "Mode not implemented yet! Try --unsorted_input flag." << std::endl;
			exit(1);
		}
		for (auto &bucket : anchors) {
			const string &query = bucket.first;
			vector<GAFHit> &anchorlist = bucket.second;
			//std::cerr << "Chaining for query " << query << "..." << std::endl;
			anchorlist.push_back(dummy_end(anchorlist[0], graph));
			std::sort(anchorlist.begin(), anchorlist.end(), [](GAFHit a1, GAFHit a2)
					{
					return (a1.get_query_start() < a2.get_query_start());
					});

#ifdef CHAINX_DEBUG
			std::cerr << "DEBUG Sorted anchors for query " << query << " are" << std::endl;
			for (auto &a : anchorlist)
				cerr << "DEBUG " << a.to_string(graph, query) << std::endl;
			std::cerr << std::endl;
#endif

			vector<GAFHit> solution;
			int initial_guess;
			if (params.initialguesscov == 0) {
				initial_guess = params.initialguess;
			} else {
				initial_guess = anchorlist.at(0).get_query_length() - (compute_coverage_greedy(anchorlist) * params.initialguesscov);
			}
			if (params.global) {
				if (params.alternativealignments == 0) {
					solution = chain_global_eds(anchorlist, graph, initial_guess, params.rampupfactor, stats);
				} else {
					vector<GAFHit> chain;
					for (int i = params.alternativealignments; i > 0; i--) {
						chain = chain_global_eds(anchorlist, graph, initial_guess, params.rampupfactor, stats, true);
						solution.reserve(solution.size() + chain.size());
						move(chain.begin(), chain.end(), back_inserter(solution));
						chain.clear();
					}
					chain = chain_global_eds(anchorlist, graph, initial_guess, params.rampupfactor, stats);
					solution.reserve(solution.size() + chain.size());
					move(chain.begin(), chain.end(), back_inserter(solution));
					//chain.clear();
				}
			} else if (params.semiglobal) {
				if (params.alternativealignments == 0) {
					solution = chain_semiglobal_eds(anchorlist, graph, initial_guess, params.rampupfactor, stats);
				} else {
					vector<GAFHit> chain;
					for (int i = params.alternativealignments; i > 0; i--) {
						chain = chain_semiglobal_eds(anchorlist, graph, initial_guess, params.rampupfactor, stats, true);
						solution.reserve(solution.size() + chain.size());
						move(chain.begin(), chain.end(), back_inserter(solution));
						chain.clear();
					}
					chain = chain_semiglobal_eds(anchorlist, graph, initial_guess, params.rampupfactor, stats);
					solution.reserve(solution.size() + chain.size());
					move(chain.begin(), chain.end(), back_inserter(solution));
					//chain.clear();
				}
			}
			//std::cerr << "...done." << std::endl;

			bool reverse = false;
			if (query.find("rev_") != std::string::npos) {
				for (auto &a : solution)
					a.reverse();
				reverse = true;
			}
			if (params.nosplit) {
				for (auto &a : solution)
					params.outputfs << a.to_string(graph, (reverse) ? query.substr(4) : query) << std::endl;
			} else {
				for (auto &a : solution) {
					if (params.splitgraphaligner) {
						for (auto &b : a.split_single_graphaligner(graph))
							params.outputfs << b.to_string(graph, (reverse) ? query.substr(4) : query) << std::endl;
					} else {
						for (auto &b : a.split(graph))
							params.outputfs << b.to_string(graph, (reverse) ? query.substr(4) : query) << std::endl;
					}
				}
			}
		}
	} else {
		std::atomic<bool> workers_done = false;
		std::atomic<bool> input_done = params.unsorted_anchors ? true : false;
		std::thread inputworker;
		if (!params.unsorted_anchors) {
			inputworker = std::thread(reader_worker, std::ref(graph), std::ref(input_done), std::ref(params));
		}
		std::thread outputworker(writer_worker, std::ref(input_done), std::ref(workers_done), std::ref(params));
		vector<std::thread> workers;
		vector<Stats> workerstats(params.threads);
		for (int i = 0; i < params.threads; i++) {
			workers.push_back(std::thread(chain_worker, std::ref(graph), std::ref(anchors), std::ref(input_done), std::ref(params), std::ref(workerstats), i));
		}
		for (int i = 0; i < workers.size(); i++)
			workers[i].join();
		workers_done = true;

		for (const Stats &s : workerstats) {
			stats = mergestats(stats, s);
		}

		if (!params.unsorted_anchors)
			inputworker.join();
		outputworker.join();
		// sanity check?
		outputworker = std::thread(writer_worker, std::ref(input_done), std::ref(workers_done), std::ref(params));
		outputworker.join();
	}

	std::cerr << "chained " << stats.seeds << " seeds for " << stats.reads << " reads\n";
	std::cerr << "min number of revisions is " << stats.miniterations << "\n";
	std::cerr << "max number of revisions is " << stats.maxiterations << "\n";
	std::cerr << "average number of revisions is " << (double)stats.totaliterations / stats.reads << "\n";
	std::cerr << "min chaining cost is " << stats.mincost << "\n";
	std::cerr << "max chaining cost is " << stats.maxcost << "\n";
	std::cerr << "average chaining cost is " << (double)stats.totalcost / stats.reads << "\n";
	std::cerr << "min relative chaining cost is " << stats.minrelativecost << "\n";
	std::cerr << "max relative chaining cost is " << stats.maxrelativecost << "\n";
	std::cerr << "average relative chaining cost is " << (double)stats.totalrelativecost / stats.reads << "\n";
}
