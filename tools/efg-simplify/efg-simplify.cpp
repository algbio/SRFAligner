/**
 * Program to check whether a given graph is repeat-free/semi-repeat-free.
 * See output of ./efg-simplify --help or command-line-parsing/config.ggo
**/

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <filesystem>
#include <algorithm> // std::find
#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <iterator>

#include "command-line-parsing/cmdline.h" // gengetopt-generated parser

//#define DEBUG

using std::vector, std::map, std::set, std::pair, std::get, std::string;

class Elasticfoundergraph {
	private:
		int m = 0, n = 0; // rows, cols
		vector<int> cuts, heights, cumulative_height;
		vector<string> ordered_node_ids;
		vector<string> ordered_node_labels;
		map<int,vector<int>> edges; // adjacency lists
		vector<string> walk_ids;
		vector<vector<int>> walks;
		vector<vector<bool>> orientations; // true for +, false for -


		bool is_to_simplify(int blocki, bool ignore_only, string &ignorechars) // check if the next block can be merged
		{
			// check if edges between block i and i+1 are a bijection

			if (heights[blocki] != heights[blocki+1])
				return false;

			if (ignore_only) {
				for (int n = cumulative_height[blocki]; n < cumulative_height[blocki] + heights[blocki]; n++) {
					for (char c : ordered_node_labels[n]) {
						if (ignorechars.find(c) == std::string::npos) {
							return false;
						}
					}
				}
			}

			vector<bool> hitsnextblock(heights[blocki+1], false);
			for (int n = cumulative_height[blocki]; n < cumulative_height[blocki] + heights[blocki]; n++) {
				if (edges[n].size() != 1) {
					return false;
				}
				hitsnextblock[edges[n][0] - cumulative_height[blocki] - heights[blocki]] = true;
			}

			for (int n = 0; n < hitsnextblock.size(); n++) {
				if (!hitsnextblock[n]) {
					return false;
				}
			}

			return true;
		}
	public:
		Elasticfoundergraph(std::ifstream *graphstream)
		{
			int nodes = 0;
			map<string,int> node_index;

			for (std::string line; std::getline(*graphstream, line); ) {
				if (line[0] == 'M') {
					assert(m == 0 && n == 0);

					char c;
					std::istringstream(line) >> c >> m >> n;
				} else if (line[0] == 'X') {
					assert(cuts.size() == 0);

					std::istringstream buffer(line.substr(2)); // ignore starting "X\t"
					vector<int> ccuts({
						std::istream_iterator<int>(buffer), std::istream_iterator<int>()
					});
					std::istringstream().swap(buffer); // clear buffer
					std::swap(ccuts, cuts);
				} else if (line[0] == 'B') {
					assert(heights.size() == 0);

					std::istringstream buffer(line.substr(2)); // ignore starting "X\t"
					vector<int> hheights({
						std::istream_iterator<int>(buffer), std::istream_iterator<int>()
					});
					std::istringstream().swap(buffer); // clear buffer
					std::swap(hheights, heights);

					vector<int> cumulative_h(heights.size() + 1, 0);
					for (int i = 1; i <= heights.size(); i++) {
						cumulative_h[i] = cumulative_h[i-1] + heights[i-1];
					}
					std::swap(cumulative_h, cumulative_height);
				} else if (line[0] == 'S') {
					char c;
					string id, label;

					std::istringstream(line) >> c >> id >> label;
					node_index[id] = nodes;
					ordered_node_ids.push_back(id);
					ordered_node_labels.push_back(label);

					nodes += 1;
				} else if (line[0] == 'L') {
					char c;
					std::string id1, id2;
					std::istringstream(line) >> c >> id1 >> c >> id2;

					edges[node_index[id1]].push_back(node_index[id2]);
				} else if (line[0] == 'P') {
					char c;
					string walkid, node;
					vector<int> walk;
					vector<bool> orientation;
					std::istringstream stream(line);

					stream >> c >> walkid >> std::ws; // read 'P', id, and whitespace
					walk_ids.push_back(walkid);

					{
						string walkstring;
						std::getline(stream, walkstring, '\t');
						std::istringstream walkstream(walkstring);
						std::swap(stream, walkstream);
					}
					while (std::getline(stream, node, ',')) {
						assert(node_index.find(node.substr(0, node.size() - 1)) != node_index.end());
						walk.push_back(node_index[node.substr(0, node.size() - 1)]);
						if (node[node.size()-1] == '+')
							orientation.push_back(true);
						else if (node[node.size()-1] == '-')
							orientation.push_back(false);
						else
							assert(false);
					}

					walks.push_back(walk);
					orientations.push_back(orientation);
				} else {
					std::cerr << "Unrecognized line " << line[0] << ": skipping..." << std::endl;
				}
			}
			if (!check()) {
				exit(1);
			}
		}

		bool check()
		{
			if (heights.size() != cuts.size()) {
				std::cerr << "Error: number of cuts and blocks mismatch!" << std::endl;
				return false;
			}

			int block_sum = 0;
			for (auto h : heights) block_sum += h;
			if (block_sum != ordered_node_ids.size()) {
				std::cerr << "Error: sum of block heights does not correspond to node number!" << std::endl;
				return false;
			}

			return true;
		}

		void simplify(bool ignore_only, string &ignorechars)
		{
			// assumption: node indices (0-indexed) respect block description
			vector<pair<int,int>> unary_ranges;
			for (int i = 0; i < heights.size() - 1; i++) {
				if (is_to_simplify(i, ignore_only, ignorechars)) {
					int j = i + 1;
					while (j < heights.size() && is_to_simplify(j, ignore_only, ignorechars)) {
						j += 1;
					}

					unary_ranges.push_back({ i, j });
					i = j;
				}
			}

#ifdef DEBUG
			std::cout << "unary_ranges: ";
			for (auto r : unary_ranges) {
				std::cout << get<0>(r) << ".." << get<1>(r) << ", ";
			}
			std::cout << std::endl;
#endif
			// merge cuts and block heights
			for (auto r : unary_ranges) {
				for (int i = get<0>(r) + 1; i <= get<1>(r); i++) {
					cuts[i] = -1;
					heights[i] = 0;
				}
			}

			// (do not!) erase the edges
			for (auto r : unary_ranges) {
				// block range [i..j]
				// before, unary path that was to be compressed was a,a+1,a+2,...,b
				// now, the paths to be compressed are partitioned in chunks of size heights[i]
				int i = get<0>(r), j = get<1>(r);

				vector<int> perm(heights[i]);
				for (int n = 0; n < heights[i]; n++) {
					// perm[i] is the node that i reaches in block k
					perm[n] = n;
				}

				// propagate perm
				for (int k = i; k < j; k++) {
					vector<int> newperm(heights[i]);

					for (int node = cumulative_height[k]; node < cumulative_height[k+1]; node++) {
						newperm[node - cumulative_height[k]] = perm[edges[node][0] - cumulative_height[k+1]];
					}

					std::swap(perm, newperm);

					for (int l = 0; l < heights[i]; l++) {
						ordered_node_labels[cumulative_height[i] + l] += ordered_node_labels[cumulative_height[k+1] + perm[l]];
						ordered_node_labels[cumulative_height[k+1] + perm[l]].clear();
					}
				}

				for (int n = 0; n < heights[i]; n++) {
					std::swap(edges[n + cumulative_height[i]], edges[perm[n] + cumulative_height[j]]);
				}
				for (int n = cumulative_height[i+1]; n < cumulative_height[j+1]; n++) {
					//edges.erase(n);
				}
			}
		}

		void simplify_tunnels(bool ignore_only, string &ignorechars)
		{
			// assumption: node indices (0-indexed) respect block description
			vector<pair<int,int>> unary_ranges;
			for (int i = 0; i < heights.size() - 1; i++) {
				if (is_to_simplify(i, ignore_only, ignorechars)) {
					int j = i + 1;
					while (j < heights.size() && is_to_simplify(j, ignore_only, ignorechars)) {
						j += 1;
					}

					if (j - i + 1 >= 4) {
						unary_ranges.push_back({ i+1, j-1 });
					}
					i = j;
				}
			}

#ifdef DEBUG
			std::cout << "unary_ranges: ";
			for (auto r : unary_ranges) {
				std::cout << get<0>(r) << ".." << get<1>(r) << ", ";
			}
			std::cout << std::endl;
#endif
			// merge cuts and block heights
			for (auto r : unary_ranges) {
				for (int i = get<0>(r) + 1; i <= get<1>(r); i++) {
					cuts[i] = -1;
					heights[i] = 0;
				}
			}

			// (do not!) erase the edges
			for (auto r : unary_ranges) {
				// block range [i..j]
				// before, unary path that was to be compressed was a,a+1,a+2,...,b
				// now, the paths to be compressed are partitioned in chunks of size heights[i]
				int i = get<0>(r), j = get<1>(r);

				vector<int> perm(heights[i]);
				for (int n = 0; n < heights[i]; n++) {
					// perm[i] is the node that i reaches in block k
					perm[n] = n;
				}

				// propagate perm
				for (int k = i; k < j; k++) {
					vector<int> newperm(heights[i]);

					for (int node = cumulative_height[k]; node < cumulative_height[k+1]; node++) {
						newperm[node - cumulative_height[k]] = perm[edges[node][0] - cumulative_height[k+1]];
					}

					std::swap(perm, newperm);

					for (int l = 0; l < heights[i]; l++) {
						ordered_node_labels[cumulative_height[i] + l] += ordered_node_labels[cumulative_height[k+1] + perm[l]];
						ordered_node_labels[cumulative_height[k+1] + perm[l]].clear();
					}
				}

				for (int n = 0; n < heights[i]; n++) {
					std::swap(edges[n + cumulative_height[i]], edges[perm[n] + cumulative_height[j]]);
				}
				for (int n = cumulative_height[i+1]; n < cumulative_height[j+1]; n++) {
					//edges.erase(n);
				}
			}
		}

		void to_stream(std::ostream *out)
		{
			*out << "M\t" << m << "\t" << n << std::endl;

			*out << "X\t";
			for (auto x : cuts) {
				if (x == -1)
					continue;

				*out << x << "\t";
			}
			*out << std::endl;

			*out << "B\t";
			for (auto h : heights) {
				if (h == 0)
					continue;

				*out << h << "\t";
			}
			*out << std::endl;

			for (int i = 0; i < ordered_node_ids.size(); i++) {
				if (ordered_node_labels[i] == "")
					continue;

				*out << "S\t" << ordered_node_ids[i] << "\t" << ordered_node_labels[i] << std::endl;
			}

			for (int i = 0; i < ordered_node_ids.size(); i++) {
				if (ordered_node_labels[i] == "")
					continue;

				for (auto j : edges[i]) {
					*out << "L\t" << ordered_node_ids[i] << "\t+\t" << ordered_node_ids[j] << "\t+\t0M" << std::endl;
				}
			}

			for (int i = 0; i < walks.size(); i++) {
				*out << "P\t" << walk_ids[i] << "\t";

				bool modified = false;
				int startingnode = walks[i][0];
				if (ordered_node_labels[startingnode] != "") {
					*out << ordered_node_ids[startingnode] << ((orientations[i][0]) ? "+" : "-");
				} else {
					modified = true;
					// go backwards until first non-simplified node that reaches this one
					for (int j = startingnode - 1; j >= 0; j -= 1) {
						// we assume reachability
						if (edges[j][0] = startingnode) {
							// update?
							if (ordered_node_labels[j] != "") {
								*out << ordered_node_ids[j] << ((orientations[i][0]) ? "+" : "-");
								break;
							} else {
								startingnode = j;
							}
						}
					}
				}

				for (int j = 1; j < walks[i].size(); j++) {
					if (j == walks[i].size() - 1 && ordered_node_labels[walks[i][j]] == "")
						modified = true;
					if (ordered_node_labels[walks[i][j]] == "")
						continue;

					*out << "," << ordered_node_ids[walks[i][j]] << ((orientations[i][j]) ? "+" : "-");
				}
				*out << "\t*" << std::endl;

				if (modified) {
					std::cerr << "Warning: string spelled by path " << walk_ids[i] << " gets longer." << std::endl;
				}
			}
		}

		void to_stream_rename(std::ostream *out)
		{
			vector<int> remapped(ordered_node_ids.size(), 0);
			for (int i = 0, n = 0; i < ordered_node_labels.size(); i++) {
				if (ordered_node_labels[i] != "")
					remapped[i] = n++;
			}

			*out << "M\t" << m << "\t" << n << std::endl;

			*out << "X\t";
			for (auto x : cuts) {
				if (x == -1)
					continue;

				*out << x << "\t";
			}
			*out << std::endl;

			*out << "B\t";
			for (auto h : heights) {
				if (h == 0)
					continue;

				*out << h << "\t";
			}
			*out << std::endl;

			for (int i = 0; i < ordered_node_ids.size(); i++) {
				if (ordered_node_labels[i] == "")
					continue;

				*out << "S\t" << remapped[i] << "\t" << ordered_node_labels[i] << std::endl;
			}

			for (int i = 0; i < ordered_node_ids.size(); i++) {
				if (ordered_node_labels[i] == "")
					continue;

				for (auto j : edges[i]) {
					*out << "L\t" << remapped[i] << "\t+\t" << remapped[j] << "\t+\t0M" << std::endl;
				}
			}

			for (int i = 0; i < walks.size(); i++) {
				*out << "P\t" << walk_ids[i] << "\t";

				bool modified = false;
				int startingnode = walks[i][0];
				if (ordered_node_labels[startingnode] != "") {
					*out << remapped[startingnode] << ((orientations[i][0]) ? "+" : "-");
				} else {
					modified = true;
					// go backwards until first non-simplified node that reaches this one
					for (int j = startingnode - 1; j >= 0; j -= 1) {
						// we assume reachability
						if (edges[j][0] = startingnode) {
							// update?
							if (ordered_node_labels[j] != "") {
								*out << remapped[j] << ((orientations[i][0]) ? "+" : "-");
								break;
							} else {
								startingnode = j;
							}
						}
					}
				}

				for (int j = 1; j < walks[i].size(); j++) {
					if (j == walks[i].size() - 1 && ordered_node_labels[walks[i][j]] == "")
						modified = true;
					if (ordered_node_labels[walks[i][j]] == "")
						continue;

					*out << "," << remapped[walks[i][j]] << ((orientations[i][j]) ? "+" : "-");
				}
				*out << "\t*" << std::endl;

				if (modified) {
					std::cerr << "Warning: string spelled by path " << walk_ids[i] << " gets longer." << std::endl;
				}
			}
		}
};

int main(int argc, char* argv[])
{
	gengetopt_args_info argsinfo;
	if (cmdline_parser(argc, argv, &argsinfo) != 0) exit(1);

	if (argsinfo.inputs_num == 0)
		{std::cerr << argv[0] << ": missing input and output file paths" << std::endl; exit(1);};
	if (argsinfo.inputs_num == 1)
		{std::cerr << argv[0] << ": missing output file path" << std::endl; exit(1);};
	if (argsinfo.inputs_num > 2)
		{std::cerr << argv[0] << ": too many arguments" << std::endl; exit(1);};

	// open files
	std::filesystem::path graphpath {argsinfo.inputs[0]};
	std::ifstream graphfs {graphpath};
	if (!graphfs) {std::cerr << "Error opening file " << graphpath << "." << std::endl; exit(1);};

	// check output file
	std::filesystem::path outputpath {argsinfo.inputs[1]};
	std::ofstream outputfs;
	if (std::filesystem::exists(outputpath)) {
		if (argsinfo.overwrite_flag) {
			outputfs = std::ofstream(outputpath, std::ios::out | std::ios::trunc);
		} else {
			std::cerr << "Error: output file already exists." << std::endl;
			exit(1);
		}
	} else {
		outputfs.exceptions(std::fstream::failbit);
		outputfs.open(outputpath, std::ios_base::out);
	}

	std::string ignorechars;
	if (argsinfo.ignore_chars_arg != NULL)
		ignorechars = std::string(argsinfo.ignore_chars_arg);
	else
		ignorechars = "";

	std::cerr << "Reading the graph..." << std::flush;
	Elasticfoundergraph graph(&graphfs);
	std::cerr << " done." << std::endl;

	std::cerr << "Simplifying...";
	if (argsinfo.simplify_tunnels_flag)
		graph.simplify_tunnels(argsinfo.ignore_only_flag, ignorechars);
	else
		graph.simplify(argsinfo.ignore_only_flag, ignorechars);

	std::cerr << " done." << std::endl;

	if (argsinfo.rename_nodes_flag)
		graph.to_stream_rename(&outputfs);
	else
		graph.to_stream(&outputfs);
}
