#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>

#include "efg.hpp"

using namespace std;

int main(int argc, char * argv[])
{
	if (argc < 2)
	{
		cout << "usage: " << argv[0] << " graph.gfa" << std::endl;
		return 1;
	}

	// open graph file
	std::filesystem::path graphpath {argv[1]};
	std::ifstream graphfs = std::ifstream {graphpath};
	if (!graphfs) {std::cerr << "Error opening graph file " << graphpath << "." << std::endl; exit(1);};

	Elasticfoundergraph graph(graphfs);
	vector<bool> is_source(graph.ordered_node_ids.size() + 1, true);
	vector<bool> is_sink(graph.ordered_node_ids.size() + 1, true);
	for (int i = 0; i < graph.ordered_node_ids.size(); i++) {
		for (int j : graph.edges[i]) {
			is_sink[i] = false;
			is_source[j] = false;
		}
	}

	for (int i = 0; i < graph.ordered_node_ids.size(); i++) {
		if (!is_source[i] and !is_sink[i]) {
			std::cout << graph.ordered_node_labels[i] << "\n";
			std::cerr << graph.ordered_node_ids[i] << "\n";
		}
	}
}
