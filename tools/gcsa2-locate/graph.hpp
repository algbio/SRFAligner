#ifndef GRAPH_HPP_LOCATE
#define GRAPH_HPP_LOCATE

#include<string>
#include<fstream>
#include<unordered_map>

//#define GRAPH_HPP_DEBUG

using std::vector, std::istringstream, std::unordered_map, std::stol;

namespace graph {

class SimpleGraph {
	private:
		long int nodes = 0;
		unordered_map<unsigned long,string> ids_to_labels;

	public:
		SimpleGraph(std::ifstream &graphstream)
		{
			long int nodes = 0;

			for (string line; std::getline(graphstream, line); ) {
				if (line[0] == 'S') {
					char c;
					string id, label;

					istringstream(line) >> c >> id >> label;
					ids_to_labels[stol(id)] =label;
					nodes += 1;
				} else if (line[0] == 'L') {
					// do nothing
				} else if (line[0] == 'P') {
					// do nothing
				} else {
					std::cerr << "Unrecognized line " << line[0] << ": skipping..." << std::endl;
				}
			}
			this->nodes = nodes;
		}

		long int get_label_length(unsigned long n) const {
			assert(n < ids_to_labels.size());
			return ids_to_labels.at(n).length();
		}
};

} // namespace graph
#endif // GRAPH_HPP_LOCATE
