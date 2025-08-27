#ifndef GAFANCHOR_HPP
#define GAFANCHOR_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <filesystem>
#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <iterator>
#include <unordered_map>
#include <utility> // std::tie, <
#include <sdsl/bit_vectors.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/util.hpp> // init_support for rank queries
#include <sdsl/config.hpp> // util things

#include "graph.hpp"

namespace gafanchor {

	using std::string, std::vector, graph::SimpleGraph;

class GAFAnchor {
	private:
		// query info
		string qname;
		long unsigned int qlength, qstart, qend;
		// path info
		bool pstrand; // true for +, false for -
		vector<long unsigned int> path; // path ids
		vector<bool> orientations; // true for +, false for -
		long unsigned int plength, pstart, pend;
		// we ignore residue matches, alignment block length, mapping quality


	public:
		GAFAnchor(
				string qname,
				long unsigned int qlength, 
				long unsigned int qstart, 
				long unsigned int qend, 
				vector<long unsigned int> path, 
				long unsigned int plength,
				long unsigned int pstart,
				long unsigned int pend,
				bool reverse = false
			 )
		{
			this->qname = qname;
			this->qlength = qlength;
			this->qstart = qstart;
			this->qend = qend;
			this->pstrand = true;
			this->path = path;
			if (reverse)
				this->orientations = vector<bool>(path.size(), false);
			else
				this->orientations = vector<bool>(path.size(), true);
			this->plength = plength;
			this->pstart = pstart;
			this->pend = pend;

			assert(pend <= plength);
		}

		// empty constructor
		GAFAnchor()
		{
		}

		long unsigned int get_query_start() const
		{
			return qstart;
		}

		string get_query_id() const
		{
			return qname;
		}

		long unsigned int get_path_length() const
		{
			return path.size();
		}

		long unsigned int get_length() const
		{
			return qend - qstart;
		}

		void reverse()
		{
			long unsigned int newqstart = qlength - qend;
			long unsigned int newqend = qlength - qstart;

			qstart = newqstart;
			qend = newqend;

			assert(plength >= pend);
			long unsigned int newpstart = plength - pend;
			assert(plength >= pstart);
			long unsigned int newpend = plength - pstart;

			pstart = newpstart;
			pend = newpend;

			std::reverse(path.begin(), path.end());
			for (long unsigned int i = 0; i < orientations.size(); i++)
				orientations[i] = !orientations[i];
		}

		string to_string()
		{
			string out;
			// GAF format
			out += qname + "\t";                                  // query name
			out += std::to_string(qlength) + "\t";                // query length
			out += std::to_string(qstart) + "\t";                 // query start (0-based, closed)
			out += std::to_string(qend) + "\t";                   // query end (open)
			out += std::string() + ((pstrand) ? "+" : "-") + "\t"; // strand

			for (long unsigned int i = 0; i < path.size(); i++) { // path
				out += std::string() + ((orientations[i]) ? ">" : "<"); 
				out += std::to_string(path[i]);
			}
			out += "\t";

			out += std::to_string(plength) + "\t"; // path length
			out += std::to_string(pstart) + "\t";  // start position on the path
			out += std::to_string(pend) + "\t";    // end position on the path
			out += "0\t0\t255";    // FIXME

			return out;
		}
};

void anchors_to_stream(std::ostream *out, const vector<GAFAnchor> &matches)
{
	for (GAFAnchor m : matches) {
		*out << m.to_string() << std::endl;
	}
}

} // namespace gafanchor
#endif // GAFANCHOR_HPP
