
#include "gspan.hpp"
#include "graph_ops.hpp"

#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <string>
#include <boost/graph/graphviz.hpp>
#include <boost/ptr_container/ptr_vector.hpp>


using namespace graph_alg;

typedef GraphOpsDefault<char, char,
			boost::directedS
			//boost::bidirectionalS
			> GraphOps;
typedef GraphOps::DFSCode DFSCode;
typedef GraphOps::vertex_index_t VI;

GraphOps ops;

// ------------------ READ stdin ---------------------------------------

template <class T, class Iterator>
void tokenize (char *str, Iterator iterator)
{
    char* p = str;
    while ( (p = strtok(p, " \t")))
	*iterator++ = std::string(p);
}


std::istream& contruct_dfsc(DFSCode& dfsc, int& tr_number, std::istream& is)
{
    std::map<VI, char> vlabels;
    //BR;
    char line[1024];
    while (true)
    {
	std::streampos pos = is.tellg();
	if (! is.getline(line, 1024))
	    break;

	std::vector <std::string> result;	
	char* p = strtok(line, " \t");
	while (p)
	{
	    result.push_back(std::string(p));
	    p = strtok(0, " \t");
	}
	
	if (result.empty())
	    continue;

	if (result[0] == "t")
	{
	    if (!dfsc.empty())
	    {
		is.seekg(pos, std::ios_base::beg);
		return is;
	    }
	    
	    tr_number = atoi(result[2].c_str());
	}
	else if (result[0] == "v")
	{
	    assert(result.size() == 3);
	    vlabels[atoi(result[1].c_str())] = result[2][0];
	}
	else if (result[0] == "e")
	{
	    assert(result.size() == 4);
	    VI from   = atoi(result[1].c_str());
	    VI to     = atoi(result[2].c_str());
	    char elabel = result[3][0];	    
	    ops.dfsc_add(dfsc, from, to, vlabels[from], elabel, vlabels[to]);
	}
    }
    return is;
}

// ------------------ RESULT -------------------------------------------
std::map<const GraphOps::graph_t*, int> tr_numbers;
int result;

// ---------------------------------------------------------------------


int main(int argc, char** argv)
{
    boost::ptr_vector<GraphOps::graph_t> gr;

    while (true)
    {
	DFSCode dfsc;
	int tr_number = -1;
	contruct_dfsc(dfsc, tr_number, std::cin);
	gr.push_back(ops.create_graph(dfsc));
	tr_numbers[&gr.back()] = tr_number;
	if (dfsc.empty())
	    break;
    }

    gspan(gr.begin(), gr.end(), atoi(argv[1]), ops, result);
}
