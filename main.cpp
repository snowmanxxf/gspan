
#include "gspan.hpp"
#include "graph_ops.hpp"

#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <string>
#include <set>
#include <boost/graph/graphviz.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

using namespace graph_alg;

typedef GraphOpsDefault<std::string, std::string, boost::directedS >      GraphOpsDir;
typedef GraphOpsDefault<std::string, std::string, boost::bidirectionalS > GraphOpsBidir;

// ------------------ READ stdin ---------------------------------------
template<class GraphOps>
std::istream& contruct_dfsc(typename GraphOps::DFSCode& dfsc, int& tr_number, std::istream& is)
{
    typedef typename GraphOps::EdgeCode EdgeCode;
    typedef typename GraphOps::DFSCode DFSCode;
    typedef typename GraphOps::vertex_label_t VL;
    typedef typename GraphOps::edge_label_t   EL;
    typedef typename GraphOps::vertex_index_t VI;

    std::map<VI, VL> vlabels;

    char line[1024];
    while (true)
    {
        std::streampos pos = is.tellg();
        if (! is.getline(line, 1024))
            break;

        std::vector<std::string> result;       
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
            vlabels[atoi(result[1].c_str())] = result[2];
        }
        else if (result[0] == "e")
        {
            assert(result.size() == 4);
            VI from   = atoi(result[1].c_str());
            VI to     = atoi(result[2].c_str());
            EL elabel = result[3];
            dfsc.push(EdgeCode(from, to, vlabels[from], elabel, vlabels[to]));
        }
    }
    return is;
}


// ------------------ RESULT -------------------------------------------
template<class GraphOps>
class Result
{
public:
    typedef typename GraphOps::graph_t Graph;
    typedef typename GraphOps::EdgeCode EdgeCode;
    typedef typename GraphOps::DFSCode DFSCode;
    typedef typename GraphOps::vertex_label_t VL;
    typedef typename GraphOps::edge_label_t   EL;
    typedef typename GraphOps::vertex_index_t VI;
    
    Result(std::map<const Graph*, int>& tr_numbers)
	: tr_numbers(tr_numbers) {}

    void operator() (const Projected_<GraphOps>& projected, const DFSCode& dfsc)
	{
	    std::cout << "result: " << dfsc << "\tat ";
	    std::set<const Graph*> n;
	    BOOST_FOREACH(const typename Traits3<GraphOps>::SBG* sbg, projected)
		n.insert(sbg->get_graph_p());
	    BOOST_FOREACH(const Graph* g, n)
		std::cout << tr_numbers[g] << ",";
	    std::cout << std::endl;
	}
private:
    std::map<const Graph*, int>& tr_numbers;
};


// ---------------------------------------------------------------------

template<class GraphOps>
void run(std::istream& is, int minsup)
{
    typedef typename GraphOps::graph_t Graph;
    typedef typename GraphOps::EdgeCode EdgeCode;
    typedef typename GraphOps::DFSCode DFSCode;
    typedef typename GraphOps::vertex_label_t VL;
    typedef typename GraphOps::edge_label_t   EL;
    typedef typename GraphOps::vertex_index_t VI;

    GraphOps ops;
    boost::ptr_vector<Graph> gr;
    std::map<const Graph*, int> tr_numbers;
    Result<GraphOps> result(tr_numbers);
    unsigned int skipped = 0;

    while (true)
    {
	int tr_number = -1;
	DFSCode dfsc;
	contruct_dfsc<GraphOps>(dfsc, tr_number, is);
	if (dfsc.empty())
	    break;
	try {
	    Graph* graph = ops.create_graph(dfsc);
	    gr.push_back(graph);
	    tr_numbers[graph] = tr_number;

	    std::cerr << "INFO:    Graph " << tr_number << " was created"
		      << " at address " << graph << std::endl;
	}
	catch (typename GraphOps::VertexNotLabeledException e)
	{
	    ++skipped;
	    std::cerr << "WARNING: Graph " << tr_number << " not created, vertex "
		      << e.vertex_index << " not labeled" << std::endl;
	}
    }

    std::cerr << "Transactional graphs: " << gr.size() << " created, "
	      << skipped << " skipped" << std::endl;

    gspan(gr.begin(), gr.end(), minsup, ops, result);
}

std::ostream& usage(std::ostream& ostr)
{
    return ostr << "Usage: gspan <minsup> [-dir]" << std::endl << std::endl;
}

int main(int argc, char** argv)
{
    if (argc < 2)
    {
	usage(std::cerr);
	return 1;
    }

    int minsup = atoi(argv[1]);
    bool directed = (argc == 3 && std::string(argv[2]) == "-dir") ? true : false;

    if (directed)
    {
	std::cerr << "directed" << std::endl;
	run<GraphOpsDir>(std::cin, minsup);
    }
    else
    {
	std::cerr << "undirected" << std::endl;
	run<GraphOpsBidir>(std::cin, minsup);
    }
}
