
#include "graph_ops.hpp"

#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <string>
#include <set>
#include <boost/graph/graphviz.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

using namespace graph_alg;

typedef GraphOpsDefault<std::string, std::string, boost::bidirectionalS > GraphOpsBidir;

bool verbose = false;
bool print_dfscode = false;

// ------------------ READ stdin ---------------------------------------
template<class GraphOps>
std::istream& contruct_dfsc(typename GraphOps::DFSCode& dfsc, std::string& tr_name, std::istream& is)
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
            
            tr_name = result[2];
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
	    EdgeCode ec(from, to, vlabels[from], elabel, vlabels[to]);
            dfsc.push(ec);
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
    typedef typename gSpan::Traits<GraphOps>::Projected Projected;
    typedef typename GraphOps::vertex_label_t VL;
    typedef typename GraphOps::edge_label_t   EL;
    typedef typename GraphOps::vertex_index_t VI;
    
    Result(std::ostream& ostr, std::map<const Graph*, std::string>& tr_names, const GraphOps& ops)
	: ostr(ostr), tr_names(tr_names), ops(ops), ngraph(0), num_patterns(0) {}

    void print_dfsc(const DFSCode& dfsc, const Projected& projected);
    void print_tgf(const DFSCode& dfsc, const Projected& projected);

    void operator() (const Projected& projected, const DFSCode& dfsc)
	{
	    ++num_patterns;
	    if (print_dfscode)
		print_dfsc(dfsc, projected);
	    else
		print_tgf(dfsc, projected);
	}
private:
    std::ostream& ostr;
    std::map<const Graph*, std::string>& tr_names;
    const GraphOps& ops;
    int ngraph;

    int num_patterns;
};

template<class GraphOps>
void Result<GraphOps>::print_dfsc(const DFSCode& dfsc, const Projected& projected)
{
    ostr << std::setw(2) << num_patterns << ": " << dfsc << std::endl;

    if (verbose)
    {
	BOOST_FOREACH(const typename Projected::value_type& sbg, projected)
	    ostr << "\t" << sbg << std::endl;
	//ostr << "\t" << projected.front() << std::endl;
    }
}

template<class GraphOps>
void Result<GraphOps>::print_tgf(const DFSCode& dfsc, const Projected& projected)
{
    std::map<VI,VL> vlabels;
    
    BOOST_FOREACH(const EdgeCode& ec, dfsc)
    {
	if (!ops.void_vlabel(ec.vl_from))vlabels[ec.vi_from] = ec.vl_from;
	if (!ops.void_vlabel(ec.vl_to))  vlabels[ec.vi_to] = ec.vl_to;
    }
    
    ostr << "t # " << ++ngraph << std::endl;

    for (typename std::map<VI,VL>::const_iterator i = vlabels.begin(); i != vlabels.end(); ++i)
	ostr << "v " << i->first << " " << i->second << std::endl;

    BOOST_FOREACH(const EdgeCode& ec, dfsc)
	ostr << "e " << ec.vi_from << " " << ec.vi_to << " " << ec.el << std::endl;


    std::set<const Graph*> gg;
    BOOST_FOREACH(const typename gSpan::Traits<GraphOps>::SBG& sbg, projected)
	gg.insert(sbg.get_graph());
    ostr << "#found_in: ";
    for (typename std::set<const Graph*>::const_iterator i = gg.begin(); i != gg.end(); ++i)
    {
	ostr << tr_names[*i];
	typename std::set<const Graph*>::const_iterator tmpi = i;
	if (++tmpi != gg.end()) ostr << ", ";
	else                    ostr << std::endl;
    }
    ostr << std::endl;
}


// ---------------------------------------------------------------------

std::ostream& usage(std::ostream& ostr)
{ 
    return ostr << "Usage: gspan <minsup> [-dfsc] -v" << std::endl << std::endl;
}

int main(int argc, char** argv)
{
    // ------------------------------------------
    // parse arguments
    // ------------------------------------------
    if (argc < 2)
    {
	usage(std::cerr);
	return 1;
    }

    int minsup = atoi(argv[1]);

    print_dfscode =
	(argc > 1 && std::string(argv[1]) == "-dfsc") ||
	(argc > 2 && std::string(argv[2]) == "-dfsc") ||
	(argc > 3 && std::string(argv[3]) == "-dfsc");

    verbose =
	(argc > 1 && std::string(argv[1]) == "-v") ||
	(argc > 2 && std::string(argv[2]) == "-v") ||
	(argc > 3 && std::string(argv[3]) == "-v");

    // ------------------------------------------
    // prepare input (transactional) graphs
    // ------------------------------------------
    GraphOpsBidir ops;
    boost::ptr_vector<GraphOpsBidir::graph_t> gr;
    std::map<const GraphOpsBidir::graph_t*, std::string> tr_names;

    {
	unsigned int skipped = 0;
	while (true)
	{
	    std::string tr_name;
	    GraphOpsBidir::DFSCode dfsc;
	    contruct_dfsc<GraphOpsBidir>(dfsc, tr_name, std::cin);
	    if (dfsc.empty())
		break;
	    try {
		GraphOpsBidir::graph_t* graph = ops.create_graph(dfsc);
		gr.push_back(graph);
		tr_names[graph] = tr_name;

		if (verbose)
		    std::cerr << "INFO:    Graph " << tr_name << " was created"
			      << " at address " << graph << std::endl;
	    }
	    catch (GraphOpsBidir::VertexNotLabeledException e)
	    {
		++skipped;
		if (verbose)
		    std::cerr << "WARNING: Graph " << tr_name << " not created, vertex "
			      << e.vertex_index << " not labeled" << std::endl;
	    }
	}

	if (verbose)
	    std::cerr << "Transactional graphs: " << gr.size() << " created, "
		      << skipped << " skipped" << std::endl;
    }

    // ------------------------------------------
    // run
    // ------------------------------------------
    Result<GraphOpsBidir> result(std::cout, tr_names, ops);
    gSpan::gspan(gr.begin(), gr.end(), minsup, ops, result);

    CloseGraph::closegraph(gr.begin(), gr.end(), minsup, ops, result);
}
