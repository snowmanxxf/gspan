
#include "graph_policy.hpp"
#include "gspan.hpp"

#include <iostream>
#include <iomanip>
#include <boost/ptr_container/ptr_vector.hpp>

//using namespace gSpan;

static bool verbose = false;
static bool print_dfscode = false;

typedef GraphPolicy Policy;
typedef Policy::graph_t Graph;
typedef Policy::vertex_label_t VL;
typedef Policy::edge_label_t   EL;
typedef Policy::vertex_index_t VI;

typedef gSpan::EdgeCode<Policy> EdgeCode;
typedef gSpan::DFSCode<Policy> DFSCode;
typedef gSpan::SBG<Policy> SBG;
typedef gSpan::Projected<Policy> Projected;

std::vector<std::string> result;
bool f;

std::istream& contruct_dfsc(DFSCode& dfsc, std::string& tr_name, std::istream& is)
{
    std::map<VI, VL> vlabels;
    
    char line[1024];
    while (true)
    {
        std::streampos pos = is.tellg();

	if (!f)
	{
	    if (! is.getline(line, 1024))
		break;

	    result.clear();
	    char* p = strtok(line, " \t");
	    while (p)
	    {
		result.push_back(std::string(p));
		p = strtok(0, " \t");
	    }
	}        
	else
	    f = false;

        if (result.empty())
            continue;

        if (result[0] == "t")
        {
            if (!dfsc.empty())
            {
                //is.seekg(pos, std::ios_base::beg);
		f = true;
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
            dfsc.push_back(ec);
        }
    }
    return is;
}


class Result
{
public:
    Result(std::ostream& ostr, std::map<const Graph*, std::string>& tr_names, const Policy& ops)
	: ostr(ostr), tr_names(tr_names), ops(ops), ngraph(0), num_patterns(0) {}

    void operator() (const DFSCode& dfsc, const Projected& projected)
	{
	    ++num_patterns;
	    /*if (num_patterns == 1)
	      return;
	      if (num_patterns == 4)
	      exit(0);*/
	    if (print_dfscode)
		print_dfsc(dfsc, projected);
	    else
		print_tgf(dfsc, projected);
	}
private:
    std::ostream& ostr;
    std::map<const Graph*, std::string>& tr_names;
    const Policy& ops;
    int ngraph;
    int num_patterns;
    
    void print_dfsc(const DFSCode& dfsc, const Projected& projected);
    void print_tgf(const DFSCode& dfsc, const Projected& projected);
};

void Result::print_dfsc(const DFSCode& dfsc, const Projected& projected)
{
    ostr << std::setw(2) << num_patterns << ": supp=" << projected.mgsbg_size() << ": "
	 << dfsc << std::endl;
    if (verbose)
    {
	BOOST_FOREACH(const SBG& sbg, projected)
	    ostr << "\t" << sbg << std::endl;
    }
}

void Result::print_tgf(const DFSCode& dfsc, const Projected& projected)
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
    BOOST_FOREACH(const SBG& sbg, projected)
	gg.insert(sbg.get_graph());
    ostr << "#found_in " << projected.mgsbg_size() << ": ";
    for (typename std::set<const Graph*>::const_iterator i = gg.begin(); i != gg.end(); ++i)
    {
	ostr << tr_names[*i];
	typename std::set<const Graph*>::const_iterator tmpi = i;
	if (++tmpi != gg.end()) ostr << ", ";
	else                    ostr << std::endl;
    }

    if (verbose)
    {
	BOOST_FOREACH(const SBG& sbg, projected)
	    ostr << "#\t" << sbg << std::endl;
    }
    ostr << std::endl;
}


std::ostream& usage(std::ostream& ostr)
{ 
    return ostr << "Usage: CMD <minsup> [-dfsc] -v" << std::endl << std::endl;
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
    Policy pl;
    boost::ptr_vector<GraphPolicy::graph_t> gr_trans;
    std::map<const Policy::graph_t*, std::string> tr_names;

    {
#ifdef DEBUG_CHECK_GRAPH_LABEL
	unsigned int skipped = 0;
#endif
	while (true)
	{
	    std::string tr_name;
	    DFSCode dfsc;
	    contruct_dfsc(dfsc, tr_name, std::cin);
	    if (dfsc.empty())
		break;
#ifdef DEBUG_CHECK_GRAPH_LABEL
	    try
#endif
	    {
		Policy::graph_t* graph = pl.create_graph(dfsc);
		gr_trans.push_back(graph);
		tr_names[graph] = tr_name;

		if (verbose)
		    std::cerr << "INFO:    Graph " << tr_name << " was created"
			      << " at address " << graph << std::endl;
	    }
#ifdef DEBUG_CHECK_GRAPH_LABEL
	    catch (Policy::VertexNotLabeledException e)
	    {
		++skipped;
		if (verbose)
		    std::cerr << "WARNING: Graph " << tr_name << " not created, vertex "
			      << e.vertex_index << " not labeled" << std::endl;
	    }
#endif
	}
    }
	
    Result result(std::cout, tr_names, pl);
    gSpan::GSPAN_FUNCTION(gr_trans.begin(), gr_trans.end(), minsup, pl, result);
}
