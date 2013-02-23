
#if defined(CLOSEGRAPH_ST)
#include "closegraph_st.hpp"
#define GSPAN_FUNCTION closegraph_st

#elif defined(CLOSEGRAPH_MT)
#include "closegraph_mt.hpp"
#define GSPAN_FUNCTION closegraph_mt

#else
#error gspan not yet implemented
#endif

#if defined(GRAPH_ADJL)
#include "graph_bgl_adjl_policy.hpp"
#else
#include "graph_bgl_csr_policy.hpp"
#endif

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <boost/ptr_container/ptr_vector.hpp>

unsigned long proj_calls;

using namespace bgl_adaptor;

typedef GraphBGLPolicy Policy;
typedef Policy::graph_t Graph;
typedef Policy::vertex_label_t VL;
typedef Policy::edge_label_t   EL;
typedef Policy::vertex_index_t VI;
typedef gSpan::EdgeCode<Policy> EdgeCode;
typedef gSpan::DFSCode<Policy> DFSCode;
typedef gSpan::SBG<Policy> SBG;
typedef gSpan::SubgraphsOfOneGraph<Policy> SOG;
typedef gSpan::SubgraphsOfManyGraph<Policy> SMG;


struct InputGraph
{
    std::string name;
    std::map<int, std::string> vl;

    struct E
    {
	int from, to;
	std::string el;
    };

    std::vector<E> edges;
};


void read_input(std::list<InputGraph>& igl, std::istream& is)
{
    InputGraph* g = 0;

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
	    igl.push_back(InputGraph());
	    g = &igl.back();
	    g->name = result[2];
        }
	else if (result[0] == "v")
        {
            assert(result.size() == 3);
	    assert(g);
	    g->vl[atoi(result[1].c_str())] = result[2];
        }
	else if (result[0] == "e")
	{
	    g->edges.push_back(InputGraph::E());
	    InputGraph::E& e = g->edges.back();
	    e.from = atoi(result[1].c_str());
	    e.to   = atoi(result[2].c_str());
	    if (result.size() == 4)
		e.el = result[3];
	}

    }
}

void set_labels_from_file(std::vector<std::string>& i_to_str, std::map<std::string, int>& str_to_i, const char* fname)
{
    //BR;
    int count = 0;
    std::ifstream file;
    file.open(fname);
    if (file.is_open())
    {
	std::string line;
	while (file.good())
	{
	    getline(file, line);
	    if (! line.empty())
		str_to_i[line] = count++;
	}
    }
    else
    {
	std::cerr << "ERROR: can not open file: " << fname << std::endl;
	file.close();
	exit(1);
    }
    file.close();

    i_to_str.resize(str_to_i.size());
    for (std::map<std::string, int>::const_iterator i = str_to_i.begin(); i != str_to_i.end(); ++i)
	i_to_str[i->second] = i->first;   
}


bool check_lebels(const std::map<std::string, int>& str_to_i, const std::map<std::string, int>& lab_counts)
{
    for (std::map<std::string, int>::const_iterator i = str_to_i.begin(); i != str_to_i.end(); ++i)
	if (lab_counts.count(i->first) == 0)
	{
	    return false;
	}
    return true;
}

void relabel_by_frequency(std::vector<std::string>& i_to_str,
			  std::map<std::string, int>& str_to_i,
			  const std::map<std::string, int>& lab_counts,
			  bool sort_label_asc)
{
    std::multimap<int, std::string> frq;
    for (std::map<std::string, int>::const_iterator i = lab_counts.begin(); i != lab_counts.end(); ++i)
	frq.insert(std::pair<int, std::string>(i->second, i->first));
    
    int count = 0;
    if (sort_label_asc)
	for (std::multimap<int, std::string>::const_reverse_iterator i = frq.rbegin(); i != frq.rend(); ++i)
	    str_to_i[i->second] = count++;
    else
	for (std::multimap<int, std::string>::const_iterator i = frq.begin(); i != frq.end(); ++i)
	    str_to_i[i->second] = count++;
    
    i_to_str.resize(str_to_i.size());
    for (std::map<std::string, int>::const_iterator i = str_to_i.begin(); i != str_to_i.end(); ++i)
	i_to_str[i->second] = i->first;
}


struct WorkingGraphs
{
    boost::ptr_vector<Policy::graph_t> graphs;
    std::map<const Policy::graph_t*, std::string> names;
};

void create_working_graphs(WorkingGraphs& wrk_graphs,
			   const std::map<std::string, int>& vls_vl,
			   const std::map<std::string, int>& els_el,
			   const std::list<InputGraph>& input_graphs)
{
    BOOST_FOREACH(const InputGraph& ig, input_graphs)
    {
	int corr = ig.vl.begin()->first == 1;
	DFSCode dfsc;
	for (std::vector<InputGraph::E>::const_iterator i = ig.edges.begin(); i != ig.edges.end(); ++i)
	{
	    typedef std::map<std::string, int>::const_iterator I;

	    I iter = vls_vl.find(ig.vl.find(i->from)->second);
	    assert(iter != vls_vl.end());
	    int vl_from_i = iter->second;

	    iter = vls_vl.find(ig.vl.find(i->to)->second);
	    assert(iter != vls_vl.end());
	    int vl_to_i   = iter->second;

	    iter = els_el.find(i->el);
	    assert(iter != els_el.end());
	    int el_i = iter->second;
	    
	    EdgeCode ec(i->from - corr, i->to - corr, vl_from_i, el_i, vl_to_i);
	    dfsc.push_back(ec);
	}
	
	Graph* graph = Policy::create_graph(dfsc);
	wrk_graphs.names[graph] = ig.name;
	wrk_graphs.graphs.push_back(graph);
    }
}


// *****************************************************************************
//                          Result
// *****************************************************************************

enum PatternViewMode { PV_LG, PV_DFSC, PV_EDGE} pattern_view_mode;
enum MappingViewMode { MV_NONE, MV_TABLE, MV_MAP } mapping_view_mode;

bool inline_view;

class Result
{
public:
    Result(std::ostream& ostr, WorkingGraphs& wrk_graphs,
	   std::vector<std::string>& vlabs, std::vector<std::string>& elabs)
	: ostr(ostr),
	  wrk_graphs(wrk_graphs),
	  vlabs(vlabs), elabs(elabs), num_patterns(0) {}

    void operator() (const DFSCode& dfsc, const SMG& smg);
    void operator() (const DFSCode& dfsc, const SOG& sog);
private:
    std::ostream& ostr;

    WorkingGraphs& wrk_graphs;
    std::vector<std::string>& vlabs;
    std::vector<std::string>& elabs;

    int num_patterns;

    void print_pattern(const DFSCode& dfsc) const;

    // print:
    // #support: N
    // #foundin: N, ... 
    void print_info(const SMG& smg) const;

    // print:
    // #support: N
    void print_info(const SOG& sog) const;

    void print_mapping(const DFSCode& dfsc, const SOG& sog) const;
};

void Result::operator() (const DFSCode& dfsc, const SMG& smg)
{
    ++num_patterns;
    print_pattern(dfsc);
    print_info(smg);

    if (mapping_view_mode != MV_NONE)
    {
	for (SMG::const_iterator i = smg.begin(); i != smg.end(); ++i)
	{
	    ostr << "#tograph: " << wrk_graphs.names[i->first] << std::endl;
	    print_mapping(dfsc, i->second);
	}
    }

    ostr << std::endl;
}


void Result::operator() (const DFSCode& dfsc, const SOG& sog)
{
    ++num_patterns;
    print_pattern(dfsc);
    print_info(sog);
    print_mapping(dfsc, sog);

    ostr << std::endl;
}


void Result::print_pattern(const DFSCode& dfsc) const
{
    std::map<VI,std::string> vlm;
    BOOST_FOREACH(const EdgeCode& ec, dfsc)
    {
	if (ec.vl_from != LABEL_NULL) vlm[ec.vi_from] = vlabs[ec.vl_from];
	if (ec.vl_to != LABEL_NULL)   vlm[ec.vi_to]   = vlabs[ec.vl_to];
    }

    switch (pattern_view_mode)
    {
    case PV_LG:
	ostr << "t # " << num_patterns << std::endl;
	for (typename std::map<VI,std::string>::const_iterator i = vlm.begin(); i != vlm.end(); ++i)
	    ostr << "v " << i->first << " " << i->second << std::endl;
	BOOST_FOREACH(const EdgeCode& ec, dfsc)
	    ostr << "e " << ec.vi_from << " " << ec.vi_to << " " << elabs[ec.el] << std::endl;
	break;
	
    case PV_DFSC:
	ostr << "#pattern: " << num_patterns << std::endl;
	BOOST_FOREACH(const EdgeCode& ec, dfsc)
	{
	    ostr << '\t';
	    ostr << "(" << ec.vi_from << "," << ec.vi_to << ", "
		 << vlm[ec.vi_from] << "," << ec.el << "," << vlm[ec.vi_to] << ")" << std::endl;
	}
	break;
	
    case PV_EDGE:
	ostr << "#pattern: " << num_patterns << " : ";
	BOOST_FOREACH(const EdgeCode& ec, dfsc)
	{
	    ostr << "(" << ec.vi_from << "," << ec.vi_to << ", "
		 << vlm[ec.vi_from] << "," << elabs[ec.el] << "," << vlm[ec.vi_to] << ") ";
	}
	if (!inline_view)
	    ostr << std::endl;
	break;

    default:
	break;
    }
}


void Result::print_info(const SMG& smg) const
{
    const int supp = smg.support();

    if (!inline_view)
    {
	ostr << "#support: " << supp;
	if (mapping_view_mode == MV_NONE)
	{
	    ostr << std::endl << "#foundin:";
	    for (SMG::const_iterator it = smg.begin(); it != smg.end(); ++it)
		ostr << " " << wrk_graphs.names[it->first];
	}
    }
    else
    {
	ostr << "support=" << supp << "; ";
	ostr << "foundin:";
	for (SMG::const_iterator it = smg.begin(); it != smg.end(); ++it)
	    ostr << " " << wrk_graphs.names[it->first];
    }
    ostr << std::endl;
}

void Result::print_info(const SOG& sog) const
{
    const int supp = sog.support();
    if (!inline_view)
    {
	ostr << "#support: " << supp << std::endl;
	ostr << "#num subgraphs: " << sog.num_sbgs() << std::endl;
    }
    else
    {
	ostr << "support=" << supp << "; ";
	ostr << "num subgraphs=" << sog.num_sbgs() << "; ";
    }
}


void Result::print_mapping(const DFSCode& dfsc, const SOG& sog) const
{
    const int NUM_EDGES = dfsc.size();

    switch (mapping_view_mode)
    {
    case MV_MAP:
        BOOST_FOREACH(const SBG* s, sog)
        {
            do
            {
                std::map<VI, VI> vvm;
                for (int i = 0; i < NUM_EDGES; ++i)
                {
                    const EdgeCode& ec = dfsc[i];
                    vvm[ec.vi_from] = (*s)[i]->vi_from;
                    vvm[ec.vi_to]   = (*s)[i]->vi_to;
                }
                ostr << "#mv: ";
                for (std::map<VI, VI>::const_iterator i = vvm.begin(); i != vvm.end(); ++i)
                    ostr << i->first << "->" << i->second << " ";
                ostr << std::endl;
                ostr << "#me: ";
                for (int i = 0; i < NUM_EDGES; ++i)
                    ostr << i << "->" << (*s)[i]->ei << " ";
                ostr << std::endl;

            } while ( (s = s->automorph_list));
        }
        break;

    case MV_TABLE:
        BOOST_FOREACH(const SBG* s, sog)
        {
            do
            {
                for (int i = 0; i < NUM_EDGES; ++i)
                    ostr << " (" << (*s)[i]->vi_from << "," << (*s)[i]->vi_to << ")" << (*s)[i]->ei;
                ostr << std::endl;
            } while ( (s = s->automorph_list));
        }
        break;

    default:
        break;
    }
}



// *****************************************************************************
//                          main
// *****************************************************************************
std::ostream& usage(std::ostream& ostr)
{ 
    return ostr << "Usage: CMD <minsup> [{-dfsc -edge}] [{ -m -mt }] [-asc] [-vlfile=FILE]" << std::endl << std::endl;
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

    pattern_view_mode = PV_LG;
    for (int i = 0; i < argc; ++i)
    {
	if (! strcmp(argv[i], "-dfsc"))
	    pattern_view_mode = PV_DFSC;
	else if (! strcmp(argv[i], "-edge"))
	    pattern_view_mode = PV_EDGE;
    }

    
    mapping_view_mode = MV_NONE;
    for (int i = 0; i < argc; ++i)
    {
	if (! strcmp(argv[i], "-m"))
	    mapping_view_mode = MV_MAP;
	else if (! strcmp(argv[i], "-mt"))
	    mapping_view_mode = MV_TABLE;
    }

    bool sort_label_asc = false;
    for (int i = 0; i < argc; ++i)
	if (! strcmp(argv[i], "-asc"))
	    sort_label_asc = true;


    const char* vertex_label_file = 0;
    for (int i = 0; i < argc; ++i)
	if (! strncmp(argv[i], "-vlfile=", 8))
	{
	    vertex_label_file = argv[i]+8;
	    std::cerr << "vertex label file: " << vertex_label_file << std::endl;
	    break;
	}

    const char* edge_label_file = 0;
    for (int i = 0; i < argc; ++i)
	if (! strncmp(argv[i], "-elfile=", 8))
	{
	    edge_label_file = argv[i]+8;
	    std::cerr << "edge label file: " << edge_label_file << std::endl;
	    break;
	}
   

    inline_view = pattern_view_mode == PV_EDGE && mapping_view_mode == MV_NONE;    

    // ------------------------------------------
    // prepare input (transactional) graphs
    // ------------------------------------------
    std::list<InputGraph> input_graphs;
    read_input(input_graphs, std::cin);


    // ------------------------------------------
    // relabel
    // ------------------------------------------
    std::vector<std::string> vlabs;
    std::vector<std::string> elabs;
    std::map<std::string, int> vls_vl;
    std::map<std::string, int> els_el;
    
    std::map<std::string, int> vlab_counts;
    std::map<std::string, int> elab_counts;
    BOOST_FOREACH(const InputGraph& ig, input_graphs)
    {
	for (std::map<int,std::string>::const_iterator i = ig.vl.begin(); i != ig.vl.end(); ++i)
	    ++vlab_counts[i->second];
	for (std::vector<InputGraph::E>::const_iterator i = ig.edges.begin(); i != ig.edges.end(); ++i)
	    ++elab_counts[i->el];
    }

    if (! vertex_label_file)
	relabel_by_frequency(vlabs, vls_vl, vlab_counts, sort_label_asc);
    else
    {
	set_labels_from_file(vlabs, vls_vl, vertex_label_file);
	if (!check_lebels(vls_vl, vlab_counts))
	{
	    std::cerr << "ERROR: vertex labels from file " << vertex_label_file
		      << " not consistent to input graphs" << std::endl;
	    exit(1);
	}
    }

    if (! edge_label_file)
	relabel_by_frequency(elabs, els_el, elab_counts, sort_label_asc);
    else
    {
	set_labels_from_file(elabs, els_el, edge_label_file);
	if (!check_lebels(els_el, elab_counts))
	{
	    std::cerr << "ERROR: vertex labels from file " << edge_label_file
		      << " not consistent to input graphs" << std::endl;
	    exit(1);
	}
    }

    std::cerr << "vertex labels: \n ";
    for (unsigned int i = 0; i < vlabs.size(); ++i)
	std::cerr << "[" << i << "]=" << vlabs[i] << " ";
    std::cerr << "\nedge labels: \n ";
    for (unsigned int i = 0; i < elabs.size(); ++i)
	std::cerr << "[" << i << "]=" << elabs[i] << " ";
    std::cerr << std::endl;

    // ------------------------------------------
    // create working graphs
    // ------------------------------------------
    WorkingGraphs wrk_graphs;
    create_working_graphs(wrk_graphs, vls_vl, els_el, input_graphs);
    

    Result result(std::cout, wrk_graphs, vlabs, elabs);

    if (wrk_graphs.graphs.size() == 1)
	gSpan::GSPAN_FUNCTION<Policy>(wrk_graphs.graphs.back(), minsup, result);
    else
	gSpan::GSPAN_FUNCTION<Policy>(wrk_graphs.graphs.begin(), wrk_graphs.graphs.end(), minsup, result);
}
