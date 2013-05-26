#include "misc.hpp"
#include <fstream>
enum PatternViewMode { PV_LG, PV_DFSC, PV_EDGE} pattern_view_mode;
enum MappingViewMode { MV_NONE, MV_TABLE, MV_MAP } mapping_view_mode;
bool inline_view;

using namespace gSpan;

#if 1
class Result : public GspanResult
{
public:
    Result(std::ostream& ostr, WorkingGraphs& wrk_graphs,
	   std::vector<std::string>& vlabs, std::vector<std::string>& elabs)
	: ostr(ostr),
	  wrk_graphs(wrk_graphs),
	  vlabs(vlabs), elabs(elabs), num_patterns(0) {}

    virtual void operator() (const DFSCode& dfsc, const SubgraphsOfOneGraph& sog);
    virtual void operator() (const DFSCode& dfsc, const SubgraphsOfManyGraph& smg);
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
    void print_info(const SubgraphsOfManyGraph& smg) const;

    // print:
    // #support: N
    void print_info(const SubgraphsOfOneGraph& sog) const;

    void print_mapping(const DFSCode& dfsc, const SubgraphsOfOneGraph& sog) const;
};


void Result::operator() (const DFSCode& dfsc, const SubgraphsOfOneGraph& sog)
{
    ++num_patterns;
    print_pattern(dfsc);
    print_info(sog);
    print_mapping(dfsc, sog);

    ostr << std::endl;
}


void Result::operator() (const DFSCode& dfsc, const SubgraphsOfManyGraph& smg)
{
    ++num_patterns;
    print_pattern(dfsc);
    print_info(smg);

    if (mapping_view_mode != MV_NONE)
    {
	for (SubgraphsOfManyGraph::const_iterator i = smg.begin(); i != smg.end(); ++i)
	{
	    ostr << "#tograph: " << wrk_graphs.names[i->get_graph()] << std::endl;
	    print_mapping(dfsc, *i);
	}
    }

    ostr << std::endl;
}


void Result::print_pattern(const DFSCode& dfsc) const
{
    std::map<DfscVI,std::string> vlm;
    BOOST_FOREACH(const EdgeCode& ec, dfsc)
    {
	if (ec.vl_src() != VL_NULL) vlm[ec.vi_src()] = vlabs[ec.vl_src()];
	if (ec.vl_dst() != VL_NULL) vlm[ec.vi_dst()] = vlabs[ec.vl_dst()];
    }

    switch (pattern_view_mode)
    {
    case PV_LG:
	ostr << "t # " << num_patterns << std::endl;
	for (typename std::map<DfscVI,std::string>::const_iterator i = vlm.begin(); i != vlm.end(); ++i)
	    ostr << "v " << i->first << " " << i->second << std::endl;
	BOOST_FOREACH(const EdgeCode& ec, dfsc)
	    ostr << "e " << ec.vi_src() << " " << ec.vi_dst() << " " << elabs[ec.el()] << std::endl;
	break;
	
    case PV_DFSC:
	ostr << "#pattern: " << num_patterns << std::endl;
	BOOST_FOREACH(const EdgeCode& ec, dfsc)
	{
	    ostr << '\t';
	    ostr << "(" << ec.vi_src() << "," << ec.vi_dst() << ", "
		 << vlm[ec.vi_src()] << "," << ec.el() << "," << vlm[ec.vi_dst()] << ")" << std::endl;
	}
	break;
	
    case PV_EDGE:
	ostr << "#pattern: " << num_patterns << " : ";
	BOOST_FOREACH(const EdgeCode& ec, dfsc)
	{
	    ostr << "(" << ec.vi_src() << "," << ec.vi_dst() << ", "
		 << vlm[ec.vi_src()] << "," << elabs[ec.el()] << "," << vlm[ec.vi_dst()] << ") ";
	}
	if (!inline_view)
	    ostr << std::endl;
	break;

    default:
	break;
    }
}


void Result::print_info(const SubgraphsOfManyGraph& smg) const
{
    const int supp = smg.support();

    if (!inline_view)
    {
	ostr << "#support: " << supp;
	if (mapping_view_mode == MV_NONE)
	{
	    ostr << std::endl << "#foundin:";
	    for (SubgraphsOfManyGraph::const_iterator it = smg.begin(); it != smg.end(); ++it)
		ostr << " " << wrk_graphs.names[it->get_graph()];
	}
    }
    else
    {
	ostr << "support=" << supp << "; ";
	ostr << "foundin:";
	for (SubgraphsOfManyGraph::const_iterator it = smg.begin(); it != smg.end(); ++it)
	    ostr << " " << wrk_graphs.names[it->get_graph()];
    }
    ostr << std::endl;
}

void Result::print_info(const SubgraphsOfOneGraph& sog) const
{
    const int supp = sog.support();
    if (!inline_view)
    {
	ostr << "#support: " << supp << std::endl;
	ostr << "#num subgraphs: " << sog.size_list_all() << std::endl;
	ostr << "#num autgroups: " << sog.size_list_autgroup() << std::endl;
    }
    else
    {
	ostr << "support=" << supp << "; ";
	ostr << "num subgraphs=" << sog.size_list_all() << "; ";
	ostr << "num autgroups=" << sog.size_list_autgroup() << "; ";
    }
}


void Result::print_mapping(const DFSCode& dfsc, const SubgraphsOfOneGraph& sog) const
{
    const int NUM_EDGES = dfsc.size();

    switch (mapping_view_mode)
    {
    case MV_MAP:
        for (const SBG* s = sog.first_all(); s; s = s->next_embedding())
        {
            std::vector<const SBG*> chain;
            get_chain(chain, s);

            std::map<DfscVI, GraphVI> vvm;
            for (int i = 0; i < NUM_EDGES; ++i)
            {
                const EdgeCode& ec = dfsc[i];
                vvm[ec.vi_src()] = chain[i]->edge()->vi_src();
                vvm[ec.vi_dst()] = chain[i]->edge()->vi_dst();
            }
            ostr << "#mv: ";
            for (std::map<DfscVI, GraphVI>::const_iterator i = vvm.begin(); i != vvm.end(); ++i)
                ostr << i->first << "->" << i->second << " ";
            ostr << std::endl;
            ostr << "#me: ";
            for (int i = 0; i < NUM_EDGES; ++i)
                ostr << i << "->" << chain[i]->edge()->eid() << " ";
            ostr << std::endl;
        }
        break;

    case MV_TABLE:
        for (const SBG* s = sog.first_all(); s; s = s->next_embedding())
        {
            std::vector<const SBG*> chain;
            get_chain(chain, s);

            for (int i = 0; i < NUM_EDGES; ++i)
                ostr << " (" << chain[i]->edge()->vi_src() << "," << chain[i]->edge()->vi_dst() << ")"
                     << chain[i]->edge()->eid();
            ostr << std::endl;
        }
        break;

    default:
        break;
    }
}
#endif

std::ostream& usage(std::ostream& ostr)
{ 
    return ostr << "Usage: CMD <minsup> [{-dfsc -edge}] [{ -m -mt }] "
		<< "[-asc] [-vlfile=FILE] [-elfile=FILE] [-trace N] [-if=file]" << std::endl << std::endl;
}


int main(int argc, char** argv)
{
    using namespace gSpan;
    using namespace std;

/*
    MemAllocator ma; int S = 38;
    Bitset bs (&ma, S);

    print(std::cerr, bs, S); std::cerr << std::endl;

    bs.set(0, S);
    print(std::cerr, bs, S);
    std::cerr << " " << bs.bit_.field << std::endl;

    bs.set(1, S); print(std::cerr, bs, S); std::cerr << " " << bs.bit_.field << std::endl;
    bs.set(6, S); print(std::cerr, bs, S); std::cerr << " " << bs.bit_.field << std::endl;
    bs.set(35, S); print(std::cerr, bs, S); std::cerr << " " << bs.bit_.field << std::endl;
    bs.set(8, S); print(std::cerr, bs, S); std::cerr << " " << bs.bit_.field << std::endl;
    bs.set(15, S); print(std::cerr, bs, S); std::cerr << " " << bs.bit_.field << std::endl;
    bs.set(36, S); print(std::cerr, bs, S); std::cerr << " " << bs.bit_.field << std::endl;
    bs.set(31, S); print(std::cerr, bs, S); std::cerr << " " << bs.bit_.field << std::endl;
    bs.set(11, S); print(std::cerr, bs, S); std::cerr << " " << bs.bit_.field << std::endl;

    return 0;


    std::cerr << "SBG size=" << sizeof(SBG) << std::endl;

    EdgeCode ec1(1, 2, 0, 2, 2, true);
    EdgeCode ec2(1, 3, 0, 2, 0, true);
    EdgeCodeCmpLex lex;
    if (lex(ec1,ec2))
        cout << "less\n";
    else
        cout << "NOT less\n";
    
    return 0;
*/
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
   
    int tracedepth = 0;
    for (int i = 0; i < argc; ++i)
	if (! strncmp(argv[i], "-trace=", 7))
	{
	    tracedepth = atoi(argv[i]+7);
	    break;
	}


    const char* inputfile = 0;
    for (int i = 0; i < argc; ++i)
	if (! strncmp(argv[i], "-if=", 4))
	{
	    inputfile = argv[i]+4;
	    break;
	}


    inline_view = pattern_view_mode == PV_EDGE && mapping_view_mode == MV_NONE;    

    // ------------------------------------------
    // prepare input (transactional) graphs
    // ------------------------------------------
    std::list<InputGraph> input_graphs;
    if (inputfile)
    {
	std::ifstream ifile;
	ifile.open(inputfile);
	if (ifile.is_open())
	    read_input(input_graphs, ifile);
	else
	{
	    std::cerr << "ERROR: can not open file: " << inputfile << std::endl;
	    return 1;
	}
    }
    else
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

    // vertex labels
    if (! vertex_label_file)
	relabel(vlabs, vls_vl, vlab_counts, sort_label_asc);
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

    // edge labels
    if (! edge_label_file)
	relabel(elabs, els_el, elab_counts, sort_label_asc);
    else
    {
	set_labels_from_file(elabs, els_el, edge_label_file);
	if (!check_lebels(els_el, elab_counts))
	{
	    std::cerr << "ERROR: edge labels from file " << edge_label_file
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

    minsup=minsup;

    Result result(std::cout, wrk_graphs, vlabs, elabs);

    if (wrk_graphs.graphs.size() == 1)
	closegraph(*wrk_graphs.graphs.back(), minsup, &result, tracedepth);
    else
	closegraph(wrk_graphs.graphs, minsup, &result, tracedepth);

}
