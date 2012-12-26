
#include "graph_policy.hpp"
#include "gspan.hpp"

#include <iostream>
#include <iomanip>
#include <cstring>
#include <boost/ptr_container/ptr_vector.hpp>

typedef GraphPolicy Policy;
typedef Policy::graph_t Graph;
typedef Policy::vertex_label_t VL;
typedef Policy::edge_label_t   EL;
typedef Policy::vertex_index_t VI;

typedef gSpan::EdgeCode<Policy> EdgeCode;
typedef gSpan::DFSCode<Policy> DFSCode;
typedef gSpan::SBG<Policy> SBG;
typedef gSpan::ProjectedManyGraph_<Policy> ProjectedManyGraph;
typedef gSpan::ProjectedOneGraph_<Policy>  ProjectedOneGraph;

std::vector<std::string> result;
bool f;

std::istream& contruct_dfsc(DFSCode& dfsc, std::string& tr_name, std::istream& is)
{
    std::map<VI, VL> vlabels;
    int vcorr = 0;
    bool first_vertex = true;
    
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

	    VI vi = atoi(result[1].c_str());
	    if (first_vertex && vi == 1)
		vcorr = 1;
	    vi -= vcorr;
            vlabels[vi] = result[2];
	    first_vertex = false;
        }
        else if (result[0] == "e")
        {
            assert(result.size() == 4 || result.size() == 3);
	    
            VI from   = atoi(result[1].c_str()) - vcorr;
            VI to     = atoi(result[2].c_str()) - vcorr;

	    EL elabel = "_";
	    if (result.size() == 4)
		elabel = result[3];

	    EdgeCode ec(from, to, vlabels[from], elabel, vlabels[to]);
            dfsc.push_back(ec);
        }
    }
    return is;
}

enum PatternViewMode { PV_LG, PV_DFSC, PV_EDGE} pattern_view_mode;
enum MappingViewMode { MV_NONE, MV_TABLE, MV_MAP } mapping_view_mode;

// *****************************************************************************
//                          ResultOneGraph
// *****************************************************************************
class ResultOneGraph
{
public:
    ResultOneGraph(std::ostream& ostr, const Policy& pl)
	: ostr(ostr), pl(pl), num_patterns(0) {}
    void operator() (const DFSCode& dfsc, const ProjectedOneGraph& projected);
private:
    std::ostream& ostr;
    const Policy& pl;
    int num_patterns;
    void print_info(const ProjectedOneGraph& projected, bool onnewline, bool ismap_view) const;
};

void ResultOneGraph::operator() (const DFSCode& dfsc, const ProjectedOneGraph& projected)
{
    ++num_patterns;
    
    const int NUM_EDGES = dfsc.size();

    std::vector<int> col_width(NUM_EDGES, 8);
    std::map<VI,VL> vlm;

    bool inline_view = false;

    BOOST_FOREACH(const EdgeCode& ec, dfsc)
    {
	if (!pl.void_vlabel(ec.vl_from)) vlm[ec.vi_from] = ec.vl_from;
	if (!pl.void_vlabel(ec.vl_to))   vlm[ec.vi_to]   = ec.vl_to;
    }

    switch (pattern_view_mode)
    {
    case PV_LG:
	ostr << "t # " << num_patterns << std::endl;
	for (typename std::map<VI,VL>::const_iterator i = vlm.begin(); i != vlm.end(); ++i)
	    ostr << "v " << i->first << " " << i->second << std::endl;
	BOOST_FOREACH(const EdgeCode& ec, dfsc)
	    ostr << "e " << ec.vi_from << " " << ec.vi_to << " " << ec.el << std::endl;
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
	inline_view = (mapping_view_mode == MV_NONE);

	ostr << "#pattern: " << num_patterns << " : ";
	BOOST_FOREACH(const EdgeCode& ec, dfsc)
	{
	    ostr << "(" << ec.vi_from << "," << ec.vi_to << ", "
		 << vlm[ec.vi_from] << "," << ec.el << "," << vlm[ec.vi_to] << ") ";
	}
	if (!inline_view)
	    ostr << std::endl;
	break;

    default:
	break;
    }


    print_info(projected, inline_view, mapping_view_mode != MV_NONE);

    const ProjectedOneGraph::SBGSconst_iterator it_p_end = projected.sbgs_end();

    switch (mapping_view_mode)
    {
    case MV_MAP:
	for (ProjectedOneGraph::SBGSconst_iterator it_p = projected.sbgs_begin();
	     it_p != it_p_end; ++it_p)
	{
	    const SBG* s = *it_p;
	    do
	    {
		std::map<VI, VI> vvm;
		for (int i = 0; i < NUM_EDGES; ++i)
		{
		    const EdgeCode& ec = dfsc[i];
		    vvm[ec.vi_from] = (*s)[i].vi_from;
		    vvm[ec.vi_to]   = (*s)[i].vi_to;
		}
		ostr << "#mv: ";
		for (std::map<VI, VI>::const_iterator i = vvm.begin(); i != vvm.end(); ++i)
		    ostr << i->first << "->" << i->second << " ";
		ostr << std::endl;
		ostr << "#me: ";
		for (int i = 0; i < NUM_EDGES; ++i)
		    ostr << i << "->" << (*s)[i].ei << " ";
		ostr << std::endl;

	    } while ( (s = s->automorph_list));
	}
	break;

    case MV_TABLE:
	for (ProjectedOneGraph::SBGSconst_iterator it_p = projected.sbgs_begin();
	     it_p != it_p_end; ++it_p)
	{
	    const SBG* s = *it_p;
	    do
	    {
		for (int i = 0; i < NUM_EDGES; ++i)
		    ostr << " (" << (*s)[i].vi_from << "," << (*s)[i].vi_to << ")" << (*s)[i].ei;
		ostr << std::endl;
	    } while ( (s = s->automorph_list));
	}
	break;

    default:
	break;
    }

    ostr << std::endl;

}


void ResultOneGraph::print_info(const ProjectedOneGraph& projected, bool onnewline, bool ismap_view) const
{
    const int supp = projected.support();

    if (!onnewline)
    {
	ostr << "#support: " << supp;
	ostr << std::endl;
    }
    else
    {
	ostr << "support=" << supp << "; ";
    }
}


// *****************************************************************************
//                          ResultManyGraph
// *****************************************************************************
class ResultManyGraph
{
public:
    ResultManyGraph(std::ostream& ostr, std::map<const Graph*, std::string>& tr_names, const Policy& pl)
	: ostr(ostr), tr_names(tr_names), pl(pl), num_patterns(0) {}
    void operator() (const DFSCode& dfsc, const ProjectedManyGraph& projected);
private:
    std::ostream& ostr;
    std::map<const Graph*, std::string>& tr_names;
    const Policy& pl;
    int num_patterns;

    void print_info(const ProjectedManyGraph& projected, bool onnewline, bool ismap_view) const;
};

void ResultManyGraph::operator() (const DFSCode& dfsc, const ProjectedManyGraph& projected)
{
    ++num_patterns;
    
    const int NUM_EDGES = dfsc.size();
    const ProjectedManyGraph::MGSBGconst_iterator it_p_end = projected.mgsbg_end();

    std::vector<int> col_width(NUM_EDGES, 8);
    std::map<VI,VL> vlm;

    bool inline_view = false;

    BOOST_FOREACH(const EdgeCode& ec, dfsc)
    {
	if (!pl.void_vlabel(ec.vl_from)) vlm[ec.vi_from] = ec.vl_from;
	if (!pl.void_vlabel(ec.vl_to))   vlm[ec.vi_to]   = ec.vl_to;
    }

    switch (pattern_view_mode)
    {
    case PV_LG:
	ostr << "t # " << num_patterns << std::endl;
	for (typename std::map<VI,VL>::const_iterator i = vlm.begin(); i != vlm.end(); ++i)
	    ostr << "v " << i->first << " " << i->second << std::endl;
	BOOST_FOREACH(const EdgeCode& ec, dfsc)
	    ostr << "e " << ec.vi_from << " " << ec.vi_to << " " << ec.el << std::endl;
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
	inline_view = (mapping_view_mode == MV_NONE);

	ostr << "#pattern: " << num_patterns << " : ";
	BOOST_FOREACH(const EdgeCode& ec, dfsc)
	{
	    ostr << "(" << ec.vi_from << "," << ec.vi_to << ", "
		 << vlm[ec.vi_from] << "," << ec.el << "," << vlm[ec.vi_to] << ") ";
	}
	if (!inline_view)
	    ostr << std::endl;
	break;

    default:
	break;
    }

    print_info(projected, inline_view, mapping_view_mode != MV_NONE);

    switch (mapping_view_mode)
    {
    case MV_MAP:
	for (ProjectedManyGraph::MGSBGconst_iterator it_p = projected.mgsbg_begin(); it_p != it_p_end; ++it_p)
	    BOOST_FOREACH(const SBG* sbg, it_p->second)
	    {
		const SBG* s = sbg;
		do
		{
		    ostr << "#tograph: " << tr_names[it_p->first] << std::endl;
		    std::map<VI, VI> vvm;
		    for (int i = 0; i < NUM_EDGES; ++i)
		    {
			const EdgeCode& ec = dfsc[i];
			vvm[ec.vi_from] = (*s)[i].vi_from;
			vvm[ec.vi_to]   = (*s)[i].vi_to;
		    }
		    ostr << "#mv: ";
		    for (std::map<VI, VI>::const_iterator i = vvm.begin(); i != vvm.end(); ++i)
			ostr << i->first << "->" << i->second << " ";
		    ostr << std::endl;
		    ostr << "#me: ";
		    for (int i = 0; i < NUM_EDGES; ++i)
			ostr << i << "->" << (*s)[i].ei << " ";
		    ostr << std::endl;

		} while ( (s = s->automorph_list));
	    }
	break;

    case MV_TABLE:
	for (ProjectedManyGraph::MGSBGconst_iterator it_p = projected.mgsbg_begin(); it_p != it_p_end; ++it_p)
	    BOOST_FOREACH(const SBG* sbg, it_p->second)
	    {
		const SBG* s = sbg;
		do
		{
		    ostr << "#tograph: " << tr_names[it_p->first] << ": ";
		    for (int i = 0; i < NUM_EDGES; ++i)
			ostr << " (" << (*s)[i].vi_from << "," << (*s)[i].vi_to << ")" << (*s)[i].ei;
		    ostr << std::endl;
		} while ( (s = s->automorph_list));
	    }
	break;

    default:
	break;
    }

    ostr << std::endl;
}


void ResultManyGraph::print_info(const ProjectedManyGraph& projected, bool onnewline, bool ismap_view) const
{
    const int supp = projected.support();
    const ProjectedManyGraph::MGSBGconst_iterator it_p_end = projected.mgsbg_end();
    if (!onnewline)
    {
	ostr << "#support: " << supp;
	if (!ismap_view)
	{
	    ostr << std::endl << "#foundin:";
	    for (ProjectedManyGraph::MGSBGconst_iterator it_p = projected.mgsbg_begin();
		 it_p != it_p_end; ++it_p)
		ostr << " " << tr_names[it_p->first];
	}
	ostr << std::endl;
    }
    else
    {
	ostr << "support=" << supp << "; ";
	ostr << "foundin:";
	for (ProjectedManyGraph::MGSBGconst_iterator it_p = projected.mgsbg_begin(); it_p != it_p_end; ++it_p)
	    ostr << " " << tr_names[it_p->first];
    }
}


// *****************************************************************************
//                          main
// *****************************************************************************
std::ostream& usage(std::ostream& ostr)
{ 
    return ostr << "Usage: CMD <minsup> [{-dfsc -edge}] [{ -m -mt }] " << std::endl << std::endl;
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

#ifdef DEBUG_PRINT
		    std::cerr << "INFO:    Graph " << tr_name << " was created"
			      << " at address " << graph << std::endl;
#endif
	    }
#ifdef DEBUG_CHECK_GRAPH_LABEL
	    catch (Policy::VertexNotLabeledException e)
	    {
		++skipped;
		std::cerr << "WARNING: Graph " << tr_name << " not created, vertex "
			  << e.vertex_index << " not labeled" << std::endl;
	    }
#endif
	}
    }


    if (gr_trans.size() == 1)
    {
	ResultOneGraph result(std::cout, pl);
	gSpan::GSPAN_FUNCTION(gr_trans.back(), minsup, pl, result);
    }
    else
    {
	ResultManyGraph result(std::cout, tr_names, pl);
	gSpan::GSPAN_FUNCTION(gr_trans.begin(), gr_trans.end(), minsup, pl, result);
    }
    
}
