
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

enum PatternViewMode { PV_LG, PV_DFSC, PV_EDGE} pattern_view_mode;
enum MappingViewMode { MV_NONE, MV_TABLE, MV_MAP } mapping_view_mode;

class Result
{
public:
    Result(std::ostream& ostr, std::map<const Graph*, std::string>& tr_names, const Policy& pl)
	: ostr(ostr), tr_names(tr_names), pl(pl), num_patterns(0) {}
    void operator() (const DFSCode& dfsc, const Projected& projected);
private:
    std::ostream& ostr;
    std::map<const Graph*, std::string>& tr_names;
    const Policy& pl;
    int num_patterns;

    void print_info(const Projected& projected, bool onnewline, bool ismap_view) const;
};

void Result::print_info(const Projected& projected, bool onnewline, bool ismap_view) const
{
    const int supp = projected.mgsbg_size();
    const Projected::MGSBGconst_iterator it_p_end = projected.mgsbg_end();
    if (!onnewline)
    {
	ostr << "#support: " << supp;
	if (!ismap_view)
	{
	    ostr << std::endl << "#foundin:";
	    for (Projected::MGSBGconst_iterator it_p = projected.mgsbg_begin();
		 it_p != it_p_end; ++it_p)
		ostr << " " << tr_names[it_p->first];
	}
	ostr << std::endl;
    }
    else
    {
	ostr << "support=" << supp << "; ";
	ostr << "foundin:";
	for (Projected::MGSBGconst_iterator it_p = projected.mgsbg_begin(); it_p != it_p_end; ++it_p)
	    ostr << " " << tr_names[it_p->first];
    }
}


void Result::operator() (const DFSCode& dfsc, const Projected& projected)
{
    ++num_patterns;
    
    const int NUM_EDGES = dfsc.size();
    const Projected::MGSBGconst_iterator it_p_end = projected.mgsbg_end();

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
	for (Projected::MGSBGconst_iterator it_p = projected.mgsbg_begin(); it_p != it_p_end; ++it_p)
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
	for (Projected::MGSBGconst_iterator it_p = projected.mgsbg_begin(); it_p != it_p_end; ++it_p)
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

/*
void Result::print_tgf(const DFSCode& dfsc, const Projected& projected)
{
    std::map<VI,VL> vlabels;
    BOOST_FOREACH(const EdgeCode& ec, dfsc)
    {
	if (!pl.void_vlabel(ec.vl_from))vlabels[ec.vi_from] = ec.vl_from;
	if (!pl.void_vlabel(ec.vl_to))  vlabels[ec.vi_to] = ec.vl_to;
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
    ostr << std::endl;
}

*/

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
	
    Result result(std::cout, tr_names, pl);
    gSpan::GSPAN_FUNCTION(gr_trans.begin(), gr_trans.end(), minsup, pl, result);
}
