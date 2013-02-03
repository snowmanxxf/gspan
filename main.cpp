
#include "graph_bgl_policy.hpp"
#include "gspan.hpp"

#include <iostream>
#include <iomanip>
#include <cstring>
#include <boost/ptr_container/ptr_vector.hpp>

double total_call;
double cache_hit;

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

	    EL elabel = LABEL_NULL;
	    if (result.size() == 4)
		elabel = result[3];

	    EdgeCode ec(from, to, vlabels[from], elabel, vlabels[to]);
            dfsc.push_back(ec);
        }
    }
    return is;
}

// *****************************************************************************
//                          Result
// *****************************************************************************

enum PatternViewMode { PV_LG, PV_DFSC, PV_EDGE} pattern_view_mode;
enum MappingViewMode { MV_NONE, MV_TABLE, MV_MAP } mapping_view_mode;

bool inline_view;
std::map<const Policy::graph_t*, std::string> tr_names;

class Result
{
public:
    Result(std::ostream& ostr)
	: ostr(ostr), num_patterns(0) {}

    void operator() (const DFSCode& dfsc, const SMG& smg);
    void operator() (const DFSCode& dfsc, const SOG& sog);
private:
    std::ostream& ostr;
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
	    ostr << "#tograph: " << tr_names[i->first] << std::endl;
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
    std::map<VI,VL> vlm;
    BOOST_FOREACH(const EdgeCode& ec, dfsc)
    {
	if (ec.vl_from != LABEL_NULL) vlm[ec.vi_from] = ec.vl_from;
	if (ec.vl_to != LABEL_NULL)   vlm[ec.vi_to]   = ec.vl_to;
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
		ostr << " " << tr_names[it->first];
	}
    }
    else
    {
	ostr << "support=" << supp << "; ";
	ostr << "foundin:";
	for (SMG::const_iterator it = smg.begin(); it != smg.end(); ++it)
	    ostr << " " << tr_names[it->first];
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
    boost::ptr_vector<Policy::graph_t> gr_trans;

    {
	while (true)
	{
	    std::string tr_name;
	    DFSCode dfsc;
	    contruct_dfsc(dfsc, tr_name, std::cin);
	    if (dfsc.empty())
		break;
	    {
		Policy::graph_t* graph = Policy::create_graph(dfsc);
		gr_trans.push_back(graph);
		tr_names[graph] = tr_name;
		
#ifdef DEBUG_PRINT
		std::cerr << "INFO:    Graph " << tr_name << " was created"
			  << " at address " << graph << std::endl;
#endif
	    }
	}
    }

    inline_view = pattern_view_mode == PV_EDGE && mapping_view_mode == MV_NONE;

    Result result(std::cout);

    if (gr_trans.size() == 1)
	gSpan::GSPAN_FUNCTION<Policy>(gr_trans.back(), minsup, result);
    else
	gSpan::GSPAN_FUNCTION<Policy>(gr_trans.begin(), gr_trans.end(), minsup, result);

    std::cerr << "total_call=" << total_call << " cache_hit" << cache_hit
	      << " cache_hit%% = " << (cache_hit/total_call)*100.0 << std::endl;
}
