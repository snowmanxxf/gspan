
#include "misc.hpp"

#include <fstream>
#include <numeric>
#include <cassert>

void set_labels_from_file(std::vector<std::string>& i_to_str,
			  std::map<std::string, int>& str_to_i,
			  const char* fname)
{
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

template<class M>
void prt_map_key(const M& m)
{
    for (typename M::const_iterator j = m.begin(); j != m.end(); ++j)
	std::cerr << j->first << " ";
    std::cerr << std::endl;
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

void relabel(std::vector<std::string>& i_to_str,
	     std::map<std::string, int>& str_to_i,
	     const std::map<std::string, int>& lab_counts,
	     bool sort_label_asc)
{
    std::multimap<int, std::string> frq;
    for (typename std::map<std::string, int>::const_iterator i = lab_counts.begin(); i != lab_counts.end(); ++i)
	frq.insert(std::pair<int, std::string>(i->second, i->first));
    
    int count = 0;
    if (sort_label_asc)
	for (typename std::multimap<int, std::string>::const_reverse_iterator i = frq.rbegin(); i != frq.rend(); ++i)
	    str_to_i[i->second] = count++;
    else
	for (typename std::multimap<int, std::string>::const_iterator i = frq.begin(); i != frq.end(); ++i)
	    str_to_i[i->second] = count++;
    
    i_to_str.resize(str_to_i.size());
    for (std::map<std::string, int>::const_iterator i = str_to_i.begin(); i != str_to_i.end(); ++i)
	i_to_str[i->second] = i->first;
}


void create_working_graphs(WorkingGraphs& wrk_graphs,
			   const std::map<std::string, int>& vls_vl,
			   const std::map<std::string, int>& els_el,
			   const std::list<InputGraph>& input_graphs)
{
    for (std::list<InputGraph>::const_iterator iter = input_graphs.begin(); iter != input_graphs.end(); ++iter)
    {
        const InputGraph& ig = *iter;

	int corr = ig.vl.begin()->first == 1;
	
	gSpan::Graph* graph = new gSpan::Graph(ig.vl.size(), ig.edges.size());

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
	    
	    bool fwd = (i->from - corr) < (i->to - corr);
	    gSpan::EdgeCode ec(i->from - corr, i->to - corr, vl_from_i, el_i, vl_to_i, fwd);
	    
            graph->push_edge(ec);
	}
        
        
	wrk_graphs.names[graph] = ig.name;
	wrk_graphs.graphs.push_back(graph);
    }
}

