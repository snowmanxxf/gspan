#ifndef MISC_H_
#define MISC_H_

#include "gspan2.hpp"
#include <vector>
#include <map>
#include <boost/foreach.hpp>

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

void read_input(std::list<InputGraph>& igl, std::istream& is);

void relabel(std::vector<std::string>& i_to_str,
	     std::map<std::string, int>& str_to_i,
	     const std::map<std::string, int>& lab_counts,
	     bool sort_label_asc);

void set_labels_from_file(std::vector<std::string>& i_to_str,
			  std::map<std::string, int>& str_to_i,
			  const char* fname);

bool check_lebels(const std::map<std::string, int>& str_to_i, const std::map<std::string, int>& lab_counts);

struct WorkingGraphs
{
    std::vector<const gSpan2::Graph*> graphs;
    std::map<const gSpan2::Graph*, std::string> names;
};

void create_working_graphs(WorkingGraphs& wrk_graphs,
			   const std::map<std::string, int>& vls_vl,
			   const std::map<std::string, int>& els_el,
			   const std::list<InputGraph>& input_graphs);


#endif
