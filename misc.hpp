#ifndef MISC_H_
#define MISC_H_

#include "read_input.hpp"

#include "gspan.hpp"
#include <vector>
#include <map>
#include <list>
#include <cstring>


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
    std::vector<const gSpan::Graph*> graphs;
    std::map<const gSpan::Graph*, std::string> names;

    ~WorkingGraphs()
        {
            for (unsigned int i = 0; i < graphs.size(); ++i)
                delete graphs[i];
        }
};

void create_working_graphs(WorkingGraphs& wrk_graphs,
			   const std::map<std::string, int>& vls_vl,
			   const std::map<std::string, int>& els_el,
			   const std::list<InputGraph>& input_graphs);


#endif
