
#include "misc.hpp"

#include <fstream>
#include <numeric>
#include <cassert>


void read_input(std::list<InputGraph>& igl, std::istream& is)
{
    InputGraph* g = 0;

    char line[1024];
    while (true)
    {
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
	    int n = atoi(result[1].c_str());
	    g->vl[n] = result[2];
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
