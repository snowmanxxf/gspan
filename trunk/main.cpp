
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <boost/graph/graphviz.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "gspan.hpp"
#include "graph_ops.hpp"

using namespace graph_alg;

typedef GraphOpsDefault<char, char, boost::bidirectionalS> GraphOps;
typedef GraphOps::DFSCode DFSCode;
GraphOps ops;

int result;

int main(int argc, char** argv)
{
    DFSCode dfsc;
    ops.dfsc_add(dfsc, 0, 1, 'X', 'a', 'X');
    ops.dfsc_add(dfsc, 1, 2, 'X', 'a', 'Y');
    ops.dfsc_add(dfsc, 2, 0, 'Y', 'b', 'X');
    ops.dfsc_add(dfsc, 2, 3, 'Y', 'b', 'Z');
    ops.dfsc_add(dfsc, 3, 0, 'Z', 'c', 'X');
    ops.dfsc_add(dfsc, 2, 4, 'Y', 'd', 'Z');

    std::cout << "is_min: " << is_min(dfsc, ops) << std::endl;

    //boost::ptr_vector<GraphOps::graph_t> gr;    
    //gspan(gr.begin(), gr.end(), 2, ops, result);
}
