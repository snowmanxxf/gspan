

#include "gspan.hpp"
#include <iostream>
#include <boost/graph/graphviz.hpp>

#include <sstream>
#include <fstream>

#include <stdlib.h>

using namespace graph_alg;
//typedef WGraph<boost::directedS> WG;
typedef WGraph<boost::bidirectionalS> WG;

typedef typename boost::graph_traits<WG>::edge_descriptor   edge_descriptor;
typedef typename boost::graph_traits<WG>::vertex_descriptor vertex_descriptor;

std::map<int, char> mlv;
std::map<int, char> mle;

struct DSFCodeGraphViz
{
    class my_label_writer {
    public:
	my_label_writer(const WG& g) : g_(g) {}
	void operator()(std::ostream& out, edge_descriptor e) const {
	    out << "[label=\"" << mle[g_[e]] << "\"]";
	}

	void operator()(std::ostream& out, const vertex_descriptor v) const {
	    out << "[label=\"" << v << " "<< mlv[g_[v]] << "\"]";
	}
    private:
	const WG& g_;
    };

    DSFCodeGraphViz() { ::system("rm -rf grviz/*"); }

    void operator() (const Projected<WG>& projected, const DFSCode& dfsc)
	{
	    static int N;
	    WG g(dfsc);

	    std::ostringstream fname1;
	    std::ostringstream fname2;
	    std::filebuf f;

	    fname1 << "grviz/subgraph_" << N << ".dot";
	    fname2 << "grviz/subgraph_" << N << ".png";
	    ++N;

	    f.open(fname1.str().c_str(), std::ios::out);
	    std::ostream os(&f);
	    write_graphviz(os, g, my_label_writer(g), my_label_writer(g));
	    f.close();
	    
	    ::system((std::string("dot -Tpng ") + fname1.str() + " > " + fname2.str()).c_str());
	}
};



template<class VLM, class ELM, class VL, class EL>
void dfsc_add(DFSCode& dfsc, Vertex from, Vertex to,
	      VL fromlabel, EL elabel, VL tolabel, VLM& vlm, ELM& elm)
{
    dfsc.push(EdgeCode(from, to,
		       dfsc.empty() ? vlm[fromlabel] : -1,
		       elm[elabel],
		       from<to ? vlm[tolabel] : -1));
}

int main(int argc, char** argv)
{
    std::map<char, int> vlm;
    vlm['X'] = 0; mlv[0]='X';
    vlm['Y'] = 1; mlv[1]='Y';
    vlm['Z'] = 2; mlv[2]='Z';

    std::map<char, int> elm;
    elm['a'] = 0; mle[0]='a';
    elm['b'] = 1; mle[1]='b';
    elm['c'] = 2; mle[2]='c';
    elm['d'] = 3; mle[3]='d';

    //----------------------------------------------------------

    DFSCode dfsc;
    dfsc_add(dfsc, 0, 1, 'X', 'a', 'X', vlm, elm);
    dfsc_add(dfsc, 1, 2, 'X', 'a', 'Y', vlm, elm);
    dfsc_add(dfsc, 2, 0, 'Y', 'b', 'X', vlm, elm);
    dfsc_add(dfsc, 2, 3, 'Y', 'b', 'Z', vlm, elm);
    dfsc_add(dfsc, 3, 0, 'Z', 'c', 'X', vlm, elm);
    dfsc_add(dfsc, 2, 4, 'Y', 'd', 'Z', vlm, elm);
    std::cout << "is_min: " << is_min<WG>(dfsc) << std::endl;

    //----------------------------------------------------------
    DSFCodeGraphViz grviz_result;

    DFSCode dfsc_1;
    dfsc_add(dfsc_1, 0, 1, 'X', 'a', 'Y', vlm, elm);
    dfsc_add(dfsc_1, 1, 2, 'Y', 'b', 'X', vlm, elm);
    dfsc_add(dfsc_1, 2, 0, 'X', 'a', 'X', vlm, elm);
    dfsc_add(dfsc_1, 2, 3, 'X', 'c', 'Z', vlm, elm);
    dfsc_add(dfsc_1, 3, 1, 'Z', 'b', 'Y', vlm, elm);
    dfsc_add(dfsc_1, 1, 4, 'Y', 'd', 'Z', vlm, elm);

    DFSCode dfsc_2;
    dfsc_add(dfsc_2, 0, 1, 'Y', 'a', 'X', vlm, elm);
    dfsc_add(dfsc_2, 1, 2, 'X', 'a', 'X', vlm, elm);
    dfsc_add(dfsc_2, 2, 0, 'X', 'b', 'Y', vlm, elm);
    dfsc_add(dfsc_2, 2, 3, 'X', 'c', 'Z', vlm, elm);
    dfsc_add(dfsc_2, 3, 0, 'Z', 'b', 'Y', vlm, elm);
    dfsc_add(dfsc_2, 0, 4, 'Y', 'd', 'Z', vlm, elm);

    WG tr[2];
    tr[0] = WG(dfsc_1);
    tr[1] = WG(dfsc_2);
    gspan(tr, tr+2, 2, grviz_result);
}
