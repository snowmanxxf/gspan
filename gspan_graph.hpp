#ifndef GSPAN_GRAPH_H_
#define GSPAN_GRAPH_H_

#include <vector>
#include <limits>
#include <iterator>		// std::distance()
#include <iostream>

namespace gSpan2
{
    typedef unsigned int VI;
    typedef unsigned int EI;
    typedef int VL;
    typedef int EL;

    const VI VI_NULL = std::numeric_limits<VI>::max();
    const EI EI_NULL = std::numeric_limits<EI>::max();
    const VL VL_NULL = -1;
    const EL EL_NULL = -1;

    class Graph
    {
    public:
	template<class EdgeCodeIterator>
	Graph(EdgeCodeIterator it, const EdgeCodeIterator it_end, int num_vertices);
	
	class Edge
	{
	    EI eid_;
	    VI vi_src_;
	    VI vi_dst_;
	    VL vl_src_;
	    VL vl_dst_;
	    EL el_;
	public:
	    Edge()
		:eid_(EI_NULL),
		 vi_src_(VI_NULL),
		 vi_dst_(VI_NULL),
		 vl_src_(VL_NULL),
		 vl_dst_(VL_NULL),
		 el_(EL_NULL)
		{}

	    Edge(EI eid, VI src, VI dst, VL srclab, VL dstlab, EL elab)
		:eid_(eid),
		 vi_src_(src),
		 vi_dst_(dst),
		 vl_src_(srclab),
		 vl_dst_(dstlab),
		 el_(elab)
		{}

	    VI vi_src() const	{ return vi_src_; }
	    VI vi_dst() const	{ return vi_dst_; }
	    VL vl_src() const	{ return vl_src_; }
	    VL vl_dst() const	{ return vl_dst_; }
	    EI eid() const	{ return eid_; }
	    EL el() const	{ return el_; }
	    
	    void chgdir()
		{
		    std::swap(vi_src_, vi_dst_);
		    std::swap(vl_src_, vl_dst_);
		}

	    Edge operator- () const { Edge e(*this); e.chgdir(); return e; }
	};

	typedef std::vector<const Edge*> IncidentEdges;
	typedef IncidentEdges::const_iterator IncidentEdgesIterator;
	const IncidentEdges& incident(VI vi) const	{ return vertices_[vi]; }

	typedef std::vector<Edge> Edges;
	typedef Edges::const_iterator EdgesIterator;
	const Edges& edges() const			{ return edges_; }

	int num_vertices() const			{ return vertices_.size(); }
	int num_edges() const				{ return num_edges_; }
    private:
	std::vector<IncidentEdges> vertices_;
	unsigned int num_edges_;
	Edges edges_;
    };


    template<class EdgeCodeIterator>
    Graph::Graph(EdgeCodeIterator it, const EdgeCodeIterator it_end, int num_vertices)
	:vertices_(num_vertices),
	 num_edges_(std::distance(it, it_end)),
	 edges_(num_edges_)
    {
	edges_.resize(num_edges_ * 2U);
	EI e_idx1 = 0;
	EI e_idx2 = num_edges_;
	for (; it != it_end; ++it)
	{
	    edges_[e_idx1] = Edge(e_idx1,
				  it->vi_src(), it->vi_dst(),
				  it->vl_src(), it->vl_dst(),
				  it->el());
	    edges_[e_idx2] = Edge(e_idx1,
				  it->vi_dst(), it->vi_src(),
				  it->vl_dst(), it->vl_src(),
				  it->el());
	    
	    vertices_[it->vi_src()].push_back(&edges_[e_idx1]);
	    vertices_[it->vi_dst()].push_back(&edges_[e_idx2]);

	    ++e_idx1;
	    ++e_idx2;
	}
    }

    inline std::ostream& operator<<(std::ostream& out, const Graph::Edge& e)
    {
	return out << e.eid() << "(" << e.vi_src() << "," << e.vi_dst()<<")";
    }
}
#endif
