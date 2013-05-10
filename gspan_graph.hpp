#ifndef GSPAN_GRAPH_H_
#define GSPAN_GRAPH_H_

#include "gspan_allocator.hpp"

#include <vector>
#include <limits>
#include <iterator>		// std::distance()
#include <iostream>

namespace gSpan2
{

    namespace detail
    {
	typedef unsigned short VI;
	typedef unsigned short EI;
    }

    typedef short VL;
    typedef short EL;

    const detail::VI VI_NULL = std::numeric_limits<detail::VI>::max();
    const detail::EI EI_NULL = std::numeric_limits<detail::EI>::max();
    const VL VL_NULL = -1;
    const EL EL_NULL = -1;

#ifdef TYPE_CHECK
    template<class T, int TypeID>
    class Value
    {
	T x_;
    public:
	Value() :x_() {}
	Value(const T x) :x_(x) {}
	operator T() const { return x_; }

	const T operator= (const T& x) { x_ = x; return x_; }
	const T operator++ () { return ++x_; }
	const T operator++ (int) { return x_++; }
	const T operator-- () { return --x_; }
	const T operator-- (int) { return x_--; }
    };

    typedef Value<detail::VI, 1>	GraphVI;
    typedef Value<detail::EI, 2>	GraphEI;
    typedef Value<detail::VI, 3>	DfscVI;
    typedef Value<detail::EI, 4>	DfscEI;
#else
    typedef detail::VI	GraphVI;
    typedef detail::EI	GraphEI;
    typedef detail::VI	DfscVI;
    typedef detail::EI	DfscEI;
#endif
    
    // *****************************************************************************
    //                          GraphSimple
    // *****************************************************************************
    class GraphSimple
    {
    public:
	template<class EdgeCodeIterator>
	GraphSimple(EdgeCodeIterator it, const EdgeCodeIterator it_end, int num_vertices);
	
	class Edge
	{
	    GraphEI eid_;
	    GraphVI vi_src_;
	    GraphVI vi_dst_;
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

	    Edge(GraphEI eid, GraphVI src, GraphVI dst, VL srclab, VL dstlab, EL elab)
		:eid_(eid),
		 vi_src_(src),
		 vi_dst_(dst),
		 vl_src_(srclab),
		 vl_dst_(dstlab),
		 el_(elab)
		{}

	    GraphVI vi_src() const	{ return vi_src_; }
	    GraphVI vi_dst() const	{ return vi_dst_; }
	    VL vl_src() const	{ return vl_src_; }
	    VL vl_dst() const	{ return vl_dst_; }
	    GraphEI eid() const	{ return eid_; }
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
	const IncidentEdges& incident(GraphVI vi) const	{ return vertices_[vi]; }

	typedef std::vector<Edge> Edges;
	typedef Edges::const_iterator EdgesIterator;
	const Edges& edges() const			{ return edges_; }

	std::size_t num_vertices() const			{ return num_vertices_; }
	std::size_t num_edges() const				{ return num_edges_; }
    private:
	std::vector<IncidentEdges> vertices_;
	std::size_t num_vertices_;
	std::size_t num_edges_;
	Edges edges_;	
    };


    template<class EdgeCodeIterator>
    GraphSimple::GraphSimple(EdgeCodeIterator it, const EdgeCodeIterator it_end, int num_vertices)
	:vertices_(num_vertices),
	 num_vertices_(num_vertices),
	 num_edges_(std::distance(it, it_end)),
	 edges_(num_edges_)	 
    {
	edges_.resize(num_edges_ * 2U);
	GraphEI e_idx1 = 0;
	GraphEI e_idx2 = num_edges_;
	for (; it != it_end; ++it)
	{
	    GraphVI vi_src(it->vi_src());
	    GraphVI vi_dst(it->vi_dst());

	    edges_[e_idx1] = Edge(e_idx1,
				  vi_src, vi_dst,
				  it->vl_src(), it->vl_dst(),
				  it->el());
	    edges_[e_idx2] = Edge(e_idx1,
				  vi_dst, vi_src,
				  it->vl_dst(), it->vl_src(),
				  it->el());
	    
	    vertices_[it->vi_src()].push_back(&edges_[e_idx1]);
	    vertices_[it->vi_dst()].push_back(&edges_[e_idx2]);

	    ++e_idx1;
	    ++e_idx2;
	}
    }


    // *****************************************************************************
    //                          Graph
    // *****************************************************************************
    class Graph : public GraphSimple
    {
    public:
	template<class EdgeCodeIterator>
	Graph(EdgeCodeIterator it, const EdgeCodeIterator it_end, int num_vertices);
    };

    template<class EdgeCodeIterator>
    Graph::Graph(EdgeCodeIterator it, const EdgeCodeIterator it_end, int num_vertices)
	:GraphSimple(it, it_end, num_vertices)
    {
    }



    inline std::ostream& operator<<(std::ostream& out, const Graph::Edge& e)
    {
	return out << e.eid() << "(" << e.vi_src() << "," << e.vi_dst()<<")";
    }
}
#endif
