#ifndef GRAPH_OPS_H_
#define GRAPH_OPS_H_

#include <limits>

#include "gspan.hpp"
#include <boost/graph/compressed_sparse_row_graph.hpp>

namespace graph_alg
{
    template<class T>
    struct VoidValue
    {
	static T void_value() { return std::numeric_limits<T>::max(); }
    };

    template<>
    struct VoidValue<std::string>
    {
	static std::string void_value() { return std::string(); }
    };

	    
    template<class VL, class EL, class DirS>
    struct GraphOpsDefault
    {
	// ----------------------------
	// Label types:
	// ----------------------------
	typedef VL  vertex_label_t;
	typedef EL  edge_label_t;
	typedef const VL& vertex_label_ref_t;
	typedef const EL& edge_label_ref_t;


	// ----------------------------
	// Graph types
	// ----------------------------
	typedef unsigned int vertex_index_t;
	typedef unsigned int edge_index_t;
	typedef boost::compressed_sparse_row_graph<DirS,
						   vertex_label_t,
						   edge_label_t,
						   boost::no_property,
						   vertex_index_t,
						   edge_index_t> graph_t;

	typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_descriptor;
	typedef typename boost::graph_traits<graph_t>::edge_descriptor edge_descriptor;

	// ----------------------------
	// Graph operations
	// ----------------------------
	typedef typename Traits<GraphOpsDefault>::DFSCode DFSCode;
	typedef typename Traits<GraphOpsDefault>::EdgeCode EdgeCode;

	struct EdgePropIterator : public DFSCode::const_iterator
	{
	    explicit EdgePropIterator(typename DFSCode::const_iterator iter)
		:DFSCode::const_iterator(iter) {}
	    const edge_label_ref_t operator* () const { return (*this)->el; };
	};

	void dfsc_add(DFSCode& dfsc, vertex_index_t from, vertex_index_t to,
		      VL fromlabel, EL elabel, VL tolabel)
	    {
		typedef typename Traits<GraphOpsDefault>::EdgeCode EdgeCode;
		dfsc.push(EdgeCode(from, to,
				   dfsc.empty() ? fromlabel : void_vlabel(),
				   elabel,
				   from < to ? tolabel : void_vlabel()));
	    }

	struct VertexNotLabeledException
	{
	    vertex_index_t vertex_index;
	    explicit VertexNotLabeledException(vertex_index_t vi) :vertex_index(vi) {}
	};


	graph_t* create_graph(const DFSCode& dfsc) const;

	vertex_index_t void_vertex_index() const
	    { return std::numeric_limits<vertex_index_t>::max(); }

	edge_index_t   void_edge_index() const
	    { return std::numeric_limits<edge_index_t>::max(); }

	vertex_index_t    get_vertex_index(vertex_descriptor vd, const graph_t& g) const
	    { return vd; }

	vertex_descriptor get_vertex_descriptor(vertex_index_t vi, const graph_t& g) const
	    { return vi; }

	edge_index_t      get_edge_index(edge_descriptor ed, const graph_t& g) const
	    { return get(boost::edge_index_t(), g, ed); }
	edge_descriptor   get_edge_descriptor(edge_index_t ei, const graph_t& g) const
	    { return edge_from_index(ei,g); }

	// ----------------------------
	// Label operations:
	// ----------------------------
	vertex_label_ref_t vlabel (vertex_descriptor vd, const graph_t& g) const { return g[vd]; }
	vertex_label_ref_t vilabel(vertex_index_t vi, const graph_t& g) const    { return g[vi]; }
	edge_label_ref_t   elabel (edge_descriptor ed, const graph_t& g) const   { return g[ed]; }
	edge_label_ref_t   eilabel(edge_index_t ei, const graph_t& g) const
	    { return elabel(get_edge_descriptor(ei,g), g); }

	bool vlabel_equal         (vertex_label_ref_t vl1, vertex_label_ref_t vl2) const
	    { return vl1 == vl2; }

	bool vlabel_less          (vertex_label_ref_t vl1, vertex_label_ref_t vl2) const
	    { return vl1 < vl2; }

	bool vlabel_less_or_equal (vertex_label_ref_t vl1, vertex_label_ref_t vl2) const
	    { return vl1 <= vl2; }

	bool elabel_equal         (edge_label_ref_t el1, edge_label_ref_t el2) const
	    { return el1 == el2; }

	bool elabel_less          (edge_label_ref_t el1, edge_label_ref_t el2) const
	    { return el1 < el2; }

	bool elabel_less_or_equal (edge_label_ref_t el1, edge_label_ref_t el2) const
	    { return el1 <= el2; }

	vertex_label_ref_t void_vlabel() const
	    {
		static VoidValue<vertex_label_t> nv;
		static vertex_label_t vl = nv.void_value();
		return vl;
	    }

	edge_label_ref_t void_elabel() const
	    {
		static VoidValue<edge_label_t> nv;
		static edge_label_t el = nv.void_value();
		return el;
	    }

	bool void_vlabel(vertex_label_ref_t vl) const { return vlabel_equal(vl, void_vlabel()); }
	bool void_elabel(edge_label_ref_t el) const   { return elabel_equal(el, void_elabel()); }
    };

    template<class VL, class EL, class DirS>
    typename GraphOpsDefault<VL,EL,DirS>::graph_t*
    GraphOpsDefault<VL,EL,DirS>::create_graph(const DFSCode& dfsc) const
    {
	graph_t* ptr = new graph_t(boost::edges_are_unsorted_multi_pass_t(),
				   dfsc.begin(), dfsc.end(),
				   EdgePropIterator(dfsc.begin()),
				   dfsc.max_vertex() + 1);
	graph_t& g = *ptr;
	
	// set labels
	typename DFSCode::const_iterator it = dfsc.begin();
	typename DFSCode::const_iterator it_end = dfsc.end();
	while (it != it_end)
	{
	    if (! void_vlabel(it->vl_from))
		g[it->vi_from] = it->vl_from;
	    
	    if (! void_vlabel(it->vl_to))
		g[it->vi_to] = it->vl_to;
	    ++it;
	}

	typename boost::graph_traits<graph_t>::vertex_iterator vi, viend;
	for (boost::tie(vi,viend) = vertices(g); vi != viend; ++vi)
	{
	    if (void_vlabel(g[*vi]))
	    {
		throw VertexNotLabeledException(*vi);
	    }
	}

	return ptr;
    }
}

#endif

