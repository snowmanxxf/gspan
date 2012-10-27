#ifndef GRAPH_POLICY_H_
#define GRAPH_POLICY_H_

#include <string>
#include <boost/graph/compressed_sparse_row_graph.hpp>

struct GraphPolicy
{
    typedef std::string		vertex_label_t;
    typedef std::string		edge_label_t;
    typedef unsigned int        vertex_index_t;
    typedef unsigned int	edge_index_t;
    typedef boost::compressed_sparse_row_graph<boost::bidirectionalS,
					       vertex_label_t,
					       edge_label_t,
					       boost::no_property,
					       vertex_index_t,
					       edge_index_t> graph_t;

    typedef const vertex_label_t& vertex_label_ref_t;
    typedef const edge_label_t&   edge_label_ref_t;

    typedef boost::graph_traits<graph_t>::vertex_descriptor vertex_descriptor;
    typedef boost::graph_traits<graph_t>::edge_descriptor edge_descriptor;

    // creation graph from dfs code
    template<class DFSCode>
    struct EdgePropIterator : public DFSCode::const_iterator
    {
	explicit EdgePropIterator(typename DFSCode::const_iterator iter)
	    :DFSCode::const_iterator(iter) {}
	const edge_label_ref_t operator* () const { return (*this)->el; };
    };

#ifdef DEBUG_CHECK_GRAPH_LABEL
    struct VertexNotLabeledException
    {
	vertex_index_t vertex_index;
	VertexNotLabeledException(vertex_index_t vi) : vertex_index(vi) {}
    };
#endif

    template<class DFSCode>
    graph_t* create_graph(const DFSCode&) const;


    // vertex_descriptor <--> vertex_index_t
    vertex_descriptor get_vertex_descriptor(vertex_index_t vi, const graph_t& g) const { return vi; }
    vertex_index_t get_vertex_index(vertex_descriptor vd, const graph_t& g) const { return vd; }

    // edge_descriptor <--> edge_index_t
    edge_descriptor get_edge_descriptor(edge_index_t ei, const graph_t& g) const { return edge_from_index(ei,g); }
    edge_index_t get_edge_index(edge_descriptor ed, const graph_t& g) const { return get(boost::edge_index_t(), g, ed);}

    vertex_index_t void_vertex_index() const { return std::numeric_limits<vertex_index_t>::max(); }
    edge_index_t void_edge_index() const { return std::numeric_limits<edge_index_t>::max(); }

    // label access
    vertex_label_ref_t vdlabel(vertex_descriptor vd, const graph_t& g) const { return g[vd]; }
    vertex_label_ref_t vilabel(vertex_index_t vi, const graph_t& g) const
	{ return vdlabel(get_vertex_descriptor(vi,g),g); }

    edge_label_ref_t edlabel(edge_descriptor ed, const graph_t& g) const { return g[ed]; }
    edge_label_ref_t eilabel(edge_index_t ei, const graph_t& g) const
	{ return edlabel(get_edge_descriptor(ei,g),g); }
    
    vertex_label_ref_t void_vlabel() const { static std::string s; return s; }
    edge_label_ref_t   void_elabel() const { static std::string s; return s; }

    // label compare
    bool vlabel_equal(vertex_label_ref_t vl1, vertex_label_ref_t vl2) const { return vl1 == vl2; }
    bool vlabel_less(vertex_label_ref_t vl1, vertex_label_ref_t vl2) const { return vl1 < vl2; }
    bool vlabel_less_or_equal(vertex_label_ref_t vl1, vertex_label_ref_t vl2) const
	{ return vl1 <= vl2; }

    bool elabel_equal(edge_label_ref_t el1, edge_label_ref_t el2) const { return el1 == el2; }
    bool elabel_less(edge_label_ref_t el1, edge_label_ref_t el2) const { return el1 < el2; }

    bool void_vlabel(vertex_label_ref_t vl) const { return vlabel_equal(vl, void_vlabel()); };
    bool void_elabel(edge_label_ref_t el) const { return elabel_equal(el, void_elabel()); };
};


// creation graph from dfs code
template<class DFSCode>
GraphPolicy::graph_t* GraphPolicy::create_graph(const DFSCode& dfsc) const
{
    graph_t* ptr = new graph_t(boost::edges_are_unsorted_multi_pass_t(),
			       dfsc.begin(), dfsc.end(),
			       EdgePropIterator<DFSCode>(dfsc.begin()),
			       max_vertex(dfsc) + 1);
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

#ifdef DEBUG_CHECK_GRAPH_LABEL
    typename boost::graph_traits<graph_t>::vertex_iterator vi, viend;
    for (boost::tie(vi,viend) = vertices(g); vi != viend; ++vi)
    {
	if (void_vlabel(g[*vi]))
	    throw VertexNotLabeledException(*vi);
    }
#endif
    return ptr;
}

#endif
