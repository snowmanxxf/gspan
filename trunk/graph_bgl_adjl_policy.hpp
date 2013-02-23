#ifndef GRAPH_BGL_ADJL_POLICY_H_
#define GRAPH_BGL_ADJL_POLICY_H_
// it is NOT a part 'gspan' library

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <limits>
#include <iostream>

#ifndef BR
#define BR asm volatile ("int3;")
#endif

namespace bgl_adaptor
{
    namespace detail
    {
	typedef int	vertex_label_t;
	typedef int	edge_label_t;
	typedef unsigned int    vertex_index_t;
	typedef unsigned int	edge_index_t;

	typedef vertex_label_t	vertex_label_ref_t;
	typedef edge_label_t	edge_label_ref_t;

	const vertex_index_t INVALID_VINDEX = std::numeric_limits<vertex_index_t>::max();
	const edge_index_t   INVALID_EINDEX = std::numeric_limits<edge_index_t>::max();
	const edge_label_t LABEL_NULL = -1;

	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
				      vertex_label_t,
				      boost::property<boost::edge_index_t, std::size_t> > G_;

	typedef boost::graph_traits<G_>::vertex_descriptor vertex_descriptor;
	typedef boost::graph_traits<G_>::edge_descriptor edge_descriptor;
	typedef boost::graph_traits<G_>::vertex_iterator vertex_iterator;
	typedef boost::graph_traits<G_>::edge_iterator edge_iterator;
	typedef boost::graph_traits<G_>::out_edge_iterator out_edge_iterator;
	typedef boost::graph_traits<G_>::in_edge_iterator in_edge_iterator;

	class GraphBGL : public G_
	{
	public:
	    template<class DFSCode>
	    explicit GraphBGL(const DFSCode&);

	    vertex_index_t source_index(edge_descriptor ed) const { return source(ed, *this); }
	    vertex_index_t target_index(edge_descriptor ed) const { return target(ed, *this); }
	    edge_index_t   edge_index(edge_descriptor ed) const
		{ return get(boost::edge_index_t(), static_cast<const G_&>(*this), ed); }

	    edge_descriptor edge_descr(edge_index_t ei) const	{ return eds_[ei]; }
	    edge_label_ref_t elabel(edge_index_t ei) const	{ return elabs_[ei]; }
	private:
	    std::vector<edge_label_t> elabs_;
	    std::vector<edge_descriptor> eds_;
	};
	
	template<class DFSCode>
	GraphBGL::GraphBGL(const DFSCode& dfsc)
	    :G_(max_vertex(dfsc) + 1),
	     elabs_(dfsc.size()),
	     eds_(dfsc.size())
	{
	    edge_index_t ei = 0;
	    typename DFSCode::const_iterator it = dfsc.begin();
	    typename DFSCode::const_iterator it_end = dfsc.end();
	    while (it != it_end)
	    {
		std::pair<edge_descriptor,bool> ep = add_edge(it->vi_from, it->vi_to, ei, *this);
		assert(ep.second);

		elabs_[ei] = it->el;
		eds_[ei] = ep.first;

		if (it->vl_from != LABEL_NULL)
		    (*this)[it->vi_from] = it->vl_from;
	    
		if (it->vl_to != LABEL_NULL)
		    (*this)[it->vi_to] = it->vl_to;
		++it;
		++ei;
	    }
	}
    } // end: namespace detail

    using detail::INVALID_VINDEX;
    using detail::INVALID_EINDEX;
    using detail::LABEL_NULL;

    // it is used with 'gspan' library
    // all methods should be static
    struct GraphBGLPolicy
    {
	typedef detail::GraphBGL	graph_t;
	typedef detail::vertex_index_t  vertex_index_t;
	typedef detail::edge_index_t	edge_index_t;	    
	typedef detail::vertex_label_t	vertex_label_t;
	typedef detail::edge_label_t	edge_label_t;
	typedef detail::vertex_label_ref_t	vertex_label_ref_t;
	typedef detail::edge_label_ref_t	edge_label_ref_t;

	template<class DFSCode>
	static graph_t* create_graph(const DFSCode& dfsc) { return new graph_t(dfsc); }

	static vertex_label_ref_t vlabel(vertex_index_t v, const graph_t& g) { return g[v]; }
	static vertex_label_ref_t vlabel_null() { return LABEL_NULL; }
	static edge_label_ref_t elabel(edge_index_t ei, const graph_t& g) { return g.elabel(ei); }
	static edge_label_ref_t elabel_null() { return LABEL_NULL; }
	static int nvertices(const graph_t& g) { return num_vertices(g); }
	static int nedges(const graph_t& g) { return num_edges(g); }
	static vertex_index_t void_vindex() { return INVALID_VINDEX; }

	struct Edge
	{
	    vertex_index_t vi_from, vi_to;
	    edge_index_t ei;
	    Edge(vertex_index_t from = INVALID_VINDEX,
		 vertex_index_t to = INVALID_VINDEX,
		 edge_index_t e = INVALID_EINDEX)
		:vi_from(from), vi_to(to), ei(e) {}
	    void chgdir() { std::swap(vi_from, vi_to); }
	};
	
	friend std::ostream& operator<<(std::ostream& out, const Edge& e)
	    {
		return out<<e.ei<<"("<<e.vi_from<<","<<e.vi_to<<") ";
	    }

    private:
	
	template<class It>
	class EdgeIt_
	{
	    std::pair<It,It> ep_;
	public:
	    EdgeIt_(std::pair<It,It> ep) :ep_(ep) {}
	    bool valid() const { return ep_.first != ep_.second; }
	    Edge edge(const graph_t& g) const
		{
		    return Edge(g.source_index(*ep_.first),
				g.target_index(*ep_.first),
				g.edge_index(*ep_.first));
		}

	    Edge edge_rev(const graph_t& g) const
		{
		    return Edge(g.target_index(*ep_.first),
				g.source_index(*ep_.first),
				g.edge_index(*ep_.first));
		}
	    void increment() { ++ep_.first; }
	};
    public:
	
	class AllVertexIt
	{
	    std::pair<detail::vertex_iterator,detail::vertex_iterator> ip_;
	    const graph_t& g_;
	public:
	    AllVertexIt(const graph_t& g) :ip_(vertices(g)), g_(g) {}
	    bool valid() const { return ip_.first != ip_.second; }
	    vertex_index_t vertex() const { return *ip_.first; }
	    void increment() { ++ip_.first; }
	};

	class AllEdgeIt : public EdgeIt_<detail::edge_iterator>
	{
	    const graph_t& g_;
	public:
	    AllEdgeIt(const graph_t& g) :EdgeIt_<detail::edge_iterator>(edges(g)), g_(g) {}
	    Edge edge() const { return EdgeIt_<detail::edge_iterator>::edge(g_); }
	};

	class IncidEdgeIt
	{
	    EdgeIt_<detail::out_edge_iterator> ou_it_;
	    const graph_t& g_;
	    Edge e_;
	public:
	    IncidEdgeIt(vertex_index_t v, const graph_t& g)
		:ou_it_(out_edges(v,g)),
		 g_(g)
		{
		    if (ou_it_.valid())
			e_ = ou_it_.edge(g_);
		}
	    bool valid() const { return ou_it_.valid(); }
	    const Edge& edge() const { return e_; }
	    void increment()
		{
		    ou_it_.increment();
		    if (ou_it_.valid())
			e_ = ou_it_.edge(g_);		    
		}
	}; // end: class IncidEdgeIt
    }; // end: class GraphBGLPolicy
}
#endif
