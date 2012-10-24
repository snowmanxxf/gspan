#ifndef GSPAN_H_
#define GSPAN_H_

#include <cassert>
#include <vector>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <memory>
#include <iostream>

#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>

#define BR asm("int $3;")

namespace graph_alg
{
    // =============================================================================
    // template parameter GraphOps
    // -----------------------------------------------------------------------------
    //
    // ----------------------------
    // Graph types
    // ----------------------------
    //  graph_t
    //  vertex_index_t
    //  edge_index_t
    //
    // ----------------------------
    // Graph operations
    // ----------------------------
    //  graph_t* create_graph(const DFSCode& dfsc);
    //  vertex_index_t void_vertex_index();
    //  edge_index_t   void_edge_index();
    //  vertex_index_t    get_vertex_index(vertex_descriptor vd, const graph_t& g);
    //  vertex_descriptor get_vertex_descriptor(vertex_index vi, const graph_t& g);
    //  edge_index_t      get_edge_index(edge_descriptor ed, const graph_t& g);
    //  edge_descriptor   get_edge_descriptor(edge_index ei, const graph_t& g);
    //
    // ----------------------------
    // Label types:
    // ----------------------------
    //  vertex_label_t
    //  edge_label_t;
    //  vertex_label_ref_t
    //  edge_label_ref_t
    //
    // ----------------------------
    // Label operations:
    // ----------------------------
    //  vertex_label_ref_t vlabel (vertex_descriptor vd, const graph_t& g);
    //  vertex_label_ref_t vilabel(vertex_index_t vi, const graph_t& g);
    //  edge_label_ref_t   elabel (edge_descriptor ed, const graph_t& g);
    //  edge_label_ref_t   eilabel(edge_index_t ei, const graph_t& g);
    //  bool vlabel_equal         (vertex_label_ref_t, vertex_label_ref_t);
    //  bool vlabel_less          (vertex_label_ref_t, vertex_label_ref_t);
    //  bool vlabel_less_or_equal (vertex_label_ref_t, vertex_label_ref_t);
    //  bool elabel_equal         (edge_label_ref_t, edge_label_ref_t);
    //  bool elabel_less          (edge_label_ref_t, edge_label_ref_t);
    //  bool elabel_less_or_equal (edge_label_ref_t, edge_label_ref_t);
    //  bool               void_vlabel(vertex_label_ref_t);
    //  vertex_label_ref_t void_vlabel();
    //
    // =============================================================================

    template<class G>
    std::ostream& print_graph(std::ostream& out, const G& g)
    {
	typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
	typedef typename boost::graph_traits<G>::edge_descriptor edge_descriptor;
	typedef typename boost::graph_traits<G>::edge_iterator edge_iterator;
	out << "graph: ";
	for (std::pair<edge_iterator,edge_iterator> p = edges(g); p.first != p.second; ++p.first)
	{
	    edge_descriptor e = *p.first;
	    vertex_descriptor vs = source(e,g); 
	    vertex_descriptor vt = target(e,g);
	    out << "("<<vs<<","<<vt<<" " <<g[vs]<<" "<<g[e]<<" "<<g[vt]<< ") ";
	}
	return out << std::endl;
    }


    // *****************************************************************************
    //                          Edge
    // *****************************************************************************
    template<class VI_, class EI_>
    struct Edge_
    {
	typedef VI_ VI;
	typedef EI_ EI;	
	VI vi_from, vi_to;
	EI ei;
    };

    template<class VI, class EI>
    inline std::ostream& operator<<(std::ostream& out, const Edge_<VI,EI>& e)
    {
	return out<<"("<<e.vi_from<<","<<e.vi_to<<")";
    }


    // *****************************************************************************
    //                          EdgeCode
    // *****************************************************************************
    template<class VI_, class EI_, class VL, class EL, class VLR, class ELR>
    struct EdgeCode_
    {
	typedef VI_ VI;
	typedef EI_ EI;
	VI vi_from, vi_to;
	VL vl_from, vl_to;
	EL el;
	EdgeCode_(VI vifrom, VI vito, VLR vlfrom, ELR el, VLR vlto)
	    :vi_from(vifrom), vi_to(vito), vl_from(vlfrom), vl_to(vlto), el(el) {}

	bool is_forward() const { return vi_from < vi_to; }
	operator std::pair<VI,VI> () const { return std::pair<VI,VI>(vi_from,vi_to); }
    };

    template<class EdgeCode, class GraphOps>
    bool edgecode_equal(const EdgeCode& ec1, const EdgeCode& ec2, const GraphOps& ops)
    {
	return
	    ec1.vi_from == ec2.vi_from && ec1.vi_to == ec2.vi_to &&
	    ops.vlabel_equal(ec1.vl_from, ec2.vl_from) &&
	    ops.vlabel_equal(ec1.vl_to, ec2.vl_to) &&
	    ops.elabel_equal(ec1.el, ec2.el);
    }

    template<class VI, class EI, class VL, class EL, class VLR, class ELR>
    inline std::ostream& operator<<(std::ostream& out, const EdgeCode_<VI,EI,VL,EL,VLR,ELR>& ec)
    {
	return out<<"("<<ec.vi_from<<","<<ec.vi_to<<", "
		  <<ec.vl_from<<","<<ec.el<<","<<ec.vl_to<<")";
    }

    // *****************************************************************************
    //                          DFSCode
    // *****************************************************************************
    template<class EdgeCode>
    class DFSCode_
    {
	std::vector<EdgeCode> dfsc_;
    public:
	const EdgeCode& operator[] (int i) const { return dfsc_[i]; }
	bool empty() const { return dfsc_.empty(); }
	unsigned int size() const { return dfsc_.size(); }
	typedef typename std::vector<EdgeCode>::const_iterator const_iterator;
	const_iterator begin() const { return dfsc_.begin(); }
	const_iterator end() const { return dfsc_.end(); }
	typename EdgeCode::VI max_vertex() const;
	void push(const EdgeCode& ec) { dfsc_.push_back(ec); }
	void pop() { dfsc_.pop_back(); }
    };


    template<class EdgeCode>
    typename EdgeCode::VI DFSCode_<EdgeCode>::max_vertex() const
    {
	typename EdgeCode::VI m = 0;
	for (const_iterator i = begin(); i != end(); ++i)
	    m = std::max(m, std::max(i->vi_from, i->vi_to));
	return m;
    }


    template<class EdgeCode>
    std::ostream& operator<<(std::ostream& out, const DFSCode_<EdgeCode>& dfsc)
    {
	std::copy(dfsc.begin(), dfsc.end(), std::ostream_iterator<EdgeCode>(out, " "));
	return out;
    }


    // *****************************************************************************
    //                          RMPath
    // *****************************************************************************
    class RMPath
    {
	std::deque<unsigned int> rmp_;
    public:
	template<class EdgeCode>
	explicit RMPath(const DFSCode_<EdgeCode>& dfsc);
	unsigned int operator[] (int i) const { return rmp_[i]; }
	unsigned int size() const { return rmp_.size(); }
	unsigned int rightmost() const { return rmp_.back(); }
	friend std::ostream& operator<<(std::ostream& out, const RMPath& rpm);
    };

    template<class EdgeCode>
    RMPath::RMPath(const DFSCode_<EdgeCode>& dfsc)
    {
	typename EdgeCode::VI old_from = 0;
	for (int i = dfsc.size()-1; i >= 0; --i)
	    if (dfsc[i].is_forward() && (rmp_.empty() || old_from == dfsc[i].vi_to))
	    {
		rmp_.push_front(i);
		old_from = dfsc[i].vi_from;
	    }
    }

    inline std::ostream& operator<<(std::ostream& out, const RMPath& r)
    {
	std::copy(r.rmp_.begin(), r.rmp_.end(), std::ostream_iterator<unsigned int>(out, " "));
	return out;
    }


    // *****************************************************************************
    //                          SBG_
    // *****************************************************************************

    template<class GraphOps>
    class SBG_
    {
    public:
	typedef typename GraphOps::graph_t            G;
	typedef typename GraphOps::vertex_index_t     VI;
	typedef typename GraphOps::edge_index_t       EI;
	typedef Edge_<VI,EI> Edge;

	SBG_(const G* g, const Edge& e)    :prev_(0), edge_(e), rec_(0), graph_(g), depth_(1) {}
	SBG_(const SBG_* s, const Edge& e) :prev_(s), edge_(e), rec_(0), graph_(s->graph_), depth_(s->depth_+1) {}
	SBG_(const SBG_& s)                :prev_(s.prev_), edge_(s.edge_), rec_(0), graph_(s.graph_), depth_(s.depth_) {}
	~SBG_() { delete rec_; }

	int size() const { return depth_; }
	const Edge& operator[] (int i) const		{ if (!rec_) init_rec_(); return *rec_->e_.at(i); }
        bool has_vertex(typename Edge::VI vi) const	{ if (!rec_) init_rec_(); return rec_->vv_[vi]; }
        bool has_edge(const Edge& e) const		{ if (!rec_) init_rec_(); return rec_->ee_[e.ei]; }
	const G* get_graph() const { return graph_; }
    private:
	const SBG_* prev_;
	Edge        edge_;

	struct R
	{
	    std::vector<const Edge*> e_;
	    std::vector<short> vv_;
	    std::vector<bool> ee_;
	    R(const Edge& e, const G& g);
	    R(const Edge& e, const R* prev);
	};
	mutable R* rec_;
	void init_rec_() const;
	const G* graph_;
	int depth_;

	SBG_& operator= (const SBG_& s);
    };

    template<class GraphOps>
    SBG_<GraphOps>::R::R(const Edge& e, const G& g)
		      : e_(1, &e), vv_(num_vertices(g), 0), ee_(num_edges(g), false)
    {
	++vv_[e.vi_from];
	++vv_[e.vi_to];
	ee_[e.ei] = true;
    }

    template<class GraphOps>    
    SBG_<GraphOps>::R::R(const Edge& e, const R* prev)
	:e_(prev->e_), vv_(prev->vv_), ee_(prev->ee_)
    {
	e_.push_back(&e);
	++vv_[e.vi_from];
	++vv_[e.vi_to];
	ee_[e.ei] = true;
    }
    
    template<class GraphOps>
    void SBG_<GraphOps>::init_rec_() const
    {
	if (prev_)
	    rec_ = new R(edge_, prev_->rec_);
	else
	    rec_ = new R(edge_, *graph_);
    }

    template<class GraphOps>
    std::ostream& operator<<(std::ostream& out, const SBG_<GraphOps>& sbg)
    {
	out << "sbg:";
	for (int i = 0; i < sbg.size(); ++i)
	    out << " " << sbg[i];
	out << " at address: " << &sbg << " of the graph: " << sbg.get_graph();
	return out;
    }


    // *****************************************************************************
    //                          Traits
    // *****************************************************************************
    template<class GraphOps>
    struct Traits
    {
	typedef typename GraphOps::graph_t            G;
	typedef typename GraphOps::vertex_index_t     VI;
	typedef typename GraphOps::edge_index_t       EI;
	typedef typename GraphOps::vertex_label_t     VL;
	typedef typename GraphOps::vertex_label_ref_t VLR;
	typedef typename GraphOps::edge_label_t       EL;
	typedef typename GraphOps::edge_label_ref_t   ELR;
	
	typedef Edge_<VI,EI> Edge;
	typedef EdgeCode_<VI,EI,VL,EL,VLR,ELR> EdgeCode;
	typedef DFSCode_<EdgeCode> DFSCode;
	typedef SBG_<GraphOps> SBG;
	typedef std::list<SBG> Projected;
    };


    // *****************************************************************************
    //                          EdgeIterator
    // *****************************************************************************

    // -----------------------------------------------------
    //          EdgeIterator
    // -----------------------------------------------------
    template<class GraphOps, class DirTag>
    class EdgeIterator
    {
    public:
	typedef typename GraphOps::graph_t        Graph;
	typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
	typedef typename boost::graph_traits<Graph>::edge_descriptor edge_descriptor;	
	typedef typename Traits<GraphOps>::Edge Edge;
	EdgeIterator(typename Traits<GraphOps>::VI vi, const Graph& graph, const GraphOps& ops);
	bool is_end() const { return ou_iters_.first == ou_iters_.second; }
	const Edge& get_edge() const { return curr_edge_; }
	void increment();
    private:
	typedef typename boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;
	std::pair<out_edge_iterator,out_edge_iterator> ou_iters_;
	const Graph& graph_;
	void set_edge_ou()
	    {
		curr_edge_.vi_from = source(*ou_iters_.first, graph_);
		curr_edge_.vi_to   = target(*ou_iters_.first, graph_);
		curr_edge_.ei      = ops_.get_edge_index(*ou_iters_.first, graph_);
	    }
    protected:
	Edge curr_edge_;
	const GraphOps& ops_;
    };
    

    template<class GraphOps, class DirTag>
    EdgeIterator<GraphOps,DirTag>::EdgeIterator(typename Traits<GraphOps>::VI vi,
						const Graph& graph,
						const GraphOps& ops)
	:ou_iters_(out_edges(ops.get_vertex_descriptor(vi,graph), graph)),
	 graph_(graph), ops_(ops)
    {
	if (ou_iters_.first != ou_iters_.second)
	    set_edge_ou();
	else
	{
	    curr_edge_.vi_from = ops_.void_vertex_index();
	    curr_edge_.vi_to   = ops_.void_vertex_index();
	    curr_edge_.ei      = ops_.void_edge_index();
	}
    }

    template<class GraphOps, class DirTag>
    void EdgeIterator<GraphOps,DirTag>::increment()
    {
	if (++ou_iters_.first != ou_iters_.second)
	    set_edge_ou();
    }


    // -----------------------------------------------------
    //		EdgeIterator for boost::bidirectionalS
    // -----------------------------------------------------
    template<class GraphOps>
    class EdgeIterator<GraphOps, boost::bidirectional_tag>
    {
    public:
	typedef typename GraphOps::graph_t        Graph;
	typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
	typedef typename boost::graph_traits<Graph>::edge_descriptor edge_descriptor;
	typedef typename Traits<GraphOps>::Edge Edge;

	EdgeIterator(vertex_descriptor vd, const Graph& graph, const GraphOps& ops);
	bool is_end() const { return state_ == END; }
	const Edge& get_edge() const { return curr_edge_; }
	void increment();
	const Graph& get_graph() const { }
    private:
	typedef typename boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;
	typedef typename boost::graph_traits<Graph>::in_edge_iterator  in_edge_iterator;
	std::pair<out_edge_iterator,out_edge_iterator> ou_iters_;
	std::pair<in_edge_iterator,in_edge_iterator>   in_iters_;
	const Graph& graph_;
	enum State { OUT_ITER, IN_ITER, END } state_;
	void set_edge_ou()
	    {
		curr_edge_.vi_from = source(*ou_iters_.first, graph_);
		curr_edge_.vi_to   = target(*ou_iters_.first, graph_);
		curr_edge_.ei      = ops_.get_edge_index(*ou_iters_.first, graph_);
	    }
	void set_edge_in()
	    {
		curr_edge_.vi_from = target(*in_iters_.first, graph_);
		curr_edge_.vi_to   = source(*in_iters_.first, graph_);
		curr_edge_.ei      = ops_.get_edge_index(*in_iters_.first, graph_);
	    }
    protected:
	Edge curr_edge_;
	const GraphOps& ops_;
    };


    template<class GraphOps>
    EdgeIterator<GraphOps, boost::bidirectional_tag>::EdgeIterator(vertex_descriptor vd,
								   const Graph& graph,
								   const GraphOps& ops)
	:ou_iters_(out_edges(vd, graph)),
	 in_iters_(in_edges(vd, graph)),
	 graph_(graph), ops_(ops)
    {
	if (ou_iters_.first != ou_iters_.second)
	{
	    set_edge_ou();
	    state_ = OUT_ITER;
	}
	else if (in_iters_.first != in_iters_.second)
	{
	    set_edge_in();
	    state_ = IN_ITER;
	}
	else
	{
	    curr_edge_.vi_from = ops_.void_vertex_index();
	    curr_edge_.vi_to   = ops_.void_vertex_index();
	    curr_edge_.ei      = ops_.void_edge_index();
	    state_ = END;
	}
    }

    template<class GraphOps>
    void EdgeIterator<GraphOps, boost::bidirectional_tag>::increment()
    {
	switch (state_)
	{
	case OUT_ITER:
	    if (++ou_iters_.first != ou_iters_.second)
		set_edge_ou();
	    else
		if (in_iters_.first != in_iters_.second)
		{
		    set_edge_in();
		    state_ = IN_ITER;
		}
		else state_ = END;
	    break;

	case IN_ITER: 
	    if (++in_iters_.first != in_iters_.second)
		set_edge_in();
	    else
		state_ = END;
	    break;

	case END: assert(0); break;
	}
    }

    template<class GraphOps>
    struct EdgeIteratorSelector
    {
	typedef typename GraphOps::graph_t G;
	typedef typename boost::graph_traits<G>::directed_category DirTag;
	typedef EdgeIterator<GraphOps,DirTag> Type;
    };


    // -----------------------------------------------------
    //		EdgeIterator forward pure
    // -----------------------------------------------------

    template<class GraphOps>
    class ForwardPureEdgeIterator : public EdgeIteratorSelector<GraphOps>::Type
    {
	typedef typename EdgeIteratorSelector<GraphOps>::Type Base;
	using Base::curr_edge_;
	using Base::ops_;

	typename Traits<GraphOps>::VLR minlabel_;
	const typename Traits<GraphOps>::SBG* sbg_;
	void next_();
    public:
	ForwardPureEdgeIterator(const typename Traits<GraphOps>::Edge& e,
				typename Traits<GraphOps>::VLR minlabel,
				const typename Traits<GraphOps>::SBG* sbg,
				const GraphOps& ops)
	    :Base(e.vi_to, *sbg->get_graph(), ops),
	     minlabel_(minlabel),
	     sbg_(sbg)
	    { next_(); }
	void next() { this->increment(); next_(); }
    };

    template<class GraphOps>
    void ForwardPureEdgeIterator<GraphOps>::next_()
    {
	const typename Traits<GraphOps>::G& g = *sbg_->get_graph();
	for (; !this->is_end(); this->increment())
	{
	    if (ops_.vlabel_less_or_equal(minlabel_, ops_.vilabel(curr_edge_.vi_to, g)) &&
		!sbg_->has_vertex(curr_edge_.vi_to))
		break;
	}
    }

    // -----------------------------------------------------
    //		EdgeIterator forward rmpath
    // -----------------------------------------------------
    template<class GraphOps>
    class ForwardRMPathEdgeIterator : public EdgeIteratorSelector<GraphOps>::Type
    {
	typedef typename EdgeIteratorSelector<GraphOps>::Type Base;
	using Base::curr_edge_;
	using Base::ops_;

	typename Traits<GraphOps>::VLR minlabel_;
	const typename Traits<GraphOps>::SBG* sbg_;
	typename Traits<GraphOps>::Edge edge_;
	void next_();
    public:
	ForwardRMPathEdgeIterator(const typename Traits<GraphOps>::Edge& e,
				  typename Traits<GraphOps>::VLR minlabel,
				  const typename Traits<GraphOps>::SBG* sbg,
				  const GraphOps& ops)
	    :Base(e.vi_from, *sbg->get_graph(), ops),
	     minlabel_(minlabel),
	     sbg_(sbg), edge_(e)
	    { next_(); }
	void next() { this->increment(); next_(); }
    };

    template<class GraphOps>
    void ForwardRMPathEdgeIterator<GraphOps>::next_()
    {
	const typename Traits<GraphOps>::G& g = *sbg_->get_graph();
	for (; !this->is_end(); this->increment())
	{
	    const typename Traits<GraphOps>::Edge& e = curr_edge_;

	    typename Traits<GraphOps>::VLR tolabel  = ops_.vilabel(edge_.vi_to, g);
	    typename Traits<GraphOps>::VLR tolabel2 = ops_.vilabel(e.vi_to, g);

	    if (sbg_->has_edge(e) || sbg_->has_vertex(e.vi_to) || ops_.vlabel_less(tolabel2, minlabel_))
		continue;

	    typename Traits<GraphOps>::ELR elabel  = ops_.eilabel(edge_.ei, g);
	    typename Traits<GraphOps>::ELR elabel2 = ops_.eilabel(e.ei, g);

	    if (ops_.elabel_less(elabel, elabel2) ||
		(ops_.elabel_equal(elabel, elabel2) && ops_.vlabel_less_or_equal(tolabel, tolabel2) ) )
		break;
	}
    }


    // -----------------------------------------------------
    //		EdgeIterator backward
    // -----------------------------------------------------
    template<class GraphOps>
    class BackwardEdgeIterator : public EdgeIteratorSelector<GraphOps>::Type
    {
	typedef typename EdgeIteratorSelector<GraphOps>::Type Base;
	using Base::curr_edge_;
	using Base::ops_;
    public:
	BackwardEdgeIterator(const typename Traits<GraphOps>::Edge& e1,
			     const typename Traits<GraphOps>::Edge& e2,
			     const typename Traits<GraphOps>::SBG* sbg,
			     const GraphOps& ops);
    };

    template<class GraphOps>
    BackwardEdgeIterator<GraphOps>::BackwardEdgeIterator(const typename Traits<GraphOps>::Edge& e1,
							 const typename Traits<GraphOps>::Edge& e2,
							 const typename Traits<GraphOps>::SBG* sbg,
							 const GraphOps& ops)
	:Base(e2.vi_to, *sbg->get_graph(), ops)
    {
	//assert(e1.ei != e2.ei); // ???????????????
	const typename Traits<GraphOps>::G& g = *sbg->get_graph();
	typename Traits<GraphOps>::ELR e1label = ops_.eilabel(e1.ei, *sbg->get_graph());
	bool b = ops.vlabel_less_or_equal(ops.vilabel(e1.vi_to, g),
					  ops.vilabel(e2.vi_to, g));

	for (; !this->is_end(); this->increment())
	{
	    if (sbg->has_edge(curr_edge_))
		continue;

	    if (curr_edge_.vi_to == e1.vi_from)
	    {
		typename Traits<GraphOps>::ELR elabel = ops.eilabel(curr_edge_.ei,g);
		if (ops.elabel_less(e1label,elabel) || (ops.elabel_equal(e1label,elabel) && b))
		    break;
	    }
	}
    }

    // *****************************************************************************
    //                          Maps
    // *****************************************************************************

    template<class M1>
    inline typename M1::mapped_type&
    mapped_val(M1& mp1key, const typename M1::key_type& k1)
    {
	typedef typename M1::value_type value_type;
	return mp1key.insert(value_type(k1, typename M1::mapped_type())).first->second;
    }

    template<class M2>
    inline typename M2::mapped_type::mapped_type&
    mapped_val(M2& mp2key,
	   const typename M2::key_type& k1,
	   const typename M2::mapped_type::key_type& k2,
	   const typename M2::mapped_type::key_compare& c2)
    {
	typedef typename M2::value_type value_type;
	return mapped_val(mp2key.insert(value_type(k1, typename M2::mapped_type(c2))).first->second, k2);
    }


    template<class M3>
    inline typename M3::mapped_type::mapped_type::mapped_type&
    mapped_val(M3& mp3key,
	   const typename M3::key_type& k1,
	   const typename M3::mapped_type::key_type& k2,
	   const typename M3::mapped_type::mapped_type::key_type& k3,
	   const typename M3::mapped_type::key_compare& c2,
	   const typename M3::mapped_type::mapped_type::key_compare& c3)
    {
	typedef typename M3::value_type value_type;
	return mapped_val(mp3key.insert(value_type(k1, typename M3::mapped_type(c2))).first->second, k2, k3, c3);
    }

    
    template<class GraphOps>
    class Vlabel_less_
    {
	const GraphOps& ops_;
    public:
	explicit Vlabel_less_(const GraphOps& ops) :ops_(ops) {}
	bool operator() (typename Traits<GraphOps>::VLR vl1, typename Traits<GraphOps>::VLR vl2) const
	    { return ops_.vlabel_less(vl1, vl2); }
    };

    template<class GraphOps>    
    class Elabel_less_
    {
	const GraphOps& ops_;
    public:
	explicit Elabel_less_(const GraphOps& ops) :ops_(ops) {}
	bool operator() (typename Traits<GraphOps>::ELR el1, typename Traits<GraphOps>::ELR el2) const
	    { return ops_.elabel_less(el1, el2); }
    };


    // *****************************************************************************
    //                          MapTraits
    // *****************************************************************************
    template<class GraphOps>
    class MapTraits
    {
	// std::map wrapper without default constructor
	template<class K, class C, class T>
	struct map_ : public std::map<K,T,C> { explicit map_(const C& c) :std::map<K,T,C>(c) {} };
    public:
	typedef typename Traits<GraphOps>::VI VI;
	typedef typename Traits<GraphOps>::VL VL;
	typedef typename Traits<GraphOps>::EL EL;
	typedef typename Traits<GraphOps>::Projected Projected;

	typedef std::less<VI> VI_less;
	typedef Vlabel_less_<GraphOps> Vlabel_less;
	typedef Elabel_less_<GraphOps> Elabel_less;
	
	typedef map_<EL, Elabel_less, Projected>	Map_EL_P;
	typedef map_<VL, Vlabel_less, Projected>	Map_VL_P;
	typedef map_<EL, Elabel_less, Map_VL_P>		Map_EL_VL_P;
	typedef map_<VI, VI_less,     Map_EL_P>		Map_VI_EL_P;
	typedef map_<VI, VI_less,     Map_EL_VL_P>	Map_VI_EL_VL_P;
	typedef map_<VL, Vlabel_less, Map_EL_VL_P>	Map_VL_EL_VL_P;
    };


    // *****************************************************************************
    //                          gspan functions
    // *****************************************************************************

    template<class Projected>
    int support(const Projected& projected)
    {
	typedef typename Projected::value_type SBG;
        std::set<const typename SBG::G*> n;
        BOOST_FOREACH(const SBG& sbg, projected) n.insert(sbg.get_graph());
        return n.size();
    }


    template<class Projected, class DFSCode, class Output>
    void report(const Projected& projected, const DFSCode& dfsc, Output& result) { result(projected, dfsc); }


    template<class GraphOps>
    void one_edges(typename MapTraits<GraphOps>::Map_VL_EL_VL_P& m,
		   const typename Traits<GraphOps>::G& g,
		   const GraphOps& ops)
    {
        typename MapTraits<GraphOps>::Vlabel_less vl_less(ops);
        typename MapTraits<GraphOps>::Elabel_less el_less(ops);
        typedef typename Traits<GraphOps>::SBG SBG;
        typedef typename Traits<GraphOps>::G   G;
        typename boost::graph_traits<G>::vertex_iterator vi, viend;
	typedef typename EdgeIteratorSelector<GraphOps>::Type EdgeIter;

        for (boost::tie(vi,viend) = vertices(g); vi != viend; ++vi)
            for (EdgeIter iter(*vi, g, ops); !iter.is_end(); iter.increment())
            {
                typename Traits<GraphOps>::Edge e = iter.get_edge();
                typename Traits<GraphOps>::VLR vl_from = ops.vilabel(e.vi_from, g);
                typename Traits<GraphOps>::VLR vl_to   = ops.vilabel(e.vi_to, g);
                typename Traits<GraphOps>::ELR el      = ops.eilabel(e.ei, g);

                assert(! ops.void_vlabel(vl_from));
                assert(! ops.void_vlabel(vl_to));

                mapped_val(m, vl_from, el, vl_to, el_less, vl_less).push_back(SBG(&g, e));
            }
    }


    template<class GraphOps>
    bool project_is_min(const typename Traits<GraphOps>::Projected& projected,
                        typename Traits<GraphOps>::DFSCode& dfsc_min,
                        const typename Traits<GraphOps>::DFSCode& dfsc_tested,
                        const GraphOps& ops)
    {
        typedef typename Traits<GraphOps>::VI  VI;
        typedef typename Traits<GraphOps>::VLR VLR;
        typedef typename Traits<GraphOps>::ELR ELR;
        typedef typename Traits<GraphOps>::SBG SBG;
        typedef typename Traits<GraphOps>::Edge Edge;
        typedef typename Traits<GraphOps>::EdgeCode EdgeCode;

        if (! edgecode_equal(dfsc_min[dfsc_min.size()-1], dfsc_tested[dfsc_min.size()-1], ops))
            return false;

        // --------------------------------------------------------------
        // enumerate
        typedef typename MapTraits<GraphOps>::Map_EL_P    BckEdges;
        typedef typename MapTraits<GraphOps>::Map_EL_VL_P FwdEdges;

        RMPath rmpath(dfsc_min);
        VI maxtoc = dfsc_min[rmpath.rightmost()].vi_to;

        // backward
        {
            typename MapTraits<GraphOps>::Elabel_less el_less(ops);
            BckEdges bck_edges(el_less);
            VI newto = ops.void_vertex_index();
            bool flg = false;
            for (unsigned int i = 0; !flg && i < rmpath.size() - 1; ++i)
            {
                BOOST_FOREACH(const SBG& sbg, projected)
                {
                    const typename Traits<GraphOps>::G& g = *sbg.get_graph();
                    BackwardEdgeIterator<GraphOps> iter(sbg[rmpath[i]], sbg[rmpath.rightmost()], &sbg, ops);
                    if (!iter.is_end())
                    {
                        const Edge& e = iter.get_edge();
                        ELR el = ops.eilabel(e.ei,g);
                        mapped_val(bck_edges, el).push_back(SBG(&sbg, e));
                        newto = dfsc_min[rmpath[i]].vi_from;
                        flg = true;
                    }
                }
            }

            if (flg)
            {
                typename BckEdges::const_iterator i1 = bck_edges.begin();
                dfsc_min.push(EdgeCode(maxtoc, newto, ops.void_vlabel(), i1->first, ops.void_vlabel()));
                return project_is_min(i1->second, dfsc_min, dfsc_tested, ops);
            }
        }

        // forward
        {
            typename MapTraits<GraphOps>::Elabel_less el_less(ops);
            typename MapTraits<GraphOps>::Vlabel_less vl_less(ops);
            FwdEdges fwd_edges(el_less);
            VLR minlabel = dfsc_min[0].vl_from;
            VI newfrom = ops.void_vertex_index();
            bool flg = false;
            
            // forward pure
            BOOST_FOREACH(const SBG& sbg, projected)
            {
                const typename Traits<GraphOps>::G& g = *sbg.get_graph();
                for (ForwardPureEdgeIterator<GraphOps> iter(sbg[rmpath.rightmost()], minlabel, &sbg, ops);
                     !iter.is_end(); iter.next())
                {
                    const Edge& e = iter.get_edge();
                    ELR el = ops.eilabel(e.ei,g);
                    VLR vl = ops.vilabel(e.vi_to,g);
                    mapped_val(fwd_edges, el, vl, vl_less).push_back(SBG(&sbg, e));
                    newfrom = maxtoc;
                    flg = true;
                }
            }

            // forward rmpath
            for (int i = rmpath.size()-1; !flg && i >= 0; --i)
            {
                BOOST_FOREACH(const SBG& sbg, projected)
                {
                    const typename Traits<GraphOps>::G& g = *sbg.get_graph();
                    for (ForwardRMPathEdgeIterator<GraphOps> iter(sbg[rmpath[i]], minlabel, &sbg, ops);
                         !iter.is_end(); iter.next())
                    {
                        const Edge& e = iter.get_edge();
                        ELR el = ops.eilabel(e.ei,g);
                        VLR vl = ops.vilabel(e.vi_to,g);
                        mapped_val(fwd_edges, el, vl, vl_less).push_back(SBG(&sbg, e));
                        newfrom = dfsc_min[rmpath[i]].vi_from;
                        flg = true;
                    }
                }
            }

            if (flg)
            {
                typename FwdEdges::const_iterator i1 = fwd_edges.begin();
                typename FwdEdges::mapped_type::const_iterator i2 = i1->second.begin();
                dfsc_min.push(EdgeCode(newfrom, maxtoc+1, ops.void_vlabel(), i1->first, i2->first));
                return project_is_min(i2->second, dfsc_min, dfsc_tested, ops);
            }
        }

	return true;
    }

    template<class GraphOps>
    bool is_min(const typename Traits<GraphOps>::DFSCode& dfsc_tested, const GraphOps& ops)
    {
        typedef typename MapTraits<GraphOps>::Map_VL_EL_VL_P M3;
        typename MapTraits<GraphOps>::Vlabel_less vl_less(ops);

        std::auto_ptr<typename Traits<GraphOps>::G> graph(ops.create_graph(dfsc_tested));
        M3 root(vl_less);

        one_edges(root, *graph, ops);
        
        typename Traits<GraphOps>::DFSCode dfsc_min;
        typename M3::const_iterator i1 = root.begin();
        typename M3::mapped_type::const_iterator i2 = i1->second.begin();
        typename M3::mapped_type::mapped_type::const_iterator i3 = i2->second.begin();
        typename Traits<GraphOps>::EdgeCode ec(0, 1, i1->first, i2->first, i3->first);
        dfsc_min.push(ec);
        return project_is_min(i3->second, dfsc_min, dfsc_tested, ops);
    }
    
    enum RETS { NOTSUP, NOTMIN, HASCHILD, NOCHILD };

    template<class GraphOps, class Output>
    RETS project(const typename Traits<GraphOps>::Projected& projected,
                 typename Traits<GraphOps>::DFSCode& dfsc, int minsup,
                 const GraphOps& ops, Output& result)
    {
        typedef typename Traits<GraphOps>::VI  VI;
        typedef typename Traits<GraphOps>::VLR VLR;
        typedef typename Traits<GraphOps>::ELR ELR;
        typedef typename Traits<GraphOps>::SBG SBG;
        typedef typename Traits<GraphOps>::Edge Edge;
        typedef typename Traits<GraphOps>::EdgeCode EdgeCode;

        int sup = support(projected);
        if (sup < minsup)
            return NOTSUP;

        if (! is_min(dfsc, ops))
            return NOTMIN;

	report(projected, dfsc, result); 

        // --------------------------------------------------------------
        // enumerate
        typedef typename MapTraits<GraphOps>::Map_VI_EL_P    BckEdges;
        typedef typename MapTraits<GraphOps>::Map_VI_EL_VL_P FwdEdges;

	typename MapTraits<GraphOps>::VI_less vi_less;
        typename MapTraits<GraphOps>::Elabel_less el_less(ops);
        typename MapTraits<GraphOps>::Vlabel_less vl_less(ops);

        BckEdges bck_edges(vi_less);
        FwdEdges fwd_edges(vi_less);

        RMPath rmpath(dfsc);

        VI maxtoc   = dfsc[rmpath.rightmost()].vi_to;
        VLR minlabel = dfsc[0].vl_from;

	BOOST_FOREACH(const SBG& sbg, projected)
        {
            const typename Traits<GraphOps>::G& g = *sbg.get_graph();

            // backward
            for (unsigned int i = 0; i < rmpath.size() - 1; ++i)
            {
                BackwardEdgeIterator<GraphOps> iter(sbg[rmpath[i]], sbg[rmpath.rightmost()], &sbg, ops);
                if (!iter.is_end())
                {
                    const Edge& e = iter.get_edge();
                    VI vi  = dfsc[rmpath[i]].vi_from;
                    ELR el = ops.eilabel(e.ei,g);
                    mapped_val(bck_edges, vi, el, el_less).push_back(SBG(&sbg, e));

		    // virtual forward pure
		    VLR vl_to = ops.vilabel(e.vi_to,g);
		    mapped_val(fwd_edges, maxtoc, el, vl_to, el_less, vl_less).push_back(SBG(&sbg, e));

		    // virtual forward rmpath
		    VLR vl_from = ops.vilabel(e.vi_from,g);
		    mapped_val(fwd_edges, vi, el, vl_from, el_less, vl_less).push_back(SBG(&sbg, e));
                }
            }

            // forward
            for (ForwardPureEdgeIterator<GraphOps> iter(sbg[rmpath.rightmost()], minlabel, &sbg, ops);
                 !iter.is_end(); iter.next())
            {
                const Edge& e = iter.get_edge();
                ELR el = ops.eilabel(e.ei,g);
                VLR vl = ops.vilabel(e.vi_to,g);
                mapped_val(fwd_edges, maxtoc, el, vl, el_less, vl_less).push_back(SBG(&sbg, e));
            }

            for (int i = rmpath.size()-1; i >= 0; --i)
                for (ForwardRMPathEdgeIterator<GraphOps> iter(sbg[rmpath[i]], minlabel, &sbg, ops);
                     !iter.is_end(); iter.next())
                {
                    const Edge& e = iter.get_edge();
                    VI vi  = dfsc[rmpath[i]].vi_from;
                    ELR el = ops.eilabel(e.ei,g);
                    VLR vl = ops.vilabel(e.vi_to,g);
                    mapped_val(fwd_edges, vi, el, vl, el_less, vl_less).push_back(SBG(&sbg, e));
                }
	}

	RETS ret = NOCHILD;

        // --------------------------------------------------------------
        // recursive process SBG children
        // backward
        for (typename BckEdges::const_iterator it1 = bck_edges.begin(); it1 != bck_edges.end(); ++it1)
            for (typename BckEdges::mapped_type::const_iterator it2 = it1->second.begin();
                 it2 != it1->second.end(); ++it2)
            {
                EdgeCode ec(maxtoc, it1->first, ops.void_vlabel(), it2->first, ops.void_vlabel());
                dfsc.push(ec);
                project(it2->second, dfsc, minsup, ops, result);
                dfsc.pop();
		ret = HASCHILD;
            }
        // forward
	typedef typename FwdEdges::const_reverse_iterator FwdI1;
	typedef typename FwdEdges::mapped_type::const_iterator FwdI2;
	typedef typename FwdEdges::mapped_type::mapped_type::const_iterator FwdI3;
        for (FwdI1 i1 = fwd_edges.rbegin(); i1 != fwd_edges.rend(); ++i1)
            for (FwdI2 i2 = i1->second.begin(); i2 != i1->second.end(); ++i2)
                for (FwdI3 i3 = i2->second.begin(); i3 != i2->second.end(); ++i3)
                {
                    EdgeCode ec(i1->first, maxtoc+1, ops.void_vlabel(), i2->first, i3->first);
                    dfsc.push(ec);
                    project(i3->second, dfsc, minsup, ops, result);
                    dfsc.pop();
		    ret = HASCHILD;
                }

	if (ret == NOCHILD)
	    ;//report(projected, dfsc, result);
	return ret;
    }

    template<class TGraphIterator, class GraphOps, class Output>
    void gspan(TGraphIterator tg_begin, TGraphIterator tg_end, int minsup,
               const GraphOps& ops, Output& result)
    {
        typedef typename Traits<GraphOps>::G   G;
        typedef typename MapTraits<GraphOps>::Map_VL_EL_VL_P M3;

        typename MapTraits<GraphOps>::Vlabel_less vl_less(ops);
        M3 root(vl_less);

        for (; tg_begin != tg_end; ++tg_begin)
            one_edges(root, *tg_begin, ops);
	
        typename Traits<GraphOps>::DFSCode dfsc;

	typedef typename M3::const_iterator I1;
	typedef typename M3::mapped_type::const_iterator I2;
	typedef typename M3::mapped_type::mapped_type::const_iterator I3;
        for (I1 i1 = root.begin(); i1 != root.end(); ++i1)
            for (I2 i2 = i1->second.begin(); i2 != i1->second.end(); ++i2)
                for (I3 i3 = i2->second.begin(); i3 != i2->second.end(); ++i3)
                {
                    typename Traits<GraphOps>::EdgeCode ec(0, 1, i1->first, i2->first, i3->first);
#ifdef DEBUG_PRINT
                    std::cerr << "gspan(): top level iteration with edgecode: " << ec << std::endl;
#endif
                    dfsc.push(ec);
                    project(i3->second, dfsc, minsup, ops, result);
                    dfsc.pop();
                }
    }

}
#endif
