#ifndef CLOSEGRAPH_H_
#define CLOSEGRAPH_H_

#include <cassert>
#include <iostream>
#include <vector>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <iterator>
#include <memory>

#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>

namespace gSpan
{
    // *****************************************************************************
    //                          Edge
    // *****************************************************************************
    template<class Policy>
    struct Edge
    {
	typename Policy::vertex_index_t vi_from, vi_to;
	typename Policy::edge_index_t ei;
    };

    template<class Policy>
    inline std::ostream& operator<<(std::ostream& out, const Edge<Policy>& e)
    {
	return out<<"("<<e.vi_from<<","<<e.vi_to<<")";
    }

    // *****************************************************************************
    //                          EdgeCode
    // *****************************************************************************
    template<class Policy>
    class EdgeCode
    {
	typedef typename Policy::vertex_index_t VI;
    public:
	typename Policy::vertex_index_t vi_from, vi_to;
	typename Policy::vertex_label_t vl_from, vl_to;
	typename Policy::edge_label_t el;

	EdgeCode(typename Policy::vertex_index_t vi_from,
		 typename Policy::vertex_index_t vi_to,
		 typename Policy::vertex_label_ref_t vl_from,
		 typename Policy::edge_label_ref_t el,
		 typename Policy::vertex_label_ref_t vl_to)
	    :vi_from(vi_from), vi_to(vi_to), vl_from(vl_from), vl_to(vl_to), el(el) {}
	EdgeCode() {}
	bool is_forward() const { return vi_from < vi_to; }
	operator std::pair<VI,VI> () const { return std::pair<VI,VI>(vi_from,vi_to); }
    };

    template<class Policy, class Lcmp>
    bool equal(const EdgeCode<Policy>& ec1, const EdgeCode<Policy>& ec2, const Lcmp& lcmp)
    {
	return
	    ec1.vi_from == ec2.vi_from && ec1.vi_to == ec2.vi_to &&
	    lcmp.vlabel_equal(ec1.vl_from, ec2.vl_from) &&
	    lcmp.vlabel_equal(ec1.vl_to, ec2.vl_to) &&
	    lcmp.elabel_equal(ec1.el, ec2.el);
    }

    template<class Policy>
    inline std::ostream& operator<<(std::ostream& out, const EdgeCode<Policy>& ec)
    {
	return out<<"("<<ec.vi_from<<","<<ec.vi_to<<", "<<ec.vl_from<<","<<ec.el<<","<<ec.vl_to<<")";
    }


    // *****************************************************************************
    //                          DFSCode
    // *****************************************************************************
    template<class Policy> class DFSCode : public std::vector<EdgeCode<Policy> > {};
	
    template<class Policy>
    typename Policy::vertex_index_t max_vertex(const DFSCode<Policy>& dfsc)
    {
	typename Policy::vertex_index_t m = 0;
	for (typename DFSCode<Policy>::const_iterator i = dfsc.begin(); i != dfsc.end(); ++i)
	    m = std::max(m, std::max(i->vi_from, i->vi_to));
	return m;
    }

    template<class Policy>
    std::ostream& operator<<(std::ostream& out, const DFSCode<Policy>& dfsc)
    {
	std::copy(dfsc.begin(), dfsc.end(), std::ostream_iterator<EdgeCode<Policy> >(out, " "));
	return out;
    }

    // *****************************************************************************
    //                          RMPath
    // *****************************************************************************
    class RMPath
    {
	std::deque<unsigned int> rmp_;
    public:
	template<class Policy>
	explicit RMPath(const DFSCode<Policy>& dfsc);
	unsigned int operator[] (int i) const { return rmp_[i]; }
	unsigned int size() const { return rmp_.size(); }
	unsigned int rightmost() const { return rmp_.back(); }
	friend std::ostream& operator<<(std::ostream& out, const RMPath& r)
	    {
		std::copy(r.rmp_.begin(), r.rmp_.end(),
			  std::ostream_iterator<unsigned int>(out, " "));
		return out;
	    }
    };

    template<class Policy>
    RMPath::RMPath(const DFSCode<Policy>& dfsc)
    {
	typename Policy::vertex_index_t old_from = 0;
	for (int i = dfsc.size()-1; i >= 0; --i)
	    if (dfsc[i].is_forward() && (rmp_.empty() || old_from == dfsc[i].vi_to))
	    {
		rmp_.push_front(i);
		old_from = dfsc[i].vi_from;
	    }
    }

    // *****************************************************************************
    //                          SBG
    // *****************************************************************************
    template<class Policy>
    class SBG
    {
    public:
	typedef typename Policy::graph_t            G;
	typedef typename Policy::vertex_index_t     VI;
	typedef typename Policy::edge_index_t       EI;

	SBG(const G* g, const Edge<Policy>& e)   :prev_(0), edge_(e), rec_(0), graph_(g), depth_(1) {}
	SBG(const SBG* s, const Edge<Policy>& e) :prev_(s), edge_(e), rec_(0), graph_(s->graph_), depth_(s->depth_+1) {}
	SBG(const SBG& s)                :prev_(s.prev_), edge_(s.edge_), rec_(0), graph_(s.graph_), depth_(s.depth_) {}
	~SBG() { delete rec_; }

	int size() const { return depth_; }
	const Edge<Policy>& operator[] (int i) const	{ if (!rec_) init_rec_(); return *rec_->e_[i]; }
	bool has_vertex(VI vi) const		{ if (!rec_) init_rec_(); return rec_->vv_[vi]; }
	bool has_edge(const Edge<Policy>& e) const	{ if (!rec_) init_rec_(); return rec_->ee_[e.ei]; }
	const G* get_graph() const { return graph_; }
    private:
	const SBG* prev_;
	Edge<Policy> edge_;

	struct R
	{
	    std::vector<const Edge<Policy>*> e_;
	    std::vector<short> vv_;
	    std::vector<bool> ee_;
	    R(const Edge<Policy>& e, const G& g);
	    R(const Edge<Policy>& e, const R* prev);
	};
	mutable R* rec_;
	void init_rec_() const;
	const G* graph_;
	int depth_;

	SBG& operator= (const SBG& s);
    };


    template<class Policy>
    SBG<Policy>::R::R(const Edge<Policy>& e, const G& g)
	: e_(1, &e), vv_(num_vertices(g), 0), ee_(num_edges(g), false)
    {
	++vv_[e.vi_from];
	++vv_[e.vi_to];
	ee_[e.ei] = true;
    }

    template<class Policy>    
    SBG<Policy>::R::R(const Edge<Policy>& e, const R* prev)
	:e_(prev->e_), vv_(prev->vv_), ee_(prev->ee_)
    {
	e_.push_back(&e);
	++vv_[e.vi_from];
	++vv_[e.vi_to];
	ee_[e.ei] = true;
    }
    
    template<class Policy>
    void SBG<Policy>::init_rec_() const
    {
	if (prev_)
	    rec_ = new R(edge_, prev_->rec_);
	else
	    rec_ = new R(edge_, *graph_);
    }

    template<class Policy>
    std::ostream& operator<<(std::ostream& out, const SBG<Policy>& sbg)
    {
	out << "sbg:";
	for (int i = 0; i < sbg.size(); ++i)
	    out << " " << sbg[i];
	out << " at address: " << &sbg << " of the graph: " << sbg.get_graph();
	return out;
    }


    template<class Policy>
    class Projected : public std::list<SBG<Policy> > {};

    // *****************************************************************************
    //                          EdgeIterator
    // *****************************************************************************

    // -----------------------------------------------------
    //          EdgeIterator
    // -----------------------------------------------------
    template<class Policy, class DirTag>
    class EdgeIterator;
	
    // -----------------------------------------------------
    //		EdgeIterator for boost::undirectedS
    // -----------------------------------------------------


    // -----------------------------------------------------
    //		EdgeIterator for boost::bidirectionalS
    // -----------------------------------------------------
    template<class Policy>
    class EdgeIterator<Policy, boost::bidirectional_tag>
    {
	typedef typename Policy::graph_t G;
	typedef typename boost::graph_traits<G>::out_edge_iterator out_edge_iterator;
	typedef typename boost::graph_traits<G>::in_edge_iterator  in_edge_iterator;
	std::pair<out_edge_iterator,out_edge_iterator> ou_iters_;
	std::pair<in_edge_iterator,in_edge_iterator>   in_iters_;
	const G& graph_;
	enum State { OUT_ITER, IN_ITER, END } state_;
	void set_edge_ou();
	void set_edge_in();
    public:
	EdgeIterator(typename Policy::vertex_index_t vi, const G& g, const Policy& pl);
	bool is_end() const { return state_ == END; }
	const Edge<Policy>& get_edge() const { return curr_edge_; }
	void increment();
    protected:
	Edge<Policy> curr_edge_;
	const Policy& pl_;
	const G& get_graph() const { return graph_; }
    };

    template<class Policy>
    void EdgeIterator<Policy, boost::bidirectional_tag>::set_edge_ou()
    {
	curr_edge_.vi_from = source(*ou_iters_.first, graph_);
	curr_edge_.vi_to   = target(*ou_iters_.first, graph_);
	curr_edge_.ei      = pl_.get_edge_index(*ou_iters_.first, graph_);
    }

    template<class Policy>
    void EdgeIterator<Policy, boost::bidirectional_tag>::set_edge_in()
    {
	curr_edge_.vi_from = target(*in_iters_.first, graph_);
	curr_edge_.vi_to   = source(*in_iters_.first, graph_);
	curr_edge_.ei      = pl_.get_edge_index(*in_iters_.first, graph_);
    }

    template<class Policy>
    EdgeIterator<Policy, boost::bidirectional_tag>::EdgeIterator(typename Policy::vertex_index_t vi,
								 const G& graph,
								 const Policy& pl)
	:ou_iters_(out_edges(pl.get_vertex_descriptor(vi,graph), graph)),
	 in_iters_(in_edges(pl.get_vertex_descriptor(vi,graph), graph)),
	 graph_(graph), pl_(pl)
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
	    state_ = END;
    }

    template<class Policy>
    void EdgeIterator<Policy, boost::bidirectional_tag>::increment()
    {
	switch (state_)
	{
	case OUT_ITER:
	    if (++ou_iters_.first != ou_iters_.second)
		set_edge_ou();
	    else if (in_iters_.first != in_iters_.second)
	    {
		set_edge_in();
		state_ = IN_ITER;
	    }
	    else
		state_ = END;
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

    template<class Policy>
    struct EdgeIteratorSelector
    {
	typedef typename Policy::graph_t G;
	typedef typename boost::graph_traits<G>::directed_category DirTag;
	typedef EdgeIterator<Policy,DirTag> Type;
    };

    // *****************************************************************************
    //                         EdgeIterator forward pure
    // *****************************************************************************
    template<class Policy>
    class ForwardPureEdgeIterator : private EdgeIteratorSelector<Policy>::Type
    {
	typedef typename EdgeIteratorSelector<Policy>::Type Base;
	using Base::pl_;
	using Base::curr_edge_;
	typename Policy::vertex_label_ref_t minlabel_;
	const SBG<Policy>& sbg_;
	void next_();
    public:
	ForwardPureEdgeIterator(const Edge<Policy>& e,
				typename Policy::vertex_label_ref_t minlabel,
				const SBG<Policy>& sbg,
				const Policy& pl)
	    :Base(e.vi_to, *sbg.get_graph(), pl),
	     minlabel_(minlabel),
	     sbg_(sbg) { next_(); }
	bool is_end() const			{ return Base::is_end(); }
	const Edge<Policy>& get_edge() const	{ return Base::get_edge(); }
	void next() { this->increment(); next_(); }
    };

    template<class Policy>
    void ForwardPureEdgeIterator<Policy>::next_()
    {
	const typename Policy::graph_t& g = this->get_graph();
	for (; !this->is_end(); this->increment())
	{
	    if (pl_.vlabel_less_or_equal(minlabel_, pl_.vilabel(curr_edge_.vi_to, g)) &&
		!sbg_.has_vertex(curr_edge_.vi_to))
		break;
	}
    }

    // *****************************************************************************
    //                         EdgeIterator forward rmpath
    // *****************************************************************************
    template<class Policy>
    class ForwardRMPathEdgeIterator : private EdgeIteratorSelector<Policy>::Type
    {
	typedef typename EdgeIteratorSelector<Policy>::Type Base;
	using Base::pl_;
	using Base::curr_edge_;
	typename Policy::vertex_label_ref_t minlabel_;
	const SBG<Policy>& sbg_;
	Edge<Policy> edge_;
	void next_();
    public:
	ForwardRMPathEdgeIterator(const Edge<Policy>& e,
				  typename Policy::vertex_label_ref_t minlabel,
				  const SBG<Policy>& sbg,
				  const Policy& pl)
	    :Base(e.vi_from, *sbg.get_graph(), pl),
	     minlabel_(minlabel),
	     sbg_(sbg),
	     edge_(e) { next_(); }
	bool is_end() const		     { return Base::is_end(); }
	const Edge<Policy>& get_edge() const { return Base::get_edge(); }
	void next() { this->increment(); next_(); }
    };

    template<class Policy>
    void ForwardRMPathEdgeIterator<Policy>::next_()
    {
	const typename Policy::graph_t& g = this->get_graph();
	for (; !this->is_end(); this->increment())
	{
	    typename Policy::vertex_label_ref_t tolabel  = pl_.vilabel(edge_.vi_to, g);
	    typename Policy::vertex_label_ref_t tolabel2 = pl_.vilabel(curr_edge_.vi_to, g);

	    if (sbg_.has_edge(curr_edge_) ||
		sbg_.has_vertex(curr_edge_.vi_to) || pl_.vlabel_less(tolabel2, minlabel_))
		continue;

	    typename Policy::edge_label_ref_t elabel  = pl_.eilabel(edge_.ei, g);
	    typename Policy::edge_label_ref_t elabel2 = pl_.eilabel(curr_edge_.ei, g);

	    if (pl_.elabel_less(elabel, elabel2) ||
		(pl_.elabel_equal(elabel, elabel2) && pl_.vlabel_less_or_equal(tolabel, tolabel2) ) )
		break;
	}
    }


    // *****************************************************************************
    //                         EdgeIterator backward
    // *****************************************************************************
    template<class Policy>
    class BackwardEdgeIterator : private EdgeIteratorSelector<Policy>::Type
    {
	typedef typename EdgeIteratorSelector<Policy>::Type Base;
	using Base::pl_;
	using Base::curr_edge_;
    public:
	BackwardEdgeIterator(const Edge<Policy>& e1,
			     const Edge<Policy>& e2,
			     const SBG<Policy>& sbg,
			     const Policy& pl);
	bool is_end() const			{ return Base::is_end(); }
	const Edge<Policy>& get_edge() const	{ return Base::get_edge(); }
    };

    template<class Policy>
    BackwardEdgeIterator<Policy>::BackwardEdgeIterator(const Edge<Policy>& e1,
						       const Edge<Policy>& e2,
						       const SBG<Policy>& sbg,
						       const Policy& pl)
	:Base(e2.vi_to, *sbg.get_graph(), pl)
    {
	const typename Policy::graph_t& g = this->get_graph();
	typename Policy::edge_label_ref_t e1label = pl_.eilabel(e1.ei, g);
	bool b = pl.vlabel_less_or_equal(pl.vilabel(e1.vi_to, g),
					 pl.vilabel(e2.vi_to, g));

	for (; !this->is_end(); this->increment())
	{
	    if (sbg.has_edge(curr_edge_))
		continue;

	    if (curr_edge_.vi_to == e1.vi_from)
	    {
		typename Policy::edge_label_ref_t elabel = pl.eilabel(curr_edge_.ei,g);
		if (pl.elabel_less(e1label,elabel) || (pl.elabel_equal(e1label,elabel) && b))
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

    template<class Policy>
    class Vlabel_less
    {
	const Policy& pl_;
    public:
	explicit Vlabel_less(const Policy& pl) :pl_(pl) {}
	bool operator() (typename Policy::vertex_label_ref_t vl1,
			 typename Policy::vertex_label_ref_t vl2) const
	    { return pl_.vlabel_less(vl1, vl2); }
    };

    template<class Policy>
    class Elabel_less
    {
	const Policy& pl_;
    public:
	explicit Elabel_less(const Policy& pl) :pl_(pl) {}
	bool operator() (typename Policy::edge_label_ref_t el1,
			 typename Policy::edge_label_ref_t el2) const
	    { return pl_.elabel_less(el1, el2); }
    };

    template<class Policy>
    class MapTraits
    {
	// std::map wrapper without default constructor
	template<class K, class C, class T>
	struct map_ : public std::map<K,T,C> { explicit map_(const C& c) :std::map<K,T,C>(c) {} };
	typedef typename Policy::vertex_index_t VI;
	typedef typename Policy::vertex_label_t VL;
	typedef typename Policy::edge_label_t EL;
	typedef Projected<Policy> P;
    public:
	typedef map_<EL, Elabel_less<Policy>, P>		Map_EL_P;
	typedef map_<VL, Vlabel_less<Policy>, P>		Map_VL_P;
	typedef map_<EL, Elabel_less<Policy>, Map_VL_P>		Map_EL_VL_P;
	typedef map_<VI, std::less<VI>,       Map_EL_P>		Map_VI_EL_P;
	typedef map_<VI, std::less<VI>,       Map_EL_VL_P>	Map_VI_EL_VL_P;
	typedef map_<VL, Vlabel_less<Policy>, Map_EL_VL_P>	Map_VL_EL_VL_P;
    };

    // *****************************************************************************
    //                         functions
    // *****************************************************************************
    class MapValueCountDefault
    {
    public:
	template<class K>
	void operator() (const K&) { /* do nothing */ }
    };

    template<class M>
    class MapValueCount
    {
	M& m_;
    public:
	explicit MapValueCount(M& m) :m_(m) {}
	void operator() (const typename M::key_type& k) const { ++m_[k]; }
    };
    

    template<class Policy, class VlMapCount>
    void enum_one_edges(typename MapTraits<Policy>::Map_VL_EL_VL_P& m,
			const typename Policy::graph_t& g,
			const Policy& pl,
			VlMapCount& vlab_count)
    {
	Vlabel_less<Policy> vl_less(pl);
	Elabel_less<Policy> el_less(pl);
		
	typename boost::graph_traits<typename Policy::graph_t>::vertex_iterator v_it, v_itend;
	for (boost::tie(v_it,v_itend) = vertices(g); v_it != v_itend; ++v_it)
	{
	    vlab_count(pl.vdlabel(*v_it,g));
	    for (typename EdgeIteratorSelector<Policy>::Type iter(*v_it, g, pl);
		 !iter.is_end(); iter.increment())
	    {
		Edge<Policy> e = iter.get_edge();
		typename Policy::vertex_label_ref_t vl_from = pl.vilabel(e.vi_from, g);
		typename Policy::vertex_label_ref_t vl_to   = pl.vilabel(e.vi_to, g);
		typename Policy::edge_label_ref_t el        = pl.eilabel(e.ei, g);

		assert(! pl.void_vlabel(vl_from));
		assert(! pl.void_vlabel(vl_to));

		mapped_val(m, vl_from, el, vl_to, el_less, vl_less).push_back(SBG<Policy>(&g, e));
	    }
	}
    }


    template<class Policy>
    bool project_is_min(const Projected<Policy>& projected,
			DFSCode<Policy>& dfsc_min,
			const DFSCode<Policy>& dfsc_tested,
			const Policy& pl)
    {
	typedef typename Policy::vertex_index_t VI;
	typedef typename Policy::vertex_label_ref_t VLR;
	typedef typename Policy::edge_label_ref_t ELR;
	typedef SBG<Policy> SBG;
	typedef Edge<Policy> Edge;
	typedef EdgeCode<Policy> EdgeCode;

	if (! equal(dfsc_min[dfsc_min.size()-1], dfsc_tested[dfsc_min.size()-1], pl))
	    return false;

	// --------------------------------------------------------------
	// enumerate
	typedef typename MapTraits<Policy>::Map_EL_P    BckEdges;
	typedef typename MapTraits<Policy>::Map_EL_VL_P FwdEdges;

	RMPath rmpath(dfsc_min);
	VI maxtoc = dfsc_min[rmpath.rightmost()].vi_to;

	// backward
	{
	    Elabel_less<Policy> el_less(pl);
	    BckEdges bck_edges(el_less);
	    VI newto = pl.void_vertex_index();
	    bool flg = false;
	    for (unsigned int i = 0; !flg && i < rmpath.size() - 1; ++i)
	    {
		BOOST_FOREACH(const SBG& sbg, projected)
		{
		    const typename Policy::graph_t& g = *sbg.get_graph();
		    BackwardEdgeIterator<Policy> iter(sbg[rmpath[i]], sbg[rmpath.rightmost()], sbg, pl);
		    if (!iter.is_end())
		    {
			const Edge& e = iter.get_edge();
			ELR el = pl.eilabel(e.ei,g);
			mapped_val(bck_edges, el).push_back(SBG(&sbg, e));
			newto = dfsc_min[rmpath[i]].vi_from;
			flg = true;
		    }
		}
	    }

	    if (flg)
	    {
		typename BckEdges::const_iterator i1 = bck_edges.begin();
		dfsc_min.push_back(EdgeCode(maxtoc, newto, pl.void_vlabel(), i1->first, pl.void_vlabel()));
		return project_is_min(i1->second, dfsc_min, dfsc_tested, pl);
	    }
	}

	// forward
	{
	    Elabel_less<Policy> el_less(pl);
	    Vlabel_less<Policy> vl_less(pl);
	    FwdEdges fwd_edges(el_less);
	    VLR minlabel = dfsc_min[0].vl_from;
	    VI newfrom = pl.void_vertex_index();
	    bool flg = false;
            
	    // forward pure
	    BOOST_FOREACH(const SBG& sbg, projected)
	    {
		const typename Policy::graph_t& g = *sbg.get_graph();
		for (ForwardPureEdgeIterator<Policy> iter(sbg[rmpath.rightmost()], minlabel, sbg, pl);
		     !iter.is_end(); iter.next())
		{
		    const Edge& e = iter.get_edge();
		    ELR el = pl.eilabel(e.ei,g);
		    VLR vl = pl.vilabel(e.vi_to,g);
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
		    const typename Policy::graph_t& g = *sbg.get_graph();
		    for (ForwardRMPathEdgeIterator<Policy> iter(sbg[rmpath[i]], minlabel, sbg, pl);
			 !iter.is_end(); iter.next())
		    {
			const Edge& e = iter.get_edge();
			ELR el = pl.eilabel(e.ei,g);
			VLR vl = pl.vilabel(e.vi_to,g);
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
		dfsc_min.push_back(EdgeCode(newfrom, maxtoc+1, pl.void_vlabel(), i1->first, i2->first));
		return project_is_min(i2->second, dfsc_min, dfsc_tested, pl);
	    }
	}

	return true;
    }


    template<class Policy>
    bool is_min(const DFSCode<Policy>& dfsc_tested, const Policy& pl)
    {
	typedef typename MapTraits<Policy>::Map_VL_EL_VL_P M3;
	Vlabel_less<Policy> vl_less(pl);
	std::auto_ptr<typename Policy::graph_t> graph(pl.create_graph(dfsc_tested));

	M3 root(vl_less);
	MapValueCountDefault foo;
	enum_one_edges(root, *graph, pl, foo);

	DFSCode<Policy> dfsc_min;
	typename M3::const_iterator i1 = root.begin();
	typename M3::mapped_type::const_iterator i2 = i1->second.begin();
	typename M3::mapped_type::mapped_type::const_iterator i3 = i2->second.begin();
	dfsc_min.push_back(EdgeCode<Policy>(0, 1, i1->first, i2->first, i3->first));
	return project_is_min(i3->second, dfsc_min, dfsc_tested, pl);
    }


    namespace gspan_detail
    {
	template<class Policy>
	int support(const Projected<Policy>& projected, const Policy& pl)
	{
	    std::set<const typename Policy::graph_t*> n;
	    BOOST_FOREACH(const SBG<Policy>& sbg, projected) n.insert(sbg.get_graph());
	    return n.size();
	}

	template<class Policy, class Output>
	void project(const Projected<Policy>& projected,
		     DFSCode<Policy>& dfsc,
		     int minsup,
		     const Policy& pl,
		     Output& result)
	{
	    typedef typename Policy::vertex_index_t VI;
	    typedef typename Policy::vertex_label_ref_t VLR;
	    typedef typename Policy::edge_label_ref_t ELR;
	    typedef SBG<Policy> SBG;
	    typedef Edge<Policy> Edge;
	    typedef EdgeCode<Policy> EdgeCode;

	    int sup = support(projected, pl);
	    if (sup < minsup)
		return;

	    if (! is_min(dfsc, pl))
		return;

	    result(dfsc, projected);

	    // --------------------------------------------------------------
	    // enumerate
	    typedef typename MapTraits<Policy>::Map_VI_EL_P    BckEdges;
	    typedef typename MapTraits<Policy>::Map_VI_EL_VL_P FwdEdges;

	    std::less<VI> vi_less;
	    Elabel_less<Policy> el_less(pl);
	    Vlabel_less<Policy> vl_less(pl);

	    BckEdges bck_edges(vi_less);
	    FwdEdges fwd_edges(vi_less);

	    RMPath rmpath(dfsc);

	    VI maxtoc   = dfsc[rmpath.rightmost()].vi_to;
	    VLR minlabel = dfsc[0].vl_from;

	    BOOST_FOREACH(const SBG& sbg, projected)
	    {
		const typename Policy::graph_t& g = *sbg.get_graph();

		// backward
		for (unsigned int i = 0; i < rmpath.size() - 1; ++i)
		{
		    BackwardEdgeIterator<Policy> iter(sbg[rmpath[i]], sbg[rmpath.rightmost()], sbg, pl);
		    if (!iter.is_end())
		    {
			const Edge& e = iter.get_edge();
			VI vi  = dfsc[rmpath[i]].vi_from;
			ELR el = pl.eilabel(e.ei,g);
			mapped_val(bck_edges, vi, el, el_less).push_back(SBG(&sbg, e));

			// virtual forward pure
			VLR vl_to = pl.vilabel(e.vi_to,g);
			mapped_val(fwd_edges, maxtoc, el, vl_to, el_less, vl_less).push_back(SBG(&sbg, e));

			// virtual forward rmpath
			VLR vl_from = pl.vilabel(e.vi_from,g);
			mapped_val(fwd_edges, vi, el, vl_from, el_less, vl_less).push_back(SBG(&sbg, e));
		    }
		}

		// forward
		for (ForwardPureEdgeIterator<Policy> iter(sbg[rmpath.rightmost()], minlabel, sbg, pl);
		     !iter.is_end(); iter.next())
		{
		    const Edge& e = iter.get_edge();
		    ELR el = pl.eilabel(e.ei,g);
		    VLR vl = pl.vilabel(e.vi_to,g);
		    mapped_val(fwd_edges, maxtoc, el, vl, el_less, vl_less).push_back(SBG(&sbg, e));
		}

		for (int i = rmpath.size()-1; i >= 0; --i)
		    for (ForwardRMPathEdgeIterator<Policy> iter(sbg[rmpath[i]], minlabel, sbg, pl);
			 !iter.is_end(); iter.next())
		    {
			const Edge& e = iter.get_edge();
			VI vi  = dfsc[rmpath[i]].vi_from;
			ELR el = pl.eilabel(e.ei,g);
			VLR vl = pl.vilabel(e.vi_to,g);
			mapped_val(fwd_edges, vi, el, vl, el_less, vl_less).push_back(SBG(&sbg, e));
		    }
	    }

	    // --------------------------------------------------------------
	    // recursive process SBG children
	    // backward
	    for (typename BckEdges::const_iterator it1 = bck_edges.begin(); it1 != bck_edges.end(); ++it1)
		for (typename BckEdges::mapped_type::const_iterator it2 = it1->second.begin();
		     it2 != it1->second.end(); ++it2)
		{
		    EdgeCode ec(maxtoc, it1->first, pl.void_vlabel(), it2->first, pl.void_vlabel());
		    dfsc.push_back(ec);
		    project(it2->second, dfsc, minsup, pl, result);
		    dfsc.pop_back();
		}

	    // forward
	    typedef typename FwdEdges::const_reverse_iterator FwdI1;
	    typedef typename FwdEdges::mapped_type::const_iterator FwdI2;
	    typedef typename FwdEdges::mapped_type::mapped_type::const_iterator FwdI3;
	    for (FwdI1 i1 = fwd_edges.rbegin(); i1 != fwd_edges.rend(); ++i1)
		for (FwdI2 i2 = i1->second.begin(); i2 != i1->second.end(); ++i2)
		    for (FwdI3 i3 = i2->second.begin(); i3 != i2->second.end(); ++i3)
		    {
			EdgeCode ec(i1->first, maxtoc+1, pl.void_vlabel(), i2->first, i3->first);
			dfsc.push_back(ec);
			project(i3->second, dfsc, minsup, pl, result);
			dfsc.pop_back();
		    }
	    return;
	}

    } // end: namespace gspan_detail

    namespace closegraph_detail
    {
	template<class Policy>
	int support(const Projected<Policy>& projected, const Policy& pl)
	{
	    std::set<const typename Policy::graph_t*> n;
	    BOOST_FOREACH(const SBG<Policy>& sbg, projected) n.insert(sbg.get_graph());
	    return n.size();
	}

	template<class Policy, class Output>
	void project(const Projected<Policy>& projected,
		     DFSCode<Policy>& dfsc,
		     int minsup,
		     const Policy& pl,
		     Output& result)
	{
	    typedef typename Policy::vertex_index_t VI;
	    typedef typename Policy::vertex_label_ref_t VLR;
	    typedef typename Policy::edge_label_ref_t ELR;
	    typedef SBG<Policy> SBG;
	    typedef Edge<Policy> Edge;
	    typedef EdgeCode<Policy> EdgeCode;

	    int sup = support(projected, pl);
	    if (sup < minsup)
		return;

	    if (! is_min(dfsc, pl))
		return;

	    //report(projected, dfsc, result); 

	    // --------------------------------------------------------------
	    // enumerate
	    typedef typename MapTraits<Policy>::Map_VI_EL_P    BckEdges;
	    typedef typename MapTraits<Policy>::Map_VI_EL_VL_P FwdEdges;

	    std::less<VI> vi_less;
	    Elabel_less<Policy> el_less(pl);
	    Vlabel_less<Policy> vl_less(pl);

	    BckEdges bck_edges(vi_less);
	    FwdEdges fwd_edges(vi_less);

	    RMPath rmpath(dfsc);

	    VI maxtoc   = dfsc[rmpath.rightmost()].vi_to;
	    VLR minlabel = dfsc[0].vl_from;

	    BOOST_FOREACH(const SBG& sbg, projected)
	    {
		const typename Policy::graph_t& g = *sbg.get_graph();

		// backward
		for (unsigned int i = 0; i < rmpath.size() - 1; ++i)
		{
		    BackwardEdgeIterator<Policy> iter(sbg[rmpath[i]], sbg[rmpath.rightmost()], sbg, pl);
		    if (!iter.is_end())
		    {
			const Edge& e = iter.get_edge();
			VI vi  = dfsc[rmpath[i]].vi_from;
			ELR el = pl.eilabel(e.ei,g);
			mapped_val(bck_edges, vi, el, el_less).push_back(SBG(&sbg, e));

			// virtual forward pure
			VLR vl_to = pl.vilabel(e.vi_to,g);
			mapped_val(fwd_edges, maxtoc, el, vl_to, el_less, vl_less).push_back(SBG(&sbg, e));

			// virtual forward rmpath
			VLR vl_from = pl.vilabel(e.vi_from,g);
			mapped_val(fwd_edges, vi, el, vl_from, el_less, vl_less).push_back(SBG(&sbg, e));
		    }
		}

		// forward
		for (ForwardPureEdgeIterator<Policy> iter(sbg[rmpath.rightmost()], minlabel, sbg, pl);
		     !iter.is_end(); iter.next())
		{
		    const Edge& e = iter.get_edge();
		    ELR el = pl.eilabel(e.ei,g);
		    VLR vl = pl.vilabel(e.vi_to,g);
		    mapped_val(fwd_edges, maxtoc, el, vl, el_less, vl_less).push_back(SBG(&sbg, e));
		}

		for (int i = rmpath.size()-1; i >= 0; --i)
		    for (ForwardRMPathEdgeIterator<Policy> iter(sbg[rmpath[i]], minlabel, sbg, pl);
			 !iter.is_end(); iter.next())
		    {
			const Edge& e = iter.get_edge();
			VI vi  = dfsc[rmpath[i]].vi_from;
			ELR el = pl.eilabel(e.ei,g);
			VLR vl = pl.vilabel(e.vi_to,g);
			mapped_val(fwd_edges, vi, el, vl, el_less, vl_less).push_back(SBG(&sbg, e));
		    }
	    }

	    // --------------------------------------------------------------
	    // recursive process SBG children
	    // backward
	    for (typename BckEdges::const_iterator it1 = bck_edges.begin(); it1 != bck_edges.end(); ++it1)
		for (typename BckEdges::mapped_type::const_iterator it2 = it1->second.begin();
		     it2 != it1->second.end(); ++it2)
		{
		    EdgeCode ec(maxtoc, it1->first, pl.void_vlabel(), it2->first, pl.void_vlabel());
		    dfsc.push_back(ec);
		    project(it2->second, dfsc, minsup, pl, result);
		    dfsc.pop_back();
		}

	    // forward
	    typedef typename FwdEdges::const_reverse_iterator FwdI1;
	    typedef typename FwdEdges::mapped_type::const_iterator FwdI2;
	    typedef typename FwdEdges::mapped_type::mapped_type::const_iterator FwdI3;
	    for (FwdI1 i1 = fwd_edges.rbegin(); i1 != fwd_edges.rend(); ++i1)
		for (FwdI2 i2 = i1->second.begin(); i2 != i1->second.end(); ++i2)
		    for (FwdI3 i3 = i2->second.begin(); i3 != i2->second.end(); ++i3)
		    {
			EdgeCode ec(i1->first, maxtoc+1, pl.void_vlabel(), i2->first, i3->first);
			dfsc.push_back(ec);
			project(i3->second, dfsc, minsup, pl, result);
			dfsc.pop_back();
		    }
	    return;
	}
	
    } // end: namespace closegraph_detail

    template<class TGraphIterator, class Output, class Policy>
    void gspan(TGraphIterator tg_begin,
	       TGraphIterator tg_end,
	       int minsup,
	       const Policy& pl,
	       Output& result)
    {
	MapValueCountDefault mvc;

	typedef typename MapTraits<Policy>::Map_VL_EL_VL_P M3;
	Vlabel_less<Policy> vl_less(pl);
	M3 root(vl_less);

	for (; tg_begin != tg_end; ++tg_begin)
	    enum_one_edges(root, *tg_begin, pl, mvc);

	DFSCode<Policy> dfsc;
	
	typedef typename M3::const_iterator I1;
	typedef typename M3::mapped_type::const_iterator I2;
	typedef typename M3::mapped_type::mapped_type::const_iterator I3;
	for (I1 i1 = root.begin(); i1 != root.end(); ++i1)
	    for (I2 i2 = i1->second.begin(); i2 != i1->second.end(); ++i2)
		for (I3 i3 = i2->second.begin(); i3 != i2->second.end(); ++i3)
		{
		    EdgeCode<Policy> ec(0, 1, i1->first, i2->first, i3->first);
#ifdef DEBUG_PRINT
		    std::cerr << "gspan(): top level iteration with edgecode: " << ec << std::endl;
#endif
		    dfsc.push_back(ec);
		    gspan_detail::project(i3->second, dfsc, minsup, pl, result);
		    dfsc.pop_back();
		}
    }


    template<class TGraphIterator, class Output, class Policy>
    void closegraph(TGraphIterator tg_begin,
		    TGraphIterator tg_end,
		    int minsup,
		    const Policy& pl,
		    Output& result)
    {
	std::map<typename Policy::vertex_label_t, int> vlab_count;
	MapValueCount<std::map<typename Policy::vertex_label_t, int> > mvc(vlab_count);

	typedef typename MapTraits<Policy>::Map_VL_EL_VL_P M3;
	Vlabel_less<Policy> vl_less(pl);
	M3 root(vl_less);

	for (; tg_begin != tg_end; ++tg_begin)
	    enum_one_edges(root, *tg_begin, pl, mvc);

	DFSCode<Policy> dfsc;
	
	typedef typename M3::const_iterator I1;
	typedef typename M3::mapped_type::const_iterator I2;
	typedef typename M3::mapped_type::mapped_type::const_iterator I3;
	for (I1 i1 = root.begin(); i1 != root.end(); ++i1)
	    for (I2 i2 = i1->second.begin(); i2 != i1->second.end(); ++i2)
		for (I3 i3 = i2->second.begin(); i3 != i2->second.end(); ++i3)
		{
		    EdgeCode<Policy> ec(0, 1, i1->first, i2->first, i3->first);
#ifdef DEBUG_PRINT
		    std::cerr << "closegraph(): top level iteration with edgecode: " << ec << std::endl;
#endif
		    dfsc.push_back(ec);
		    closegraph_detail::project(i3->second, dfsc, minsup, pl, result);
		    dfsc.pop_back();
		}
    }

} // end: namespace gSpan

#endif
