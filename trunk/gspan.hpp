#ifndef CLOSEGRAPH_H_
#define CLOSEGRAPH_H_

#include "edge_iterator.hpp"

#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>

#include <cassert>
#include <iostream>
#include <vector>
#include <deque>
#include <queue>
#include <list>
#include <map>
#include <set>
#include <iterator>
#include <memory>
#include <algorithm>

#ifndef BR
#define BR asm volatile ("int3;")
#endif

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
	return out<<e.ei<<"("<<e.vi_from<<","<<e.vi_to<<") ";
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
	bool is_backward() const { return ! is_forward(); }
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
    template<class Policy>
    class DFSCode : public std::vector<EdgeCode<Policy> >
    {
    public:
	typename Policy::vertex_label_ref_t
	vlabel(typename Policy::vertex_index_t vi) const;
    };

    class NotFoundException {};

    template<class Policy>
    typename Policy::vertex_index_t max_vertex(const DFSCode<Policy>& dfsc)
    {
	typename Policy::vertex_index_t m = 0;
	for (typename DFSCode<Policy>::const_iterator i = dfsc.begin(); i != dfsc.end(); ++i)
	    m = std::max(m, std::max(i->vi_from, i->vi_to));
	return m;
    }

    template<class Policy>
    typename Policy::vertex_label_ref_t
    DFSCode<Policy>::vlabel(typename Policy::vertex_index_t vi) const
    {
	BOOST_FOREACH(const EdgeCode<Policy>& ec, *this)
	{
	    if (ec.vi_from == vi)
		return ec.vl_from;
	    if (ec.vi_to == vi)
		return ec.vl_to;
	}
	throw NotFoundException();
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
	std::deque<int> rmp_;
    public:
	template<class Policy>
	explicit RMPath(const DFSCode<Policy>& dfsc);
	int operator[] (int i) const { return rmp_[i]; }
	int size() const { return rmp_.size(); }
	int rightmost() const { return rmp_.back(); }
	friend std::ostream& operator<<(std::ostream& out, const RMPath& r)
	    {
		std::copy(r.rmp_.begin(), r.rmp_.end(),
			  std::ostream_iterator<int>(out, " "));
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

	SBG(const G* g, const Edge<Policy>& e)
	    :prev_(0), edge_(e), rec_(0), graph_(g), depth_(1) { automorph_list = 0; }

	SBG(const SBG* s, const Edge<Policy>& e)
	    :prev_(s), edge_(e), rec_(0), graph_(s->graph_), depth_(s->depth_+1)
	    { automorph_list = 0; }

	SBG(const SBG& s)
	    :prev_(s.prev_), edge_(s.edge_), rec_(0), graph_(s.graph_), depth_(s.depth_)
	    {
		automorph_list = s.automorph_list;
	    }

	~SBG()
	    { delete rec_; }

	int size() const { return depth_; }
	const Edge<Policy>& operator[] (int i) const	{ if (!rec_) init_rec_(); return *rec_->get_edge(i); }
	bool has_vertex(VI vi) const		        { if (!rec_) init_rec_(); return rec_->has_vertex(vi); }
	bool has_edge(EI ei) const       	        { if (!rec_) init_rec_(); return rec_->has_edge(ei); }
	const G* get_graph() const { return graph_; }
	const Edge<Policy>& last_edge() const { return edge_; }
	SBG* automorph_list;
    private:
	const SBG* prev_;
	Edge<Policy> edge_;

	class R
	{
	    std::vector<const Edge<Policy>*> e_;
	    std::vector<short> vv_;
	    std::vector<bool> ee_;
	public:
	    R(const Edge<Policy>& e, const G& g);
	    R(const Edge<Policy>& e, const R* prev);
	    const Edge<Policy>* get_edge(int i) const { return e_[i]; }
	    bool has_vertex(VI vi) const { return vv_[vi]; }
	    bool has_edge(EI ei) const { return ee_[ei]; }
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
 

    // *****************************************************************************
    //                          Projected
    // *****************************************************************************

    template<class Policy>
    int num_common_edges(const SBG<Policy>& sbg1, const SBG<Policy>& sbg2)
    {
	assert(sbg1.get_graph() == sbg2.get_graph());
	assert(sbg1.size() == sbg2.size());
	std::set<typename Policy::edge_index_t> ei1;
	std::set<typename Policy::edge_index_t> ei2;
	int num = sbg1.size();
	for (int i = 0; i < num; ++i)
	{
	    ei1.insert(sbg1[i].ei);
	    ei2.insert(sbg2[i].ei);
	}
	
	std::vector<typename Policy::edge_index_t> common;
	std::set_intersection(ei1.begin(), ei1.end(), ei2.begin(), ei2.end(), std::back_inserter(common));
	return common.size();
    }
   
    //
    // Projected part visible for user
    //

    template<class Policy>
    class SubgraphsOfOneGraph
    {	
	typedef std::vector<SBG<Policy>* > SBGS_PTR;
	SBGS_PTR sbgs_ptr_;
	int count_non_overlapped_;
    public:
	typedef typename SBGS_PTR::const_iterator const_iterator;
	const_iterator begin() const { return sbgs_ptr_.begin(); }
	const_iterator end()   const { return sbgs_ptr_.end(); }
	int support() const { return count_non_overlapped_; }
	void push(SBG<Policy>* s);
	void merge(SubgraphsOfOneGraph& sog);
    };

    template<class Policy>
    void SubgraphsOfOneGraph<Policy>::push(SBG<Policy>* s)
    {
	const int ne = s->size();
	bool has_common_edge = false;
	BOOST_FOREACH(SBG<Policy>* sbg, sbgs_ptr_)
	{
	    int nce = num_common_edges(*s, *sbg);
	    if (nce == ne)
	    {
		assert(s->automorph_list == 0);
		s->automorph_list = sbg->automorph_list;
		sbg->automorph_list = s;
		return;
	    }

	    if (nce > 0)
		has_common_edge = true;
	}
	
	sbgs_ptr_.push_back(s);

	if (! has_common_edge)
	    ++count_non_overlapped_;
    }

    template<class Policy>
    void SubgraphsOfOneGraph<Policy>::merge(SubgraphsOfOneGraph& sog)
    {
	const typename SBGS_PTR::const_iterator iend = sog.sbgs_ptr_.end();
	for (typename SBGS_PTR::const_iterator i = sog.sbgs_ptr_.begin(); i != iend; ++i)
	    push(*i);
	sog.sbgs_ptr_.clear();
	sog.count_non_overlapped_ = 0;
    }


    template<class Policy>
    class SubgraphsOfManyGraph
    {
	typedef typename Policy::graph_t Graph;
	typedef SubgraphsOfOneGraph<Policy> SOG;
	typedef std::map<const Graph* const, SOG> G_SOG;
	G_SOG g_sog_;
    public:
	typedef typename G_SOG::const_iterator const_iterator;
	const_iterator begin() const	{ return g_sog_.begin(); }
	const_iterator end()   const	{ return g_sog_.end(); }	
	int support() const		{ return g_sog_.size(); }
	void push(SBG<Policy>* s)	{ g_sog_[s->get_graph()].push(s); }
	void merge(SubgraphsOfManyGraph& sg);
    };

    template<class Policy>
    void SubgraphsOfManyGraph<Policy>::merge(SubgraphsOfManyGraph& sg)
    {
	for (typename G_SOG::iterator i = sg.g_sog_.begin(); i != sg.g_sog_.end(); ++i)
	    g_sog_[i->first].merge(i->second);
	sg.g_sog_.clear();
    }


    //
    // Projected implementation part
    //

    template<class Policy>
    class ProjectedSimple
    {
	typedef std::list<SBG<Policy> > SBGS;
	SBGS sbgs_;
    public:
	typedef typename SBGS::const_iterator const_iterator;
	const_iterator begin() const	{ return sbgs_.begin(); }
	const_iterator end()   const	{ return sbgs_.end(); }
	SBG<Policy>& push(const SBG<Policy>& s) { sbgs_.push_back(s); return sbgs_.back(); }
	SBG<Policy>& back() { return sbgs_.back(); }
	int size() const { return sbgs_.size(); }
	void merge(ProjectedSimple& p)	{ sbgs_.splice(sbgs_.begin(), p.sbgs_); }
    };


    template<class Policy, template<class> class S>
    class Projected : public ProjectedSimple<Policy>
    {
	S<Policy> sg_;
	typedef ProjectedSimple<Policy> Base;
	Base& base()			{ return *this; }
	const Base& base() const	{ return *this; }
    public:
	int support() const		{ return sg_.support(); }
	int size() const		{ return base().size(); }
	void push(const SBG<Policy>& s) { sg_.push(&base().push(s)); }
	void merge(Projected& p)	{ base().merge(p); sg_.merge(p.sg_); }
	const S<Policy>& sg() const { return sg_; }
    };

    template<class Policy>
    std::ostream& operator<<(std::ostream& out, const ProjectedSimple<Policy>& p)
    {
        for (typename ProjectedSimple<Policy>::const_iterator i = p.begin(); i != p.end(); ++i)
            out << '\t' << *i << std::endl;
        return out;
    }

    // *****************************************************************************
    //                          EdgeIterator
    // *****************************************************************************
    
    template<class Policy>
    class GetEdge
    {
	typedef typename Policy::graph_t G;
	const Policy& pl_;
	const G& g_;
    public:
	GetEdge(const Policy& pl, const G& g) :pl_(pl), g_(g) {}

	Edge<Policy> operator() (const typename boost::graph_traits<G>::out_edge_iterator& it) const
	    {
		Edge<Policy> e;
		e.vi_from = source(*it, g_);
		e.vi_to   = target(*it, g_);
		e.ei      = pl_.get_edge_index(*it, g_);
		return e;
	    }

	Edge<Policy> operator() (const typename boost::graph_traits<G>::in_edge_iterator& it) const
	    {
		Edge<Policy> e;
		e.vi_from = target(*it, g_);
		e.vi_to   = source(*it, g_);
		e.ei      = pl_.get_edge_index(*it, g_);
		return e;
	    }
    };


    template<class Policy>
    class EdgeIter : public EdgeIterator<typename Policy::graph_t,
					 Edge<Policy>,
					 GetEdge<Policy> >
    {
	typedef EdgeIterator<typename Policy::graph_t, Edge<Policy>, GetEdge<Policy> > Base;
    public:
	EdgeIter(typename Policy::vertex_index_t vi,
		 const typename Policy::graph_t& g,
		 const Policy& pl)
	    :Base(pl.get_vertex_descriptor(vi,g), g, GetEdge<Policy>(pl, g))
	    {}
    };


    template<class Policy>
    struct EdgeInSbgPred
    {
	const Policy& pl_;
	const SBG<Policy>& sbg_;
    public:
	EdgeInSbgPred(const Policy& pl, const SBG<Policy>& sbg) :pl_(pl), sbg_(sbg) {}
	bool operator() (const Edge<Policy>& e) const { return ! sbg_.has_edge(e.ei); }
    };


    template<class Policy>
    class BackwardEdgePred
    {
	const SBG<Policy>& sbg_;
	const Edge<Policy>& e1_;
	const Edge<Policy>& e2_;
	const Policy& pl_;
    public:
	BackwardEdgePred(const SBG<Policy>& sbg, const Policy& pl,
			 const Edge<Policy>& e1,
			 const Edge<Policy>& e2)
	    :sbg_(sbg), e1_(e1), e2_(e2), pl_(pl) {}
	bool operator() (const Edge<Policy>& e) const
	    {
		const typename Policy::graph_t& g = *sbg_.get_graph();
		typename Policy::edge_label_ref_t e1label = pl_.eilabel(e1_.ei, g);
		typename Policy::edge_label_ref_t elabel = pl_.eilabel(e.ei,g);
		typename Policy::vertex_label_ref_t vl1_to = pl_.vilabel(e1_.vi_to, g);
		typename Policy::vertex_label_ref_t vl2_to = pl_.vilabel(e2_.vi_to, g);
		return
		    e.vi_to == e1_.vi_from && (
			pl_.elabel_less(e1label,elabel) ||
			(pl_.elabel_equal(e1label,elabel) && pl_.vlabel_less_or_equal(vl1_to, vl2_to)));
	    }
    };


    template<class Policy>
    class ForwardPureEdgePred
    {
	const SBG<Policy>& sbg_;
	const Policy& pl_;
	typename Policy::vertex_label_ref_t minlabel_;
    public:
	ForwardPureEdgePred(const SBG<Policy>& sbg,
			    const Policy& pl,
			    typename Policy::vertex_label_ref_t minlabel)
	    :sbg_(sbg), pl_(pl), minlabel_(minlabel) {}
	bool operator() (const Edge<Policy>& e) const
	    {
		return
		    pl_.vlabel_less_or_equal(minlabel_, pl_.vilabel(e.vi_to, *sbg_.get_graph())) &&
		    !sbg_.has_vertex(e.vi_to); // skip backward edge
	    }
    };


    template<class Policy>
    class ForwardRmpathEdgePred
    {
	const SBG<Policy>& sbg_;
	const Policy& pl_;
	typename Policy::vertex_label_ref_t minlabel_;
	Edge<Policy> edge_;
    public:
	ForwardRmpathEdgePred(const SBG<Policy>& sbg,
			      const Policy& pl,
			      typename Policy::vertex_label_ref_t minlabel,
			      const Edge<Policy> edge)
	    :sbg_(sbg), pl_(pl), minlabel_(minlabel), edge_(edge) {}
	bool operator() (const Edge<Policy>& e) const;
    };

    template<class Policy>
    bool ForwardRmpathEdgePred<Policy>::operator() (const Edge<Policy>& e) const
    {
	const typename Policy::graph_t& g = *sbg_.get_graph();

	typename Policy::vertex_label_ref_t tolabel  = pl_.vilabel(edge_.vi_to, g);
	typename Policy::vertex_label_ref_t tolabel2 = pl_.vilabel(e.vi_to, g);

	if (sbg_.has_vertex(e.vi_to) || pl_.vlabel_less(tolabel2, minlabel_))
	    return false;
	
	typename Policy::edge_label_ref_t elabel  = pl_.eilabel(edge_.ei, g);
	typename Policy::edge_label_ref_t elabel2 = pl_.eilabel(e.ei, g);
	
	return
	    pl_.elabel_less(elabel, elabel2) ||
	    (pl_.elabel_equal(elabel, elabel2) && pl_.vlabel_less_or_equal(tolabel, tolabel2) );
    }


    struct DoNothing
    {
	template<class T>
	void operator() (const T&) {}
    };

    // -----------------------------------------------------
    //			EdgeExtensionIterator
    //		   filter out edges that exists in SBG
    // -----------------------------------------------------
    template<class Policy>
    class EdgeExtensionIterator
	: public EdgeConditionalIterator<EdgeIter<Policy>,
					 EdgeInSbgPred<Policy>,
					 DoNothing>
    {
	typedef EdgeIter<Policy> Iter;
	typedef EdgeConditionalIterator<Iter, EdgeInSbgPred<Policy>, DoNothing> Base;
    public:
	typedef Edge<Policy> value_type;
	EdgeExtensionIterator(typename Policy::vertex_index_t vi,
			      const SBG<Policy>& sbg,
			      const Policy& pl)
	    :Base(Iter(vi, *sbg.get_graph(), pl),
		  EdgeInSbgPred<Policy>(pl, sbg),
		  DoNothing()) {}
    };

    
    // -----------------------------------------------------
    //			EdgeExtensionPredIterator
    // -----------------------------------------------------  
    template<class Policy, class Pred, class OnFalse>
    class EdgeExtensionPredIterator
	: public EdgeConditionalIterator<EdgeExtensionIterator<Policy>, Pred, OnFalse>
    {
	typedef EdgeExtensionIterator<Policy> Iter;
	typedef EdgeConditionalIterator<Iter, Pred, OnFalse> Base;
    public:
	EdgeExtensionPredIterator(typename Policy::vertex_index_t vi,
			      const SBG<Policy>& sbg,
			      const Policy& pl,
			      const Pred& pred,
			      OnFalse& onfalse)
	    :Base(Iter(vi, sbg, pl), pred, onfalse) {}
    };

    // -----------------------------------------------------
    //			BckEdgeIterator
    // -----------------------------------------------------
    template<class Policy, class OnFalse = DoNothing>
    class BckEdgeIterator
	: public EdgeExtensionPredIterator<Policy, BackwardEdgePred<Policy>, OnFalse>
    {
	typedef BackwardEdgePred<Policy> Pred;
	typedef EdgeExtensionPredIterator<Policy, Pred, OnFalse> Base;
    public:
	BckEdgeIterator(const Edge<Policy>& e1,
			const Edge<Policy>& e2,
			const SBG<Policy>& sbg,
			const Policy& pl,
			OnFalse onfalse = DoNothing())
	    :Base(e2.vi_to, sbg, pl,
		  Pred(sbg, pl, e1, e2),
		  onfalse) {}
    };

    // -----------------------------------------------------
    //			FwdPureEdgeIterator
    // -----------------------------------------------------
    template<class Policy, class OnFalse = DoNothing>
    class FwdPureEdgeIterator
	: public EdgeExtensionPredIterator<Policy, ForwardPureEdgePred<Policy>, OnFalse>
    {
	typedef ForwardPureEdgePred<Policy> Pred;
	typedef EdgeExtensionPredIterator<Policy, Pred, OnFalse> Base;
    public:
	FwdPureEdgeIterator(const Edge<Policy>& e,
			    typename Policy::vertex_label_ref_t minlabel,
			    const SBG<Policy>& sbg,
			    const Policy& pl,
			    OnFalse onfalse = DoNothing())
	    :Base(e.vi_to, sbg, pl,
		  Pred(sbg, pl, minlabel),
		  onfalse) {}
    };
    

    // -----------------------------------------------------
    //			FwdRmpathEdgeIterator
    // -----------------------------------------------------
    template<class Policy, class OnFalse = DoNothing>
    class FwdRmpathEdgeIterator
	: public EdgeExtensionPredIterator<Policy, ForwardRmpathEdgePred<Policy>, OnFalse>
    {
	typedef ForwardRmpathEdgePred<Policy> Pred;
	typedef EdgeExtensionPredIterator<Policy, Pred, OnFalse> Base;
    public:
	FwdRmpathEdgeIterator(const Edge<Policy>& e,
			      typename Policy::vertex_label_ref_t minlabel,
			      const SBG<Policy>& sbg,
			      const Policy& pl,
			      OnFalse onfalse = DoNothing())
	    :Base(e.vi_from, sbg, pl,
		  Pred(sbg, pl, minlabel, e),
		  onfalse) {}
    };


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


    template<class M4>
    inline typename M4::mapped_type::mapped_type::mapped_type::mapped_type&
    mapped_val(M4& mp4key,
	       const typename M4::key_type& k1,
	       const typename M4::mapped_type::key_type& k2,
	       const typename M4::mapped_type::mapped_type::key_type& k3,
	       const typename M4::mapped_type::mapped_type::mapped_type::key_type& k4,
	       const typename M4::mapped_type::key_compare& c2,
	       const typename M4::mapped_type::mapped_type::key_compare& c3,
	       const typename M4::mapped_type::mapped_type::mapped_type::key_compare& c4)
    {
	typedef typename M4::value_type value_type;
	return mapped_val(mp4key.insert(value_type(k1, typename M4::mapped_type(c2))).first->second,
			  k2, k3, k4, c3, c4);
    }

    template<class M5>
    inline typename M5::mapped_type::mapped_type::mapped_type::mapped_type::mapped_type&
    mapped_val(M5& mp5key,
	       const typename M5::key_type& k1,
	       const typename M5::mapped_type::key_type& k2,
	       const typename M5::mapped_type::mapped_type::key_type& k3,
	       const typename M5::mapped_type::mapped_type::mapped_type::key_type& k4,
	       const typename M5::mapped_type::mapped_type::mapped_type::mapped_type::key_type& k5,
	       const typename M5::mapped_type::key_compare& c2,
	       const typename M5::mapped_type::mapped_type::key_compare& c3,
	       const typename M5::mapped_type::mapped_type::mapped_type::key_compare& c4,
	       const typename M5::mapped_type::mapped_type::mapped_type::mapped_type::key_compare& c5)
    {
	typedef typename M5::value_type value_type;
	return mapped_val(mp5key.insert(value_type(k1, typename M5::mapped_type(c2))).first->second,
			  k2, k3, k4, k5, c3, c4, c5);
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

    template<class Policy, class P>
    class MapTraits
    {
	// std::map wrapper without default constructor
	template<class K, class C, class T>
	struct map_ : public std::map<K,T,C> { explicit map_(const C& c) :std::map<K,T,C>(c) {} };
	typedef typename Policy::vertex_index_t VI;
	typedef typename Policy::vertex_label_t VL;
	typedef typename Policy::edge_label_t EL;
    public:
	typedef map_<EL, Elabel_less<Policy>, P>		Map_EL_P;
	typedef map_<VL, Vlabel_less<Policy>, P>		Map_VL_P;
	typedef map_<EL, Elabel_less<Policy>, Map_VL_P>		Map_EL_VL_P;
	typedef map_<VI, std::less<VI>,       Map_EL_P>		Map_VI_EL_P;
	typedef map_<VI, std::less<VI>,       Map_EL_VL_P>	Map_VI_EL_VL_P;
	typedef map_<VL, Vlabel_less<Policy>, Map_EL_VL_P>	Map_VL_EL_VL_P;

	typedef map_<VI, std::less<VI>,       Map_VL_EL_VL_P>		Map_VI_VL_EL_VL_P;
	typedef map_<VI, std::less<VI>,       Map_VI_VL_EL_VL_P>	Map_VI_VI_VL_EL_VL_P;
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
    

    template<class Policy, class Map_VL_EL_VL_P, class VlMapCount>
    void enum_one_edges(Map_VL_EL_VL_P& m,
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
	    typedef EdgeIterator<typename Policy::graph_t,
				 Edge<Policy>,
				 GetEdge<Policy> > EIter;

	    for (EIter iter(*v_it, g, GetEdge<Policy>(pl,g));
		 !iter.is_end(); iter.increment())
	    {
		Edge<Policy> e = iter.edge();
		typename Policy::vertex_label_ref_t vl_from = pl.vilabel(e.vi_from, g);
		typename Policy::vertex_label_ref_t vl_to   = pl.vilabel(e.vi_to, g);
		typename Policy::edge_label_ref_t el        = pl.eilabel(e.ei, g);

		assert(! pl.void_vlabel(vl_from));
		assert(! pl.void_vlabel(vl_to));

		mapped_val(m, vl_from, el, vl_to, el_less, vl_less).push(SBG<Policy>(&g, e));
	    }
	}
    }


    template<class Policy, class P>
    bool project_is_min(const P& projected,
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
	typedef typename MapTraits<Policy, P>::Map_EL_P    BEdges;
	typedef typename MapTraits<Policy, P>::Map_EL_VL_P FEdges;

	RMPath rmpath(dfsc_min);
	VI maxtoc = dfsc_min[rmpath.rightmost()].vi_to;

	// backward
	{
	    Elabel_less<Policy> el_less(pl);
	    BEdges b_edges(el_less);
	    VI newto = pl.void_vertex_index();
	    bool flg = false;
	    for (int i = 0; !flg && i < rmpath.size() - 1; ++i)
	    {
		BOOST_FOREACH(const SBG& sbg, projected)
		{
		    const typename Policy::graph_t& g = *sbg.get_graph();
		    BckEdgeIterator<Policy> iter(sbg[rmpath[i]],
						 sbg[rmpath.rightmost()], sbg, pl);
		    if (!iter.is_end())
		    {
			const Edge& e = iter.edge();
			ELR el = pl.eilabel(e.ei,g);
			mapped_val(b_edges, el).push(SBG(&sbg, e));
			newto = dfsc_min[rmpath[i]].vi_from;
			flg = true;
		    }
		}
	    }

	    if (flg)
	    {
		typename BEdges::const_iterator i1 = b_edges.begin();
		dfsc_min.push_back(EdgeCode(maxtoc, newto, pl.void_vlabel(), i1->first, pl.void_vlabel()));
		return project_is_min(i1->second, dfsc_min, dfsc_tested, pl);
	    }
	}

	// forward
	{
	    Elabel_less<Policy> el_less(pl);
	    Vlabel_less<Policy> vl_less(pl);
	    FEdges f_edges(el_less);
	    VLR minlabel = dfsc_min[0].vl_from;
	    VI newfrom = pl.void_vertex_index();
	    bool flg = false;
            
	    // forward pure
	    BOOST_FOREACH(const SBG& sbg, projected)
	    {
		const typename Policy::graph_t& g = *sbg.get_graph();

		for (FwdPureEdgeIterator<Policy> iter(sbg[rmpath.rightmost()], minlabel, sbg, pl);
		     !iter.is_end(); iter.increment())
		{
		    const Edge& e = iter.edge();
		    ELR el = pl.eilabel(e.ei,g);
		    VLR vl = pl.vilabel(e.vi_to,g);
		    mapped_val(f_edges, el, vl, vl_less).push(SBG(&sbg, e));
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
		    for (FwdRmpathEdgeIterator<Policy> iter(sbg[rmpath[i]], minlabel, sbg, pl);
			 !iter.is_end(); iter.increment())
		    {
			const Edge& e = iter.edge();
			ELR el = pl.eilabel(e.ei,g);
			VLR vl = pl.vilabel(e.vi_to,g);
			mapped_val(f_edges, el, vl, vl_less).push(SBG(&sbg, e));
			newfrom = dfsc_min[rmpath[i]].vi_from;
			flg = true;
		    }
		}
	    }

	    if (flg)
	    {
		typename FEdges::const_iterator i1 = f_edges.begin();
		typename FEdges::mapped_type::const_iterator i2 = i1->second.begin();
		dfsc_min.push_back(EdgeCode(newfrom, maxtoc+1, pl.void_vlabel(), i1->first, i2->first));
		return project_is_min(i2->second, dfsc_min, dfsc_tested, pl);
	    }
	}

	return true;
    }


    template<class Policy>
    bool is_min(const DFSCode<Policy>& dfsc_tested, const Policy& pl)
    {
	typedef ProjectedSimple<Policy> P;
	typedef typename MapTraits<Policy, P>::Map_VL_EL_VL_P M3;
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


    // ================================================================
    //                   gspan_detail
    // ================================================================
    namespace gspan_detail
    {
	template<class Policy, class P, class Output>
	void project(const P& projected,
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

	    int sup = projected.support();
	    if (sup < minsup)
		return;

	    if (! is_min(dfsc, pl))
		return;

	    result(dfsc, projected.sg());

	    // --------------------------------------------------------------
	    // enumerate
	    typedef typename MapTraits<Policy, P>::Map_VI_EL_P    BEdges;
	    typedef typename MapTraits<Policy, P>::Map_VI_EL_VL_P FEdges;

	    std::less<VI> vi_less;
	    Elabel_less<Policy> el_less(pl);
	    Vlabel_less<Policy> vl_less(pl);

	    BEdges b_edges(vi_less);
	    FEdges f_edges(vi_less);

	    RMPath rmpath(dfsc);

	    VI maxtoc   = dfsc[rmpath.rightmost()].vi_to;
	    VLR minlabel = dfsc[0].vl_from;

	    BOOST_FOREACH(const SBG& sbg, projected)
	    {
		const typename Policy::graph_t& g = *sbg.get_graph();

		// backward
		for (int i = 0; i < rmpath.size() - 1; ++i)
		{
		    BckEdgeIterator<Policy> iter(sbg[rmpath[i]], sbg[rmpath.rightmost()], sbg, pl);
		    if (!iter.is_end())
		    {
			const Edge& e = iter.edge();
			VI vi  = dfsc[rmpath[i]].vi_from;
			ELR el = pl.eilabel(e.ei,g);
			mapped_val(b_edges, vi, el, el_less).push(SBG(&sbg, e));
		    }
		}

		// forward
		for (FwdPureEdgeIterator<Policy> iter(sbg[rmpath.rightmost()], minlabel, sbg, pl);
		     !iter.is_end(); iter.increment())
		{
		    const Edge& e = iter.edge();
		    ELR el = pl.eilabel(e.ei,g);
		    VLR vl = pl.vilabel(e.vi_to,g);
		    mapped_val(f_edges, maxtoc, el, vl, el_less, vl_less).push(SBG(&sbg, e));
		}

		for (int i = rmpath.size()-1; i >= 0; --i)
		    for (FwdRmpathEdgeIterator<Policy> iter(sbg[rmpath[i]], minlabel, sbg, pl);
			 !iter.is_end(); iter.increment())
		    {
			const Edge& e = iter.edge();
			VI vi  = dfsc[rmpath[i]].vi_from;
			ELR el = pl.eilabel(e.ei,g);
			VLR vl = pl.vilabel(e.vi_to,g);
			mapped_val(f_edges, vi, el, vl, el_less, vl_less).push(SBG(&sbg, e));
		    }
	    }

	    // --------------------------------------------------------------
	    // recursive process SBG children
	    // backward
	    for (typename BEdges::const_iterator it1 = b_edges.begin(); it1 != b_edges.end(); ++it1)
		for (typename BEdges::mapped_type::const_iterator it2 = it1->second.begin();
		     it2 != it1->second.end(); ++it2)
		{
		    EdgeCode ec(maxtoc, it1->first, pl.void_vlabel(), it2->first, pl.void_vlabel());
		    dfsc.push_back(ec);
		    project(it2->second, dfsc, minsup, pl, result);
		    dfsc.pop_back();
		}

	    // forward
	    typedef typename FEdges::const_reverse_iterator FwdI1;
	    typedef typename FEdges::mapped_type::const_iterator FwdI2;
	    typedef typename FEdges::mapped_type::mapped_type::const_iterator FwdI3;
	    for (FwdI1 i1 = f_edges.rbegin(); i1 != f_edges.rend(); ++i1)
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

    // ================================================================
    //                   closegraph_detail
    // ================================================================
    namespace closegraph_detail
    {
	template<class Policy>
	typename Policy::vertex_index_t dfsc_vindex(typename Policy::vertex_index_t sbg_vi,
						    const DFSCode<Policy>& dfsc,
						    const SBG<Policy>& s,
						    const Policy& pl)
	{
	    for (int i = s.size()-1; i >= 0; --i)
	    {
		if (sbg_vi == s[i].vi_from)
		    return dfsc[i].vi_from;
		if (sbg_vi == s[i].vi_to)
		    return dfsc[i].vi_to;
	    }
	    return pl.void_vertex_index();
	}


	template<class Policy, class Output, class P>
	class Closegraph_alg
	{
	    typedef typename Policy::graph_t Graph;
	    typedef typename Policy::vertex_index_t VI;
	    typedef typename Policy::edge_index_t EI;
	    typedef typename Policy::vertex_label_ref_t VLR;
	    typedef typename Policy::edge_label_ref_t ELR;

	    typedef typename MapTraits<Policy, P>::Map_VI_EL_P    BEdges;
	    typedef typename MapTraits<Policy, P>::Map_VI_EL_VL_P FEdges;
	    typedef typename MapTraits<Policy, P>::Map_VI_VI_VL_EL_VL_P XEdges;

	    const Policy& pl_;
	    const int minsup_;
	    Output& result_;

	    DFSCode<Policy> dfsc_;
	    std::map<typename Policy::vertex_label_t, int> vlab_count_;
	    std::map<typename Policy::vertex_label_t, bool> early_termin_1v_;
	    std::vector<bool> early_termin_;

	    const Vlabel_less<Policy> vl_less_;
	    const Elabel_less<Policy> el_less_;

	    void run(typename MapTraits<Policy, P>::Map_VL_EL_VL_P& root);
	    void find_xedges_vito(const XEdges& x_edges, VI vi, VLR vl, ELR el, std::set<VI>& vii) const;
	    void fail_et(const EdgeCode<Policy>& ec, bool rightmost,
			 const BEdges& b_edges, const XEdges& x_edges,
			 std::set<int>& vidx) const;
	    void project(const P* projected, const P* prev_projected);
	public:
	    template<class TGraphIterator>
	    Closegraph_alg(TGraphIterator tg_begin,
			   TGraphIterator tg_end,
			   int minsup,
			   const Policy& pl,
			   Output& result);

	    Closegraph_alg(const Graph& graph, int minsup, const Policy& pl, Output& result);
	};

	template<class Policy, class Output, class P>
	template<class TGraphIterator>
	Closegraph_alg<Policy,Output,P>::Closegraph_alg(TGraphIterator tg_begin,
							TGraphIterator tg_end,
							int minsup,
							const Policy& pl,
							Output& result)
	    :pl_(pl), minsup_(minsup), result_(result), vl_less_(pl), el_less_(pl)
	{
	    MapValueCount<std::map<typename Policy::vertex_label_t, int> > mvc(vlab_count_);
	    typename MapTraits<Policy, P>::Map_VL_EL_VL_P root(vl_less_);
	    for (; tg_begin != tg_end; ++tg_begin)
		enum_one_edges(root, *tg_begin, pl, mvc);
	    run(root);
	}

	template<class Policy, class Output, class P>
	Closegraph_alg<Policy,Output,P>::Closegraph_alg(const Graph& graph,
							int minsup,
							const Policy& pl,
							Output& result)
	    :pl_(pl), minsup_(minsup), result_(result), vl_less_(pl), el_less_(pl)
	{
	    MapValueCount<std::map<typename Policy::vertex_label_t, int> > mvc(vlab_count_);
	    typename MapTraits<Policy, P>::Map_VL_EL_VL_P root(vl_less_);
	    enum_one_edges(root, graph, pl, mvc);
	    run(root);
	}

	template<class Policy, class Output, class P>
	void Closegraph_alg<Policy,Output,P>::run(typename MapTraits<Policy, P>::Map_VL_EL_VL_P& root)
	{
	    typedef typename MapTraits<Policy, P>::Map_VL_EL_VL_P M3;
	    typedef typename M3::const_iterator I1;
	    typedef typename M3::mapped_type::const_iterator I2;
	    typedef typename M3::mapped_type::mapped_type::const_iterator I3;
	    for (I1 i1 = root.begin(); i1 != root.end(); ++i1)
	    {
		bool break_loop = false;
		if (early_termin_1v_[i1->first])
		{
#ifdef DEBUG_PRINT
		    std::cerr << "Early Termination for edgecode beginning with " << i1->first << std::endl;
#endif
		    break_loop = true;
		}

		for (I2 i2 = i1->second.begin(); i2 != i1->second.end() && !break_loop; ++i2)
		    for (I3 i3 = i2->second.begin(); i3 != i2->second.end() && !break_loop; ++i3)
		    {
			if (early_termin_1v_[i1->first]) // goto to first loop
			{
			    break_loop = true;
			    break;
			}

			EdgeCode<Policy> ec(0, 1, i1->first, i2->first, i3->first);
			dfsc_.push_back(ec);

#ifdef DEBUG_PRINT
			std::cerr << "================ "
				  << "TOP LEVEL iteration with EDGECODE: " << ec
				  << " ================" << std::endl;
#endif
			const P& projected = i3->second;
			if (projected.support() >=  minsup_)
			{
			    if (projected.size() == vlab_count_[ec.vl_from])
				early_termin_1v_[ec.vl_from] = true;
			    if (projected.size() == vlab_count_[ec.vl_to])
				early_termin_1v_[ec.vl_to] = true;
			    
			    project(&projected, 0);
			}

			dfsc_.pop_back();
		    }
	    }
	}


	template<class Policy, class Output, class P>
	void Closegraph_alg<Policy,Output,P>::find_xedges_vito(const XEdges& x_edges,
							       VI vi, VLR vl, ELR el,
							       std::set<VI>& vii) const
	{

	    typedef typename XEdges::const_iterator XI1;
	    typedef typename XEdges::mapped_type::const_iterator XI2;
	    typedef typename XEdges::mapped_type::mapped_type::const_iterator XI3;
	    typedef typename XEdges::mapped_type::mapped_type::mapped_type::const_iterator XI4;
	    typedef typename XEdges::mapped_type::mapped_type::mapped_type::mapped_type::const_iterator XI5;
	    
	    const VI VI_VOID = pl_.void_vertex_index();

	    for (XI1 i1 = x_edges.begin(); i1 != x_edges.end(); ++i1)
		if (i1->first != VI_VOID && i1->first == vi)
		    for (XI2 i2 = i1->second.begin(); i2 != i1->second.end(); ++i2)
			if (i2->first != VI_VOID)
			    for (XI3 i3 = i2->second.begin(); i3 != i2->second.end(); ++i3)
				for (XI4 i4 = i3->second.begin(); i4 != i3->second.end(); ++i4)
				    if (i4->first == el)
					for (XI5 i5 = i4->second.begin(); i5 != i4->second.end(); ++i5)
					    if (i5->first == vl)
						vii.insert(i2->first);

	    for (XI1 i1 = x_edges.begin(); i1 != x_edges.end(); ++i1)
		if (i1->first != VI_VOID)
		    for (XI2 i2 = i1->second.begin(); i2 != i1->second.end(); ++i2)
			if (i2->first != VI_VOID && i2->first == vi)
			    for (XI3 i3 = i2->second.begin(); i3 != i2->second.end(); ++i3)
				if (i3->first == vl)
				    for (XI4 i4 = i3->second.begin(); i4 != i3->second.end(); ++i4)
					if (i4->first == el)
					    for (XI5 i5 = i4->second.begin(); i5 != i4->second.end(); ++i5)
						vii.insert(i1->first);
	}
 

	template<class Policy, class Output, class P>
	void Closegraph_alg<Policy,Output,P>::fail_et(const EdgeCode<Policy>& ec, bool rightmost,
						      const BEdges& b_edges, const XEdges& x_edges,
						      std::set<int>& vidx) const
	{
	    std::set<VI> vii_1;

	    // find Edge in BEdges
	    if (rightmost)
	    {
		typedef typename BEdges::const_iterator BckCI1;
		typedef typename BEdges::mapped_type::const_iterator BckCI2;
		for (BckCI1 i1 = b_edges.begin(); i1 != b_edges.end(); ++i1)
		    if (dfsc_.vlabel(i1->first) == ec.vl_to)
			for (BckCI2 i2 = i1->second.begin(); i2 != i1->second.end(); ++i2)
			    if (i2->first == ec.el)
				vii_1.insert(i1->first);
	    }
	    
	    // find Edge in XEdges
	    find_xedges_vito(x_edges, ec.vi_from, ec.vl_to, ec.el, vii_1);


	    std::set<VI> vii_2;
	    for (typename std::set<VI>::const_iterator i = vii_1.begin(); i != vii_1.end(); ++i)
		for (typename DFSCode<Policy>::const_iterator d = dfsc_.begin(); d != dfsc_.end(); ++d)
		{
		    if (d->vi_from == *i)
			vii_2.insert(d->vi_to);
		    if (d->vi_to == *i)
			vii_2.insert(d->vi_from);
		}


	    for (typename std::set<VI>::const_iterator i = vii_2.begin(); i != vii_2.end(); ++i)
		for (unsigned int idx = 0; idx != dfsc_.size(); ++idx)
		{
		    if (dfsc_[idx].vi_from == *i)
			vidx.insert(idx);
		    if (dfsc_[idx].vi_to == *i)
			vidx.insert(idx);
		}
	}


	template<class Policy, class Output, class P>
	void Closegraph_alg<Policy,Output,P>::project(const P* projected, const P* prev_projected)
	{
#ifdef DEBUG_PRINT
	    std::cerr << "=========== project() ===========\n"
		      << "DFSC: " << dfsc_ << std::endl
		      << "P supp=" << projected->support() << "\n"<< *projected
		      << std::endl;
#endif
	    early_termin_.push_back(false);
	    assert(early_termin_.size() == dfsc_.size());

	    typedef SBG<Policy> SBG;
	    typedef Edge<Policy> Edge;
	    typedef EdgeCode<Policy> EdgeCode;

	    // --------------------------------------------------------------
	    // enumerate
	    
	    std::less<VI> vi_less;

	    BEdges b_edges(vi_less);
	    FEdges f_edges(vi_less);
	    XEdges x_edges(vi_less);

	    RMPath rmpath(dfsc_);
	    VI maxtoc   = dfsc_[rmpath.rightmost()].vi_to;
	    VLR minlabel = dfsc_[0].vl_from;
	    const int NUM_EDGES = dfsc_.size(); 

	    BOOST_FOREACH(const SBG& sbg, *projected)
	    {
		assert(NUM_EDGES == sbg.size());
		const typename Policy::graph_t& g = *sbg.get_graph();
		
		// -----------------------------------
		// discover RMPath extension edges
		// -----------------------------------
		std::set<EI> r_extension;

		// backward
		for (int i = 0; i < rmpath.size() - 1; ++i)
		{
		    BckEdgeIterator<Policy> iter(sbg[rmpath[i]], sbg[rmpath.rightmost()], sbg, pl_);
		    if (!iter.is_end())
		    {
			const Edge& e = iter.edge();
			VI vi  = dfsc_[rmpath[i]].vi_from;
			ELR el = pl_.eilabel(e.ei,g);
			mapped_val(b_edges, vi, el, el_less_).push(SBG(&sbg, e));
			r_extension.insert(e.ei);
		    }
		}

		// forward
		for (FwdPureEdgeIterator<Policy> iter(sbg[rmpath.rightmost()], minlabel, sbg, pl_);
		     !iter.is_end(); iter.increment())
		{
		    const Edge& e = iter.edge();
		    ELR el = pl_.eilabel(e.ei,g);
		    VLR vl = pl_.vilabel(e.vi_to,g);
		    mapped_val(f_edges, maxtoc, el, vl, el_less_, vl_less_).push(SBG(&sbg, e));
		    r_extension.insert(e.ei);
		}

		for (int i = rmpath.size()-1; i >= 0; --i)
		    for (FwdRmpathEdgeIterator<Policy> iter(sbg[rmpath[i]], minlabel, sbg, pl_);
			 !iter.is_end(); iter.increment())
		    {
			const Edge& e = iter.edge();
			VI vi  = dfsc_[rmpath[i]].vi_from;
			ELR el = pl_.eilabel(e.ei,g);
			VLR vl = pl_.vilabel(e.vi_to,g);
			mapped_val(f_edges, vi, el, vl, el_less_, vl_less_).push(SBG(&sbg, e));
			r_extension.insert(e.ei);
		    }

		// -----------------------------------
		// discover NOT RMPath extension edges
		// -----------------------------------
		std::set<EI> x_extension;
		for (int i = 0; i < NUM_EDGES; ++i)
		{
		    for (EdgeExtensionIterator<Policy> iter(sbg[i].vi_from, sbg, pl_);
			 !iter.is_end(); iter.increment())
		    {
			const Edge& e = iter.edge();
			if (r_extension.count(e.ei) == 0 && x_extension.count(e.ei) == 0)
			{
			    typename Policy::vertex_label_ref_t vl_from = pl_.vilabel(e.vi_from, g);
			    typename Policy::vertex_label_ref_t vl_to   = pl_.vilabel(e.vi_to, g);
			    typename Policy::edge_label_ref_t el        = pl_.eilabel(e.ei, g);
			    mapped_val(x_edges,
				       dfsc_vindex(e.vi_from, dfsc_, sbg, pl_), dfsc_vindex(e.vi_to, dfsc_, sbg, pl_),
				       vl_from, el, vl_to,
				       vi_less, vl_less_, el_less_, vl_less_).push(SBG(&sbg, e));
			    x_extension.insert(e.ei);
			}
		    }
		    
		    for (EdgeExtensionIterator<Policy> iter(sbg[i].vi_to, sbg, pl_); !iter.is_end(); iter.increment())
		    {
			const Edge& e = iter.edge();
			if (r_extension.count(e.ei) == 0 && x_extension.count(e.ei) == 0)
			{
			    typename Policy::vertex_label_ref_t vl_from = pl_.vilabel(e.vi_from, g);
			    typename Policy::vertex_label_ref_t vl_to   = pl_.vilabel(e.vi_to, g);
			    typename Policy::edge_label_ref_t el        = pl_.eilabel(e.ei, g);
			    mapped_val(x_edges,
				       dfsc_vindex(e.vi_from, dfsc_, sbg, pl_), dfsc_vindex(e.vi_to, dfsc_, sbg, pl_),
				       vl_from, el, vl_to,
				       vi_less, vl_less_, el_less_, vl_less_).push(SBG(&sbg, e));
			    x_extension.insert(e.ei);
			}
		    }
		}
	    }


	    // --------------------------------------------------------------
	    // recursive process SBG children
	    bool r_closed = true;

	    // backward
	    typedef typename BEdges::iterator BckI1;
	    typedef typename BEdges::mapped_type::iterator BckI2;
	    for (BckI1 i1 = b_edges.begin(); i1 != b_edges.end() && !early_termin_.back(); ++i1)
		for (BckI2 i2 = i1->second.begin(); i2 != i1->second.end() && !early_termin_.back(); ++i2)
		{
		    EdgeCode ec(maxtoc, i1->first, pl_.void_vlabel(), i2->first, pl_.void_vlabel());
		    dfsc_.push_back(ec);
		    P& new_prj = i2->second;

		    if (is_min(dfsc_, pl_))
		    {
			// detect early termination
			if (new_prj.size() == projected->size())
			    early_termin_.back() = true;

			if (new_prj.support() >= minsup_)
			{
			    if (new_prj.support() == projected->support())
				r_closed = false; // projected is not closed

			    project(&new_prj, projected);
			}
		    }
		    else
		    {
			mapped_val(x_edges, ec.vi_from, ec.vi_to, ec.vl_from, ec.el, ec.vl_to,
				   vi_less, vl_less_, el_less_, vl_less_).merge(new_prj);
		    }

		    dfsc_.pop_back();
		}

	    // forward
	    typedef typename FEdges::reverse_iterator FwdRevI1;
	    typedef typename FEdges::mapped_type::iterator FwdI2;
	    typedef typename FEdges::mapped_type::mapped_type::iterator FwdI3;
	    for (FwdRevI1 i1 = f_edges.rbegin(); i1 != f_edges.rend() && !early_termin_.back(); ++i1)
		for (FwdI2 i2 = i1->second.begin(); i2 != i1->second.end() && !early_termin_.back(); ++i2)
		    for (FwdI3 i3 = i2->second.begin(); i3 != i2->second.end() && !early_termin_.back(); ++i3)
		    {
			EdgeCode ec(i1->first, maxtoc+1, pl_.void_vlabel(), i2->first, i3->first);
			dfsc_.push_back(ec);
			P& new_prj = i3->second;

			if (is_min(dfsc_, pl_))
			{
			    // detect early termination
			    if (new_prj.size() == projected->size())
				early_termin_.back() = true;
			    
			    if (new_prj.support() >= minsup_)
			    {
				if (new_prj.support() == projected->support())
				    r_closed = false; // projected is not closed

				project(&new_prj, projected);
			    }
			    else
			    {
				std::set<int> idxs;
				fail_et(ec, ec.vi_from == maxtoc, b_edges, x_edges, idxs);
				for (std::set<int>::const_iterator ii = idxs.begin(); ii != idxs.end(); ++ii)
				    early_termin_[*ii] = false;
			    }
			}
			else
			{
			    mapped_val(x_edges, ec.vi_from, ec.vi_to, ec.vl_from, ec.el, ec.vl_to,
				       vi_less, vl_less_, el_less_, vl_less_).merge(new_prj);
			}

			dfsc_.pop_back();
		    }
	    

	    // if current dfsc may be extented by x_edges to any supergraph with the same support,
	    // then dfsc is not closed
	    bool x_closed = true;
	    typedef typename XEdges::const_iterator XI1;
	    typedef typename XEdges::mapped_type::const_iterator XI2;
	    typedef typename XEdges::mapped_type::mapped_type::const_iterator XI3;
	    typedef typename XEdges::mapped_type::mapped_type::mapped_type::const_iterator XI4;
	    typedef typename XEdges::mapped_type::mapped_type::mapped_type::mapped_type::const_iterator XI5;
	    for (XI1 i1 = x_edges.begin(); i1 != x_edges.end() && x_closed; ++i1)
		for (XI2 i2 = i1->second.begin(); i2 != i1->second.end() && x_closed; ++i2)
		    for (XI3 i3 = i2->second.begin(); i3 != i2->second.end() && x_closed; ++i3)
			for (XI4 i4 = i3->second.begin(); i4 != i3->second.end() && x_closed; ++i4)
			    for (XI5 i5 = i4->second.begin(); i5 != i4->second.end() && x_closed; ++i5)
			    {
				if (i5->second.support() == projected->support())
				{
				    x_closed = false;
				    break;
				}
			    }

	    if (r_closed && x_closed)
		result_(dfsc_, projected->sg());

	    early_termin_.pop_back();

#ifdef DEBUG_PRINT
	    std::cerr << "RET\n";
#endif
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

	typedef Projected<Policy, SubgraphsOfManyGraph> P;
	typedef typename MapTraits<Policy, P>::Map_VL_EL_VL_P M3;
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
		    std::cerr << "================ "
			      << "TOP LEVEL iteration with EDGECODE: " << ec
			      << " ================" << std::endl;
#endif
		    dfsc.push_back(ec);
		    gspan_detail::project(i3->second, dfsc, minsup, pl, result);
		    dfsc.pop_back();
		}
    }



    template<class Graph, class Output, class Policy>
    void gspan(const Graph& graph, int minsup, const Policy& pl, Output& result)
    {
	MapValueCountDefault mvc;

	typedef Projected<Policy, SubgraphsOfOneGraph> P;
	typedef typename MapTraits<Policy, P>::Map_VL_EL_VL_P M3;
	Vlabel_less<Policy> vl_less(pl);
	M3 root(vl_less);

	enum_one_edges(root, graph, pl, mvc);

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
		    std::cerr << "================ "
			      << "TOP LEVEL iteration with EDGECODE: " << ec
			      << " ================" << std::endl;
#endif
		    dfsc.push_back(ec);
		    gspan_detail::project(i3->second, dfsc, minsup, pl, result);
		    dfsc.pop_back();
		}
    }


    template<class Graph, class Output, class Policy>
    void closegraph(const Graph& graph, int minsup, const Policy& pl, Output& result)
    {
	using namespace closegraph_detail;
	Closegraph_alg<Policy, Output, Projected<Policy, SubgraphsOfOneGraph> >
	    (graph, minsup, pl, result);
    }


    template<class TGraphIterator, class Output, class Policy>
    void closegraph(TGraphIterator tg_begin,
		       TGraphIterator tg_end,
		       int minsup,
		       const Policy& pl,
		       Output& result)
    {
	using namespace closegraph_detail;
	Closegraph_alg<Policy, Output, Projected<Policy, SubgraphsOfManyGraph> >
	    (tg_begin, tg_end, minsup, pl, result);
    }

} // end: namespace gSpan

//#ifdef BR
//#undef BR
//#endif

#endif
