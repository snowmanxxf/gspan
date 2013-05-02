#ifndef GSPAN2_H_
#define GSPAN2_H_

#include <iostream>
#include <vector>
#include <boost/unordered_set.hpp>

#include "gspan_graph.hpp"
#include <boost/graph/adjacency_list.hpp>

#ifndef BR
#define BR asm volatile ("int3;")
#endif

#include <cstring>

namespace gSpan2
{

    class EdgeCode_
    {
    public:
	VI vi_src() const;
	VI vi_dst() const;
	VL vl_src() const;
	VL vl_dst() const;
	EL el() const;

	bool is_fwd() const;
	bool is_bck() const;
	void chgdir();
    };

    // *****************************************************************************
    //                          EdgeCode
    // *****************************************************************************
    class EdgeCode
    {
	VI vi_src_, vi_dst_;
	VL vl_src_, vl_dst_;
	EL el_;
	bool is_fwd_;
    public:
	EdgeCode()
	    :vi_src_(VI_NULL), vi_dst_(VI_NULL),
	     vl_src_(VL_NULL), vl_dst_(VL_NULL), el_(EL_NULL),
	     is_fwd_(false)
	    {}

	EdgeCode(VI vi_from, VI vi_to, VL vl_from, EL el, VL vl_to, bool fwd)
	    :vi_src_(vi_from), vi_dst_(vi_to),
	     vl_src_(vl_from), vl_dst_(vl_to), el_(el),
	     is_fwd_(fwd) {}

	VI vi_src() const	{ return vi_src_; }
	VI vi_dst() const	{ return vi_dst_; }
	VL vl_src() const	{ return vl_src_; }
	VL vl_dst() const	{ return vl_dst_; }
	EL el() const		{ return el_; }

	bool is_forward() const		{ return is_fwd_; }
	bool is_backward() const	{ return !is_fwd_; }
	void chgdir() { std::swap(vi_src_, vi_dst_); std::swap(vl_src_, vl_dst_); }
	EdgeCode operator- () const { EdgeCode ec(*this); ec.chgdir(); return ec; }
	bool operator!= (const EdgeCode& ec) const	{ return ! (*this == ec); }
	bool operator== (const EdgeCode& ec) const;
    };

    struct EdgeCodeCmpDfs
    {
	bool operator() (const EdgeCode& ec1, const EdgeCode& ec2) const;
    };

    struct EdgeCodeCmpLex
    {
	bool operator() (const EdgeCode& ec1, const EdgeCode& ec2) const;
    };


    std::ostream& operator<<(std::ostream& out, const EdgeCode& ec);
	
    // *****************************************************************************
    //                          DFSCode
    // *****************************************************************************
    typedef std::vector<EdgeCode> DFSCode;
    
    VI max_vertex(const DFSCode& dfsc);
    std::ostream& operator<<(std::ostream& out, const DFSCode& dfsc);

    // *****************************************************************************
    //                          SBGSimple
    // *****************************************************************************
    class SBGSimple
    {
	const SBGSimple* prev_;
	Graph::Edge edge_;
	
	// Embeddings
	mutable std::vector<const Graph::Edge*>* e_;
	mutable std::vector<short>* vv_;
	mutable std::vector<char>* ee_;

	int depth_;
	SBGSimple& operator= (const SBGSimple& s);
    public:
	
	SBGSimple(const Graph::Edge& edge)
	    :prev_(0),
	     edge_(edge),
	     e_(0),
	     vv_(0),
	     ee_(0),
	     depth_(1)
	    {
	    }

	SBGSimple(const SBGSimple* s, const Graph::Edge& e)
	    :prev_(s),
	     edge_(e),
	     e_(0),
	     vv_(0),
	     ee_(0),
	     depth_(s->depth_+1)
	    {
	    }

	SBGSimple(const SBGSimple& s)
	    :prev_(s.prev_),
	     edge_(s.edge_),
	     e_(0),
	     vv_(0),
	     ee_(0),
	     depth_(s.depth_)
	    {
	    }

	~SBGSimple()
	    {
		delete e_;
		delete vv_;
		delete ee_;
	    }

	int size() const { return depth_; }
	const Graph::Edge* operator[] (int i) const	{ return (*e_)[i]; }
	bool has_vertex(VI vi) const		        { return (*vv_)[vi]; }
	bool has_edge(EI ei) const       	        { return (*ee_)[ei]; }
	const Graph::Edge& edge() const { return edge_; }
	const SBGSimple* parent() const { return prev_; }
    };

    // *****************************************************************************
    //                          SBG
    // *****************************************************************************
    typedef std::vector<unsigned char> VecEE;

    class SBG
    {
	const SBG* prev_;
	Graph::Edge edge_;


	// Embeddings
	mutable std::vector<const Graph::Edge*>* e_;
	mutable std::vector<short>* vv_;
	mutable VecEE* ee_;
	mutable unsigned int sum_;
	mutable std::vector<VI>* s2g_vv_;
	
	void init_e_() const;
	void init_vv_() const;
	void init_ee_() const;
	void init_s2g_vv_() const;

	const Graph* graph_;
	int depth_;

	VI vi_src_dfsc_;
	VI vi_dst_dfsc_;

	SBG& operator= (const SBG& s);
    public:
	SBG(const Graph* g, const Graph::Edge& e)
	    :prev_(0),
	     edge_(e),
	     e_(0),
	     vv_(0),
	     ee_(0),
	     s2g_vv_(0),
	     graph_(g),
	     depth_(1),
	     vi_src_dfsc_(0),
	     vi_dst_dfsc_(1)
	    {
		automorph_next = automorph_prev = this;
	    }

	SBG(const SBG* s, const Graph::Edge& e, VI vi_src_dfsc, VI vi_dst_dfsc)
	    :prev_(s),
	     edge_(e),
	     e_(0),
	     vv_(0),
	     ee_(0),
	     s2g_vv_(0),
	     graph_(s->graph_),
	     depth_(s->depth_+1),
	     vi_src_dfsc_(vi_src_dfsc),
	     vi_dst_dfsc_(vi_dst_dfsc)
	    {
		automorph_next = automorph_prev = this;
	    }

	SBG(const SBG& s)
	    :prev_(s.prev_),
	     edge_(s.edge_),
	     e_(0),
	     vv_(0),
	     ee_(0),
	     s2g_vv_(0),
	     graph_(s.graph_),
	     depth_(s.depth_),
	     vi_src_dfsc_(s.vi_src_dfsc_),
	     vi_dst_dfsc_(s.vi_dst_dfsc_)
	    {
		automorph_next = automorph_prev = this;
	    }


	~SBG()
	    {
		delete e_;
		delete vv_;
		delete ee_;
		delete s2g_vv_;
	    }

	int size() const { return depth_; }

	void init_e() const		{ if (!e_)  init_e_(); }
	void init_vv() const		{ if (!vv_) init_vv_(); }
	void init_ee() const		{ if (!ee_) init_ee_(); }
	void init_s2g_vv() const	{ if (!s2g_vv_) init_s2g_vv_(); }

	const std::vector<const Graph::Edge*>& get_e() const	{ init_e(); return *e_; }
	const VecEE& get_ee() const		{ init_ee(); return *ee_; }
	const std::vector<short>& get_vv() const	{ init_vv(); return *vv_; }

	//
	// return array of the graph VI, indexed by the dfsc VI
	// so, sbg_vertices[dfsc_vi] == graph vi
	//
	const std::vector<VI> get_s2g_vv() const	{ init_s2g_vv(); return *s2g_vv_; }

	//
	// make array of the dfsc VI, indexed by the graph VI
	// so, sbg_vertices[graph_vi] == dfsc vi
	//
	void graph_to_dfsc_v(std::vector<VI>& vv, VI vi_default = VI_NULL) const;


	int num_vertices() const			{ return get_s2g_vv().size(); }
	const Graph::Edge* operator[] (int i) const	{ return (*e_)[i]; }
	bool has_vertex(VI vi) const		        { return (*vv_)[vi]; }
	bool has_edge(EI ei) const       	        { return (*ee_)[ei]; }

	const Graph* get_graph() const { return graph_; }
	const Graph::Edge& edge() const { return edge_; }
	const SBG* parent() const { return prev_; }

	VI dfsc_vindex(VI sbg_vi, VI vi_default = VI_NULL) const;
	const SBG* find_sbg_by_depth(int depth) const;
	void insert_to_automorph_list(SBG*);
	SBG* automorph_next;
	SBG* automorph_prev;

	bool equal(const SBG& s) const
	    {
		if (sum_ != s.sum_)
		    return false;

		const VecEE& ee1 = *ee_;
		const VecEE& ee2 = *s.ee_;
		const VecEE::value_type* p1 = ee1.data();
		const VecEE::value_type* p2 = ee2.data();
		int n = ee1.size();
		for (int i = n; i >= 0; --i)
		    if (((p1[i]^p2[i]) | (p1[n-i]^p2[n-i]))) return false;
		return true;
	    }

	bool operator== (const SBG& s) const { return equal(s); }
    };


    // *****************************************************************************
    //                          Projected
    // *****************************************************************************
    //
    // Projected part visible for user
    //

    class SubgraphsOfOneGraph
    {
	typedef std::vector<SBG*> SBGS_PTR;
	SBGS_PTR sbgs_uniq_ptr_;
	static void insert_sbgs_uniq_ptr(SBGS_PTR& cont, SBG* s);
	
	std::set<const SBG*> sbg_parents_;
	mutable bool support_valid_;
	mutable int support_;
	void calc_support_v() const;
    public:
	SubgraphsOfOneGraph() :support_valid_(false) {}
	typedef typename SBGS_PTR::const_iterator const_iterator;
	const_iterator begin() const { return sbgs_uniq_ptr_.begin(); }
	const_iterator end()   const { return sbgs_uniq_ptr_.end(); }
	int support() const { if (! support_valid_) calc_support_v(); return support_; }
	int size() const { return sbgs_uniq_ptr_.size(); }
	int num_sbgs_uniq() const { return size(); }
	void push(SBG* s);
	bool is_equal_occurrence(const SubgraphsOfOneGraph& prj) const;
    };

    
    class SubgraphsOfManyGraph
    {
	typedef SubgraphsOfOneGraph SOG;
	typedef std::map<const Graph* const, SOG> G_SOG;
	G_SOG g_sog_;
	int num_sbgs_uniq_;
    public:
	SubgraphsOfManyGraph() :num_sbgs_uniq_(0) {}
	typedef typename G_SOG::const_iterator const_iterator;
	const_iterator begin() const	{ return g_sog_.begin(); }
	const_iterator end()   const	{ return g_sog_.end(); }
	int support() const		{ return g_sog_.size(); }
	int num_sbgs_uniq() const	{ return num_sbgs_uniq_; }
	void push(SBG* s);
	bool is_equal_occurrence(const SubgraphsOfManyGraph& prj) const;
    };

    // *****************************************************************************
    //                          GspanResult
    // *****************************************************************************
    class GspanResult
    {
    public:
	virtual ~GspanResult() {}
	virtual void operator() (const DFSCode& dfsc, const SubgraphsOfOneGraph& sog) {}
	virtual void operator() (const DFSCode& dfsc, const SubgraphsOfManyGraph& smg) {}
    };

    
    void closegraph(const Graph& graph, int minsup, GspanResult* result, int max_trace_depth = 0);
    void closegraph(const std::vector<const Graph*>&, int minsup, GspanResult* result, int max_trace_depth = 0);

}

#endif
