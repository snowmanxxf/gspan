#ifndef GSPAN2_H_
#define GSPAN2_H_

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <cassert>

#ifdef USE_ASM
#include <stdint.h>
#endif

#include <boost/noncopyable.hpp>

#include "gspan_allocator.hpp"
#include "gspan_graph.hpp"


#ifndef BR
#define BR asm volatile ("int3;")
#endif


extern unsigned int nctor;
extern unsigned int ndtor;

namespace gSpan2
{

    // *****************************************************************************
    //                          EdgeCode
    // *****************************************************************************
    class EdgeCode
    {
	friend struct EdgeCodeCmpLex;
#ifdef USE_ASM
	// TODO: Compatibility check for DfscVI, VL and EL for uint16_t

	enum { VI_SRC_I, VI_DST_I, VL_SRC_I, EL_I, VL_DST_I };

	uint16_t x_[5];
#else
	DfscVI vi_src_, vi_dst_;
	VL vl_src_, vl_dst_;
	EL el_;
#endif
	bool is_fwd_;
    public:
#ifdef USE_ASM
	EdgeCode()
	    :is_fwd_(false)
	    {
		x_[VI_SRC_I] = x_[VI_DST_I] = VI_NULL;
		x_[VL_SRC_I] = x_[VL_DST_I] = VL_NULL;
		x_[EL_I] = EL_NULL;
	    }

	EdgeCode(DfscVI vi_src, DfscVI vi_dst, VL vl_src, EL el, VL vl_dst, bool fwd)
	    :is_fwd_(fwd)
	    {
		x_[VI_SRC_I] = vi_src; x_[VI_DST_I] = vi_dst;
		x_[VL_SRC_I] = vl_src; x_[VL_DST_I] = vl_dst;
		x_[EL_I] = el;
	    }
	
	DfscVI vi_src() const	{ return x_[VI_SRC_I]; }
	DfscVI vi_dst() const	{ return x_[VI_DST_I]; }
	VL vl_src() const	{ return x_[VL_SRC_I]; }
	VL vl_dst() const	{ return x_[VL_DST_I]; }
	EL el() const		{ return x_[EL_I]; }

	void chgdir()
	    {
		std::swap(x_[VI_SRC_I], x_[VI_DST_I]);
		std::swap(x_[VL_SRC_I], x_[VL_DST_I]);
	    }
#else
	EdgeCode()
	    :vi_src_(VI_NULL), vi_dst_(VI_NULL),
	     vl_src_(VL_NULL), vl_dst_(VL_NULL), el_(EL_NULL),
	     is_fwd_(false)
	    {}

	EdgeCode(DfscVI vi_from, DfscVI vi_to, VL vl_from, EL el, VL vl_to, bool fwd)
	    :vi_src_(vi_from), vi_dst_(vi_to),
	     vl_src_(vl_from), vl_dst_(vl_to), el_(el),
	     is_fwd_(fwd) {}

	DfscVI vi_src() const	{ return vi_src_; }
	DfscVI vi_dst() const	{ return vi_dst_; }
	VL vl_src() const	{ return vl_src_; }
	VL vl_dst() const	{ return vl_dst_; }
	EL el() const		{ return el_; }
	void chgdir() { std::swap(vi_src_, vi_dst_); std::swap(vl_src_, vl_dst_); }
#endif

	bool is_forward() const		{ return is_fwd_; }
	bool is_backward() const	{ return !is_fwd_; }
	EdgeCode operator- () const { EdgeCode ec(*this); ec.chgdir(); return ec; }
	bool operator!= (const EdgeCode& ec) const	{ return ! (*this == ec); }
	bool operator== (const EdgeCode& ec) const;
    };

    struct EdgeCodeCmpDfs
    {
	bool operator() (const EdgeCode& ec1, const EdgeCode& ec2) const  __attribute__((noinline));
    };

    struct EdgeCodeCmpLex
    {
	bool operator() (const EdgeCode& ec1, const EdgeCode& ec2) const  __attribute__((noinline));
    };


    std::ostream& operator<<(std::ostream& out, const EdgeCode& ec);
	
    // *****************************************************************************
    //                          DFSCode
    // *****************************************************************************
    typedef std::vector<EdgeCode> DFSCode;
    
    DfscVI max_vertex(const DFSCode& dfsc);
    std::ostream& operator<<(std::ostream& out, const DFSCode& dfsc);

    // *****************************************************************************
    //                          SBGBase
    // *****************************************************************************

    struct EVBool
    {
	bool e_in_sbg : 1;
	bool v_in_sbg : 1;
	bool v_no_exts : 1;
    };

    inline std::size_t graph_size(const Graph* g)
    { return std::max(g->num_vertices(), g->num_edges()); }


    class SBG_Creator;

    template<class SBG_Derived>
    class SBGBase : private boost::noncopyable
    {
	friend class SBG_Creator;

	Graph::Edge edge_;

	// single linked list of the subgraph edges
	const SBG_Derived* parent_;

	unsigned short depth_;

	EVBool* ev_array_;

	const Graph* graph_;
    protected:
	SBGBase(const Graph::Edge& e, const Graph* g)
	    :edge_(e), parent_(0), depth_(1), ev_array_(0), graph_(g) {}
	
	SBGBase(const Graph::Edge& e, const SBG_Derived* s)
	    :edge_(e), parent_(s), depth_(s->depth_ + 1), ev_array_(0), graph_(s->graph_) {}
    public:
	const Graph::Edge& edge() const		{ return edge_; }
	const Graph* get_graph() const		{ return graph_; }
	const SBG_Derived* parent() const	{ return parent_; }
	std::size_t num_edges() const		{ return depth_; }
	unsigned short depth() const		{ return depth_; }

	bool has_vertex(GraphVI v) const	{ return ev_array_[v].v_in_sbg; }
	bool has_edge(GraphEI e) const		{ return ev_array_[e].e_in_sbg; }
	
	bool has_no_extension(GraphVI v)	{ return ev_array_[v].v_no_exts; }
	void set_no_has_extension(GraphVI v)	{ ev_array_[v].v_no_exts = true; }
	void set_has_extension(GraphVI v)	{ ev_array_[v].v_no_exts = false; }

	const EVBool* ev_array() const		{ return ev_array_; }
    };

    template<class S>
    const S* parent(const S* sbg, unsigned short depth)
    {
	const S* s = sbg;
	do
	{
	    if (depth == s->depth())
		return s;
	    s = s->parent();
	} while (s);
	return 0;
    }
    

    // *****************************************************************************
    //                          SBGSimple
    // *****************************************************************************
    class SBGSimple : public SBGBase<SBGSimple>
    {
	friend class SBG_Creator;

	Graph::Edge** edge_ptr_array_;

	SBGSimple(const Graph::Edge& e, const Graph* g)
	    :SBGBase<SBGSimple>(e, g),
	     edge_ptr_array_(0)
	    {}

	SBGSimple(const Graph::Edge& e, const SBGSimple* s)
	    :SBGBase<SBGSimple>(e, s),
	     edge_ptr_array_(0)
	    {}
    public:
	const Graph::Edge* operator[] (int i) const	{ return edge_ptr_array_[i]; }
    };


    // *****************************************************************************
    //                          SBG
    // *****************************************************************************
    class SBG : public SBGBase<SBG>
    {
	friend class SBG_Creator;

	DfscVI vi_src_dfsc_;
	DfscVI vi_dst_dfsc_;

	// size vi_dfsc_to_graph_ array
	unsigned short num_vertices_;

	GraphVI* vi_dfsc_to_graph_;

	// double linked list of the automorphic subgraphs
	SBG* automorph_next_;
	SBG* automorph_prev_;

	unsigned short sum_;

	// private ctor and dtor
	SBG(const Graph::Edge& e, const Graph* g)
	    :SBGBase<SBG>(e, g),
	     vi_src_dfsc_(0),
	     vi_dst_dfsc_(1),
	     num_vertices_(2),
	     vi_dfsc_to_graph_(0),
	     automorph_next_(this),
	     automorph_prev_(this),
	     sum_(e.eid())
	    {}

	SBG(const Graph::Edge& e, const SBG* s, const EdgeCode& ec)
	    :SBGBase<SBG>(e, s),
	     vi_src_dfsc_(ec.vi_src()),
	     vi_dst_dfsc_(ec.vi_dst()),
	     num_vertices_(s->num_vertices() + ec.is_forward()),
	     vi_dfsc_to_graph_(0),
	     automorph_next_(this),
	     automorph_prev_(this),
	     sum_(s->sum_ + e.eid())
	    {}
    public:

	/*
	 * return array of the graph VI, indexed by the dfsc VI
	 * so, sbg_vertices[dfsc_vi] == graph vi
	 * size is number of the vertices, that formed subgraph
	 */
	const GraphVI* get_dfsc_to_graph_v() const { return vi_dfsc_to_graph_; }

	/*
	 * vertices and edge number formed subgraph
	 */
	std::size_t num_vertices() const	{ return num_vertices_; }


	static void insert_to_automorph_list(SBG* pos, SBG* s);
	const SBG* next_automorph() const { return automorph_next_; }

	/*
	 * true in case of the automorphism
	 */
	friend bool operator== (const SBG& s1, const SBG& s2);
	
	friend std::ostream& operator<<(std::ostream& out, const SBG& sbg);
    };

    void get_chain(std::vector<const SBG*>, const SBG* sbg);

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
