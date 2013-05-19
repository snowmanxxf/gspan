#ifndef GSPAN_H_
#define GSPAN_H_

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <cassert>
#include <cstring>	// memset(), memcpy()


#ifdef USE_ASM
#include <stdint.h>
#endif

#include "gspan_allocator.hpp"
#include "gspan_graph.hpp"

#ifndef BR
#define BR asm volatile ("int3;")
#endif

#ifdef USE_ASM
#define PREFETCH(addr)  asm("prefetcht0 %0\n" : :"m"((addr)))
#else
#define PREFETCH(addr)
#endif


namespace gSpan
{

    // *****************************************************************************
    //                          EdgeCode
    //				DFSCode
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

    // *****************************************************************************
    //                          DFSCode
    // *****************************************************************************
    class DFSCode
    {
	std::vector<EdgeCode> dfsc_;
	Graph graph_;
    public:
        DFSCode(std::size_t max_num_vertices) :graph_(max_num_vertices) {}

	// dfscode is a stack
	typedef std::vector<EdgeCode>::const_iterator const_iterator;
	const EdgeCode& top() const { return dfsc_.back(); }
	const_iterator begin() const { return dfsc_.begin(); }
	const_iterator end() const { return dfsc_.end(); }
	std::size_t size() const { return dfsc_.size(); }
	void push(const EdgeCode& ec)	{ dfsc_.push_back(ec); graph_.push_edge(ec); }
	void pop()			{ dfsc_.pop_back(); graph_.pop_edge(); }
	const EdgeCode& operator[] (int i) const { return dfsc_[i]; }
	void reserve(std::size_t s) { dfsc_.reserve(s); }

	const Graph& get_graph() const { return graph_; }
    };


    //typedef std::vector<EdgeCode> DFSCode;
    DfscVI max_vertex(const DFSCode& dfsc);
    

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

    template<class S> class SBGCreator;
    template<class S> class SBG_List;

    template<class SBG_Derived>
    class SBGBase : private boost::noncopyable
    {
	friend class SBGCreator<SBG_Derived>;
	friend class SBG_List<SBG_Derived>;

	const Graph::Edge* edge_;

	// single linked list of the subgraph edges
	const SBG_Derived* parent_;

	unsigned short depth_;

	EVBool* ev_array_;

	const Graph* graph_;

	// single linked list all SBG in given embedding
	SBG_Derived* next_embedding_;

	void init_ev_array1(MemAllocator* mem_alloc);
	void init_ev_array2(MemAllocator* mem_alloc);
	void free_ev_array(MemAllocator* mem_alloc)
	    { mem_alloc->dealloc_array(ev_array_, graph_size(get_graph())); }

    protected:
	SBGBase(const Graph::Edge* e, const Graph* g)
	    :edge_(e), parent_(0), depth_(1), ev_array_(0), graph_(g), next_embedding_(0) {}
	
	SBGBase(const Graph::Edge* e, const SBG_Derived* s)
	    :edge_(e),
	     parent_(s),
	     depth_(s->depth_ + 1),
	     ev_array_(0),
	     graph_(s->graph_),
	     next_embedding_(0)
	    {}

    public:
	const Graph::Edge* edge() const		{ return edge_; }
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

	SBG_Derived* next_embedding() const	{ return next_embedding_; }
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
    
    template<class SBG_Derived>
    void SBGBase<SBG_Derived>::init_ev_array1(MemAllocator* mem_alloc)
    {
	std::size_t ev_size = graph_size(graph_);
	std::size_t alloc_size = ev_size * sizeof(*ev_array_);
	ev_array_ = mem_alloc->alloc_array<EVBool>(ev_size);
	::memset(ev_array_, 0, alloc_size);

        ev_array_[edge_->eid()].e_in_sbg = true;
        ev_array_[edge_->vi_src()].v_in_sbg = true;
        ev_array_[edge_->vi_dst()].v_in_sbg = true;
    }

    template<class SBG_Derived>
    void SBGBase<SBG_Derived>::init_ev_array2(MemAllocator* mem_alloc)
    {
	PREFETCH(*parent_->ev_array_);
	std::size_t ev_size = graph_size(graph_);
	std::size_t alloc_size = ev_size * sizeof(*ev_array_);
	ev_array_ = mem_alloc->alloc_array<EVBool>(ev_size);
	::memcpy(ev_array_, parent_->ev_array_, alloc_size);
        ev_array_[edge_->eid()].e_in_sbg = true;
        ev_array_[edge_->vi_src()].v_in_sbg = true;
        ev_array_[edge_->vi_dst()].v_in_sbg = true;
    }


    // *****************************************************************************
    //                          SBG
    // *****************************************************************************
    class SBG : public SBGBase<SBG>
    {
	friend class SBGCreator<SBG>;

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
	SBG(const Graph::Edge* e, const Graph* g)
	    :SBGBase<SBG>(e, g),
	     vi_src_dfsc_(0),
	     vi_dst_dfsc_(1),
	     num_vertices_(2),
	     vi_dfsc_to_graph_(0),
	     automorph_next_(this),
	     automorph_prev_(this),
	     sum_(e->eid())
	    {}

	SBG(const Graph::Edge* e, const SBG* s, const EdgeCode& ec)
	    :SBGBase<SBG>(e, s),
	     vi_src_dfsc_(ec.vi_src()),
	     vi_dst_dfsc_(ec.vi_dst()),
	     num_vertices_(s->num_vertices() + ec.is_forward()),
	     vi_dfsc_to_graph_(0),
	     automorph_next_(this),
	     automorph_prev_(this),
	     sum_(s->sum_ + e->eid())
	    {}

	void init_dfsc_to_graph_array1(MemAllocator* mem_alloc);
	void init_dfsc_to_graph_array2(MemAllocator* mem_alloc);
	void free_dfsc_to_graph_array(MemAllocator* mem_alloc)
	    { mem_alloc->dealloc_array(vi_dfsc_to_graph_, num_vertices()); }
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

	//
	// array of the dfsc VI, indexed by the graph VI
	// so, sbg_vertices[graph_vi] == dfsc vi
	//
	DfscVI* create_graph_to_dfsc_v(MemAllocator* mem_alloc, DfscVI vi_default);
	void free_graph_to_dfsc_v(DfscVI*, MemAllocator* mem_alloc);

	static void insert_to_automorph_list(SBG* pos, SBG* s);
	const SBG* next_automorph() const { return automorph_next_; }


	/*
	 * true in case of the automorphism
	 */
	friend bool operator== (const SBG& s1, const SBG& s2);
	
	friend std::ostream& operator<<(std::ostream& out, const SBG& sbg);
    };

    void get_chain(std::vector<const SBG*>&, const SBG* sbg);


    // *****************************************************************************
    //                          SubgraphsOfGraph
    // *****************************************************************************

    class SubgraphsOfOneGraph
    {
	typedef std::vector<SBG*, STL_Allocator<SBG*> > SBGS_PTR;
	typedef std::set<const SBG*,
			 std::less<const SBG*>,
			 STL_Allocator<const SBG*> > SBG_Parents;

	SBGS_PTR sbgs_uniq_ptr_;	
	SBG_Parents sbg_parents_;
	mutable bool support_valid_;
	mutable int support_;

	void calc_support_v(const SBG* slist) const;
	static void insert_sbgs_uniq_ptr(SBGS_PTR& cont, SBG* s);
    public:
	explicit
	SubgraphsOfOneGraph(MemAllocator* mem_alloc)
	    :sbgs_uniq_ptr_(STL_Allocator<SBG*>(mem_alloc)),
	     sbg_parents_(std::less<const SBG*>(),
			  STL_Allocator<const SBG*>(mem_alloc)),
	     support_valid_(false)
	    {}

	typedef typename SBGS_PTR::const_iterator const_iterator;
	const_iterator begin() const	{ return sbgs_uniq_ptr_.begin(); }
	const_iterator end()   const	{ return sbgs_uniq_ptr_.end(); }
	std::size_t size() const	{ return sbgs_uniq_ptr_.size(); }
	int support() const
	    {
		assert(support_valid_);
		assert(support_ > 0);
		return support_;
	    }

	// ----------------------------
	// for internal use
	// ----------------------------
	int support(const SBG* slist) const
	    { if (!support_valid_) calc_support_v(slist); return support_; }
	
	bool is_equal_occurrence(const SubgraphsOfOneGraph& prj) const;

	void insert(SBG* s);
    };

    
    
    class SubgraphsOfManyGraph
    {
	typedef const Graph* const Key;
	typedef SubgraphsOfOneGraph SOG;
	typedef STL_Allocator<std::pair<Key, SOG> > G_SOG_Allocator;
	typedef std::map<Key, SOG, std::less<Key>, G_SOG_Allocator> G_SOG;
	G_SOG g_sog_;
	std::size_t size_;
    public:
	explicit
	SubgraphsOfManyGraph(MemAllocator* mem_alloc)
	    :g_sog_(std::less<Key>(), G_SOG_Allocator(mem_alloc)),
	     size_(0)
	    {}

	typedef typename G_SOG::iterator iterator;
	typedef typename G_SOG::const_iterator const_iterator;
	const_iterator begin() const	{ return g_sog_.begin(); }
	const_iterator end()   const	{ return g_sog_.end(); }
	int support() const		{ return g_sog_.size(); }
	int support(const SBG*) const	{ return g_sog_.size(); }
	std::size_t size() const	{ return size_; }
	bool is_equal_occurrence(const SubgraphsOfManyGraph& prj) const;
	void insert(SBG* s);
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
