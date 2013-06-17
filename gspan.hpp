#ifndef GSPAN_H_
#define GSPAN_H_

#include "gspan_allocator.hpp"
#include "gspan_graph.hpp"

#include <iostream>
#include <vector>
#include <cassert>

#include <boost/intrusive/list.hpp>

#ifdef USE_ASM
#include <stdint.h>
#endif

#ifndef BR
#define BR asm volatile ("int3;")
#endif

#ifndef SERIALIZE
#define SERIALIZE asm volatile ("xor %%eax, %%eax\n\t cpuid\n\t" : : : "eax", "ebx", "ecx", "edx")
#endif

#ifndef NOINLINE
#define NOINLINE __attribute__((noinline))
#endif

#if defined(USE_ASM)
#define PREFETCH(addr)  asm("prefetcht0 %0\n\t" : :"m"((addr)))
#else
#define PREFETCH(addr)
#endif

namespace gSpan
{
    typedef unsigned int support_type;

    // *****************************************************************************
    //                          EdgeCode
    // *****************************************************************************
    class EdgeCode
    {
        friend struct EdgeCodeCmpLex;
#if defined(USE_ASM)
        /* for little-endian
         * +-----------------+-----------------+-----------------+
         * |                 |                 |                 |
         * 0        2        4        6        8        10       |
         * +--------+--------+--------+--------+--------+--------+
         * |        |        |        |        |        |        |
         * |  ZERO  |VL_DST_I|  EL_I  |VL_SRC_I|VI_DST_I|VI_SRC_I|
         * +--------+--------+--------+--------+--------+--------+
         *
         * for big-endian (logical sense of that)
         * +-----------------+-----------------+-----------------+
         * |                 |                 |                 |
         * 0        2        4        6        8        10       |
         * +--------+--------+--------+--------+--------+--------+
         * |        |        |        |        |        |        |
         * |VI_SRC_I|VI_DST_I|VL_SRC_I|  EL_I  |VL_DST_I|  ZERO  |
         * +--------+--------+--------+--------+--------+--------+
         */
        enum { ZERO, VL_DST_I, EL_I, VL_SRC_I, VI_DST_I, VI_SRC_I };
        uint16_t x_[6];
#else
        DfscVI vi_src_, vi_dst_;
        VL vl_src_, vl_dst_;
        EL el_;
#endif
        bool is_fwd_;

    public:

#if defined(USE_ASM)
        EdgeCode()
            :is_fwd_(false)
            {
                assert(is_aligned(this));
                x_[VI_SRC_I] = x_[VI_DST_I] = VI_NULL;
                x_[VL_SRC_I] = x_[VL_DST_I] = VL_NULL;
                x_[EL_I] = EL_NULL;
                x_[ZERO] = 0U;
            }

        EdgeCode(DfscVI vi_src, DfscVI vi_dst, VL vl_src, EL el, VL vl_dst, bool fwd)
            :is_fwd_(fwd)
            {
                assert(is_aligned(this));
                x_[VI_SRC_I] = vi_src; x_[VI_DST_I] = vi_dst;
                x_[VL_SRC_I] = vl_src; x_[VL_DST_I] = vl_dst;
                x_[EL_I] = el;
                x_[ZERO] = 0U;
            }

        DfscVI vi_src() const   { return x_[VI_SRC_I]; }
        DfscVI vi_dst() const   { return x_[VI_DST_I]; }
        VL vl_src() const       { return x_[VL_SRC_I]; }
        VL vl_dst() const       { return x_[VL_DST_I]; }
        EL el() const           { return x_[EL_I]; }
        void chgdir() { std::swap(x_[VI_SRC_I], x_[VI_DST_I]); std::swap(x_[VL_SRC_I], x_[VL_DST_I]); }
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

        DfscVI vi_src() const   { return vi_src_; }
        DfscVI vi_dst() const   { return vi_dst_; }
        VL vl_src() const       { return vl_src_; }
        VL vl_dst() const       { return vl_dst_; }
        EL el() const           { return el_; }
        void chgdir() { std::swap(vi_src_, vi_dst_); std::swap(vl_src_, vl_dst_); }
#endif

        bool is_forward() const         { return is_fwd_; }
        bool is_backward() const        { return !is_fwd_; }
        EdgeCode operator- () const { EdgeCode ec(*this); ec.chgdir(); return ec; }
        bool operator!= (const EdgeCode& ec) const      { return ! (*this == ec); }
        bool operator== (const EdgeCode& ec) const;
    } __attribute__ ((aligned));

    struct EdgeCodeCmpDfs
    {
        bool operator() (const EdgeCode& ec1, const EdgeCode& ec2) const;
    };

    struct EdgeCodeCmpLex
    {
        bool operator() (const EdgeCode& ec1, const EdgeCode& ec2) const;
    };


    // *****************************************************************************
    //                          DFSCode
    // *****************************************************************************

    // for dfsc check minimality, we use std::vector as DFSCode
    typedef std::vector<EdgeCode> EdgeCodeVector;

    class DFSCode
    {
        EdgeCodeVector dfsc_;
        Graph graph_;
    public:
        DFSCode(std::size_t max_num_vertices, std::size_t max_num_edges)
            :graph_(max_num_vertices, max_num_edges) { dfsc_.reserve(max_num_edges); }
        
        // dfscode is a stack
        typedef std::vector<EdgeCode>::const_iterator const_iterator;
        const EdgeCode& top() const     { return dfsc_.back(); }
        const_iterator begin() const    { return dfsc_.begin(); }
        const_iterator end() const      { return dfsc_.end(); }
        std::size_t size() const        { return dfsc_.size(); }
        void push(const EdgeCode& ec)   { dfsc_.push_back(ec); graph_.push_edge(ec); }
        void pop()                      { dfsc_.pop_back(); graph_.pop_edge(); }
        const EdgeCode& operator[] (std::size_t i) const { return dfsc_[i]; }
        const Graph& get_graph() const  { return graph_; }
        const EdgeCodeVector& get_ecvector() const { return dfsc_; }
    };

    // *****************************************************************************
    //                          RMPath
    // *****************************************************************************
    class RMPath : private boost::noncopyable
    {
        typedef std::vector<int, STL_Allocator<int> >	RMPath_I_EdgeCode;
        typedef std::vector<DfscVI, STL_Allocator<DfscVI> >	DfscVI_Vertices;

        RMPath_I_EdgeCode       rmpath_i_ec_;
        DfscVI_Vertices         dfsc_vertices_;
        std::size_t             n_rmp_;

        static std::size_t make_rmpath(RMPath_I_EdgeCode& rmpath_i_ec,
                                       DfscVI_Vertices& dfsc_vertices,
                                       const EdgeCodeVector& dfsc);
    public:

        RMPath(const EdgeCodeVector& dfsc, MemAllocator* mem_alloc)
	    :rmpath_i_ec_(STL_Allocator<int>(mem_alloc)),
	     dfsc_vertices_(STL_Allocator<DfscVI>(mem_alloc))
	    {
		rmpath_i_ec_.reserve(dfsc.size());
		n_rmp_ = make_rmpath(rmpath_i_ec_, dfsc_vertices_, dfsc);
	    }

        //
        // there are two subset of the vertices 
        // [ begin()   ... rmp_end() )    right most path vertices
        // [ rmp_end() ... end()    ]     all other vertices
        //
        typedef DfscVI_Vertices::const_iterator const_iterator;
        const_iterator begin() const            { return dfsc_vertices_.begin(); }
        const_iterator rmp_end() const          { return dfsc_vertices_.begin() + n_rmp_; }
        const_iterator end() const              { return dfsc_vertices_.end(); }

        std::size_t num_all_vertices() const    { return dfsc_vertices_.size(); }
        std::size_t num_rmp_vertices() const    { return n_rmp_; }


        int operator[] (int i) const    { return rmpath_i_ec_[i]; }
        int num_edges() const		{ return rmpath_i_ec_.size(); }
        int rightmost_edgeindex() const { return rmpath_i_ec_[0]; }
        bool is_rightmost_edgeindex(int i) const
            {
                return std::find(rmpath_i_ec_.begin(), rmpath_i_ec_.end(), i) != rmpath_i_ec_.end();
            }

        bool is_rightmost_vertex(DfscVI vi) const;
        DfscVI rightmost_vertex() const { return dfsc_vertices_[0]; }

        friend std::ostream& operator<<(std::ostream& out, const RMPath& r);
    };


    // *****************************************************************************
    //                          Bitset
    // *****************************************************************************
    class Bitset : private boost::noncopyable
    {
    public:
        typedef unsigned int size_type;
        Bitset(MemAllocator*, size_type size);
        Bitset(MemAllocator*, size_type size, const Bitset& bf);
        ~Bitset();
        void free_resource(MemAllocator*, size_type);
        bool test(size_type pos, size_type size) const;
        bool set(size_type pos, size_type size);
        bool clear(size_type pos, size_type size);
        bool all_set(size_type size) const;
        static bool is_equal(const Bitset& bf1, const Bitset& bf2, size_type size);
    private:
        typedef unsigned long long Chunk;
        static size_type num_chunkbits()                { return CHAR_BIT * sizeof(Chunk); }
        static std::size_t num_chanks(size_type size)   { return (size + num_chunkbits() - 1U) / num_chunkbits(); }
        bool use_bitfield(size_type size) const         { return size <= num_chunkbits(); }
        union
        {
            Chunk       field;
            Chunk*      ptr;            
        } bit_;

        friend void print(std::ostream& out, const Bitset& bf, Bitset::size_type size);
    };
    
    // *****************************************************************************
    //                          SBGBase
    // *****************************************************************************
    template<class S> class SBGCreator;
    template<class S> struct Make_SBG_List_Embedding;
    
    template<class SBG_Derived>
    class SBGBase : private boost::noncopyable
    {
        friend class SBGCreator<SBG_Derived>;
        
        template<class S>
        friend struct Make_SBG_List_Embedding;

        const Graph::Edge* edge_;

        // single linked list of the subgraph edges
	const SBG_Derived* parent_;

        unsigned short depth_;
        
        Bitset  v_bits_;
        Bitset  e_bits_;
        Bitset  v_no_exts_bits_;

        const Graph* graph_;
        
        // single linked list all SBG in given embedding
	SBG_Derived* embedding_next_;

        void free_bitsets(MemAllocator* ma)
            {
                v_bits_.free_resource(ma, graph_->num_vertices());
                e_bits_.free_resource(ma, graph_->num_edges());
                v_no_exts_bits_.free_resource(ma, graph_->num_vertices());
            }

    protected:
	SBGBase(MemAllocator* ma, const Graph::Edge* e, const Graph* g)
	    :edge_(e),
             parent_(0),
             depth_(1),
             v_bits_(ma, g->num_vertices()),
             e_bits_(ma, g->num_edges()),
             v_no_exts_bits_(ma, g->num_vertices()),
             graph_(g),
             embedding_next_(0)
            {
                v_bits_.set(e->vi_src(), g->num_vertices());
                v_bits_.set(e->vi_dst(), g->num_vertices());
                e_bits_.set(e->eid(), g->num_edges());

                assert(is_aligned(this));
            }
	
	SBGBase(MemAllocator* ma, const Graph::Edge* e, const SBG_Derived* s)
	    :edge_(e),
	     parent_(s),
	     depth_(s->depth_ + 1),
             v_bits_(ma, s->graph_->num_vertices(), s->v_bits_),
             e_bits_(ma, s->graph_->num_edges(), s->e_bits_),
             v_no_exts_bits_(ma, s->graph_->num_vertices(), s->v_no_exts_bits_),
	     graph_(s->graph_),
	     embedding_next_(0)
	    {
                v_bits_.set(e->vi_src(), graph_->num_vertices());
                v_bits_.set(e->vi_dst(), graph_->num_vertices());
                e_bits_.set(e->eid(), graph_->num_edges());

                assert(is_aligned(this));
            }
    public:
        const Graph::Edge* edge() const		{ return edge_; }
	const Graph* get_graph() const		{ return graph_; }
	const SBG_Derived* parent() const	{ return parent_; }
	std::size_t num_edges() const		{ return depth_; }
	unsigned short depth() const		{ return depth_; }

        // return true if vertex or edge is part of the subgraph
	bool has_vertex(GraphVI v) const	{ return v_bits_.test(v, graph_->num_vertices()); }
	bool has_edge(GraphEI e) const		{ return e_bits_.test(e, graph_->num_edges()); }
	
	bool has_no_extension(GraphVI v)	{ return v_no_exts_bits_.test(v, graph_->num_vertices()); }
	void set_no_has_extension(GraphVI v)	{ v_no_exts_bits_.set(v, graph_->num_vertices()); }
	void set_has_extension(GraphVI v)	{ v_no_exts_bits_.clear(v, graph_->num_vertices()); }

	SBG_Derived* next_embedding() const	{ return embedding_next_; }

	//
        // true in case of the automorphism
	// this two graphs consist of the same edges
        //
	friend bool operator== (const SBGBase& s1, const SBGBase& s2)
            {
                assert(s1.graph_->num_edges() == s2.graph_->num_edges());
                return Bitset::is_equal(s1.e_bits_, s2.e_bits_, s1.graph_->num_edges());
            }
    };

    // *****************************************************************************
    //                          SBG
    // *****************************************************************************
    class SBG : public SBGBase<SBG>
    {
        friend class SBGCreator<SBG>;
        friend struct Make_SBG_List_Autgroup;

        // single linked list
        SBG* automorph_group_next_;

	// double linked list of the automorphic subgraphs
	SBG* automorph_next_;
	SBG* automorph_prev_;

	GraphVI* vi_dfsc_to_graph_;

	DfscVI vi_src_dfsc_;
	DfscVI vi_dst_dfsc_;

	// size vi_dfsc_to_graph_ array
	unsigned short num_vertices_;

	// private ctor and dtor
	SBG(MemAllocator* ma, const Graph::Edge* e, const Graph* g)
	    :SBGBase<SBG>(ma, e, g),
             automorph_group_next_(0),
	     automorph_next_(this),
	     automorph_prev_(this),
	     vi_dfsc_to_graph_(0),
	     vi_src_dfsc_(0),
	     vi_dst_dfsc_(1),
	     num_vertices_(2)
	    { init_dfsc_to_graph_array1(ma); }

	SBG(MemAllocator* ma, const Graph::Edge* e, const SBG* s, const EdgeCode& ec)
	    :SBGBase<SBG>(ma, e, s),
             automorph_group_next_(0),
	     automorph_next_(this),
	     automorph_prev_(this),
	     vi_dfsc_to_graph_(0),
	     vi_src_dfsc_(ec.vi_src()),
	     vi_dst_dfsc_(ec.vi_dst()),
	     num_vertices_(s->num_vertices_ + ec.is_forward())
	    { init_dfsc_to_graph_array2(ma); }

	void init_dfsc_to_graph_array1(MemAllocator* mem_alloc);
	void init_dfsc_to_graph_array2(MemAllocator* mem_alloc);
	void free_dfsc_to_graph_array(MemAllocator* mem_alloc)
	    { mem_alloc->dealloc_array(vi_dfsc_to_graph_, num_vertices()); }
    public:

        DfscVI dfsc_src() const { return vi_src_dfsc_; }
        DfscVI dfsc_dst() const { return vi_dst_dfsc_; }

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
	DfscVI* create_graph_to_dfsc_v(MemAllocator* mem_alloc, DfscVI vi_default = VI_NULL) const;
	void free_graph_to_dfsc_v(DfscVI*, MemAllocator* mem_alloc) const;

        static void insert_to_automorph_list(SBG* pos, SBG* s);
	const SBG* next_automorph() const       { return automorph_next_; }
        const SBG* next_automorph_group() const { return automorph_group_next_; }
        SBG* next_automorph_group()             { return automorph_group_next_; }
        

	friend std::ostream& operator<<(std::ostream& out, const SBG& sbg);
    };

    void get_chain(std::vector<const SBG*>&, const SBG* sbg);

    // *****************************************************************************
    //                          specialization SBGCreator
    //                               for SBG
    // *****************************************************************************
    template<>
    class SBGCreator<SBG>
    {
	MemAllocator*	ma_;
	FixedAllocator*	sbg_alloc_;
    public:
	explicit SBGCreator(MemAllocator* mem_alloc)
            :ma_(mem_alloc), sbg_alloc_(ma_->get_fixed_allocator(sizeof(SBG))) {}        
	MemAllocator* get_mem_allocator() const { return ma_; }
	
        SBG* new_sbg(const Graph::Edge* e, const Graph* g)
            { return new (sbg_alloc_->allocate()) SBG(ma_, e, g); }

	SBG* new_sbg(const Graph::Edge* e, const SBG* s, const EdgeCode& ec)
            { return new (sbg_alloc_->allocate()) SBG(ma_, e, s, ec); }

	void delete_sbg(SBG* s)
            {
                s->free_bitsets(ma_);
                s->free_dfsc_to_graph_array(ma_);
                s->~SBG();
                sbg_alloc_->deallocate(s);
            }
    };


    // *****************************************************************************
    //                          SBG_List
    // *****************************************************************************
    namespace sbg_list_impl
    {
        template<class S>
        class DeletePolicy
        {
            SBGCreator<S>* sbg_creator_;
        public:
            DeletePolicy(SBGCreator<S>* sbg_creator) :sbg_creator_(sbg_creator) {}
            void operator() (S* s) { sbg_creator_->delete_sbg(s); }
        };

        class DoNothing {};

        template<class S, class S2, S* (S2::*mp_next), class DelPolicy>
        class SBG_List : private boost::noncopyable
        {
            S* lh_;
            std::size_t size_;
            DelPolicy del_policy_;

            void clear(DeletePolicy<S>)
                {
                    S* p = lh_;
                    while (p)
                    {
                        S* next = p->*mp_next;
                        del_policy_(p);
                        p = next;
                    }
                    lh_ = 0;
                    size_ = 0;
                }

            void clear(DoNothing)           {}
        public:
            explicit SBG_List(DelPolicy del_policy = DelPolicy()) :lh_(0), size_(0), del_policy_(del_policy)    {}
            ~SBG_List()                     { clear(); }
            void clear()                    { clear(del_policy_); }
            void insert(S* s)               { s->*mp_next = lh_; lh_ = s; ++size_; }
            const S* first() const          { return lh_; }
            S* first()                      { return lh_; }
            std::size_t size() const        { return size_; }
            bool empty() const              { return size_ == 0; }

            friend std::ostream& operator<<(std::ostream& out, const SBG_List<S,S2,mp_next,DelPolicy>& li)
                {
                    for (const S* s = li.first(); s; s = s->*mp_next) out << '\t' << *s << std::endl;
                    return out;
                }
        };
    } // end: namespace sbg_list_impl

    template<class S>
    struct Make_SBG_List_Embedding {
        typedef SBGBase<S> S2;
        typedef sbg_list_impl::SBG_List<S, S2, &S2::embedding_next_, sbg_list_impl::DeletePolicy<S> > Type;
    };

    struct Make_SBG_List_Autgroup {
        typedef sbg_list_impl::SBG_List<SBG, SBG, &SBG::automorph_group_next_, sbg_list_impl::DoNothing> Type;
    };


    // *****************************************************************************
    //                          Intrusive Containers
    // *****************************************************************************
    
    template<class Key, class SG>
    class SGSetNode : public SG
    {
        Key key_;
    public:
        SGSetNode(const Key& key, SBGCreator<SBG>* sbg_creator) : SG(sbg_creator), key_(key) {}
        typedef Key key_type;
        const Key& get_key() const { return key_; }
        typedef boost::intrusive::set_member_hook<
            boost::intrusive::link_mode<boost::intrusive::normal_link>,
            boost::intrusive::optimize_size<true>
            > Hook_;
        Hook_ hook_;
        typedef boost::intrusive::member_hook<SGSetNode, Hook_, &SGSetNode::hook_> HookOptions_;
    };


    template<class Key, class SG>
    class SGListNode : public SG
    {
        Key key_;
    public:
        SGListNode(const Key& key, SBGCreator<SBG>* sbg_creator) : SG(sbg_creator), key_(key) {}
        const Key& get_key() const { return key_; }
        typedef boost::intrusive::list_member_hook<boost::intrusive::link_mode<boost::intrusive::normal_link> > Hook_;
        Hook_ hook_;
        typedef boost::intrusive::member_hook<SGListNode, Hook_, &SGListNode::hook_> HookOptions_;
    };

    
    template<class Node, class Cmp>
    class SGNodeCompare
    {
        Cmp cmp;
    public:
        bool operator() (const Node& n1, const Node& n2) const                  { return cmp(n1.get_key(), n2.get_key()); }
        bool operator() (const typename Node::key_type& k, const Node& n) const { return cmp(k, n.get_key()); }
        bool operator() (const Node& n, const typename Node::key_type& k) const { return cmp(n.get_key(), k); }
    };


    template<class Node, bool const_time_size, class size_type, class Cmp>
    class SGSet
        : public boost::intrusive::set<Node,
                                       typename Node::HookOptions_,
                                       boost::intrusive::constant_time_size<const_time_size>,
                                       boost::intrusive::size_type<size_type>,
                                       boost::intrusive::compare<Cmp>
                                       >
    {
    };


    template<class Node, class Cmp>
    class SGSet<Node, false, void, Cmp>
        : public boost::intrusive::set<Node,
                                       typename Node::HookOptions_,
                                       boost::intrusive::constant_time_size<false>,
                                       boost::intrusive::compare<Cmp>
                                       >
    {
    };



    template<class Node, bool const_time_size, class size_type, class Cmp>
    class SGList
        : public boost::intrusive::list<Node,
                                       typename Node::HookOptions_,
                                       boost::intrusive::constant_time_size<const_time_size>,
                                       boost::intrusive::size_type<size_type>
                                       >
    {
    };

    template<class Node, class Cmp>
    class SGList<Node, false, void, Cmp>
        : public boost::intrusive::list<Node,
                                       typename Node::HookOptions_,
                                       boost::intrusive::constant_time_size<false>
                                       >
    {
    };


    // intrusive container type constructor
    enum ContSelector { SET, LIST };

    template<ContSelector, class Key, class SG, class Cmp, bool const_time_size, class size_type>
    struct ContainerType;
    
    template<class Key, class SG, class Cmp, bool const_time_size, class size_type>
    struct ContainerType<SET, Key, SG, Cmp, const_time_size, size_type>
    {
        typedef SGSetNode<Key, SG> NodeType;
        typedef SGNodeCompare<NodeType, Cmp> NodeCompare;
        typedef SGSet<NodeType, const_time_size, size_type, NodeCompare> ContType;
    };

    template<class Key, class SG, class Cmp, bool const_time_size, class size_type>
    struct ContainerType<LIST, Key, SG, Cmp, const_time_size, size_type>
    {
        typedef SGListNode<Key, SG> NodeType;
        typedef SGNodeCompare<NodeType, Cmp> NodeCompare;
        typedef SGList<NodeType, const_time_size, size_type, NodeCompare> ContType;
    };



    // *****************************************************************************
    //                          Intrusive Set
    // *****************************************************************************

    template<class Key, class SG, class Cmp = std::less<Key> >
    class SetNode : public SG
    {
        Key key_;
    public:
        SetNode(const Key& key, SBGCreator<SBG>* sbg_creator)
            :SG(sbg_creator), key_(key) {}
        
        const Key& get_key() const { return key_; }
        
        typedef boost::intrusive::set_member_hook<
            boost::intrusive::link_mode<boost::intrusive::normal_link>,
            boost::intrusive::optimize_size<true>
            > SetHook;
        SetHook set_hook;
        
        typedef boost::intrusive::member_hook<SetNode, SetHook, &SetNode::set_hook> SetOptions;

        friend bool operator< (const SetNode& n1, const SetNode& n2)
            { Cmp cmp; return cmp(n1.key_, n2.key_); }
            
        struct CompareKey
        {
            bool operator() (const Key& key, const SetNode& node) const
                { Cmp cmp; return cmp(key, node.key_); }
            
            bool operator() (const SetNode& node, const Key& key) const
                { Cmp cmp; return cmp(node.key_, key); }
        };
        friend struct CompareKey;

        typedef Cmp Compare;
    };
    
    template<class T, bool const_time_size, class size_type>
    class Set
        : public boost::intrusive::set<T,
                                       typename T::SetOptions,
                                       boost::intrusive::constant_time_size<const_time_size>,
                                       boost::intrusive::size_type<size_type>,
                                       boost::intrusive::compare<typename T::Compare>
                                       >
    {};

    template<class T>
    class Set<T,false,void>
        : public boost::intrusive::set<T,
                                       typename T::SetOptions,
                                       boost::intrusive::constant_time_size<false>,
                                       boost::intrusive::compare<typename T::Compare>
                                       >
    {};
    


    // *****************************************************************************
    //                          SubgraphsOfOneGraph
    // *****************************************************************************
    class SubgraphsOfOneGraph : private boost::noncopyable
    {
        typedef typename Make_SBG_List_Embedding<SBG>::Type     List_Embedding;
        typedef Make_SBG_List_Autgroup::Type                    List_Autgroup;
        List_Embedding  list_sbg_all_;
        List_Autgroup   list_sbg_autgroup_;
        mutable support_type support_;
        void calc_support() const;
    public:
        explicit SubgraphsOfOneGraph(SBGCreator<SBG>* sbg_creator)
            :list_sbg_all_(sbg_creator), support_(0)
            {}

        const SBG* first_all() const            { return list_sbg_all_.first(); }
        const SBG* first_autgroup() const       { return list_sbg_autgroup_.first(); }
        SBG* first_all()                        { return list_sbg_all_.first(); }
        SBG* first_autgroup()                   { return list_sbg_autgroup_.first(); }
        
        std::size_t size_list_all() const       { return list_sbg_all_.size(); }
        std::size_t size_list_autgroup() const  { return list_sbg_autgroup_.size(); }

        const Graph* get_graph() const          { return first_all()->get_graph(); }
        support_type support() const            { if (0 == support_) calc_support(); return support_; }
        bool is_equal_occurence(const SubgraphsOfOneGraph& parent) const;

        void insert(SBG* s);

        friend std::ostream& operator<<(std::ostream& out, const SubgraphsOfOneGraph& sog)
            {
                return out << sog.list_sbg_all_;
            }
    };

    // *****************************************************************************
    //                          SubgraphsOfManyGraph
    // *****************************************************************************
    class SubgraphsOfManyGraph : private boost::noncopyable
    {
        typedef ContainerType<SET, const Graph*, SubgraphsOfOneGraph, std::less<const Graph*>, true, support_type> CT;
        typedef typename CT::NodeType SOG;
        typedef typename CT::NodeCompare CompareKey;
        typedef typename CT::ContType SOGSet;
        
        SOGSet set_;
        FixedAllocator* fa_;
        SBGCreator<SBG>* sbg_creator_;
        std::size_t size_list_all_;
    public:
        explicit SubgraphsOfManyGraph(SBGCreator<SBG>* sbg_creator);
        ~SubgraphsOfManyGraph();

        void insert(SBG* s);

        typedef SOGSet::const_iterator  const_iterator;
        typedef SOGSet::iterator        iterator;
        const_iterator begin() const    { return set_.begin(); }
        const_iterator end() const      { return set_.end(); }
        iterator begin()                { return set_.begin(); }
        iterator end()                  { return set_.end(); }
        const SubgraphsOfOneGraph& find(const Graph*) const;
        
        std::size_t size_list_all() const       { return size_list_all_; }

        support_type support() const    { return set_.size(); }
        bool is_equal_occurence(const SubgraphsOfManyGraph&) const;

        friend std::ostream& operator<<(std::ostream& out, const SubgraphsOfManyGraph& smg)
            {
                for (const_iterator it = smg.begin(); it != smg.end(); ++it)
                    out << *it;
                return out;
            }
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


    void closegraph(const Graph& graph,               support_type minsup, GspanResult* result, int max_trace_depth = 0);
    void closegraph(const std::vector<const Graph*>&, support_type minsup, GspanResult* result, int max_trace_depth = 0);

}

#endif
