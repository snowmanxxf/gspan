#include "gspan.hpp"

#include <boost/foreach.hpp>

#include <list>
#include <queue>
#include <algorithm>

namespace gSpan
{
    // *****************************************************************************
    //                          EdgeCode
    // *****************************************************************************

    struct EdgeCodeCmpDfs
    {
	bool operator() (const EdgeCode& ec1, const EdgeCode& ec2) const __attribute__((noinline));
    };

    bool EdgeCodeCmpDfs::operator() (const EdgeCode& ec1, const EdgeCode& ec2) const
    {
        bool ec1_f = ec1.is_forward();
        bool ec2_f = ec2.is_forward();

        if (!ec1_f && ec2_f)
            return true;
	else if (!ec1_f && !ec2_f)
	{
	    if (ec1.vi_dst()  < ec2.vi_dst())
		return true;
	    else if (ec1.vi_dst() == ec2.vi_dst() && ec1.el() < ec2.el())
		return true;
	    else
		return false;
	}
	else if (ec1_f && ec2_f)
	{
	    if (ec1.vi_src() > ec2.vi_src())
		return true;
	    else if (ec1.vi_src() == ec2.vi_src() && ec1.vl_src()  < ec2.vl_src())
		return true;
	    else if (ec1.vi_src() == ec2.vi_src() && ec1.vl_src() == ec2.vl_src() && ec1.el() < ec2.el())
		return true;
	    else if (ec1.vi_src() == ec2.vi_src() && ec1.vl_src() == ec2.vl_src() && ec1.el() == ec2.el() && ec1.vl_dst() < ec2.vl_dst())
		return true;
	    else
		return false;
	}
	else
	    return false;
    }

    struct EdgeCodeCmpLex
    {
	bool operator() (const EdgeCode& ec1, const EdgeCode& ec2) const __attribute__((noinline));
    };

    bool EdgeCodeCmpLex::operator() (const EdgeCode& ec1, const EdgeCode& ec2) const
    {
#ifdef USE_ASM
        uint8_t res;
        asm (
            "movw 8(%1), %%ax\n\t"
            "subw 8(%2), %%ax\n\t"
            "movl 4(%1), %%eax\n\t"
            "sbbl 4(%2), %%eax\n\t"
            "movl (%1), %%eax\n\t"
            "sbbl (%2), %%eax\n\t"
            "setc %0"
            : "=r" (res)
            : "r" (ec1.x_), "r" (ec2.x_)
            : "eax"
            );
        return res;
#else
        if (ec1.vi_src() < ec2.vi_src())
            return true;
        if (ec1.vi_src() > ec2.vi_src())
            return false;
        if (ec1.vi_dst() < ec2.vi_dst())
            return true;
        if (ec1.vi_dst() > ec2.vi_dst())
            return false;
        if (ec1.vl_src() < ec2.vl_src())
            return true;
        if (ec1.vl_src() > ec2.vl_src())
            return false;
        if (ec1.el() < ec2.el())
            return true;
        if (ec1.el() > ec2.el())
            return false;
        if (ec1.vl_dst() < ec2.vl_dst())
            return true;
        else
            return false;
#endif
    }


    bool EdgeCode::operator== (const EdgeCode& ec) const
    {
        return !
            ((vi_src()^ec.vi_src()) | (vi_dst()^ec.vi_dst()) |
             (vl_src()^ec.vl_src()) | (vl_dst()^ec.vl_dst()) | (el()^ec.el()));
    }

    std::ostream& operator<<(std::ostream& out, const EdgeCode& ec)
    {
        return out<<"("<<ec.vi_src()<<","<<ec.vi_dst()<<", "<<ec.vl_src()<<","<<ec.el()<<","<<ec.vl_dst()<<")"
                  << (ec.is_forward() ? " fwd" : " bck");
    }

    // *****************************************************************************
    //                          DFSCode
    // *****************************************************************************
    DfscVI max_vertex(const DFSCode& dfsc)
    {
        DfscVI m = 0;
        for (DFSCode::const_iterator i = dfsc.begin(); i != dfsc.end(); ++i)
            m = std::max(m, std::max(i->vi_src(), i->vi_dst()));
        return m;
    }

    std::ostream& operator<<(std::ostream& out, const DFSCode& dfsc)
    {
        std::copy(dfsc.begin(), dfsc.end(), std::ostream_iterator<EdgeCode>(out, "\n"));
        return out;
    }

    // *****************************************************************************
    //                          RMPath
    // *****************************************************************************

    typedef std::vector<int, STL_Allocator<int> >	RMPath_I_EdgeCode;
    typedef std::vector<DfscVI, STL_Allocator<DfscVI> >	DfscVI_Vertices;

    // ----------------------------------------------------
    // dfsc_vertices[0        ... retvalue )                    : rmpath vertices
    // dfsc_vertices[retvalue ... dfsc_vertices.size() )        : notrmpath vertices
    //
    // return: rmpath vertices number
    //
    std::size_t
    make_rmpath(RMPath_I_EdgeCode& rmpath_i_ec, DfscVI_Vertices& dfsc_vertices, const DFSCode& dfsc)
    {
        // until first forward edge
        signed int idx = dfsc.size();
        while (dfsc[--idx].is_backward());
        
        // first forward edge is the right most edge
        rmpath_i_ec.push_back(idx);

        // right most vertex is the greatest
        // we get number of all vertices
        const std::size_t n_all = dfsc[idx].vi_dst() + 1;

	dfsc_vertices.resize(n_all, VI_NULL);
        dfsc_vertices[0] = dfsc[idx].vi_dst();
        dfsc_vertices[1] = dfsc[idx].vi_src();
        std::size_t n_rmp = 2;
        std::size_t n_not_rmp = 0;

        
        DfscVI old_src = dfsc[idx].vi_src();
        while (--idx >= 0)
        {
            if (dfsc[idx].is_forward())
            {
                if (old_src == dfsc[idx].vi_dst())
                {
                    rmpath_i_ec.push_back(idx);
                    dfsc_vertices[n_rmp++] = dfsc[idx].vi_src();
                    old_src = dfsc[idx].vi_src();
                }
                else
                {
                    dfsc_vertices[n_all - ++n_not_rmp] = dfsc[idx].vi_dst();
                }
            }
        }

        assert(n_rmp + n_not_rmp == n_all);
        
        return n_rmp;
    }

    class RMPath : private boost::noncopyable
    {
        RMPath_I_EdgeCode       rmpath_i_ec_;
        DfscVI_Vertices         dfsc_vertices_;
        std::size_t             n_rmp_;
    public:
        RMPath(const DFSCode& dfsc, MemAllocator* mem_alloc)
	    :rmpath_i_ec_(STL_Allocator<int>(mem_alloc)),
	     dfsc_vertices_(STL_Allocator<DfscVI>(mem_alloc))
	    {
		rmpath_i_ec_.reserve(dfsc.size());
		n_rmp_ = make_rmpath(rmpath_i_ec_, dfsc_vertices_, dfsc);
	    }
        
        typedef DfscVI_Vertices::const_iterator const_iterator;
        
        const_iterator begin() const            { return dfsc_vertices_.begin(); }
        const_iterator rmp_end() const          { return dfsc_vertices_.begin() + n_rmp_; }
        const_iterator end() const              { return dfsc_vertices_.end(); }

        std::size_t num_all_vertices() const    { return dfsc_vertices_.size(); }
        std::size_t num_rmp_vertices() const    { return n_rmp_; }

        int operator[] (int i) const    { return rmpath_i_ec_[i]; }
        int num_edges() const		{ return rmpath_i_ec_.size(); }
        int rightmost_edgeindex() const { return rmpath_i_ec_[0]; }

        DfscVI rightmost_vertex() const { return dfsc_vertices_[0]; }

        friend std::ostream& operator<<(std::ostream& out, const RMPath& r);
    };

    std::ostream& operator<<(std::ostream& out, const RMPath& r)
    {
        out << "rmpath edge indexes: ";
        std::copy(r.rmpath_i_ec_.begin(), r.rmpath_i_ec_.end(), std::ostream_iterator<int>(out, " "));
        out << "rmpath vertices    : ";
        std::copy(r.begin(), r.rmp_end(), std::ostream_iterator<DfscVI>(out, " "));
        out << "not rmpath vertices: ";
        std::copy(r.rmp_end(), r.end(), std::ostream_iterator<DfscVI>(out, " "));
        return out;
    }


    // *****************************************************************************
    //                          SBGSimple
    // *****************************************************************************
    class SBGSimple : public SBGBase<SBGSimple>
    {
	friend class SBGCreator<SBGSimple>;
	SBGSimple(const Graph::Edge& e, const Graph* g)
	    :SBGBase<SBGSimple>(e, g)
	    {}

	SBGSimple(const Graph::Edge& e, const SBGSimple* s)
	    :SBGBase<SBGSimple>(e, s)
	    {}
    };


    // *****************************************************************************
    //                          SBG
    // *****************************************************************************

    void SBG::init_dfsc_to_graph_array1(MemAllocator* mem_alloc)
    {
	vi_dfsc_to_graph_ = mem_alloc->alloc_array<GraphVI>(2);
	vi_dfsc_to_graph_[0] = edge().vi_src();
        vi_dfsc_to_graph_[1] = edge().vi_dst();
    }


    void SBG::init_dfsc_to_graph_array2(MemAllocator* mem_alloc)
    {
	PREFETCH(*parent()->vi_dfsc_to_graph_);
	vi_dfsc_to_graph_ = mem_alloc->alloc_array<GraphVI>(num_vertices());
	::memcpy(vi_dfsc_to_graph_,
		 parent()->vi_dfsc_to_graph_,
		 parent()->num_vertices() * sizeof(GraphVI));
	vi_dfsc_to_graph_[vi_src_dfsc_] = edge().vi_src();
        vi_dfsc_to_graph_[vi_dst_dfsc_] = edge().vi_dst();
    }


    //
    // array of the dfsc VI, indexed by the graph VI
    // so, sbg_vertices[graph_vi] == dfsc vi
    //
    DfscVI* SBG::create_graph_to_dfsc_v(MemAllocator* mem_alloc, DfscVI vi_default)
    {
        std::size_t n = get_graph()->num_vertices();
        DfscVI* p = mem_alloc->alloc_array<DfscVI>(n);
        std::fill(p, p + n, vi_default);

        const SBG* s = this;
        do
        {
            p[s->edge().vi_src()] = s->vi_src_dfsc_;
            p[s->edge().vi_dst()] = s->vi_dst_dfsc_;
            s = s->parent();
        } while (s);

        return p;
    }

    void SBG::free_graph_to_dfsc_v(DfscVI* array, MemAllocator* mem_alloc)
    {
	mem_alloc->dealloc_array(array, get_graph()->num_vertices());
    }


    void SBG::insert_to_automorph_list(SBG* pos, SBG* s)
    {
        pos->automorph_next_->automorph_prev_ = s;
        s->automorph_next_ = pos->automorph_next_;
        s->automorph_prev_ = pos;
        pos->automorph_next_ = s;
    }

    bool operator== (const SBG& s1, const SBG& s2)
    {
        assert(s1.get_graph() == s2.get_graph());

        if (s1.sum_ != s2.sum_)
            return false;
        unsigned int n = s1.get_graph()->num_vertices();
        for (unsigned int i = 0; i < n; ++i)
	{
	    if (s1.ev_array()[i].e_in_sbg != s2.ev_array()[i].e_in_sbg)
		return false;
	}
        return true;
    }

    std::ostream& operator<<(std::ostream& out, const SBG& sbg)
    {
        std::vector<const SBG*> chain;
        get_chain(chain, &sbg);

        out << "sbg:";
        for (std::size_t i = 0; i < chain.size(); ++i)
            out << " " << chain[i]->edge();	
        //out << " at address=" << &sbg << " parent=" << sbg.parent();
        return out;
    }


    void get_chain(std::vector<const SBG*>& chain, const SBG* sbg)
    {
        chain.resize(sbg->num_edges());
        int i = 0;
        const SBG* s = sbg;
        do
        {
            chain[chain.size() - ++i] = s;
            s = s->parent();
        } while (s);
    }

    // *****************************************************************************
    //                          SBGCreator
    // *****************************************************************************
    template<class S>
    class SBGCreator;

    // ---- creator SBGSimple ------------------
    template<>
    class SBGCreator<SBGSimple>
    {
	MemAllocator*	mem_alloc_;
	FixedAllocator*	sbg_alloc_;
    public:
	explicit SBGCreator(MemAllocator* mem_alloc)
	    :mem_alloc_(mem_alloc), sbg_alloc_(mem_alloc->get_fixed_allocator(sizeof(SBGSimple))) {}

	MemAllocator* mem_allocator() const { return mem_alloc_; }

	SBGSimple* new_sbg(const Graph::Edge& e, const Graph* g);
	SBGSimple* new_sbg(const Graph::Edge& e, const SBGSimple* s);
	void delete_sbg(SBGSimple*);
    };

    // root sbg
    SBGSimple* SBGCreator<SBGSimple>::new_sbg(const Graph::Edge& e, const Graph* g)
    {
	SBGSimple* p = new (sbg_alloc_->allocate()) SBGSimple(e, g);
	p->init_ev_array1(mem_alloc_);
	return p;
    }

    // other sbg
    SBGSimple* SBGCreator<SBGSimple>::new_sbg(const Graph::Edge& e, const SBGSimple* s)
    {
	SBGSimple* p = new (sbg_alloc_->allocate()) SBGSimple(e, s);
	p->init_ev_array2(mem_alloc_);
	return p;
    }

    void SBGCreator<SBGSimple>::delete_sbg(SBGSimple* p)
    {
	p->free_ev_array(mem_alloc_);
	p->~SBGSimple();
	sbg_alloc_->deallocate(p);
    }

    // ---- creator SBG ------------------------
    template<>
    class SBGCreator<SBG>
    {
	MemAllocator*	mem_alloc_;
	FixedAllocator*	sbg_alloc_;
    public:
	explicit SBGCreator(MemAllocator* mem_alloc)
	    :mem_alloc_(mem_alloc), sbg_alloc_(mem_alloc->get_fixed_allocator(sizeof(SBG))) {}

	MemAllocator*	mem_allocator() const { return mem_alloc_; }

	SBG* new_sbg(const Graph::Edge& e, const Graph* g);
	SBG* new_sbg(const Graph::Edge& e, const SBG* s, const EdgeCode& ec);
	void delete_sbg(SBG*);
    };


    SBG* SBGCreator<SBG>::new_sbg(const Graph::Edge& e, const Graph* g)
    {
	SBG* p = new (sbg_alloc_->allocate()) SBG(e, g);
	p->init_ev_array1(mem_alloc_);
	p->init_dfsc_to_graph_array1(mem_alloc_);
	return p;
    }


    SBG* SBGCreator<SBG>::new_sbg(const Graph::Edge& e, const SBG* s, const EdgeCode& ec)
    {
	//PREFETCH(*s);
	SBG* p = new (sbg_alloc_->allocate()) SBG(e, s, ec);
	p->init_ev_array2(mem_alloc_);
	p->init_dfsc_to_graph_array2(mem_alloc_);
	return p;
    }


    void SBGCreator<SBG>::delete_sbg(SBG* p)
    {
	p->free_dfsc_to_graph_array(mem_alloc_);
	p->free_ev_array(mem_alloc_);
	p->~SBG();
	sbg_alloc_->deallocate(p);
    }


    // *****************************************************************************
    //                          SBG_List
    //                          its owner
    // *****************************************************************************
    template<class S>
    class SBG_List
    {
	S* sbg_list_;
	std::size_t size_;
	SBGCreator<S>* sbg_creator_;
    public:
	typedef S SBG_Type;

	SBG_List(SBGCreator<S>* sbg_creator)
	    :sbg_list_(0), size_(0), sbg_creator_(sbg_creator) {}

	SBG_List(const SBG_List& r)
	    :sbg_list_(0), size_(0), sbg_creator_(r.sbg_creator_) { assert(r.empty()); }

	~SBG_List() { clear(); }

	void insert(S* s);
	void clear();
	bool empty() const			{ return !sbg_list_; }
	std::size_t size() const		{ return size_; }
	S* get_first() const			{ return sbg_list_; }
	SBGCreator<S>* sbg_creator() const	{ return sbg_creator_; }
    };

    template<class S>
    void SBG_List<S>::insert(S* s)
    {
	assert(! s->next_embedding_);
	s->next_embedding_ = sbg_list_;
	sbg_list_ = s;
	++size_;
    }


    template<class S>
    void SBG_List<S>::clear()
    {
	S* s = sbg_list_;
	S* next = 0;
	while (s)
	{
	    next = s->next_embedding_;
	    sbg_creator_->delete_sbg(s);
	    s = next;
	}
	sbg_list_ = 0;
	size_ = 0;
    }

    template<class S>
    std::ostream& operator<<(std::ostream& out, const SBG_List<S>& sbgs)
    {
	for (const S* s = sbgs.get_first(); s; s = s->next_embedding())
	    out << '\t' << *s << std::endl;
	return out;
    }

    // *****************************************************************************
    //                          SubgraphsOfOneGraph
    // *****************************************************************************

    void SubgraphsOfOneGraph::calc_support_v(const SBG* slist) const
    {
        int support = std::numeric_limits<int>::max();
	const GraphVI n_ver_g = slist->get_graph()->num_vertices();
	const DfscVI  n_ver_s = slist->num_vertices();

        for (DfscVI vi = 0; vi < n_ver_s; ++vi)
	{
            int n = 0;
            std::vector<bool> vvg(n_ver_g, false);
	    for (const SBG* s = slist; s; s = s->next_embedding())
	    {
		GraphVI graph_vi = s->get_dfsc_to_graph_v()[vi];
		
		assert(graph_vi >= 0 && graph_vi < n_ver_g);
		
		if (! vvg[graph_vi])
                {
                    ++n;
                    vvg[graph_vi] = true;
                }
	    }
            if (n < support)
                support = n;
	}
	
	assert(support > 0);

        support_ = support;
        support_valid_ = true;
    }

    void SubgraphsOfOneGraph::insert_sbgs_uniq_ptr(SBGS_PTR& cont, SBG* s)
    {
        BOOST_FOREACH(SBG* sbg, cont)
        {
            if (*sbg == *s)
            {
                SBG::insert_to_automorph_list(sbg, s);
                return;
            }
        }
        cont.push_back(s);
    }


    bool SubgraphsOfOneGraph::is_equal_occurrence(const SubgraphsOfOneGraph& parent) const
    {
	unsigned int n = 0;
        for (const_iterator i_gr = parent.begin(); i_gr != parent.end(); ++i_gr)
        {
            const SBG* s = *i_gr;
            const SBG* const S_END = s;
            
            do
            {
                if (sbg_parents_.count(s) > 0)
                {
                    ++n;
                    break;
                }
                s = s->next_automorph();
            } while(s != S_END);
        }

        return n == parent.size();
    }


    void SubgraphsOfOneGraph::insert(SBG* s)
    {
        sbg_parents_.insert(s->parent());
        insert_sbgs_uniq_ptr(sbgs_uniq_ptr_, s);
    }

    // *****************************************************************************
    //                          SubgraphsOfManyGraph
    // *****************************************************************************

    bool SubgraphsOfManyGraph::is_equal_occurrence(const SubgraphsOfManyGraph& embd) const
    {
	return 0;
    }

    void SubgraphsOfManyGraph::insert(SBG* s)
    {
	MemAllocator* mem_alloc = g_sog_.get_allocator().get_mem_alloc();
	G_SOG::value_type val(s->get_graph(), SOG(mem_alloc));
	std::pair<iterator,bool> pr = g_sog_.insert(val);
	SOG& sog = pr.first->second;
        std::size_t old_num_sbgs_uniq = sog.size();
        sog.insert(s);
        int new_num_sbgs_uniq = sog.size();
        size_ += new_num_sbgs_uniq - old_num_sbgs_uniq;
    }

    // *****************************************************************************
    //                          Embeddings
    // *****************************************************************************
    template<class SubgraphsOfGraph>
    class Embeddings : public SBG_List<SBG>
    {
	typedef SBG_List<SBG> Base;
	Base& base() { return *this; }
	const Base& base() const { return *this; }
	SubgraphsOfGraph sg_;
    public:
	typedef SBG SBG_Type;

	Embeddings(SBGCreator<SBG>* creator)
	    :SBG_List<SBG>(creator),
	     sg_(creator->mem_allocator())
	    {}
	
	int support() const	{ return sg_.support(base().get_first()); }
	void insert(SBG* s)	{ base().insert(s); sg_.insert(s); }
	const SubgraphsOfGraph& container() const { return sg_; }

	SBGCreator<SBG_Type>* sbg_creator() const { return base().sbg_creator(); }

	template<class T>
	friend std::ostream& operator<<(std::ostream& out, const Embeddings<T>& sbg_collection);
    };

    template<class SubgraphsOfGraph>
    std::ostream& operator<<(std::ostream& out, const Embeddings<SubgraphsOfGraph>& embd)
    {
	out << static_cast<const SBG_List<SBG>&>(embd);
	out << "\tsupport=" << embd.support()
	    << " size=" << embd.size()
	    << " num_uniq=" << embd.sg_.size()
	    << std::endl;
	return out;
    }

    template<class SG>
    inline bool is_equal_occurrence(const Embeddings<SG>& new_embd, const Embeddings<SG>& embd)
    { return new_embd.container().is_equal_occurrence(embd.container()); }


    // *****************************************************************************
    //                          VertexRMPathStatus
    // *****************************************************************************
    class VertexRMPathStatus : private boost::noncopyable
    {
        unsigned char* status_;
	std::size_t size_;
	MemAllocator* mem_alloc_;
    public:

	explicit VertexRMPathStatus(MemAllocator* mem_alloc)
	    :status_(0), size_(0), mem_alloc_(mem_alloc) {}

        VertexRMPathStatus(const RMPath& rmpath, const SBG* sbg, MemAllocator* mem_alloc)
	    :status_(0), size_(0), mem_alloc_(mem_alloc) { create(rmpath, sbg); }

        ~VertexRMPathStatus()
	    {
		mem_alloc_->dealloc_array(status_, size_);
	    }

        void create(const RMPath& rmpath, const SBG* sbg);
        bool operator[] (GraphVI vi) const      { return status_[vi]; }
    };

    void VertexRMPathStatus::create(const RMPath& rmpath, const SBG* sbg)
    {
        assert(status_ == 0);
	size_ = sbg->get_graph()->num_vertices();
	status_ = mem_alloc_->alloc_array<unsigned char>(size_);
        const GraphVI* dfsc_to_graph_v = sbg->get_dfsc_to_graph_v();
        for (RMPath::const_iterator it = rmpath.begin(); it != rmpath.end(); ++it)
        {
            assert(dfsc_to_graph_v[*it] >= 0);
            assert(dfsc_to_graph_v[*it] < size_);
            status_[dfsc_to_graph_v[*it]] = it < rmpath.rmp_end();
        }
    }

    // *****************************************************************************
    //                          Extension
    // *****************************************************************************
    template<class SG, class Cmp>
    class Extension : public std::map<EdgeCode,
				      Embeddings<SG>,
				      Cmp,
				      STL_Allocator<std::pair<EdgeCode, Embeddings<SG> > > >
    {
	typedef STL_Allocator<std::pair<EdgeCode, Embeddings<SG> > > MapAlloc;
	typedef std::map<EdgeCode, Embeddings<SG>, Cmp, MapAlloc> Base;
	Base& base() { return *this; }
	SBGCreator<SBG>* sbg_creator_;
    public:
	explicit
	Extension(MemAllocator* mem_alloc, SBGCreator<SBG>* sbg_creator)
	    :Base(Cmp(), MapAlloc(mem_alloc)),
	     sbg_creator_(sbg_creator) {}

        void insert(const EdgeCode& ec, const Graph::Edge& e, const Graph* g);
        void insert(const EdgeCode& ec, const Graph::Edge& e, const SBG* s);	
    };

    template<class SG, class Cmp>
    void Extension<SG,Cmp>::insert(const EdgeCode& ec, const Graph::Edge& e, const Graph* g)
    {
	typename Base::value_type val(ec, Embeddings<SG>(sbg_creator_));
	std::pair<typename Base::iterator,bool> pr = base().insert(val);
	pr.first->second.insert(sbg_creator_->new_sbg(e, g));
    }

    template<class SG, class Cmp>
    void Extension<SG,Cmp>::insert(const EdgeCode& ec, const Graph::Edge& e, const SBG* s)
    {
	typename Base::value_type val(ec, Embeddings<SG>(sbg_creator_));
	std::pair<typename Base::iterator,bool> pr = base().insert(val);
	pr.first->second.insert(sbg_creator_->new_sbg(e, s, ec));
    }

    template<class SG, class Cmp>
    std::ostream& operator<<(std::ostream& out, const Extension<SG,Cmp>& ext)
    {
	typedef typename Extension<SG,Cmp>::const_iterator I;
	for (I it = ext.begin(); it != ext.end(); ++it)
	    out << it->first << "\n" << it->second;
	return out;
    }

    // *****************************************************************************
    //                          MinimalExtension
    // *****************************************************************************    
    class MinimalExtension : public SBG_List<SBGSimple>
    {
	EdgeCode ec_;
	SBG_List<SBGSimple>& base() { return *this; }
	void insert(const EdgeCode& ec, SBGSimple* s);
    public:
	explicit
	MinimalExtension(SBGCreator<SBGSimple>* sbg_creator)
	    :SBG_List<SBGSimple>(sbg_creator) {}
	
        void insert(const EdgeCode& ec, const Graph::Edge& e, const Graph* g)
	    { insert(ec, base().sbg_creator()->new_sbg(e, g)); }

	void insert(const EdgeCode& ec, const Graph::Edge& e, const SBGSimple* s)
	    { insert(ec, base().sbg_creator()->new_sbg(e, s)); }

	const EdgeCode& get_ec() const { return ec_; }
    };

    void MinimalExtension::insert(const EdgeCode& ec, SBGSimple* s)
    {
	if (base().empty())
	{
	    ec_ = ec;
	    base().insert(s);
	}
	else if (ec == ec_)
	{
	    // equal
	    base().insert(s);
	}
	else
	{
	    EdgeCodeCmpDfs cmp;
	    if (cmp(ec, ec_))
	    {
		// lesser
		ec_ = ec;
		base().clear();
		base().insert(s);
	    }
	    else
	    {
		// greater
		base().sbg_creator()->delete_sbg(s);
	    }
	}
    }


#ifdef GSPAN_WITH_STATISTICS
    // *****************************************************************************
    //                         Statistics
    // *****************************************************************************
    struct Statistics
    {
        unsigned int num_project_calls;
        unsigned long num_ismin_calls;
        unsigned long num_ismin_true_ret;
        unsigned long num_detect_early_termin_x_f;
        unsigned long num_detect_early_termin_x_b;
        unsigned long num_detect_early_termin_r_f;
        unsigned long num_detect_early_termin_r_b;
        Statistics()
            :num_project_calls(0),
             num_ismin_calls(0),
             num_ismin_true_ret(0),
             num_detect_early_termin_x_f(0),
             num_detect_early_termin_x_b(0),
             num_detect_early_termin_r_f(0),
             num_detect_early_termin_r_b(0)
            {}
        friend std::ostream& operator<< (std::ostream& strm, const Statistics& s);
    };

    std::ostream& operator<< (std::ostream& strm, const Statistics& s)
    {
        float ismin_true_ret_prc = s.num_ismin_true_ret;
        ismin_true_ret_prc /= (float)s.num_ismin_calls;
        int et_total =
            s.num_detect_early_termin_x_f + s.num_detect_early_termin_x_b +
            s.num_detect_early_termin_r_f + s.num_detect_early_termin_r_b;
        strm << "--------- statistics ------------------------";
        strm << "\nnum_project_calls           = " << s.num_project_calls
             << "\nnum_ismin_calls             = " << s.num_ismin_calls
             << "\nismin_true_ret              = " << s.num_ismin_true_ret
             << "\nismin_true_ret%             = " << ismin_true_ret_prc*100.0f
             << "\nnum_detect_early_termin_x_f = " << s.num_detect_early_termin_x_f
             << "\nnum_detect_early_termin_x_b = " << s.num_detect_early_termin_x_b
             << "\nnum_detect_early_termin_r_f = " << s.num_detect_early_termin_r_f
             << "\nnum_detect_early_termin_r_b = " << s.num_detect_early_termin_r_b
             << "\nnum_detect_early_termin     = " << et_total
             << std::endl;
        strm << "---------------------------------------------" << std::endl;
        return strm;
    }
#endif


    // *****************************************************************************
    //                         SharedData
    // *****************************************************************************
    enum ExtType { EXT_X, EXT_R, EXT_NONE };
    
    struct SharedData
    {
        DFSCode dfsc;
 
	MemAllocator		mem_alloc;
	SBGCreator<SBG>		sbg_creator;
	SBGCreator<SBGSimple>	sbgsimple_creator;

	const int minsup_;
        GspanResult* result_;
 
#ifdef GSPAN_TRACE
        unsigned int max_trace_depth_;
#endif

#ifdef GSPAN_WITH_STATISTICS
        Statistics statistics_;
        void on_enter_frame()           { ++statistics_.num_project_calls; }
        unsigned int num_calls() const  { return statistics_.num_project_calls; }
        bool on_ismin(bool ismin)
            { ++statistics_.num_ismin_calls; if (ismin) ++statistics_.num_ismin_true_ret; return ismin; }
        void on_early_term_x_f()        { ++statistics_.num_detect_early_termin_x_f; }
        void on_early_term_x_b()        { ++statistics_.num_detect_early_termin_x_b; }
        void on_early_term_r_f()        { ++statistics_.num_detect_early_termin_r_f; }
        void on_early_term_r_b()        { ++statistics_.num_detect_early_termin_r_b; }
        void on_early_term_f(ExtType t) { if (t == EXT_X) on_early_term_x_f(); else on_early_term_r_f(); }
        void on_early_term_b(ExtType t) { if (t == EXT_X) on_early_term_x_b(); else on_early_term_r_b(); }
        void statistics_report() const  { std::cerr << statistics_; }
#else
        bool on_ismin(bool ismin)	{ return ismin; }
        void on_early_term_x_f()        {}
        void on_early_term_x_b()        {}
        void on_early_term_r_f()        {}
        void on_early_term_r_b()        {}
        void on_early_term_f(ExtType t) {}
        void on_early_term_b(ExtType t) {}
        void statistics_report() const  {}
#endif

#if not defined(GSPAN_WITH_STATISTICS) && (defined(GSPAN_TRACE) || defined(CHECK_MODE))
        unsigned int num_project_calls;
        void on_enter_frame()           { ++num_project_calls; }
        unsigned int num_calls() const  { return num_project_calls; }
#endif
#if not defined(GSPAN_WITH_STATISTICS) && not defined(GSPAN_TRACE) && not defined(CHECK_MODE)
        void on_enter_frame()           {}
#endif

        SharedData(int minsup, GspanResult* result, unsigned int max_trace_depth)
            : sbg_creator(&mem_alloc)
	    , sbgsimple_creator(&mem_alloc)
	    , minsup_(minsup)
            , result_(result)
#ifdef GSPAN_TRACE
            , max_trace_depth_(max_trace_depth)
#endif

#if not defined(GSPAN_WITH_STATISTICS) && (defined(GSPAN_TRACE) || defined(CHECK_MODE))
            , num_project_calls(0)
#endif
            {}

        ~SharedData() { statistics_report(); }
    };


    // *****************************************************************************
    //                         FrameState
    // *****************************************************************************
    template<class SG>
    struct FrameState
    {
	typedef Embeddings<SG> Embd;
	const Embd* embd;

        RMPath rmpath;

	typedef Extension<SG, EdgeCodeCmpDfs> REdges;
	typedef Extension<SG, EdgeCodeCmpLex> XEdges;
        typedef typename REdges::const_iterator REcIter;
        typedef typename XEdges::const_iterator XEcIter;
        typedef typename REdges::iterator REIter;
        typedef typename XEdges::iterator XEIter;
        XEdges x_edges;
        REdges r_edges;

        bool closed;
        bool early_term;

        // cause of not close
        ExtType exttype_notclose;
        const EdgeCode* ec_notclose;
        const Embd* embd_ext_notclose;

        // cause of early_termin
        ExtType exttype_early_term;
        const EdgeCode* ec_early_term;
        const Embd* embd_ext_early_term;

	typedef std::pair<const EdgeCode*, const Embd*> RChild;
        typedef std::vector<RChild, STL_Allocator<RChild> > RChildren;
        typedef typename RChildren::const_iterator RChildIter;
        RChildren children;
	
        FrameState* prev_state;
        const std::size_t depth;

#if defined(GSPAN_TRACE) || defined(CHECK_MODE)
        unsigned int id;
#endif

#if defined(CHECK_MODE)
        const bool TERMINATED_BRANCH;
        const FrameState* early_term_frame;
#endif

	FrameState(SharedData* shared, FrameState* prev, const Embd* embd);
    };


    template<class SG>
    FrameState<SG>::FrameState(SharedData* shared, FrameState* prev, const Embd* embeddings)
	: embd(embeddings)
	, rmpath(shared->dfsc, &shared->mem_alloc)
	, x_edges(&shared->mem_alloc, &shared->sbg_creator)
	, r_edges(&shared->mem_alloc, &shared->sbg_creator)
	, closed(true)
	, early_term(false)
	, exttype_notclose(EXT_NONE), ec_notclose(0), embd_ext_notclose(0)
	, exttype_early_term(EXT_NONE), ec_early_term(0), embd_ext_early_term(0)
	, children(STL_Allocator<RChild>(&shared->mem_alloc))
	, prev_state(prev)
	, depth(shared->dfsc.size())
#if defined(GSPAN_TRACE) || defined(CHECK_MODE)
	, id(shared->num_calls())
#endif
#if defined(CHECK_MODE)
        , TERMINATED_BRANCH(prev && prev->early_term_frame != 0)
        , early_term_frame(prev ? prev->early_term_frame : 0)
#endif
    {
#if defined(CHECK_MODE)
        if (prev)
            assert(TERMINATED_BRANCH ? prev->early_term_frame!=0 : prev->early_term_frame==0);
#endif
    }

    // *****************************************************************************
    //                         functions
    // *****************************************************************************
    template<class Ext>
    void enum_one_edges(Ext& ext, const Graph& g)
    {
        const Graph::Edges& g_edges = g.edges();
        for (Graph::EdgesIterator it = g_edges.begin(); it != g_edges.end(); ++it)
            ext.insert(EdgeCode(0, 1, it->vl_src(), it->el(), it->vl_dst(), true), *it, &g);
    }

    // -----------------------------------
    //	dfscode minimality check
    // -----------------------------------

    void enumerate_min_bck(MinimalExtension& ext,
			   const DFSCode& dfsc_min,
                           const SBG_List<SBGSimple>& elist,
			   const RMPath& rmpath,
			   const Graph& g)
    {
	const DfscVI vi_dfsc_rmost = dfsc_min[rmpath.rightmost_edgeindex()].vi_dst();
	const VL vl_rmost = dfsc_min[rmpath.rightmost_edgeindex()].vl_dst();

	// from first vertex toward right most vertex
	for (int i = rmpath.num_edges() - 1; ext.empty() && i > 0; --i)
	{
            const EdgeCode& ec_rmpath = dfsc_min[rmpath[i]];
            const EL el_rmpath = ec_rmpath.el();
            const bool vl_less_eq = ec_rmpath.vl_dst() <= vl_rmost;

	    for (SBGSimple* s = elist.get_first(); s; s = s->next_embedding())
	    {
		const SBGSimple* s_rmost = parent(s, rmpath.rightmost_edgeindex() + 1);
		const Graph::Edge* e_rmost = &s_rmost->edge();
		const GraphVI graph_vi = parent(s_rmost, rmpath[i] + 1)->edge().vi_src();

		// const Graph::Edge* e_rmost = (*s)[rmpath.rightmost_edgeindex()];
                // const GraphVI graph_vi = (*s)[rmpath[i]]->vi_src();

                if (s->has_no_extension(graph_vi))
                    continue;
                s->set_no_has_extension(graph_vi);

		const Graph::IncidentEdges& incid_edges = g.incident(graph_vi);
                Graph::IncidentEdgesIterator it = incid_edges.begin();
                Graph::IncidentEdgesIterator it_end = incid_edges.end();
                for (; it != it_end; ++it)
                {
                    const Graph::Edge* pe = *it;
                    if (s->has_edge(pe->eid()))
                        continue;
		    s->set_has_extension(graph_vi);
		    
		    if (pe->vi_dst() == e_rmost->vi_dst() &&
                        ((vl_less_eq && el_rmpath == pe->el()) || el_rmpath < pe->el()))
                    {
                        Graph::Edge e = *pe;
                        e.chgdir();
                        EdgeCode ec(vi_dfsc_rmost, ec_rmpath.vi_src(), vl_rmost, e.el(), ec_rmpath.vl_src(), false);
                        ext.insert(ec, e, s);
                        break;
                    }
		}
	    }
	}
    }
 
    void enumerate_min_fwd(MinimalExtension& ext,
			   const DFSCode& dfsc_min,
                           const SBG_List<SBGSimple>& elist,
			   const RMPath& rmpath,
			   const Graph& g)
    {
        const DfscVI vi_dfsc_rmost = dfsc_min[rmpath.rightmost_edgeindex()].vi_dst();
        const VL vl_rmost = dfsc_min[rmpath.rightmost_edgeindex()].vl_dst();
        const VL vl_minimum = dfsc_min[0].vl_src();

	// forward pure
	for (SBGSimple* s = elist.get_first(); s; s = s->next_embedding())
	{
	    const SBGSimple* s_rmost = parent(s, rmpath.rightmost_edgeindex() + 1);
	    const Graph::Edge* e_rmost = &s_rmost->edge();
	    GraphVI graph_vi = e_rmost->vi_dst();

            // const Graph::Edge* e_rmost = (*s)[rmpath.rightmost_edgeindex()];
            // GraphVI graph_vi = e_rmost->vi_dst();

            if (s->has_no_extension(graph_vi))
                continue;
            s->set_no_has_extension(graph_vi);

            const Graph::IncidentEdges& incid_edges = g.incident(graph_vi);
            Graph::IncidentEdgesIterator it = incid_edges.begin();
            Graph::IncidentEdgesIterator it_end = incid_edges.end();
            for (; it != it_end; ++it)
            {
                const Graph::Edge* e = *it;
                if (s->has_edge(e->eid()))
                    continue;
                s->set_has_extension(graph_vi);

                if (! s->has_vertex(e->vi_dst()) && vl_minimum <= e->vl_dst())
                {
                    EdgeCode ec(vi_dfsc_rmost, vi_dfsc_rmost+1, vl_rmost, e->el(), e->vl_dst(), true);
                    ext.insert(ec, *e, s);
                }
            }
	}

	// forward rmpath
	// from right most vertex toward first vertex
	for (int i = 0; ext.empty() && i < rmpath.num_edges(); ++i)
	{
	    const EdgeCode& ec_rmpath = dfsc_min[rmpath[i]];
	    for (SBGSimple* s = elist.get_first(); s; s = s->next_embedding())
	    {
		const GraphVI graph_vi = parent(s, rmpath[i] + 1)->edge().vi_src();
		//GraphVI graph_vi = (*s)[rmpath[i]]->vi_src();

		if (s->has_no_extension(graph_vi))
		    continue;
		s->set_no_has_extension(graph_vi);

                const Graph::IncidentEdges& incid_edges = g.incident(graph_vi);
                Graph::IncidentEdgesIterator it = incid_edges.begin();
                Graph::IncidentEdgesIterator it_end = incid_edges.end();
                for (; it != it_end; ++it)
                {
                    const Graph::Edge* e = *it;
                    if (s->has_edge(e->eid()))
                        continue;
		    s->set_has_extension(graph_vi);

                    if (! s->has_vertex(e->vi_dst()) && vl_minimum <= e->vl_dst() &&
                        ((ec_rmpath.vl_dst() <= e->vl_dst() && ec_rmpath.el() == e->el()) || ec_rmpath.el() < e->el()))
                    {
                        EdgeCode ec(ec_rmpath.vi_src(), vi_dfsc_rmost+1, ec_rmpath.vl_src(), e->el(), e->vl_dst(), true);
                        ext.insert(ec, *e, s);
                    }
                }
	    }
	}
    }
    

    bool is_min_iterative(SharedData* shared)
    {
	const DFSCode& dfsc_tested = shared->dfsc;
	Graph graph(dfsc_tested.begin(), dfsc_tested.end());

	typedef STL_Allocator<MinimalExtension> ExtAllocator;
	std::list<MinimalExtension, ExtAllocator> exts(1,
						       MinimalExtension(&shared->sbgsimple_creator),
						       ExtAllocator(&shared->mem_alloc));
	MinimalExtension* exts1 = &exts.back();
	
	enum_one_edges(*exts1, graph);

	DFSCode dfsc_min;
	dfsc_min.reserve(dfsc_tested.size());

	while (true)
	{
	    dfsc_min.push_back(exts1->get_ec());
            if (dfsc_min[dfsc_min.size()-1] != dfsc_tested[dfsc_min.size()-1])
                return false;
	    
	    RMPath rmpath(dfsc_min, &shared->mem_alloc);
	    
	    exts.push_back(MinimalExtension(&shared->sbgsimple_creator));
	    MinimalExtension* exts2 = &exts.back();

	    enumerate_min_bck(*exts2, dfsc_min, *exts1, rmpath, graph);
	    if (! exts2->empty())
            {
                exts1 = exts2;
                continue;
            }

	    enumerate_min_fwd(*exts2, dfsc_min, *exts1, rmpath, graph);
	    if (! exts2->empty())
            {
                exts1 = exts2;
                continue;
            }
	    
	    return true;
	}
    }


    bool is_min_fast(SharedData* shared, EdgeCode ec, const RMPath& rmpath,
		     bool& last_result, EdgeCode& ec_last)
    {
        EdgeCode f;
        if (ec.vi_src() == ec_last.vi_src() && ec.is_forward() && last_result && ec_last.is_forward())
            return true;
        else
        {
            last_result = is_min_iterative(shared);
            shared->on_ismin(last_result);
            return last_result;
        }
    }


    template<class P>
    bool is_min(SharedData* shared, FrameState<P>* frame,
		const EdgeCode& ec, bool& last_result, EdgeCode& ec_last)
    {
        shared->dfsc.push_back(ec);
        bool r = is_min_fast(shared, ec, frame->rmpath, last_result, ec_last);
        if (r)
            ec_last = ec;
        shared->dfsc.pop_back();
        return r;
    }


    // -----------------------------------
    //	trace, check and debug functions
    // -----------------------------------

#ifdef GSPAN_TRACE
    template<class SG>
    void print_frame_trace(const FrameState<SG>* frame, const SharedData* shared)
    {
        using namespace std;
        cerr << shared->dfsc.size() << ":";
        for (unsigned int i = 0; i < shared->dfsc.size(); ++i)
            cerr << " ";
        cerr << shared->dfsc.back()
             << " support=" << frame->embd->support()
             << " size=" << frame->embd->size()
             << " num_uniq=" << frame->embd->container().size()
             << " frame=" << frame->id;
#ifdef CHECK_MODE
        if (frame->TERMINATED_BRANCH)
        {
            assert(frame->early_term_frame && frame->early_term_frame != frame);
            cerr << "; Was Terminated at frame="
		 << frame->early_term_frame->id
		 << " with depth=" << frame->early_term_frame->depth;
        }
#endif
        if (frame->early_term)
        {
            assert(frame->exttype_early_term == EXT_X || frame->exttype_early_term == EXT_R);
            cerr << "; Detect Early Termination: "
		 << (frame->exttype_early_term == EXT_X ? "X" : "R")
		 << *frame->ec_early_term;
        }
        cerr << endl;
    }
#endif

#ifdef CHECK_MODE
    template<class SG>
    bool is_child(const FrameState<SG>* frame, const EdgeCode* ec)
    {
	typedef typename FrameState<SG>::RChildren RChildren;
        typedef typename RChildren::const_iterator RChildIter;
        for (RChildIter it = frame->children.begin(); it != frame->children.end(); ++it)
            if (it->first == ec)
                return true;
        return false;
    }
        
    template<class SG>
    void print_fail_et_report(const FrameState<SG>* frame, const SharedData* shared)
    {
	typedef typename FrameState<SG>::REdges REdges;
	typedef typename FrameState<SG>::XEdges XEdges;
        typedef typename REdges::const_iterator REcIter;
        typedef typename XEdges::const_iterator XEcIter;
        using namespace std;

        cerr << "Missed pattern detected" << endl;

        cerr << "============ MISSED FRAME: " << frame->id << " =======================================" << endl;
        cerr << "DFSC:\n";
        cerr << "--------------------------\n";
        for (unsigned int i = 0; i < frame->depth; ++i)
            cerr << i+1 << ":\t"
                 << (find(frame->rmpath.begin(),frame->rmpath.end(), i) != frame->rmpath.end() ? "*" : " ")
                 << shared->dfsc[i]
                 << endl;
        cerr << "--------------------------\n";
        cerr << "Embeddings:\n";
        cerr << *frame->embd;
        cerr << "--------------------------\n";
        cerr << "XExtentions:\n";
            
        for (XEcIter it = frame->x_edges.begin(); it != frame->x_edges.end(); ++it)
            cerr << it->first << "\n" << it->second;
            
        cerr << endl;
        const FrameState<SG>* fet_frame = frame->early_term_frame;
        cerr << "============ FAIL Eearly Termination FRAME: " << fet_frame->id << " ======================" << endl;
        cerr << "DFSC:\n";
        cerr << "--------------------------\n";
        for (unsigned int i = 0; i < fet_frame->depth; ++i)
            cerr << i+1 << ":\t"
                 << (find(fet_frame->rmpath.begin(),fet_frame->rmpath.end(), i) != fet_frame->rmpath.end() ? "*" : " ")
                 << shared->dfsc[i]
                 << endl;
        cerr << "--------------------------\n";
        cerr << "Embeddings:\n";
        cerr << *fet_frame->embd;
        cerr << "--------------------------\n";
        cerr << "Cause Early Termination:\n";
        cerr << (fet_frame->exttype_early_term == EXT_X ? "X" : "R") << *fet_frame->ec_early_term << endl;
        cerr << *fet_frame->embd_ext_early_term;
        cerr << "--------------------------\n";
        cerr << "RExtensions:\n";
        const EdgeCode& ec_ext = shared->dfsc[fet_frame->depth];
        for (XEcIter it = fet_frame->r_edges.begin(); it != fet_frame->r_edges.end(); ++it)
        {
            bool ch = it->first == ec_ext;
            cerr << it->first
                 << (is_child(fet_frame, &it->first) ? " * " : " ")
                 << (ch ? "CHILD" : " ")
                 << "\n" << it->second;
            if (ch)
                break;
        }
        cerr << endl;
    }
#endif

    template<class SG>
    void debug_print(const FrameState<SG>* frame, const SharedData* shared)
    {
        std::cerr << "===================================================================\n";
        std::cerr << "Frame: " << shared->num_calls() << std::endl;
        std::cerr << "DFSC:\n"
                  << shared->dfsc;
        std::cerr << "REDGES---------------------------\n" << frame->r_edges;
        std::cerr << "XEDGES---------------------------\n" << frame->x_edges;
    }

    template<class SG>
    void trace_frame(FrameState<SG>* frame, SharedData* shared)
    {
#ifdef GSPAN_TRACE
        if (frame->depth < shared->max_trace_depth_)
            print_frame_trace(frame, shared);
#endif
    }


    template<class SG>
    void result(FrameState<SG>* frame, SharedData* shared)
    {
        if (frame->closed)
        {
	    // calc support if nessesary
	    frame->embd->support();

            (*shared->result_)(shared->dfsc, frame->embd->container());
#ifdef CHECK_MODE
            if (frame->TERMINATED_BRANCH)
            {
                std::cerr << "# !!! MISSED !!!" << std::endl;
                print_fail_et_report(frame, shared);
                exit(1);
            }
#endif
        }
    }


    // -----------------------------------
    //	closegraph
    // -----------------------------------

    template<class SG>
    void enumerate(FrameState<SG>* frame, SharedData* shared)
    {
        typedef typename FrameState<SG>::REdges REdges;
        typedef typename FrameState<SG>::XEdges XEdges;

	MemAllocator* mem_alloc = &shared->mem_alloc;
        REdges& r_edges = frame->r_edges;
        XEdges& x_edges = frame->x_edges;
        const Embeddings<SG>& embd = *frame->embd;
        const RMPath& rmpath = frame->rmpath;
        const DfscVI vi_dfsc_rmost   = rmpath.rightmost_vertex();
        const DfscVI vi_dfsc_new = vi_dfsc_rmost + 1;

	assert(rmpath.rightmost_vertex() == shared->dfsc[rmpath.rightmost_edgeindex()].vi_dst());

	for (SBG* s = embd.get_first(); s; s = s->next_embedding())
	{
	    const Graph* g = s->get_graph();

            //
            // array of the graph VI, indexed by the dfsc VI
            // so, sbg_vertices[dfsc_vi] == graph vi
            //
            const GraphVI* dfsc_to_graph_v = s->get_dfsc_to_graph_v();

            const GraphVI vi_rmost = dfsc_to_graph_v[rmpath.rightmost_vertex()];

            //
            // array of the dfsc VI, indexed by the graph VI
            // so, sbg_vertices[graph_vi] == dfsc vi
            //
            DfscVI* graph_to_dfsc_v = s->create_graph_to_dfsc_v(mem_alloc, vi_dfsc_new);

            VertexRMPathStatus vertex_rmpath_status(rmpath, s, mem_alloc);
      
            for (RMPath::const_iterator dfsc_vi_iter = rmpath.begin();
                 dfsc_vi_iter != rmpath.end(); ++dfsc_vi_iter)
            {
                const DfscVI dfsc_vi = *dfsc_vi_iter;
                const GraphVI graph_vi = dfsc_to_graph_v[dfsc_vi];
                const bool rmpath_vertex = dfsc_vi_iter < rmpath.rmp_end();

                assert(graph_to_dfsc_v[graph_vi] == dfsc_vi);
                assert(graph_vi >= 0 && graph_vi < s->get_graph()->num_vertices());
                assert(graph_vi < graph_size(g));

                if (s->has_no_extension(graph_vi))
                    continue;
                s->set_no_has_extension(graph_vi);
                                
                const Graph::IncidentEdges& incid_edges = g->incident(graph_vi);
                Graph::IncidentEdgesIterator it = incid_edges.begin();
                const Graph::IncidentEdgesIterator it_end = incid_edges.end();
                for (; it != it_end; ++it)
                {
                    const Graph::Edge* e = *it;
                    if (s->has_edge(e->eid()))
                        continue;
                    s->set_has_extension(graph_vi);
                    
                    //
                    // R edges will be
                    // IF
                    // 1) vi_src is rmpath vertex AND
                    // 2) vi_dst is new vertex (forward) OR 
                    //     (vi_src is rmost vertex AND vi_dst is any rmpath vertex)
                    // ELSE
                    // X edges
                    // 

                    if (rmpath_vertex)
                    {
                        // vi_src is rmpath vertex

                        if (! s->has_vertex(e->vi_dst()))
                        {
                            // vi_dst is new vertex
                            // R forward
                            EdgeCode ec(dfsc_vi, vi_dfsc_new, e->vl_src(), e->el(), e->vl_dst(), true);
                            r_edges.insert(ec, *e, s);
                        }
                        else if (e->vi_src() == vi_rmost && vertex_rmpath_status[e->vi_dst()])
                        {
                            // vi_src is rmost vertex AND vi_dst is any rmpath vertex
                            // R backward
                            DfscVI dfsc_vi_dst = graph_to_dfsc_v[e->vi_dst()];
                            assert(dfsc_vi_dst != vi_dfsc_new);
                            EdgeCode ec(dfsc_vi, dfsc_vi_dst, e->vl_src(), e->el(), e->vl_dst(), false);
                            r_edges.insert(ec, *e, s);
                        }
                        else
                        {
                            if (! (e->vi_dst() == vi_rmost && vertex_rmpath_status[e->vi_src()]))
                            {
                                // X backward
                                DfscVI dfsc_vi_dst = graph_to_dfsc_v[e->vi_dst()];
                                assert(dfsc_vi_dst != vi_dfsc_new);
                                EdgeCode ec(dfsc_vi, dfsc_vi_dst, e->vl_src(), e->el(), e->vl_dst(), false);
                                x_edges.insert(ec, *e, s);
                            }
                        }
                    }
                    else
                    {
                        // from vertex not on rmpath
                           
                        DfscVI dfsc_vi_dst = graph_to_dfsc_v[e->vi_dst()];
                        EdgeCode ec(dfsc_vi, dfsc_vi_dst,
                                    e->vl_src(), e->el(), e->vl_dst(), ! s->has_vertex(e->vi_dst()));
                        x_edges.insert(ec, *e, s);
                    }
                }
            }
            
	    s->free_graph_to_dfsc_v(graph_to_dfsc_v, mem_alloc);
	}
    }

    template<class Ext>
    void remove_not_frequents(Ext& ext, int minsup)
    {
        if (minsup > 1)
	{
            for (typename Ext::iterator it = ext.begin(); it != ext.end();)
	    {
                if (it->second.support() < minsup)
                    ext.erase(it++);
                else
                    ++it;
	    }
	}
    }


    bool fail_et(const SBG* sbg, const RMPath& rmpath, MemAllocator* mem_alloc, ExtType exttype)
    {
	const Graph& g = *sbg->get_graph();
        std::vector<char> visited(g.num_vertices(), false);

        VertexRMPathStatus vertex_rmpath_status(mem_alloc);
        if (exttype == EXT_R)
            vertex_rmpath_status.create(rmpath, sbg);
                
        std::queue<GraphVI> q;

        const Graph::IncidentEdges& incid_edges = g.incident(sbg->edge().vi_dst());
        Graph::IncidentEdgesIterator it = incid_edges.begin();
        Graph::IncidentEdgesIterator it_end = incid_edges.end();
        for (; it != it_end; ++it)
        {
            const Graph::Edge* e = *it;
            if (sbg->has_edge(e->eid()))
                continue;
            if (sbg->has_vertex(e->vi_dst()))
            {
                if (exttype == EXT_X || vertex_rmpath_status[e->vi_dst()])
                    return true;
            }
            else
                q.push(e->vi_dst());
        }
        visited[sbg->edge().vi_dst()] = true;

        bool failure = false;
        while (!failure && !q.empty())
        {
            GraphVI vi = q.front();
            q.pop();
            visited[vi] = true;


            const Graph::IncidentEdges& incid_edges = g.incident(vi);
            Graph::IncidentEdgesIterator it = incid_edges.begin();
            Graph::IncidentEdgesIterator it_end = incid_edges.end();
            for (; it != it_end; ++it)
            {
                const Graph::Edge* e = *it;
                GraphVI vi2 = e->vi_dst();
                if (sbg->has_vertex(vi2))
                {
                    if (exttype == EXT_X || vertex_rmpath_status[vi2])
                    {
                        failure = true;
                        break;
                    }
                }
                else
                    if (! visited[vi2])
                        q.push(vi2);
            }
        }
        return failure;
    }


    bool fail_et(const SBG_List<SBG>& elist, const RMPath& rmpath,
		 MemAllocator* mem_alloc, ExtType exttype)
    {
	for (SBG* s = elist.get_first(); s; s = s->next_embedding())
	{
	    if (fail_et(s, rmpath, mem_alloc, exttype))
                return true;
	}
	return false;
    }
   

    template<class P, class ExtIter>
    void detect_nc(FrameState<P>* frame, SharedData* shared, const ExtIter& it, bool equiv, ExtType exttype)
    {
        if (frame->closed && (equiv || it->second.support() == frame->embd->support()))
        {
            frame->closed = false;
            assert(frame->exttype_notclose == EXT_NONE);
            frame->exttype_notclose = exttype;
            frame->ec_notclose = &it->first;
            frame->embd_ext_notclose = &it->second;
        }

    }

    template<class P, class ExtIter>
    void detect_et(FrameState<P>* frame, SharedData* shared, const ExtIter& it, bool equiv, ExtType exttype)
    {
        if (equiv && !frame->early_term
#ifdef CHECK_MODE
            && !frame->TERMINATED_BRANCH
#endif
            )
        {
            if (it->first.is_forward())
            {
                frame->early_term = ! fail_et(it->second, frame->rmpath, &shared->mem_alloc, exttype);
                if (frame->early_term)
                    shared->on_early_term_f(exttype);
            }
            else
            {
                frame->early_term = it->second.size() == frame->embd->size();
                if (exttype == EXT_R)
                    frame->early_term = true;

                if (frame->early_term)
                    shared->on_early_term_b(exttype);
            }

            if (frame->early_term)
            {
                assert(frame->exttype_early_term == EXT_NONE);
                frame->exttype_early_term = exttype;
                frame->ec_early_term = &it->first;
                frame->embd_ext_early_term = &it->second;
            }
        }
    }
    
    template<class SG>
    void project(SharedData* shared, FrameState<SG>* parent_frame, const Embeddings<SG>* embd)
    {
        typedef typename FrameState<SG>::REdges REdges;
        typedef typename FrameState<SG>::XEdges XEdges;
        typedef typename REdges::const_iterator REcIter;
        typedef typename XEdges::const_iterator XEcIter;
        typedef typename FrameState<SG>::RChild RChild;
        typedef typename FrameState<SG>::RChildren RChildren;

        shared->on_enter_frame();
        FrameState<SG> frame(shared, parent_frame, embd);

        enumerate(&frame, shared);
        remove_not_frequents(frame.r_edges, shared->minsup_);

        for (XEcIter it = frame.x_edges.begin(); it != frame.x_edges.end(); ++it)
        {
            bool equiv = is_equal_occurrence(it->second, *frame.embd);
            detect_nc(&frame, shared, it, equiv, EXT_X);
            detect_et(&frame, shared, it, equiv, EXT_X);
            if (frame.early_term && frame.closed)
                break;
        }

#ifdef CHECK_MODE
        if (! frame.TERMINATED_BRANCH && frame.early_term)
            frame.early_term_frame = &frame;
#endif

        frame.children.reserve(frame.r_edges.size());

        bool is_min_last = false;
        EdgeCode ec_last;
#ifdef CHECK_MODE
        for (REcIter it = frame.r_edges.begin(); it != frame.r_edges.end(); ++it)
        {
            bool equiv = is_equal_occurrence(it->second, *frame.embd);
            detect_nc(&frame, shared, it, equiv, EXT_R);
            if (is_min(shared, &frame, it->first, is_min_last, ec_last))
            {
                detect_et(&frame, shared, it, equiv, EXT_R);
                frame.children.push_back(RChild(&it->first, &it->second));
            }
        }
#else
        for (REcIter it = frame.r_edges.begin(); it != frame.r_edges.end(); ++it)
        {
            bool equiv = is_equal_occurrence(it->second, *frame.embd);
            detect_nc(&frame, shared, it, equiv, EXT_R);
            if (!frame.early_term)
            {
                if (is_min(shared, &frame, it->first, is_min_last, ec_last))
                {
                    detect_et(&frame, shared, it, equiv, EXT_R);
                    frame.children.push_back(RChild(&it->first, &it->second));
                }
            }
            if (frame.early_term && frame.closed)
                break;      
        }
#endif

        trace_frame(&frame, shared);
        result(&frame, shared);
            
        typedef typename RChildren::const_iterator ChIter;
        for (ChIter it = frame.children.begin(); it != frame.children.end(); ++it)
        {
            shared->dfsc.push_back(*it->first);
            project(shared, &frame, it->second);
            shared->dfsc.pop_back();

#ifdef CHECK_MODE
            if (! frame.TERMINATED_BRANCH && it->first == frame.ec_early_term)
                frame.early_term_frame = &frame;
#endif
        }

    }

    template<class SG>
    void run(SharedData* shared, Extension<SG, EdgeCodeCmpDfs>& r_edges)
    {
        for (typename Extension<SG, EdgeCodeCmpDfs>::iterator i = r_edges.begin(); i != r_edges.end(); ++i)
        {
            Embeddings<SG>& embd = i->second;
            if (embd.support() >= shared->minsup_)
            {
                shared->dfsc.push_back(i->first);
                project<SG>(shared, 0, &embd);
                shared->dfsc.pop_back();
            }
        }
    }

    void closegraph(const Graph& graph, int minsup, GspanResult* result, int max_trace_depth)
    {
/*
	std::cerr << "SIZE= " << sizeof(typename FrameState<SubgraphsOfOneGraph>::REdges) << std::endl;
	std::cerr << "SBG SIZE= " << sizeof(SBG) << std::endl;
	std::cerr << std::endl;
	return;
*/
        SharedData shared(minsup, result, max_trace_depth);
        typedef Embeddings<SubgraphsOfOneGraph> Embd;
        Extension<SubgraphsOfOneGraph, EdgeCodeCmpDfs> r_edges(&shared.mem_alloc, &shared.sbg_creator);

        enum_one_edges(r_edges, graph);
        run(&shared, r_edges);
    }


    void closegraph(const std::vector<const Graph*>& graphs, int minsup, GspanResult* result, int max_trace_depth)
    {
        SharedData shared(minsup, result, max_trace_depth);
        typedef Embeddings<SubgraphsOfManyGraph> Embd;
        Extension<SubgraphsOfManyGraph, EdgeCodeCmpDfs> r_edges(&shared.mem_alloc, &shared.sbg_creator);
        for (std::vector<const Graph*>::const_iterator i = graphs.begin(); i != graphs.end(); ++i)
            enum_one_edges(r_edges, **i);
        run(&shared, r_edges);
    }

}
