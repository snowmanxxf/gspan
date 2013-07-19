#include "gspan.hpp"

#include <list>
#include <algorithm>
#include <map>
#include <climits>
#include <cstring>

#ifdef WITH_CHECKS
#include <set>
#include <algorithm>
#endif

namespace gSpan
{
    bool br = false;

    enum ExtType { EXT_X, EXT_R, EXT_NONE };
    
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
    //                          EdgeCode
    // *****************************************************************************
    bool EdgeCode::operator== (const EdgeCode& ec) const
    {
        assert(is_aligned(this));
        assert(is_aligned(&ec));
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
    //                          EdgeCodeCmpDfs
    // *****************************************************************************
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
            else if (ec1.vi_src() == ec2.vi_src() && ec1.vl_src() == ec2.vl_src() &&
                     ec1.el() < ec2.el())
                return true;
            else if (ec1.vi_src() == ec2.vi_src() && ec1.vl_src() == ec2.vl_src() &&
                     ec1.el() == ec2.el() && ec1.vl_dst() < ec2.vl_dst())
                return true;
            else
                return false;
        }
        else
            return false;
    }

    // *****************************************************************************
    //                          EdgeCodeCmpLex
    // *****************************************************************************
    bool EdgeCodeCmpLex::operator() (const EdgeCode& ec1, const EdgeCode& ec2) const
    {
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
    }

    // *****************************************************************************
    //                          DFSCode
    // *****************************************************************************
    std::ostream& operator<<(std::ostream& out, const DFSCode& dfsc)
    {
        std::copy(dfsc.begin(), dfsc.end(), std::ostream_iterator<EdgeCode>(out, "\n"));
        return out;
    }

    // *****************************************************************************
    //                          RMPath
    // *****************************************************************************

    // ----------------------------------------------------
    // dfsc_vertices[0        ... retvalue )                    : rmpath vertices
    // dfsc_vertices[retvalue ... dfsc_vertices.size() )        : notrmpath vertices
    //
    // return: rmpath vertices number
    //
    std::size_t RMPath::make_rmpath(RMPath_I_EdgeCode& rmpath_i_ec,
                                    DfscVI_Vertices& dfsc_vertices,
                                    const EdgeCodeVector& dfsc)
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


    bool RMPath::is_rightmost_vertex(DfscVI vi) const
    {
        assert(vi != VI_NULL);

        DfscVI_Vertices::const_iterator it = std::find(dfsc_vertices_.begin(),
                                                       dfsc_vertices_.begin() + num_rmp_vertices(),
                                                       vi);
        if (it != dfsc_vertices_.begin() + num_rmp_vertices())
            return true;

        assert(dfsc_vertices_.begin() + num_all_vertices() == dfsc_vertices_.end());
        assert(std::find(dfsc_vertices_.begin() + num_rmp_vertices(),
                         dfsc_vertices_.begin() + num_all_vertices(),
                         vi) != dfsc_vertices_.end());

        return false;
    }

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
    //                          RMPathSimple
    // *****************************************************************************
    class RMPathSimple
    {
        int* p_;
        const EdgeCode** ec_;
        int num_edges_;
        MemAllocator* ma_;
        std::size_t max_size_;
    public:
        RMPathSimple(std::size_t max_size, MemAllocator* ma)
            :p_(ma->alloc_array<int>(max_size)),
             ec_(ma->alloc_array<const EdgeCode*>(max_size)),
             num_edges_(0),
             ma_(ma),
             max_size_(max_size) {}

        ~RMPathSimple()
            {
                ma_->dealloc_array(p_, max_size_);
                ma_->dealloc_array(ec_, max_size_);                
            }

        void push(const EdgeCode& ec, int i);
        int num_edges() const { return num_edges_; }
        int operator[] (int i) const { return p_[i]; }
        int rightmost_edgeindex() const { return p_[num_edges_ - 1]; }
        DfscVI rightmost_vertex() const { return ec_[num_edges_ - 1]->vi_dst(); }
    };

    void RMPathSimple::push(const EdgeCode& ec, int i)
    {
        if (ec.is_forward())
        {
            if (num_edges_ == 0)
            {
                p_[num_edges_] = i;
                ec_[num_edges_] = &ec;
                ++num_edges_;
                return;
            }

            int n = num_edges_;
            while (n > 0 && ec_[n - 1]->vi_dst() != ec.vi_src())
                --n;
            
            p_[n] = i;
            ec_[n] = &ec;            
            num_edges_ = n + 1;
        }
    }


    // *****************************************************************************
    //                          Bitset
    // *****************************************************************************
    Bitset::Bitset(MemAllocator* ma, size_type size)
    {
        if (use_bitfield(size))
            bit_.field = Chunk(0);
        else
        {
            bit_.ptr = new (ma->alloc_array<Chunk>(num_chanks(size))) Chunk;
            std::fill(bit_.ptr, bit_.ptr + num_chanks(size), Chunk(0));
        }
    }

    Bitset::Bitset(MemAllocator* ma, size_type size, const Bitset& bf)
    {
        if (use_bitfield(size))
            bit_.field = bf.bit_.field;
        else
        {
            bit_.ptr = new (ma->alloc_array<Chunk>(num_chanks(size))) Chunk;
            std::copy(bf.bit_.ptr, bf.bit_.ptr + num_chanks(size), bit_.ptr);
        }
    }

    Bitset::~Bitset()
    {
        assert(! bit_.ptr);
    }

    void Bitset::free_resource(MemAllocator* ma, size_type size)
    {
        if (! use_bitfield(size))
            ma->dealloc_array(bit_.ptr, num_chanks(size));
        bit_.ptr = 0;
    }

    bool Bitset::test(size_type pos, size_type size) const
    {
        assert(pos < size);
        if (use_bitfield(size))
        {
            Chunk mask = Chunk(1U) << pos;
            return Chunk(0U) != (bit_.field & mask);
        }

        size_type i = pos / num_chunkbits();
        Chunk mask = Chunk(1U) << (pos % num_chunkbits());
        return Chunk(0U) != (bit_.ptr[i] & mask);
    }

    bool Bitset::set(size_type pos, size_type size)
    {
        assert(pos < size);
        if (use_bitfield(size))
        {
            Chunk mask = Chunk(1U) << pos;
            bool oldval = Chunk(0U) != (bit_.field & mask);
            bit_.field |= mask;
            return oldval;
        }

        size_type i = pos / num_chunkbits();
        Chunk mask = Chunk(1U) << (pos % num_chunkbits());
        bool oldval = Chunk(0U) != (bit_.ptr[i] & mask);
        bit_.ptr[i] |= mask;
        return oldval;
    }

    bool Bitset::clear(size_type pos, size_type size)
    {
        assert(pos < size);
        if (use_bitfield(size))
        {
            Chunk mask = Chunk(1U) << pos;
            bool oldval = Chunk(0U) != (bit_.field & mask);
            bit_.field &= ~mask;
            return oldval;
        }

        size_type i = pos / num_chunkbits();
        Chunk mask = Chunk(1U) << (pos % num_chunkbits());
        bool oldval = Chunk(0U) != (bit_.ptr[i] & mask);
        bit_.ptr[i] &= ~mask;
        return oldval;
    }

    bool Bitset::all_set(size_type size) const
    {
        for (Bitset::size_type i = 0; i < size; ++i)
            if (!test(i, size))
                return false;
        return true;
    }

    bool Bitset::is_equal(const Bitset& bf1, const Bitset& bf2, size_type size)
    {
        assert(bf1.use_bitfield(size) == bf2.use_bitfield(size));
        if (bf1.use_bitfield(size))
            return bf1.bit_.field == bf2.bit_.field;
        else
        {
            size_type n = num_chanks(size);
            for (size_type i = 0; i < n; ++i)
                if (bf1.bit_.ptr[i] != bf2.bit_.ptr[i])
                    return false;
            return true;
        }
    }
    
    void print(std::ostream& out, const Bitset& bf, Bitset::size_type size)
    {
        for (Bitset::size_type i = 0; i < size; ++i)
            out << int(bf.test(i, size));
    }

    
    // *****************************************************************************
    //                          SBGSimple
    // *****************************************************************************
    class SBGSimple : public SBGBase<SBGSimple>
    {
        friend class SBGCreator<SBGSimple>;

	SBGSimple(MemAllocator* ma, const Graph::Edge* e, const Graph* g)
	    :SBGBase<SBGSimple>(ma, e, g), edge_array_(0)
	    {}

	SBGSimple(MemAllocator* ma, const Graph::Edge* e, const SBGSimple* s)
	    :SBGBase<SBGSimple>(ma, e, s), edge_array_(0)
	    {}

        mutable const Graph::Edge** edge_array_;
    public:
        const Graph::Edge* find_edge(std::size_t depth, MemAllocator*) const NOINLINE;
        void free_edge_array(MemAllocator* ma)
            { if (edge_array_) { ma->dealloc_array(edge_array_, num_edges()); edge_array_ = 0; } }
    };

    const Graph::Edge* SBGSimple::find_edge(std::size_t depth, MemAllocator* ma) const
    {
        if (!edge_array_)
        {
            int i = num_edges();
            edge_array_ = ma->alloc_array<const Graph::Edge*>(i);

            const SBGSimple* s = this;
            do
            {
                edge_array_[--i] = s->edge();
                s = s->parent();
            } while (s);
            
        }
        return edge_array_[depth - 1];
    }

    // *****************************************************************************
    //                          specialization SBGCreator
    //                               for SBGSimple
    // *****************************************************************************
    template<>
    class SBGCreator<SBGSimple>
    {
	MemAllocator*	ma_;
	FixedAllocator*	sbg_alloc_;
    public:
	explicit SBGCreator(MemAllocator* ma)
	    :ma_(ma), sbg_alloc_(ma->get_fixed_allocator(sizeof(SBGSimple))) {}
	MemAllocator* get_mem_allocator() const { return ma_; }
	SBGSimple* new_sbg(const Graph::Edge* e, const Graph* g)
            { return new (sbg_alloc_->allocate()) SBGSimple(ma_, e, g); }
	SBGSimple* new_sbg(const Graph::Edge* e, const SBGSimple* s)
            { return new (sbg_alloc_->allocate()) SBGSimple(ma_, e, s); }
	void delete_sbg(SBGSimple*);
    };

    void SBGCreator<SBGSimple>::delete_sbg(SBGSimple* s)
    {
        s->free_bitsets(ma_);
        s->free_edge_array(ma_);
	s->~SBGSimple();
	sbg_alloc_->deallocate(s);
    }

    // *****************************************************************************
    //                          SBG
    // *****************************************************************************
    void SBG::init_dfsc_to_graph_array1(MemAllocator* mem_alloc)
    {
	vi_dfsc_to_graph_ = mem_alloc->alloc_array<GraphVI>(2);
	vi_dfsc_to_graph_[0] = edge()->vi_src();
        vi_dfsc_to_graph_[1] = edge()->vi_dst();
    }

    void SBG::init_dfsc_to_graph_array2(MemAllocator* mem_alloc)
    {
	PREFETCH(*parent()->vi_dfsc_to_graph_);
	vi_dfsc_to_graph_ = mem_alloc->alloc_array<GraphVI>(num_vertices());
	::memcpy(vi_dfsc_to_graph_,
		 parent()->vi_dfsc_to_graph_,
		 parent()->num_vertices() * sizeof(GraphVI));
	vi_dfsc_to_graph_[vi_src_dfsc_] = edge()->vi_src();
        vi_dfsc_to_graph_[vi_dst_dfsc_] = edge()->vi_dst();
    }

    void SBG::insert_to_automorph_list(SBG* pos, SBG* s)
    {
        pos->automorph_next_->automorph_prev_ = s;
        s->automorph_next_ = pos->automorph_next_;
        s->automorph_prev_ = pos;
        pos->automorph_next_ = s;
    }

    //
    // array of the dfsc VI, indexed by the graph VI
    // so, sbg_vertices[graph_vi] == dfsc vi
    //
    DfscVI* SBG::create_graph_to_dfsc_v(MemAllocator* mem_alloc, DfscVI vi_default) const
    {
        std::size_t n = get_graph()->num_vertices();
        DfscVI* p = mem_alloc->alloc_array<DfscVI>(n);
        std::fill(p, p + n, vi_default);

        const SBG* s = this;
        do
        {
            p[s->edge()->vi_src()] = s->vi_src_dfsc_;
            p[s->edge()->vi_dst()] = s->vi_dst_dfsc_;
            s = s->parent();
        } while (s);

        return p;
    }

    void SBG::free_graph_to_dfsc_v(DfscVI* array, MemAllocator* mem_alloc) const
    {
	mem_alloc->dealloc_array(array, get_graph()->num_vertices());
    }


    std::ostream& operator<<(std::ostream& out, const SBG& sbg)
    {
        std::vector<const SBG*> chain;
        get_chain(chain, &sbg);
        out << "sbg:";
        for (std::size_t i = 0; i < chain.size(); ++i)
            out << " " << *chain[i]->edge();
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
    //                          SubgraphsOfOneGraph
    // *****************************************************************************

    void SubgraphsOfOneGraph::calc_support() const
    {
        int support = std::numeric_limits<int>::max();
        const GraphVI n_ver_g = first_all()->get_graph()->num_vertices();
	const DfscVI  n_ver_s = first_all()->num_vertices();

        for (DfscVI vi = 0; vi < n_ver_s; ++vi)
	{
            int n = 0;
            std::vector<bool> vvg(n_ver_g, false);
	    for (const SBG* s = first_all(); s; s = s->next_embedding())
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
    }


    void SubgraphsOfOneGraph::insert(SBG* s)
    {
        list_sbg_all_.insert(s);
        
        SBG* sbg = list_sbg_autgroup_.first();
        while (sbg)
        {
            if (*sbg == *s)
            {
                SBG::insert_to_automorph_list(sbg, s);
                break;
            }
            sbg = sbg->next_automorph_group();
        }

        if (!sbg)
            list_sbg_autgroup_.insert(s);
    }

    bool has_parent(const SubgraphsOfOneGraph& sog, const SBG* par)
    {
        for (const SBG* s = sog.first_all(); s; s = s->next_embedding())
            if (s->parent() == par)
                return true;
        return false;
    }

    bool SubgraphsOfOneGraph::is_equal_occurence(const SubgraphsOfOneGraph& parent) const
    {
        unsigned int n = 0;
        for (const SBG* s = parent.first_autgroup(); s; s = s->next_automorph_group())
        {
            const SBG* aut = s;
            do
            {
                if (has_parent(*this, aut))
                {
                    ++n;
                    break;
                }
                aut = aut->next_automorph();
            } while (aut != s);
        }
        return n == parent.size_list_autgroup();
    }


    // *****************************************************************************
    //                          SubgraphsOfManyGraph
    // *****************************************************************************
    SubgraphsOfManyGraph::SubgraphsOfManyGraph(SBGCreator<SBG>* sbg_creator)
        :fa_(sbg_creator->get_mem_allocator()->get_fixed_allocator(sizeof(SOG))),
         sbg_creator_(sbg_creator),
         size_list_all_(0)
    {
    }
    

    SubgraphsOfManyGraph::~SubgraphsOfManyGraph()
    {
        iterator it = set_.begin();
        while (it != set_.end())
        {
            SOG* p = &*it;
            set_.erase(it++);
            p->~SOG();
            fa_->deallocate(p);
        }
    }

    void SubgraphsOfManyGraph::insert(SBG* s)
    {
        ++size_list_all_;
        SOGSet::insert_commit_data insert_data;
        std::pair<iterator,bool> r = set_.insert_check(s->get_graph(), CompareKey(), insert_data);
        if(r.second)
            r.first = set_.insert_commit(*new (fa_->allocate()) SOG(s->get_graph(), sbg_creator_), insert_data);
        r.first->insert(s);
    }

    const SubgraphsOfOneGraph& SubgraphsOfManyGraph::find(const Graph* g) const
    {
        const_iterator it = set_.find(g, CompareKey());
        assert(it != set_.end());
        return *it;
    }

    bool SubgraphsOfManyGraph::is_equal_occurence(const SubgraphsOfManyGraph& parent) const
    {
        if (support() != parent.support())
            return false;
        for (const_iterator it = begin(); it != end(); ++it)
        {
            const SOG& node = *parent.set_.find(it->get_graph(), CompareKey());
            if (! it->is_equal_occurence(node))
                return false;
        }
        return true;
    }

    // *****************************************************************************
    //                          Extension
    //                          set based
    // *****************************************************************************

    // SG       : SubgraphsOfOneGraph or SubgraphsOfManyGraph
    // Compare  : EdgeCodeCmpDfs or EdgeCodeCmpLex
    template<ContSelector CS, class SG, class Compare>   
    class Extension
        : public ContainerType<CS, EdgeCode, SG, Compare, false, void>::ContType
    {
        typedef ContainerType<CS, EdgeCode, SG, Compare, false, void> CT;
        typedef typename CT::NodeType Node;
        typedef typename CT::NodeCompare CompareKey;
        typedef typename CT::ContType Base;

        FixedAllocator* fa_;
        SBGCreator<SBG>* sbg_creator_;

        typename Base::iterator last_used_iter_;
    public:
        explicit Extension(SBGCreator<SBG>* sbg_creator)
            :fa_(sbg_creator->get_mem_allocator()->get_fixed_allocator(sizeof(Node))),
             sbg_creator_(sbg_creator) { last_used_iter_ = this->end(); }
        ~Extension();
        
        void insert(const EdgeCode& ec, const Graph::Edge* e, const Graph* g);
        void insert(const EdgeCode& ec, const Graph::Edge* e, const SBG* s);
    };

    template<ContSelector CS, class SG, class Compare>
    Extension<CS, SG, Compare>::~Extension()
    {
        typename Base::iterator it = this->begin();
        while (it != this->end())
        {
            Node* p = &*it;
            this->erase(it++);
            p->~Node();
            fa_->deallocate(p);
        }
    }

    template<ContSelector CS, class SG, class Compare>        
    void Extension<CS, SG, Compare>::insert(const EdgeCode& ec, const Graph::Edge* e, const Graph* g)
    {
        typename Base::insert_commit_data insert_data;
        std::pair<typename Base::iterator,bool> r = this->insert_check(last_used_iter_,
                                                                       ec, CompareKey(), insert_data);
        if(r.second)
            r.first = this->insert_commit(*new (fa_->allocate()) Node(ec, sbg_creator_), insert_data);
        r.first->insert(sbg_creator_->new_sbg(e, g));
        last_used_iter_ = r.first;
    }

    template<ContSelector CS, class SG, class Compare>
    void Extension<CS, SG, Compare>::insert(const EdgeCode& ec, const Graph::Edge* e, const SBG* s)
    {
        typename Base::insert_commit_data insert_data;
        std::pair<typename Base::iterator,bool> r = this->insert_check(last_used_iter_,
                                                                       ec, CompareKey(), insert_data);
        if(r.second)
            r.first = this->insert_commit(*new (fa_->allocate()) Node(ec, sbg_creator_), insert_data);
        r.first->insert(sbg_creator_->new_sbg(e, s, ec));
        last_used_iter_ = r.first;
    }

    template<class SG, template <class,class> class Node>
    inline bool is_equal_occurrence(const Node<EdgeCode,SG>& ext_sg, const SG& sg_parent)
    {
        return ext_sg.is_equal_occurence(sg_parent);
    }



    // *****************************************************************************
    //                          Extension
    //                          list based
    // *****************************************************************************

    // SG       : SubgraphsOfOneGraph or SubgraphsOfManyGraph
    // Compare  : EdgeCodeCmpDfs or EdgeCodeCmpLex
    template<class SG, class Compare>
    class Extension<LIST,SG,Compare>
        : public ContainerType<LIST, EdgeCode, SG, Compare, false, void>::ContType
    {
        typedef ContainerType<LIST, EdgeCode, SG, Compare, false, void> CT;
        typedef typename CT::NodeType Node;
        typedef typename CT::NodeCompare CompareKey;
        typedef typename CT::ContType Base;

        FixedAllocator* fa_;
        SBGCreator<SBG>* sbg_creator_;

        typename Base::iterator last_used_iter_;
        Base& base() { return *this; }
    public:
        explicit Extension(SBGCreator<SBG>* sbg_creator)
            :fa_(sbg_creator->get_mem_allocator()->get_fixed_allocator(sizeof(Node))),
             sbg_creator_(sbg_creator) { last_used_iter_ = this->end(); }
        ~Extension();

        void insert(const EdgeCode& ec, const Graph::Edge* e, const SBG* s);
    };

    template<class SG, class Compare>
    Extension<LIST, SG, Compare>::~Extension()
    {
        typename Base::iterator it = this->begin();
        while (it != this->end())
        {
            Node* p = &*it;
            this->erase(it++);
            p->~Node();
            fa_->deallocate(p);
        }
    }

    template <class ForwardIterator, class Compare>
    bool is_sorted (ForwardIterator first, ForwardIterator last, Compare comp)
    {
        if (br)
            BR;

        Compare cmp;
        ForwardIterator prev = first;
        ForwardIterator next = first;
        if (next != last)
            ++next;

        while (next != last)
        {
            if (! cmp(prev->get_key(), next->get_key()))
                return false;
            ++prev;
            ++next;
        }
        return true;
    }

    template<class SG, class Compare>
    void Extension<LIST, SG, Compare>::insert(const EdgeCode& ec, const Graph::Edge* e, const SBG* s)
    {
        SBG* sbg = sbg_creator_->new_sbg(e, s, ec);
        Compare cmp;
        typename Base::iterator it = base().begin();
        typename Base::iterator it_end = base().end();

        if (last_used_iter_ != base().end())
        {
            if (! cmp(ec, last_used_iter_->get_key()))
            {
                // last_used_iter_->get_key() <= ec
                it = last_used_iter_;
            }
        }

        while (it != it_end)
        {
            if (cmp(ec, it->get_key()))
            {
                // ec < it->get_key()
                Node* node = new (fa_->allocate()) Node(ec, sbg_creator_);
                last_used_iter_ = base().insert(it, *node);
                node->insert(sbg);

#ifdef WITH_CHECKS
                assert(is_sorted(base().begin(), base().end(), Compare()));
#endif
                return;
            }
            
            if (!cmp(it->get_key(), ec))
            {
                // equal
                it->insert(sbg);
                last_used_iter_ = it;
                return;
            }
            
            ++it;
        }
        
        Node* node = new (fa_->allocate()) Node(ec, sbg_creator_);
        base().push_back(*node);
        node->insert(sbg);
        last_used_iter_ = base().begin();

#ifdef WITH_CHECKS
        assert(is_sorted(base().begin(), base().end(), Compare()));
#endif
    }

    // *****************************************************************************
    //                          MinimalExtension
    // *****************************************************************************
    class MinimalExtension
    {
    public:
        typedef typename Make_SBG_List_Embedding<SBGSimple>::Type List;

        explicit MinimalExtension(SBGCreator<SBGSimple>* sbg_creator) :list_(sbg_creator), sbg_creator_(sbg_creator)  {}
        const SBGSimple* first() const  { return list_.first(); }  
        bool empty() const              { return list_.empty(); }
        void insert(const EdgeCode& ec, const Graph::Edge* e, const Graph* g);
        void insert(const EdgeCode& ec, const Graph::Edge* e, const SBGSimple* s);
        const EdgeCode& get_ec() const  { return ec_; }
        List& get_list()                { return list_; }

        MinimalExtension* next;
    private:
        EdgeCode ec_;
        List list_;
        SBGCreator<SBGSimple>* sbg_creator_;
    };


    void MinimalExtension::insert(const EdgeCode& ec, const Graph::Edge* e, const Graph* g)
    {
        EdgeCodeCmpDfs cmp;
        if (! empty() && cmp(ec_, ec))
            return;
        SBGSimple* p = sbg_creator_->new_sbg(e, g);
        if (empty())
        {
            ec_ = ec;
            list_.insert(p);
        }
        else if (ec_ == ec)
            list_.insert(p);
        else
        {
            assert(cmp(ec, ec_));
            ec_ = ec;
            list_.clear();
            list_.insert(p);            
        }
    }

    void MinimalExtension::insert(const EdgeCode& ec, const Graph::Edge* e, const SBGSimple* s)
    {
        EdgeCodeCmpDfs cmp;
        if (! empty() && cmp(ec_, ec))
            return;
        SBGSimple* p = sbg_creator_->new_sbg(e, s);
        if (empty())
        {
            ec_ = ec;
            list_.insert(p);
        }
        else if (ec_ == ec)
            list_.insert(p);
        else
        {
            assert(cmp(ec, ec_));
            ec_ = ec;
            list_.clear();
            list_.insert(p);            
        }
    }


    // *****************************************************************************
    //                         MinimalExtensionList
    // *****************************************************************************
    class MinimalExtensionList : private boost::noncopyable
    {
        MinimalExtension* head_;
        SBGCreator<SBGSimple>* sbg_creator_;
        FixedAllocator* fa_;
    public:
        MinimalExtensionList(SBGCreator<SBGSimple>* sbg_creator)
            :head_(0),
             sbg_creator_(sbg_creator),
             fa_(sbg_creator->get_mem_allocator()->get_fixed_allocator(sizeof(*head_))) {}
        ~MinimalExtensionList()
            {
                MinimalExtension* p = head_;
                while (p)
                {
                    MinimalExtension* next = p->next;
                    p->~MinimalExtension();
                    fa_->deallocate(p);
                    p = next;
                }
            }
        MinimalExtension* push()
            {
                MinimalExtension* p = new (fa_->allocate()) MinimalExtension(sbg_creator_);
                p->next = head_;
                head_ = p;
                return p;
            }
    };



    // *****************************************************************************
    //                         SharedData
    // *****************************************************************************
    
    struct SharedData
    {
        DFSCode dfsc;
 
	MemAllocator		mem_alloc;
	SBGCreator<SBG>		sbg_creator;
	SBGCreator<SBGSimple>	sbgsimple_creator;

	const support_type minsup_;
        GspanResult* result_;
 
        bool flg;
        std::map<const Graph*, bool> whole_was_reached;

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

        SharedData(support_type minsup,
                   std::size_t max_num_vertices,
                   std::size_t max_num_edges,
                   GspanResult* result,
                   unsigned int max_trace_depth)
            : dfsc(max_num_vertices, max_num_edges)
            , sbg_creator(&mem_alloc)
	    , sbgsimple_creator(&mem_alloc)
	    , minsup_(minsup)
            , result_(result)
            , flg(false)
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
	SG* sg;

        RMPath rmpath;

	typedef Extension<SET,  SG, EdgeCodeCmpDfs> REdges;
	typedef Extension<LIST, SG, EdgeCodeCmpLex> XEdges;
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
        const SG* sg_ext_notclose;
        
        // cause of early_termin
        ExtType exttype_early_term;
        const EdgeCode* ec_early_term;
        const SG* sg_ext_early_term;

	typedef std::pair<const EdgeCode*, SG*> RChild;
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

	FrameState(SharedData* shared, FrameState* prev, SG* subgraphs);
    };


    template<class SG>
    FrameState<SG>::FrameState(SharedData* shared, FrameState* prev, SG* subgraphs)
	: sg(subgraphs)
	, rmpath(shared->dfsc.get_ecvector(), &shared->mem_alloc)
	, x_edges(&shared->sbg_creator)
	, r_edges(&shared->sbg_creator)
	, closed(true)
	, early_term(false)
	, exttype_notclose(EXT_NONE), ec_notclose(0), sg_ext_notclose(0)
	, exttype_early_term(EXT_NONE), ec_early_term(0), sg_ext_early_term(0)
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
        const Graph::EdgesSet& g_edges = g.edges();
        for (Graph::EdgesIterator it = g_edges.begin(); it != g_edges.end(); ++it)
            ext.insert(EdgeCode(0, 1, it->vl_src(), it->el(), it->vl_dst(), true), &*it, &g);
    }



    // -----------------------------------
    //	dfscode minimality check
    // -----------------------------------

    // for dfsc check minimality, we use std::vector as DFSCode
    typedef EdgeCodeVector DFSCodeLite;


    void enumerate_min_bck(MinimalExtension& ext,
			   const DFSCodeLite& dfsc_min,
                           typename Make_SBG_List_Embedding<SBGSimple>::Type& slist,
			   const RMPathSimple& rmpath,
			   const Graph& g,
                           MemAllocator* mem_allocator)
    {
        int rmpath_numedges = rmpath.num_edges();
        int rmpath_rmost_e_i = rmpath.rightmost_edgeindex();
	const DfscVI vi_dfsc_rmost = rmpath.rightmost_vertex();
        assert(rmpath.rightmost_vertex() == dfsc_min[rmpath_rmost_e_i].vi_dst());
	const VL vl_rmost = dfsc_min[rmpath_rmost_e_i].vl_dst();

	// from first vertex toward right most vertex
	for (int i = 0; ext.empty() && i < rmpath_numedges - 1; ++i)
	{
            int rmpath_e_i = rmpath[i];

            const EdgeCode& ec_rmpath = dfsc_min[rmpath_e_i];
            const EL el_rmpath = ec_rmpath.el();
            const bool vl_less_eq = ec_rmpath.vl_dst() <= vl_rmost;

	    for (SBGSimple* s = slist.first(); s; s = s->next_embedding())
	    {
		const Graph::Edge* e_rmost = s->find_edge(rmpath_rmost_e_i + 1, mem_allocator);
		const GraphVI graph_vi = s->find_edge(rmpath_e_i + 1, mem_allocator)->vi_src();

                if (s->has_no_extension(graph_vi))
                    continue;
                s->set_no_has_extension(graph_vi);

                for (const Graph::Edge* e = g.incident(graph_vi); e; e = e->next_adjacent())
                {
                    if (s->has_edge(e->eid()))
                        continue;
		    s->set_has_extension(graph_vi);
		    
		    if (e->vi_dst() == e_rmost->vi_dst() &&
                        ((vl_less_eq && el_rmpath == e->el()) || el_rmpath < e->el()))
                    {
                        const Graph::Edge* e_rev = e->reverse();
                        EdgeCode ec(vi_dfsc_rmost, ec_rmpath.vi_src(), vl_rmost, e_rev->el(), ec_rmpath.vl_src(), false);
                        ext.insert(ec, e_rev, s);
                        break;
                    }
		}
	    }
	}
    }

    void enumerate_min_fwd(MinimalExtension& ext,
			   const DFSCodeLite& dfsc_min,
                           typename Make_SBG_List_Embedding<SBGSimple>::Type& slist,
			   const RMPathSimple& rmpath,
			   const Graph& g,
                           MemAllocator* mem_allocator)
    {
	const DfscVI vi_dfsc_rmost = rmpath.rightmost_vertex();
        assert(rmpath.rightmost_vertex() == dfsc_min[rmpath.rightmost_edgeindex()].vi_dst());

        const VL vl_rmost = dfsc_min[rmpath.rightmost_edgeindex()].vl_dst();
        const VL vl_minimum = dfsc_min[0].vl_src();

	// forward pure
	for (SBGSimple* s = slist.first(); s; s = s->next_embedding())
	{
	    const Graph::Edge* e_rmost = s->find_edge(rmpath.rightmost_edgeindex() + 1, mem_allocator);
	    GraphVI graph_vi = e_rmost->vi_dst();

            if (s->has_no_extension(graph_vi))
                continue;
            s->set_no_has_extension(graph_vi);

            for (const Graph::Edge* e = g.incident(graph_vi); e; e = e->next_adjacent())
            {
                if (s->has_edge(e->eid()))
                    continue;
                s->set_has_extension(graph_vi);

                if (! s->has_vertex(e->vi_dst()) && vl_minimum <= e->vl_dst())
                {
                    EdgeCode ec(vi_dfsc_rmost, vi_dfsc_rmost+1, vl_rmost, e->el(), e->vl_dst(), true);
                    ext.insert(ec, e, s);
                }
            }
	}

	// forward rmpath
	// from right most vertex toward first vertex
	for (int i = rmpath.num_edges() - 1; ext.empty() && i >= 0; --i)
	{
	    const EdgeCode& ec_rmpath = dfsc_min[rmpath[i]];
	    for (SBGSimple* s = slist.first(); s; s = s->next_embedding())
	    {
		const GraphVI graph_vi = s->find_edge(rmpath[i] + 1, mem_allocator)->vi_src();

		if (s->has_no_extension(graph_vi))
		    continue;
		s->set_no_has_extension(graph_vi);

                for (const Graph::Edge* e = g.incident(graph_vi); e; e = e->next_adjacent())
                {
                    if (s->has_edge(e->eid()))
                        continue;
		    s->set_has_extension(graph_vi);

                    if (! s->has_vertex(e->vi_dst()) && vl_minimum <= e->vl_dst() &&
                        ((ec_rmpath.vl_dst() <= e->vl_dst() && ec_rmpath.el() == e->el()) || ec_rmpath.el() < e->el()))
                    {
                        EdgeCode ec(ec_rmpath.vi_src(), vi_dfsc_rmost+1, ec_rmpath.vl_src(), e->el(), e->vl_dst(), true);
                        ext.insert(ec, e, s);
                    }
                }
	    }
	}
    }
    
    bool is_min_iterative(const DFSCode& dfsc_tested, SBGCreator<SBGSimple>* sbgsimple_creator)
    {
        MemAllocator* mem_alloc = sbgsimple_creator->get_mem_allocator();
        DFSCodeLite dfsc_min;
	dfsc_min.reserve(dfsc_tested.size());

	const Graph& graph = dfsc_tested.get_graph();
        MinimalExtensionList ext_list(sbgsimple_creator);
        MinimalExtension* exts1 = ext_list.push();      

        enum_one_edges(*exts1, graph);

        RMPathSimple rmpath(dfsc_tested.size(), mem_alloc);

        while (true)
	{
	    dfsc_min.push_back(exts1->get_ec());
            int lasti = dfsc_min.size() - 1;
            if (dfsc_min[lasti] != dfsc_tested[lasti])
                return false;

            rmpath.push(exts1->get_ec(), lasti);

            MinimalExtension* exts2 = ext_list.push();

            enumerate_min_bck(*exts2, dfsc_min, exts1->get_list(), rmpath, graph, mem_alloc);
	    if (! exts2->empty())
            {
                exts1 = exts2;
                continue;
            }

	    enumerate_min_fwd(*exts2, dfsc_min, exts1->get_list(), rmpath, graph, mem_alloc);
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
            last_result = is_min_iterative(shared->dfsc, &shared->sbgsimple_creator);
            shared->on_ismin(last_result);
            return last_result;
        }
    }


    template<class P>
    bool is_min(SharedData* shared, FrameState<P>* frame,
		const EdgeCode& ec, bool& last_result, EdgeCode& ec_last)
    {
        shared->dfsc.push(ec);
        bool r = is_min_fast(shared, ec, frame->rmpath, last_result, ec_last);
        if (r)
            ec_last = ec;
        shared->dfsc.pop();
        return r;
    }

    // -----------------------------------
    //	trace, check and debug functions
    // -----------------------------------

#ifdef GSPAN_TRACE
    void print_msg(const FrameState<SubgraphsOfOneGraph>* frame, const EdgeCode& ec)
    {
        std::cerr << ec
                  << " support=" << frame->sg->support()
                  << " num_subgraphs=" << frame->sg->size_list_all()
                  << " num_autgroups=" << frame->sg->size_list_autgroup()
                  << " frame=" << frame->id;
    }

    void print_msg(const FrameState<SubgraphsOfManyGraph>* frame, const EdgeCode& ec)
    {
        std::cerr << ec
                  << " support=" << frame->sg->support()
                  << " num_subgraphs=" << frame->sg->size_list_all()
                  << " frame=" << frame->id;
    }

    template<class SG>
    void print_frame_trace(const FrameState<SG>* frame, const SharedData* shared)
    {
        using namespace std;
        cerr << shared->dfsc.size() << ":";
        for (unsigned int i = 0; i < shared->dfsc.size(); ++i) cerr << " ";
        print_msg(frame, shared->dfsc.top());
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
                 << (frame->rmpath.is_rightmost_edgeindex(i) ? "*" : " ")
                 << shared->dfsc[i]
                 << endl;
        cerr << "--------------------------\n";
        cerr << "Embeddings:\n";
        cerr << *frame->sg;
        cerr << "--------------------------\n";
        cerr << "XExtentions:\n";
            
        for (XEcIter it = frame->x_edges.begin(); it != frame->x_edges.end(); ++it)
            cerr << it->get_key() << "\n" << *it;
            
        cerr << endl;
        const FrameState<SG>* fet_frame = frame->early_term_frame;
        cerr << "============ FAIL Eearly Termination FRAME: " << fet_frame->id << " ======================" << endl;
        cerr << "DFSC:\n";
        cerr << "--------------------------\n";
        for (unsigned int i = 0; i < fet_frame->depth; ++i)
            cerr << i+1 << ":\t"
                 << (fet_frame->rmpath.is_rightmost_edgeindex(i) ? "*" : " ")
                 << shared->dfsc[i]
                 << endl;
        cerr << "--------------------------\n";
        cerr << "Embeddings:\n";
        cerr << *fet_frame->sg;
        cerr << "--------------------------\n";
        cerr << "Cause Early Termination:\n";
        cerr << (fet_frame->exttype_early_term == EXT_X ? "X" : "R") << *fet_frame->ec_early_term << endl;
        cerr << *fet_frame->sg_ext_early_term;
        cerr << "--------------------------\n";
        cerr << "RExtensions:\n";
        const EdgeCode& ec_ext = shared->dfsc[fet_frame->depth];
        for (REcIter it = fet_frame->r_edges.begin(); it != fet_frame->r_edges.end(); ++it)
        {
            bool ch = it->get_key() == ec_ext;
            cerr << it->get_key()
                 << (is_child(fet_frame, &it->get_key()) ? " * " : " ")
                 << (ch ? "CHILD" : " ")
                 << "\n" << *it;
            if (ch)
                break;
        }
        cerr << endl;
    }
#endif

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
	    frame->sg->support();

            (*shared->result_)(shared->dfsc, *frame->sg);
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
    void enumerate_sog(FrameState<SG>* frame, SharedData* shared, SubgraphsOfOneGraph& sg)
    {
        typedef typename FrameState<SG>::REdges REdges;
        typedef typename FrameState<SG>::XEdges XEdges;

	MemAllocator* mem_alloc = &shared->mem_alloc;
        REdges& r_edges = frame->r_edges;
        XEdges& x_edges = frame->x_edges;
        const RMPath& rmpath = frame->rmpath;
        const DfscVI vi_dfsc_rmost   = rmpath.rightmost_vertex();
        const DfscVI vi_dfsc_new = vi_dfsc_rmost + 1;

        assert(rmpath.rightmost_vertex() == shared->dfsc[rmpath.rightmost_edgeindex()].vi_dst());

        const Graph* g = sg.get_graph();

        for (SBG* s = sg.first_all(); s; s = s->next_embedding())
        {

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

            //
            // R edges will be
            // IF
            // 1) vi_src is rmpath vertex AND
            // 2) vi_dst is new vertex (forward) OR 
            //     (vi_src is rmost vertex AND vi_dst is any rmpath vertex)
            // ELSE
            // X edges
            //

#ifdef WITH_CHECKS
            std::set<const Graph::Edge*> e_found;
#endif

            for (RMPath::const_iterator dfsc_vi_iter = rmpath.begin();
                 dfsc_vi_iter != rmpath.rmp_end(); ++dfsc_vi_iter)
            {
                const DfscVI dfsc_vi = *dfsc_vi_iter;
                const GraphVI graph_vi = dfsc_to_graph_v[dfsc_vi];

                assert(graph_to_dfsc_v[graph_vi] == dfsc_vi);
                assert(graph_vi >= 0 && graph_vi < s->get_graph()->num_vertices());
                assert(graph_vi < graph_size(g));

                if (s->has_no_extension(graph_vi))
                    continue;
                s->set_no_has_extension(graph_vi);

                for (const Graph::Edge* e = g->incident(graph_vi); e; e = e->next_adjacent())
                {
                    if (s->has_edge(e->eid()))
                        continue;
                    s->set_has_extension(graph_vi);
                            
                    // vi_src is rmpath vertex

                    if (! s->has_vertex(e->vi_dst()))
                    {
                        // vi_dst is new vertex
                        // R forward
                        EdgeCode ec(dfsc_vi, vi_dfsc_new, e->vl_src(), e->el(), e->vl_dst(), true);
                        r_edges.insert(ec, e, s);
                    }
                    else if (e->vi_src() == vi_rmost && vertex_rmpath_status[e->vi_dst()])
                    {
                        // vi_src is rmost vertex AND vi_dst is any rmpath vertex
                        // R backward
                        DfscVI dfsc_vi_dst = graph_to_dfsc_v[e->vi_dst()];
                        assert(dfsc_vi_dst != vi_dfsc_new);
                        EdgeCode ec(dfsc_vi, dfsc_vi_dst, e->vl_src(), e->el(), e->vl_dst(), false);
                        r_edges.insert(ec, e, s);
                    }
                    else
                    {
                        if (! (e->vi_dst() == vi_rmost && vertex_rmpath_status[e->vi_src()]))
                        {
                            // X backward
                            DfscVI dfsc_vi_dst = graph_to_dfsc_v[e->vi_dst()];
                            assert(dfsc_vi_dst != vi_dfsc_new);
                            EdgeCode ec(dfsc_vi, dfsc_vi_dst, e->vl_src(), e->el(), e->vl_dst(), false);
                            x_edges.insert(ec, e, s);
                        }
                    }

#ifdef WITH_CHECKS
                    assert(e_found.count(e) == 0);
                    e_found.insert(e);
#endif
                }                
            }


            for (RMPath::const_iterator dfsc_vi_iter = rmpath.rmp_end();
                 dfsc_vi_iter != rmpath.end(); ++dfsc_vi_iter)
            {
                const DfscVI dfsc_vi = *dfsc_vi_iter;
                const GraphVI graph_vi = dfsc_to_graph_v[dfsc_vi];

                assert(graph_to_dfsc_v[graph_vi] == dfsc_vi);
                assert(graph_vi >= 0 && graph_vi < s->get_graph()->num_vertices());
                assert(graph_vi < graph_size(g));

                if (s->has_no_extension(graph_vi))
                    continue;
                s->set_no_has_extension(graph_vi);

                for (const Graph::Edge* e = g->incident(graph_vi); e; e = e->next_adjacent())
                {
                    if (s->has_edge(e->eid()))
                        continue;

                    s->set_has_extension(graph_vi);

                    // from vertex not on rmpath
                    
                    DfscVI dfsc_vi_dst = graph_to_dfsc_v[e->vi_dst()];
                    EdgeCode ec(dfsc_vi, dfsc_vi_dst,
                                e->vl_src(), e->el(), e->vl_dst(), ! s->has_vertex(e->vi_dst()));
                    x_edges.insert(ec, e, s);

#ifdef WITH_CHECKS
                    assert(e_found.count(e) == 0);
                    e_found.insert(e);
#endif
                }
            }
            
	    s->free_graph_to_dfsc_v(graph_to_dfsc_v, mem_alloc);
        } // end:: for all embeddings

    }



    void enumerate(FrameState<SubgraphsOfOneGraph>* frame, SharedData* shared)
    {
        enumerate_sog(frame, shared, *frame->sg);
    }



    void enumerate(FrameState<SubgraphsOfManyGraph>* frame, SharedData* shared)
    {
        SubgraphsOfManyGraph* smg = frame->sg;
        for (SubgraphsOfManyGraph::iterator it = smg->begin(); it != smg->end(); ++it)
            enumerate_sog(frame, shared, *it);
    }



    template<class Ext>
    void remove_not_frequents(Ext& ext, support_type minsup)
    {
        if (minsup > 1)
	{
            for (typename Ext::iterator it = ext.begin(); it != ext.end();)
	    {
                if (it->support() < minsup)
                    ext.erase(it++);
                else
                    ++it;
	    }
	}
    }


    // -----------------------------------
    //	early termination
    // -----------------------------------

    template<class SG>
    struct et_param
    {
        const SG* sg_ext;
        const EdgeCode* ec_ext;
        const SG* sg_ori;
        DFSCode* dfsc_ori;
        const RMPath* rmpath_ori;
        SBGCreator<SBGSimple>* sbgsimple_creator;
        ExtType exttype;
    };
    
    struct EdgeList
    {
        const Graph::Edge* edge;
        const EdgeList* prev;
        EdgeList(const Graph::Edge* e = 0, const EdgeList* p = 0) :edge(e), prev(p) {}
    };


#ifdef ET_TEST_PATH_MIN
    bool is_path_minimum(const EdgeList* path, const SBG* s_ext, const DfscVI* graph_to_dfsc_ori,
                         const et_param<SubgraphsOfOneGraph>* param)
    {        
        std::vector<const Graph::Edge*> rev_path;

        const DfscVI vi_connect = graph_to_dfsc_ori[path->edge->vi_dst()];

        // ----------------------------------------------------------
        // check path from some embeddings vertex to extention vertex
        // ----------------------------------------------------------
        DfscVI vi_first = vi_connect;
        DfscVI vi_second = param->rmpath_ori->rightmost_vertex() + 1;
        assert(vi_first != VI_NULL);
        assert(param->rmpath_ori->is_rightmost_vertex(vi_first));

        int num_edgecodes = 0;
        for (const EdgeList* e = path; e; e = e->prev)
        {
            ++num_edgecodes;

            EdgeCode ec(vi_first, vi_second, e->edge->vl_dst(), e->edge->el(), e->edge->vl_src(), true);
            param->dfsc_ori->push(ec);
            rev_path.push_back(e->edge);
            vi_first = vi_second;
            ++vi_second;
        }

        bool ism_direct = is_min_iterative(*param->dfsc_ori, param->sbgsimple_creator);

        for (int i = 0; i < num_edgecodes; ++i)
            param->dfsc_ori->pop();

        return ism_direct;
    }
#endif


    bool exist_path_minimum(const SBG* s_ext, const et_param<SubgraphsOfOneGraph>* param)
    {
        MemAllocator* ma = param->sbgsimple_creator->get_mem_allocator();
	const Graph& g = *s_ext->get_graph();
        std::vector<char> visited(g.num_vertices(), false);
        
        VertexRMPathStatus vertex_rmpath_status(ma);
        DfscVI* graph_to_dfsc_ori = 0;
        //if (param->exttype == EXT_R)
        {
            vertex_rmpath_status.create(*param->rmpath_ori, s_ext);
            graph_to_dfsc_ori = s_ext->parent()->create_graph_to_dfsc_v(ma);
        }
        
        typedef STL_Allocator<EdgeList> Alloc;
        typedef STL_Allocator<const EdgeList*> AllocPtr;
        Alloc alloc(ma);
        AllocPtr alloc_ptr(ma);
        std::list<EdgeList, Alloc> li(alloc);
        std::list<const EdgeList*, AllocPtr> qu(alloc_ptr);

        bool ret = false;

        for (const Graph::Edge* e = g.incident(s_ext->edge()->vi_dst()); !ret && e; e = e->next_adjacent())
        {
            if (s_ext->has_edge(e->eid()))
                continue;

            EdgeList pthe(e);
            if (s_ext->has_vertex(e->vi_dst()))
            {
                if (//param->exttype == EXT_X ||
                    (vertex_rmpath_status[e->vi_dst()]
#ifdef ET_TEST_PATH_MIN
                     && is_path_minimum(&pthe, s_ext, graph_to_dfsc_ori, param)
#endif
                        ))
                {
                    ret = true;
                    break;
                }
            }
            else
            {
                li.push_back(pthe);
                qu.push_back(&li.back());
            }
        }
        visited[s_ext->edge()->vi_dst()] = true;

        while (!ret && ! qu.empty())
        {
            const EdgeList& pthe1 = *qu.front();
            qu.pop_front();
            visited[pthe1.edge->vi_dst()] = true;

            for (const Graph::Edge* e = g.incident(pthe1.edge->vi_dst()); e; e = e->next_adjacent())
            {
                EdgeList pthe2(e, &pthe1);
                if (!visited[e->vi_dst()])
                {
                    if (s_ext->has_vertex(e->vi_dst()))
                    {
                        if (//param->exttype == EXT_X ||
                            (vertex_rmpath_status[e->vi_dst()]
#ifdef ET_TEST_PATH_MIN
                             && is_path_minimum(&pthe2, s_ext, graph_to_dfsc_ori, param)
#endif
                                ))
                        {
                            ret = true;
                            break;
                        }
                    }
                    else
                    {
                        li.push_back(pthe2);
                        qu.push_back(&li.back());
                    }
                }
            }
        }

        if (graph_to_dfsc_ori)
            s_ext->parent()->free_graph_to_dfsc_v(graph_to_dfsc_ori, ma);

        return ret;
    }



    bool fail_early_termination(const et_param<SubgraphsOfOneGraph>* param)
    {
        if (param->ec_ext->is_forward())
        {          
            for (const SBG* s = param->sg_ext->first_all(); s; s = s->next_embedding())
                if (exist_path_minimum(s, param))
                    return true;
            return false;
        }
        else
        {
            if (param->exttype == EXT_R)
                return false;
            return param->sg_ext->size_list_all() != param->sg_ori->size_list_all();
        }
    }



    bool fail_early_termination(const et_param<SubgraphsOfManyGraph>* param)
    {
        et_param<SubgraphsOfOneGraph> psog;
        psog.ec_ext = param->ec_ext;
        psog.dfsc_ori = param->dfsc_ori;
        psog.rmpath_ori = param->rmpath_ori;
        psog.sbgsimple_creator = param->sbgsimple_creator;
        psog.exttype = param->exttype;

        const SubgraphsOfManyGraph& smg = *param->sg_ext;
        for (SubgraphsOfManyGraph::const_iterator i = smg.begin(); i != smg.end(); ++i)
        {
            psog.sg_ext = &*i;
            psog.sg_ori = &param->sg_ori->find(i->get_graph());
            if (fail_early_termination(&psog))
                return true;
        }
        return false;
    }
    


    template<class SG, class SGNode>
    void detect_et(FrameState<SG>* frame, SharedData* shared,
                   const SGNode& ext_sg,
                   bool equiv, ExtType exttype)
    {
#ifdef CHECK_MODE
        if (frame->TERMINATED_BRANCH)
            return;
#endif
        if (!equiv)
            return;
        if (frame->early_term)
            return;

        et_param<SG> param;
        param.sg_ext = &ext_sg;
        param.ec_ext = &ext_sg.get_key();
        param.sg_ori = frame->sg;
        param.dfsc_ori = &shared->dfsc;
        param.rmpath_ori = &frame->rmpath;
        param.sbgsimple_creator = &shared->sbgsimple_creator;
        param.exttype = exttype;

        if (fail_early_termination(&param))
            return;
        
        if (param.ec_ext->is_forward())
            shared->on_early_term_f(exttype);
        else
            shared->on_early_term_b(exttype);            

        assert(frame->exttype_early_term == EXT_NONE);
        frame->early_term = true;
        frame->exttype_early_term = exttype;
        frame->ec_early_term = &ext_sg.get_key();
        frame->sg_ext_early_term = &ext_sg;
    }



    bool detect_et_2(SharedData* shared, FrameState<SubgraphsOfOneGraph>* frame)
    {
        bool& whole_was_reached = shared->whole_was_reached[frame->sg->get_graph()];
        if (! whole_was_reached)
        {
            whole_was_reached = frame->x_edges.empty() && frame->r_edges.empty();
            if (whole_was_reached)
                assert(frame->sg->support() == 1);
        }        
        return whole_was_reached && frame->sg->support() == 1;
    }



    bool detect_et_2(SharedData* shared, FrameState<SubgraphsOfManyGraph>* frame)
    {
        if (frame->sg->support() == 1)
        {
            bool& whole_was_reached = shared->whole_was_reached[frame->sg->begin()->get_graph()];
            if (whole_was_reached)
                return true;
            else
                whole_was_reached = frame->x_edges.empty() && frame->r_edges.empty();
        }
        return false;
    }

    // -----------------------------------



    template<class SG, class SGNode>
    void detect_nc(FrameState<SG>* frame, const SGNode& ext_sg, bool equiv, ExtType exttype)
    {
        if (frame->closed && (equiv || ext_sg.support() == frame->sg->support()))
        {
            frame->closed = false;
            assert(frame->exttype_notclose == EXT_NONE);
            frame->exttype_notclose = exttype;
            frame->ec_notclose = &ext_sg.get_key();
            frame->sg_ext_notclose = &ext_sg;
        }
    }


    template<class SG>
    void project(SharedData* shared, FrameState<SG>* parent_frame, SG* sg_parent)
    {
        typedef typename FrameState<SG>::REdges REdges;
        typedef typename FrameState<SG>::XEdges XEdges;
        typedef typename REdges::const_iterator REcIter;
        typedef typename XEdges::const_iterator XEcIter;
        typedef typename REdges::iterator REIter;
        typedef typename XEdges::iterator XEIter;
        typedef typename FrameState<SG>::RChild RChild;
        typedef typename FrameState<SG>::RChildren RChildren;

        shared->on_enter_frame();
        FrameState<SG> frame(shared, parent_frame, sg_parent);
        
        enumerate(&frame, shared);
        remove_not_frequents(frame.r_edges, shared->minsup_);

        br = frame.id == 2;

        for (XEIter it = frame.x_edges.begin(); it != frame.x_edges.end(); ++it)
        {
            bool equiv = is_equal_occurrence(*it, *sg_parent);
            detect_nc(&frame, *it, equiv, EXT_X);
            detect_et(&frame, shared, *it, equiv, EXT_X);
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
        for (REIter it = frame.r_edges.begin(); it != frame.r_edges.end(); ++it)
        {
            bool equiv = is_equal_occurrence(*it, *frame.sg);
            detect_nc(&frame, *it, equiv, EXT_R);
            if (is_min(shared, &frame, it->get_key(), is_min_last, ec_last))
            {
                detect_et(&frame, shared, *it, equiv, EXT_R);
                frame.children.push_back(RChild(&it->get_key(), &*it));
            }
        }
#else
        for (REIter it = frame.r_edges.begin(); it != frame.r_edges.end(); ++it)
        {
            bool equiv = is_equal_occurrence(*it, *frame.sg);
            detect_nc(&frame, *it, equiv, EXT_R);
            if (!frame.early_term)
            {
                if (is_min(shared, &frame, it->get_key(), is_min_last, ec_last))
                {
                    detect_et(&frame, shared, *it, equiv, EXT_R);
                    frame.children.push_back(RChild(&it->get_key(), &*it));
                }
            }
            if (frame.early_term && frame.closed)
                break;      
        }
#endif
        
        trace_frame(&frame, shared);
        result(&frame, shared);

        if (detect_et_2(shared, &frame))
            return;

        typedef typename RChildren::iterator ChIter;
        for (ChIter it = frame.children.begin(); it != frame.children.end(); ++it)
        {
            shared->dfsc.push(*it->first);
            project(shared, &frame, it->second);
            shared->dfsc.pop();

#ifdef CHECK_MODE
            if (! frame.TERMINATED_BRANCH && it->first == frame.ec_early_term)
                frame.early_term_frame = &frame;
#endif
        }
    }

    template<class SG>
    void run(SharedData* shared, Extension<SET, SG, EdgeCodeCmpDfs>& r_edges)
    {
        for (typename Extension<SET, SG, EdgeCodeCmpDfs>::iterator i = r_edges.begin(); i != r_edges.end(); ++i)
        {
            SG& sg = *i;
            if (sg.support() >= shared->minsup_)
            {
                shared->dfsc.push(i->get_key());
                project<SG>(shared, 0, &sg);
                shared->dfsc.pop();
            }
        }
    }

    void closegraph(const Graph& graph,
                    support_type minsup,
                    GspanResult* result,
                    int max_trace_depth)
    {
        SharedData shared(minsup, graph.num_vertices(), graph.num_edges(), result, max_trace_depth);
        Extension<SET, SubgraphsOfOneGraph, EdgeCodeCmpDfs> r_edges(&shared.sbg_creator);
        enum_one_edges(r_edges, graph);
        run(&shared, r_edges);
    }

    void closegraph(const std::vector<const Graph*>& graphs,
                    support_type minsup,
                    GspanResult* result,
                    int max_trace_depth)
    {
        std::size_t max_num_vertices = 0;
        std::size_t max_num_edges = 0;
        for (std::vector<const Graph*>::const_iterator i = graphs.begin(); i != graphs.end(); ++i)
        {
            max_num_vertices    = std::max(max_num_vertices, (*i)->num_vertices());
            max_num_edges       = std::max(max_num_edges, (*i)->num_vertices());
        }
        
        SharedData shared(minsup, max_num_vertices, max_num_edges, result, max_trace_depth);
        Extension<SET, SubgraphsOfManyGraph, EdgeCodeCmpDfs> r_edges(&shared.sbg_creator);
        for (std::vector<const Graph*>::const_iterator i = graphs.begin(); i != graphs.end(); ++i)
            enum_one_edges(r_edges, **i);
        run(&shared, r_edges);
    }




    void closegraph(const Graph* array, std::size_t size,
                    support_type minsup, GspanResult* result, int max_trace_depth)
    {
        if (size > 1)
        {

            std::size_t max_num_vertices = 0;
            std::size_t max_num_edges = 0;

            for (std::size_t i = 0; i < size; ++i)
            {
                const Graph* g = array + i;
                max_num_vertices    = std::max(max_num_vertices, g->num_vertices());
                max_num_edges       = std::max(max_num_edges, g->num_vertices());
            }
        
            SharedData shared(minsup, max_num_vertices, max_num_edges, result, max_trace_depth);
            Extension<SET, SubgraphsOfManyGraph, EdgeCodeCmpDfs> r_edges(&shared.sbg_creator);

            for (std::size_t i = 0; i < size; ++i)
                enum_one_edges(r_edges, array[i]);
            run(&shared, r_edges);
        }
        else
        {
            const Graph& graph = *array;
            SharedData shared(minsup, graph.num_vertices(), graph.num_edges(), result, max_trace_depth);
            Extension<SET, SubgraphsOfOneGraph, EdgeCodeCmpDfs> r_edges(&shared.sbg_creator);
            enum_one_edges(r_edges, graph);
            run(&shared, r_edges);
        }
    }


    void closegraph(const Graph* const * ptr_array, std::size_t size,
                    support_type minsup, GspanResult* result, int max_trace_depth)
    {
        if (size > 1)
        {

            std::size_t max_num_vertices = 0;
            std::size_t max_num_edges = 0;

            for (std::size_t i = 0; i < size; ++i)
            {
                const Graph* g = ptr_array[i];
                max_num_vertices    = std::max(max_num_vertices, g->num_vertices());
                max_num_edges       = std::max(max_num_edges, g->num_vertices());
            }
        
            SharedData shared(minsup, max_num_vertices, max_num_edges, result, max_trace_depth);
            Extension<SET, SubgraphsOfManyGraph, EdgeCodeCmpDfs> r_edges(&shared.sbg_creator);

            for (std::size_t i = 0; i < size; ++i)
                enum_one_edges(r_edges, *ptr_array[i]);
            run(&shared, r_edges);
        }
        else
        {
            const Graph& graph = **ptr_array;
            SharedData shared(minsup, graph.num_vertices(), graph.num_edges(), result, max_trace_depth);
            Extension<SET, SubgraphsOfOneGraph, EdgeCodeCmpDfs> r_edges(&shared.sbg_creator);
            enum_one_edges(r_edges, graph);
            run(&shared, r_edges);
        }
    }
}
