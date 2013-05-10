#include "gspan.hpp"

#include <boost/foreach.hpp>

#include <list>
#include <deque>	// used in RMPath
#include <memory>
#include <queue>
#include <algorithm>
#include <cstring>	// memset()



unsigned int nctor;
unsigned int ndtor;

#define PREFETCH(addr)	asm("prefetcht0 %0\n" : :"m"((addr)))


namespace gSpan2
{
    enum ExtType { EXT_X, EXT_R, EXT_NONE };


    // *****************************************************************************
    //                          EdgeCode
    // *****************************************************************************
    bool EdgeCodeCmpDfs::operator() (const EdgeCode& ec1, const EdgeCode& ec2) const
    {
	bool ec1_f = ec1.is_forward();
	bool ec2_f = ec2.is_forward();

	if (!ec1_f && ec2_f)
	    return true;
	else if (!ec1_f && !ec2_f && ec1.vi_dst()  < ec2.vi_dst())
	    return true;
	else if (!ec1_f && !ec2_f && ec1.vi_dst() == ec2.vi_dst() && ec1.el() < ec2.el() )
	    return true;
	else if (ec1_f && ec2_f && ec1.vi_src()  > ec2.vi_src())
	    return true;
	else if (ec1_f && ec2_f && ec1.vi_src() == ec2.vi_src() && ec1.vl_src()  < ec2.vl_src())
	    return true;
	else if (ec1_f && ec2_f && ec1.vi_src() == ec2.vi_src() && ec1.vl_src() == ec2.vl_src() && ec1.el() < ec2.el())
	    return true;
	else if (ec1_f && ec2_f && ec1.vi_src() == ec2.vi_src() && ec1.vl_src() == ec2.vl_src() && ec1.el() == ec2.el() && ec1.vl_dst() < ec2.vl_dst())
	    return true;
	else
	    return false;
    }

    
    inline static bool cmp_lex_c(const EdgeCode& ec1, const EdgeCode& ec2)
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

    bool EdgeCodeCmpLex::operator() (const EdgeCode& ec1, const EdgeCode& ec2) const
    {
#ifdef USE_ASM
	uint8_t res;

	/* res = 0; a = 1; b = 2 */
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
	return cmp_lex_c(ec1, ec2);
#endif
    }


    bool EdgeCode::operator== (const EdgeCode& ec) const
    {
	return !
	    ((vi_src()^ec.vi_src()) | (vi_dst()^ec.vi_dst()) |
	     (vl_src()^ec.vl_src()) | (vl_dst()^ec.vl_dst()) | (el()^ec.el()));
	
	return
	    vi_src() == ec.vi_src() && vi_dst() == ec.vi_dst() &&
	    vl_src() == ec.vl_src() && vl_dst() == ec.vl_dst() && el() == ec.el();
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

    typedef std::deque<int>		RMPath_I_EdgeCode;
    typedef std::vector<DfscVI>		DfscVI_Vertices;

    // ----------------------------------------------------
    // dfsc_vertices[0        ... retvalue )			: rmpath vertices
    // dfsc_vertices[retvalue ... dfsc_vertices.size() )	: notrmpath vertices
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
	rmpath_i_ec.push_front(idx);

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
		    rmpath_i_ec.push_front(idx);
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


    class RMPath
    {
	RMPath_I_EdgeCode	rmpath_i_ec_;
	DfscVI_Vertices		dfsc_vertices_;
	std::size_t		n_rmp_;
    public:
	RMPath(const DFSCode& dfsc) { n_rmp_ = make_rmpath(rmpath_i_ec_, dfsc_vertices_, dfsc); }
	
	typedef DfscVI_Vertices::const_iterator const_iterator;
	
	const_iterator begin() const		{ return dfsc_vertices_.begin(); }
	const_iterator rmp_end() const		{ return dfsc_vertices_.begin() + n_rmp_; }
	const_iterator end() const		{ return dfsc_vertices_.end(); }

	std::size_t num_all_vertices() const	{ return dfsc_vertices_.size(); }
	std::size_t num_rmp_vertices() const	{ return n_rmp_; }


	int operator[] (int i) const	{ return rmpath_i_ec_[i]; }
	int size() const		{ return rmpath_i_ec_.size(); }
	int rightmost_edgeindex() const	{ return rmpath_i_ec_.back(); }

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
    //                          SBG
    // *****************************************************************************
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
	    if (s1.ev_array()[i].e_in_sbg != s2.ev_array()[i].e_in_sbg)
		return false;
	return true;
    }

    std::ostream& operator<<(std::ostream& out, const SBG& sbg)
    {
	std::vector<const SBG*> chain;
	get_chain(chain, &sbg);

	out << "sbg:";
	for (std::size_t i = 0; i < chain.size(); ++i)
	    out << " " << *chain[i];
	out << " at address=" << &sbg << " parent=" << sbg.parent();
	return out;
    }


    void get_chain(std::vector<const SBG*> chain, const SBG* sbg)
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
    
    class VertexRMPathStatus
    {
	unsigned char* status_;
	VertexRMPathStatus(const VertexRMPathStatus&);
	VertexRMPathStatus& operator= (const VertexRMPathStatus&);
    public:
	VertexRMPathStatus() :status_(0) {}
	VertexRMPathStatus(const RMPath& rmpath, const SBG* sbg) :status_(0) { create(rmpath, sbg); }
	~VertexRMPathStatus() { delete[] status_; }
	void create(const RMPath& rmpath, const SBG* sbg);
	bool operator[] (GraphVI vi) const	{ return status_[vi]; }
    };

    void VertexRMPathStatus::create(const RMPath& rmpath, const SBG* sbg)
    {
	assert(status_ == 0);
	status_ = new unsigned char [sbg->get_graph()->num_vertices()];
	const GraphVI* dfsc_to_graph_v = sbg->get_dfsc_to_graph_v();
	for (RMPath::const_iterator it = rmpath.begin(); it != rmpath.end(); ++it)
	{
	    assert(dfsc_to_graph_v[*it] >= 0);
	    assert(dfsc_to_graph_v[*it] < GraphVI(sbg->get_graph()->num_vertices()));
	    status_[dfsc_to_graph_v[*it]] = it < rmpath.rmp_end();
	}
    }
    
    // *****************************************************************************
    //                          SBG_Creator
    // *****************************************************************************
    class SBG_Creator : private boost::noncopyable
    {
	MemAllocator alloc_;

    public:

	// create root sbg
	SBGSimple* new_sbg_simple(const Graph::Edge& e, const Graph* g);
	SBG* new_sbg(const Graph::Edge& e, const Graph* g);
	
	// create other sbg
	SBGSimple* new_sbg(const Graph::Edge& e, const SBGSimple* s);
	SBG* new_sbg(const Graph::Edge& e, const SBG* s, const EdgeCode& ec);

	// delete sbg
	void delete_sbg(SBGSimple*);
	void delete_sbg(SBG*);


	// allocate and make
	// array of the dfsc VI, indexed by the graph VI
	// so, sbg_vertices[graph_vi] == dfsc vi
	//
	DfscVI* create_graph_to_dfsc_v(const SBG*,
				       DfscVI vi_default = VI_NULL);
	void destroy_graph_to_dfsc_v(DfscVI*, const SBG*);
    };

    // --------------------------------------------------
    // create root sbg
    // --------------------------------------------------

    SBGSimple* SBG_Creator::new_sbg_simple(const Graph::Edge& e, const Graph* g)
    {
	SBGSimple* p = new (alloc_.allocate(sizeof(SBGSimple))) SBGSimple(e, g);

	std::size_t ev_size = graph_size(g);
	p->ev_array_ = alloc_.alloc_array<EVBool>(ev_size);
	::memset(p->ev_array_, 0, ev_size * sizeof(EVBool));
	p->ev_array_[e.eid()].e_in_sbg = true;
	p->ev_array_[e.vi_src()].v_in_sbg = true;
	p->ev_array_[e.vi_dst()].v_in_sbg = true;

	p->edge_ptr_array_ = alloc_.alloc_array<Graph::Edge*>(1);
	p->edge_ptr_array_[0] = &p->edge_;

	return p;
    }

    SBG* SBG_Creator::new_sbg(const Graph::Edge& e, const Graph* g)
    {
	SBG* p = new (alloc_.allocate(sizeof(SBG))) SBG(e, g);

	std::size_t ev_size = graph_size(g);
	p->ev_array_ = alloc_.alloc_array<EVBool>(ev_size);
	::memset(p->ev_array_, 0, ev_size * sizeof(EVBool));
	p->ev_array_[e.eid()].e_in_sbg = true;
	p->ev_array_[e.vi_src()].v_in_sbg = true;
	p->ev_array_[e.vi_dst()].v_in_sbg = true;

	p->vi_dfsc_to_graph_ = alloc_.alloc_array<GraphVI>(2);
	p->vi_dfsc_to_graph_[0] = e.vi_src();
	p->vi_dfsc_to_graph_[1] = e.vi_dst();

	return p;
    }



    // --------------------------------------------------
    // create other sbg
    // --------------------------------------------------

    SBGSimple* SBG_Creator::new_sbg(const Graph::Edge& e, const SBGSimple* s)
    {
	PREFETCH(*s);
	SBGSimple* p = new (alloc_.allocate(sizeof(SBGSimple))) SBGSimple(e, s);

	std::size_t ev_size = graph_size(s->get_graph());
	p->ev_array_ = alloc_.alloc_array<EVBool>(ev_size);
	PREFETCH(*s->ev_array_);
	::memcpy(p->ev_array_, s->ev_array_, ev_size * sizeof(EVBool));
	p->ev_array_[e.eid()].e_in_sbg = true;
	p->ev_array_[e.vi_src()].v_in_sbg = true;
	p->ev_array_[e.vi_dst()].v_in_sbg = true;

	p->edge_ptr_array_ = alloc_.alloc_array<Graph::Edge*>(p->num_edges());
	PREFETCH(*s->edge_ptr_array_);	
	::memcpy(p->edge_ptr_array_, s->edge_ptr_array_, s->num_edges() * sizeof(Graph::Edge*));
	p->edge_ptr_array_[p->num_edges() - 1] = &p->edge_;

	return p;
    }

    SBG* SBG_Creator::new_sbg(const Graph::Edge& e, const SBG* s, const EdgeCode& ec)
    {
	PREFETCH(*s);
	SBG* p = new (alloc_.allocate(sizeof(SBG))) SBG(e, s, ec);
	
	std::size_t ev_size = graph_size(s->get_graph());
	p->ev_array_ = alloc_.alloc_array<EVBool>(ev_size);
	PREFETCH(*s->ev_array_);
	::memcpy(p->ev_array_, s->ev_array_, ev_size * sizeof(EVBool));
	p->ev_array_[e.eid()].e_in_sbg = true;
	p->ev_array_[e.vi_src()].v_in_sbg = true;
	p->ev_array_[e.vi_dst()].v_in_sbg = true;

	p->vi_dfsc_to_graph_ = alloc_.alloc_array<GraphVI>(p->num_vertices());
	PREFETCH(*s->vi_dfsc_to_graph_);
	::memcpy(p->vi_dfsc_to_graph_, s->vi_dfsc_to_graph_, s->num_vertices() * sizeof(GraphVI));
	p->vi_dfsc_to_graph_[p->vi_src_dfsc_] = e.vi_src();
	p->vi_dfsc_to_graph_[p->vi_dst_dfsc_] = e.vi_dst();

	return p;
    }


    // --------------------------------------------------
    // delete sbg
    // --------------------------------------------------

    void SBG_Creator::delete_sbg(SBGSimple* s)
    {
	alloc_.dealloc_array(s->edge_ptr_array_, s->num_edges());
	alloc_.dealloc_array(s->ev_array_, graph_size(s->get_graph()));
	alloc_.deallocate(s, sizeof(*s));
    }

    void SBG_Creator::delete_sbg(SBG* s)
    {
	alloc_.dealloc_array(s->vi_dfsc_to_graph_, s->num_vertices());
	alloc_.dealloc_array(s->ev_array_, graph_size(s->get_graph()));
	alloc_.deallocate(s, sizeof(*s));
    }


    DfscVI* SBG_Creator::create_graph_to_dfsc_v(const SBG* sbg,
						DfscVI vi_default)
    {
	std::size_t n = sbg->get_graph()->num_vertices();
	DfscVI* p = alloc_.alloc_array<DfscVI>(n);
	std::fill(p, p + n, vi_default);

	const SBG* s = sbg;
	do
	{
	    PREFETCH(*s->parent());

	    p[s->edge().vi_src()] = s->vi_src_dfsc_;
	    p[s->edge().vi_dst()] = s->vi_dst_dfsc_;
	    s = s->parent();
	} while (s);

	return p;
    }

    void SBG_Creator::destroy_graph_to_dfsc_v(DfscVI* p, const SBG* sbg)
    {
	alloc_.dealloc_array(p, sbg->get_graph()->num_vertices());
    }

    // *****************************************************************************
    //                          Projected
    // *****************************************************************************

    //
    // Projected implementation part
    //

    template<class S>
    class ProjectedSimple
    {
	typedef std::vector<S*> SBGS;
	SBGS sbgs_;
	typedef typename SBGS::iterator iterator;
    public:
	ProjectedSimple() {}

	typedef typename SBGS::const_iterator const_iterator;
	const_iterator begin() const	{ return sbgs_.begin(); }
	const_iterator end()   const	{ return sbgs_.end(); }
	S* push(S* s)			{ sbgs_.push_back(s); return s; }
	S* back()			{ return sbgs_.back(); }
	const S* back() const		{ return sbgs_.back(); }
	int size() const		{ return sbgs_.size(); }

	void release_sbgs(SBG_Creator* sbg_allocator);
	void clear()			{ sbgs_.clear(); }

	bool empty() const		{ return sbgs_.empty(); }
	void swap(ProjectedSimple& r)	{ sbgs_.swap(r.sbgs_); }
    };

    template<class S>
    void ProjectedSimple<S>::release_sbgs(SBG_Creator* sbg_allocator)
    {
	for (iterator i = sbgs_.begin(); i != sbgs_.end(); ++i)
	    sbg_allocator->delete_sbg(*i);
    }


    template<class Container>
    class Projected : public ProjectedSimple<SBG>
    {
	Container sg_;
	typedef ProjectedSimple Base;
	Base& base()			{ return *this; }
	const Base& base() const	{ return *this; }
    public:
	int support() const		{ return sg_.support(); }
	int size() const		{ return base().size(); }
	int num_sbgs_uniq() const	{ return sg_.num_sbgs_uniq(); }
	void push(SBG* s)		{ sg_.push(base().push(s)); }
	const Container& container() const	{ return sg_; }
    };


    int calc_supp(const std::vector<const SBG*>& ss)
    {
	int support = std::numeric_limits<int>::max();
	const GraphVI n_ver_g = ss.front()->get_graph()->num_vertices();
	const DfscVI n_ver_s = ss.front()->num_vertices();
	const int n_ss = ss.size();
	for (DfscVI vi = 0; vi < n_ver_s; ++vi)
	{
	    int n = 0;
	    std::vector<bool> vvg(n_ver_g, false);
	    for (int i = 0; i < n_ss; ++i)
	    {
		GraphVI graph_vi = ss[i]->get_dfsc_to_graph_v()[vi];

		assert(graph_vi >= 0);
		assert(graph_vi < n_ver_g);

		if (! vvg[graph_vi])
		{
		    ++n;
		    vvg[graph_vi] = true;
		}
	    }
	    if (n < support)
		support = n;
	}
	return support;
    }


    void SubgraphsOfOneGraph::calc_support_v() const
    {
	std::vector<const SBG*> ss;
	ss.reserve(sbgs_uniq_ptr_.size() * 2);
	BOOST_FOREACH(const SBG* sbg, sbgs_uniq_ptr_)
	{
	    const SBG* s = sbg;
	    do
	    {
		ss.push_back(s);
		s = s->next_automorph();
	    } while (s != sbg);
	}

	support_ = calc_supp(ss);
	support_valid_ = true;
    }
    
    int list_size(const SBG* sbg)
    {
	int n = 0;
	const SBG* s = sbg;
	do 
	{
	    ++n;
	    s = s->next_automorph();
	} while (s != sbg);
	return n;
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

    void SubgraphsOfOneGraph::push(SBG* s)
    {
	sbg_parents_.insert(s->parent());
	insert_sbgs_uniq_ptr(sbgs_uniq_ptr_, s);
    }
    
    bool SubgraphsOfOneGraph::is_equal_occurrence(const SubgraphsOfOneGraph& parent) const
    {
	int n = 0;
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

	return n == parent.num_sbgs_uniq();
    }

    void SubgraphsOfManyGraph::push(SBG* s)
    {
	SOG& sog = g_sog_[s->get_graph()];
	int old_num_sbgs_uniq = sog.num_sbgs_uniq();
	sog.push(s);
	int new_num_sbgs_uniq = sog.num_sbgs_uniq();
	num_sbgs_uniq_ += new_num_sbgs_uniq - old_num_sbgs_uniq;
    }

    bool SubgraphsOfManyGraph::is_equal_occurrence(const SubgraphsOfManyGraph& prj) const
    {
	return false;
    }


    template<class S>
    std::ostream& operator<<(std::ostream& out, const ProjectedSimple<S>& p)
    {
	for (typename ProjectedSimple<S>::const_iterator i = p.begin(); i != p.end(); ++i)
	    out << '\t' << **i << std::endl;
	return out;
    }

    template<class Container>
    std::ostream& operator<<(std::ostream& out, const Projected<Container>& p)
    {
	const ProjectedSimple<SBG>& prj_simple = p;	
	out << prj_simple
	    << "\t support="<<p.support()
	    << " size="<<p.size()
	    << " num_uniq=" << p.num_sbgs_uniq()
	    << std::endl;
	return out;
    }


    template<class P>
    inline bool is_equal_occurrence(const P& new_prj, const P& prj)
    {
	return new_prj.container().is_equal_occurrence(prj.container());
    }


    // *****************************************************************************
    //                         ExtEdges
    // *****************************************************************************
    template<class P, class C>
    class Extension : public std::map<EdgeCode, P, C>
    {
	SBG_Creator* sbg_allocator_;
	void release_sbgs();
    public:
	explicit Extension(SBG_Creator* sbg_allocator) :sbg_allocator_(sbg_allocator) {}
	~Extension() { release_sbgs(); }

	void insert(const EdgeCode& ec, const Graph::Edge& e, const Graph* g);
	void insert(const EdgeCode& ec, const Graph::Edge& e, const SBG* s);
    };

    template<class P, class C>
    void Extension<P,C>::release_sbgs()
    {
	for (typename std::map<EdgeCode, P, C>::iterator i = this->begin(); i != this->end(); ++i)
	    i->second.release_sbgs(sbg_allocator_);
    }

    template<class P, class C>
    void Extension<P,C>::insert(const EdgeCode& ec, const Graph::Edge& e, const Graph* g)
    {
	SBG* p = sbg_allocator_->new_sbg(e, g);
	
	typedef typename std::map<EdgeCode, P, C>::iterator It;
	std::map<EdgeCode, P, C>& m = *this;
	std::pair<It,bool> r = m.insert(std::pair<EdgeCode, P>(ec, P()));
	r.first->second.push(p);
	    
	//(*this)[ec].push(p);
    }

    template<class P, class C>
    void Extension<P,C>::insert(const EdgeCode& ec, const Graph::Edge& e, const SBG* s)
    {
	SBG* p = sbg_allocator_->new_sbg(e, s, ec);
	(*this)[ec].push(p);
    }

    
    template<class P, class C>
    std::ostream& operator<<(std::ostream& out, const std::map<EdgeCode, P, C>& m)
    {
	for (typename std::map<EdgeCode, P, C>::const_iterator it = m.begin(); it != m.end(); ++it)
	    out << it->first << "\n" << it->second;
	return out;
    }

    class MostMinEdgeCode
    {
	EdgeCode ec_;
	typedef ProjectedSimple<SBGSimple> PrjSimple;
	PrjSimple prj_;
	SBG_Creator* sbg_allocator_;
	void insert(const EdgeCode& ec, SBGSimple* sbg);
    public:
	MostMinEdgeCode(SBG_Creator* sbg_allocator)
	    :sbg_allocator_(sbg_allocator) {}

	~MostMinEdgeCode() { prj_.release_sbgs(sbg_allocator_); }

	void insert(const EdgeCode& ec, const Graph::Edge& e, const Graph*);
	void insert(const EdgeCode& ec, const Graph::Edge& e, const SBGSimple* s);

	const EdgeCode& get_ec() const		{ return ec_; }
	const PrjSimple& get_prj() const	{ return prj_; }
	bool empty() const			{ return prj_.empty(); }
    };

    void MostMinEdgeCode::insert(const EdgeCode& ec, SBGSimple* sbg)
    {
	if (prj_.empty())
	{
	    ec_ = ec;
	    prj_.push(sbg);
	}
	else
	{
	    EdgeCodeCmpDfs cmp;
	    if (cmp(ec, ec_))
	    {
		ec_ = ec;
		prj_.release_sbgs(sbg_allocator_);
		prj_.clear();
		prj_.push(sbg);
	    }
	    else
	    {
		if (ec == ec_)
		{
		    // equal
		    assert(! cmp(ec_, ec));
		    prj_.push(sbg);	
		}
		else
		    sbg_allocator_->delete_sbg(sbg);
	    }
	}
    }

    void MostMinEdgeCode::insert(const EdgeCode& ec, const Graph::Edge& e, const Graph*g)
    {
	SBGSimple* p = sbg_allocator_->new_sbg_simple(e, g);
	insert(ec, p);
    }

    void MostMinEdgeCode::insert(const EdgeCode& ec, const Graph::Edge& e, const SBGSimple* s)
    {
	SBGSimple* p = sbg_allocator_->new_sbg(e, s);
	insert(ec, p);
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
    struct SharedData
    {
	DFSCode dfsc_;
	const int minsup_;
	GspanResult* result_;
	    
	SBG_Creator sbg_allocator;

#ifdef GSPAN_TRACE
	int max_trace_depth_;
#endif

#ifdef GSPAN_WITH_STATISTICS
	Statistics statistics_;
	void on_enter_frame()		{ ++statistics_.num_project_calls; }
	unsigned int num_calls() const	{ return statistics_.num_project_calls; }
	bool on_ismin(bool ismin)
	    { ++statistics_.num_ismin_calls; if (ismin) ++statistics_.num_ismin_true_ret; return ismin; }
	void on_early_term_x_f()	{ ++statistics_.num_detect_early_termin_x_f; }
	void on_early_term_x_b()	{ ++statistics_.num_detect_early_termin_x_b; }
	void on_early_term_r_f()	{ ++statistics_.num_detect_early_termin_r_f; }
	void on_early_term_r_b()	{ ++statistics_.num_detect_early_termin_r_b; }
	void on_early_term_f(ExtType t)	{ if (t == EXT_X) on_early_term_x_f(); else on_early_term_r_f(); }
	void on_early_term_b(ExtType t)	{ if (t == EXT_X) on_early_term_x_b(); else on_early_term_r_b(); }
	void statistics_report() const	{ std::cerr << statistics_; }
#else
	bool on_ismin(bool ismin)		{ return ismin; }
	void on_early_term_x_f()	{}
	void on_early_term_x_b()	{}
	void on_early_term_r_f()	{}
	void on_early_term_r_b()	{}
	void on_early_term_f(ExtType t)	{}
	void on_early_term_b(ExtType t)	{}
	void statistics_report() const	{}
#endif

#if not defined(GSPAN_WITH_STATISTICS) && (defined(GSPAN_TRACE) || defined(CHECK_MODE))
	unsigned int num_project_calls;
	void on_enter_frame()		{ ++num_project_calls; }
	unsigned int num_calls() const	{ return num_project_calls; }
#endif
#if not defined(GSPAN_WITH_STATISTICS) && not defined(GSPAN_TRACE) && not defined(CHECK_MODE)
	void on_enter_frame()		{}
#endif

	SharedData(int minsup, GspanResult* result, unsigned int max_trace_depth)
	    : minsup_(minsup)
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

    template<class P>
    struct FrameState
    {
	FrameState* prev_state;
	const int depth;

	const P* prj;
	RMPath rmpath;

	typedef Extension<P, EdgeCodeCmpDfs> REdges;
	typedef Extension<P, EdgeCodeCmpLex> XEdges;
	typedef typename REdges::const_iterator REcIter;
	typedef typename XEdges::const_iterator XEcIter;
	typedef typename REdges::iterator REIter;
	typedef typename XEdges::iterator XEIter;

	XEdges x_edges;
	REdges r_edges;

	typedef std::vector<std::pair<const EdgeCode*, const P*> > RChildren;
	typedef typename RChildren::const_iterator RChildIter;
	RChildren children;
	    
	bool closed;
	bool early_term;

	// cause of not close
	ExtType exttype_notclose;
	const EdgeCode* ec_notclose;
	const P* prj_ext_notclose;
  
	// cause of early_termin
	ExtType exttype_early_term;
	const EdgeCode* ec_early_term;
	const P* prj_ext_early_term;

#if defined(GSPAN_TRACE) || defined(CHECK_MODE)
	unsigned int id;
#endif

#if defined(CHECK_MODE)
	const bool TERMINATED_BRANCH;
	const FrameState* early_term_frame;
#endif
	FrameState(SharedData* shared, FrameState* prev, const P* prj);
    };

    template<class P>
    FrameState<P>::FrameState(SharedData* shared, FrameState* prev, const P* prj)
	: prev_state(prev)
	, depth(shared->dfsc_.size())
	, prj(prj)
	, rmpath(shared->dfsc_)
	, x_edges(&shared->sbg_allocator)
	, r_edges(&shared->sbg_allocator)
	, closed(true)
	, early_term(false)
	, exttype_notclose(EXT_NONE)
	, ec_notclose(0)
	, prj_ext_notclose(0)
	, exttype_early_term(EXT_NONE)
	, ec_early_term(0)
	, prj_ext_early_term(0)
#if defined(CHECK_MODE)
	, TERMINATED_BRANCH(prev && prev->early_term_frame != 0)
	, early_term_frame(prev ? prev->early_term_frame : 0)
#endif
    {
#if defined(CHECK_MODE)
	if (prev)
	    assert(TERMINATED_BRANCH ? prev->early_term_frame!=0 : prev->early_term_frame==0);
#endif

#if defined(GSPAN_TRACE) || defined(CHECK_MODE)
	id = shared->num_calls();
#endif
    }


    // *****************************************************************************
    //                         functions
    // *****************************************************************************
    template<class ExtEdges>
    void enum_one_edges(ExtEdges& m, const Graph& g)
    {
	const Graph::Edges& g_edges = g.edges();
	for (Graph::EdgesIterator it = g_edges.begin(); it != g_edges.end(); ++it)
	    m.insert(EdgeCode(0, 1, it->vl_src(), it->el(), it->vl_dst(), true), *it, &g);
    }

    
    void enumerate_min_bck(MostMinEdgeCode& most_min_b, const DFSCode& dfsc_min,
			   const ProjectedSimple<SBGSimple>& prj, const RMPath& rmpath, const Graph& g)
    {
	DfscVI vi_dfsc_rmost = dfsc_min[rmpath.rightmost_edgeindex()].vi_dst();
	VL vl_rmost = dfsc_min[rmpath.rightmost_edgeindex()].vl_dst();
	bool flg = false;

	const int N = rmpath.size() - 1;
	for (int i = 0; !flg && i < N; ++i)
	{
	    const EdgeCode& ec_rmpath = dfsc_min[rmpath[i]];
	    EL el_rmpath = ec_rmpath.el();
	    bool vl_less_eq = ec_rmpath.vl_dst() <= vl_rmost;

	    BOOST_FOREACH(SBGSimple* psbg, prj)
	    {
		SBGSimple& sbg = *psbg;
		const Graph::Edge* e_rmost = sbg[rmpath.rightmost_edgeindex()];
		GraphVI graph_vi = sbg[rmpath[i]]->vi_src();

		if (sbg.has_no_extension(graph_vi))
		    continue;
		sbg.set_no_has_extension(graph_vi);

		const Graph::IncidentEdges& incid_edges = g.incident(graph_vi);
		Graph::IncidentEdgesIterator it = incid_edges.begin();
		Graph::IncidentEdgesIterator it_end = incid_edges.end();
		for (; it != it_end; ++it)
		{
		    const Graph::Edge* pe = *it;
		    if (sbg.has_edge(pe->eid()))
			continue;

		    sbg.set_has_extension(graph_vi);

		    if (pe->vi_dst() == e_rmost->vi_dst() &&
			((vl_less_eq && el_rmpath == pe->el()) || el_rmpath < pe->el()))
		    {
			Graph::Edge e = *pe;
			e.chgdir();
			EdgeCode ec(vi_dfsc_rmost, ec_rmpath.vi_src(), vl_rmost, e.el(), ec_rmpath.vl_src(), false);
			most_min_b.insert(ec, e, &sbg);
			flg = true;
			break;
		    }
		}
	    }
	}
    }


    void enumerate_min_fwd(MostMinEdgeCode& most_min_f, const DFSCode& dfsc_min,
			   const ProjectedSimple<SBGSimple>& prj, const RMPath& rmpath, const Graph& g)
    {
	DfscVI vi_dfsc_rmost = dfsc_min[rmpath.rightmost_edgeindex()].vi_dst();
	VL vl_rmost = dfsc_min[rmpath.rightmost_edgeindex()].vl_dst();
	VL vl_minimum = dfsc_min[0].vl_src();
	bool flg = false;

	// forward pure
	BOOST_FOREACH(SBGSimple* psbg, prj)
	{
	    SBGSimple& sbg = *psbg;
	    const Graph::Edge* e_rmost = sbg[rmpath.rightmost_edgeindex()];
	    GraphVI graph_vi = e_rmost->vi_dst();

	    if (sbg.has_no_extension(graph_vi))
		continue;
	    sbg.set_no_has_extension(graph_vi);

	    const Graph::IncidentEdges& incid_edges = g.incident(graph_vi);
	    Graph::IncidentEdgesIterator it = incid_edges.begin();
	    Graph::IncidentEdgesIterator it_end = incid_edges.end();
	    for (; it != it_end; ++it)
	    {
		const Graph::Edge* e = *it;
		if (sbg.has_edge(e->eid()))
		    continue;

		sbg.set_has_extension(graph_vi);

		if (! sbg.has_vertex(e->vi_dst()) && vl_minimum <= e->vl_dst())
		{
		    EdgeCode ec(vi_dfsc_rmost, vi_dfsc_rmost+1, vl_rmost, e->el(), e->vl_dst(), true);
		    most_min_f.insert(ec, *e, &sbg);
		    flg = true;
		}
	    }
	}

	// forward rmpath
	for (int i = rmpath.size()-1; !flg && i >= 0; --i)
	{
	    const EdgeCode& ec_rmpath = dfsc_min[rmpath[i]];
	    BOOST_FOREACH(const SBGSimple* psbg, prj)
	    {
		const SBGSimple& sbg = *psbg;

		const Graph::IncidentEdges& incid_edges = g.incident(sbg[rmpath[i]]->vi_src());
		Graph::IncidentEdgesIterator it = incid_edges.begin();
		Graph::IncidentEdgesIterator it_end = incid_edges.end();
		for (; it != it_end; ++it)
		{
		    const Graph::Edge* e = *it;
		    if (sbg.has_edge(e->eid()))
			continue;
		    if (!sbg.has_vertex(e->vi_dst()) && vl_minimum <= e->vl_dst() &&
			((ec_rmpath.vl_dst() <= e->vl_dst() && ec_rmpath.el() == e->el()) || ec_rmpath.el() < e->el()))
		    {
			EdgeCode ec(ec_rmpath.vi_src(), vi_dfsc_rmost+1, ec_rmpath.vl_src(), e->el(), e->vl_dst(), true);
			most_min_f.insert(ec, *e, &sbg);
			flg = true;
		    }
		}
	    }
	}

    }


    bool is_min_iterative(const DFSCode& dfsc_tested, SharedData* shared, bool debug = false)
    {
	std::auto_ptr<Graph> graph(new Graph(dfsc_tested.begin(), dfsc_tested.end(), max_vertex(dfsc_tested)+1));
	std::list<MostMinEdgeCode> exts;
	
	exts.push_back(MostMinEdgeCode(&shared->sbg_allocator));
	MostMinEdgeCode* exts1 = &exts.back();

	enum_one_edges(*exts1, *graph);
	DFSCode dfsc_min;

	while (true)
	{
	    dfsc_min.push_back(exts1->get_ec());
	    if (dfsc_min[dfsc_min.size()-1] != dfsc_tested[dfsc_min.size()-1])
		return false;

	    RMPath rmpath(dfsc_min);

	    exts.push_back(MostMinEdgeCode(&shared->sbg_allocator));
	    MostMinEdgeCode* exts2 = &exts.back();
	    
	    enumerate_min_bck(*exts2, dfsc_min, exts1->get_prj(), rmpath, *graph);
	    if (! exts2->empty())
	    {
		exts1 = exts2;
		continue;
	    }
	    
	    enumerate_min_fwd(*exts2, dfsc_min, exts1->get_prj(), rmpath, *graph);
	    if (! exts2->empty())
	    {
		exts1 = exts2;
		continue;
	    }

	    return true;   
	}
    }

#ifdef GSPAN_TRACE
    template<class P>
    void print_frame_trace(const FrameState<P>* frame, const SharedData* shared)
    {
	using namespace std;
	cerr << shared->dfsc_.size() << ":";
	for (unsigned int i = 0; i < shared->dfsc_.size(); ++i)
	    cerr << " ";
	cerr << shared->dfsc_.back()
	     << " support=" << frame->prj->support()
	     << " size=" << frame->prj->size()
	     << " num_uniq=" << frame->prj->num_sbgs_uniq()
	     << " frame=" << frame->id;
#ifdef CHECK_MODE
	if (frame->TERMINATED_BRANCH)
	{
	    assert(frame->early_term_frame && frame->early_term_frame != frame);
	    cerr << "; Was Terminated at frame=" << frame->early_term_frame->id << " with depth=" << frame->early_term_frame->depth;
	}
#endif
	if (frame->early_term)
	{
	    assert(frame->exttype_early_term == EXT_X || frame->exttype_early_term == EXT_R);
	    cerr << "; Detect Early Termination: " << (frame->exttype_early_term == EXT_X ? "X" : "R") << *frame->ec_early_term;
	}
	cerr << endl;
    }
#endif

#ifdef CHECK_MODE
    template<class P>
    bool is_child(const FrameState<P>* frame, const EdgeCode* ec)
    {
	typedef std::vector<std::pair<const EdgeCode*, const P*> > RChildren;
	typedef typename RChildren::const_iterator RChildIter;
	for (RChildIter it = frame->children.begin(); it != frame->children.end(); ++it)
	    if (it->first == ec)
		return true;
	return false;
    }
	
    template<class P>
    void print_fail_et_report(const FrameState<P>* frame, const SharedData* shared)
    {
	typedef R_ExtEdges<P> REdges;
	typedef X_ExtEdges<P> XEdges;
	typedef typename REdges::const_iterator REcIter;
	typedef typename XEdges::const_iterator XEcIter;
	using namespace std;

	cerr << "Missed pattern detected" << endl;

	cerr << "============ MISSED FRAME: " << frame->id << " =======================================" << endl;
	cerr << "DFSC:\n";
	cerr << "--------------------------\n";
	for (int i = 0; i < frame->depth; ++i)
	    cerr << i+1 << ":\t"
		 << (find(frame->rmpath.begin(),frame->rmpath.end(), i) != frame->rmpath.end() ? "*" : " ")
		 << shared->dfsc_[i]
		 << endl;
	cerr << "--------------------------\n";
	cerr << "Embeddings:\n";
	cerr << *frame->prj;
	cerr << "--------------------------\n";
	cerr << "XExtentions:\n";
	    
	for (XEcIter it = frame->x_edges.begin(); it != frame->x_edges.end(); ++it)
	    cerr << it->first << "\n" << it->second;
	    
	cerr << endl;
	const FrameState<P>* fet_frame = frame->early_term_frame;
	cerr << "============ FAIL Eearly Termination FRAME: " << fet_frame->id << " ======================" << endl;
	cerr << "DFSC:\n";
	cerr << "--------------------------\n";
	for (int i = 0; i < fet_frame->depth; ++i)
	    cerr << i+1 << ":\t"
		 << (find(fet_frame->rmpath.begin(),fet_frame->rmpath.end(), i) != fet_frame->rmpath.end() ? "*" : " ")
		 << shared->dfsc_[i]
		 << endl;
	cerr << "--------------------------\n";
	cerr << "Embeddings:\n";
	cerr << *fet_frame->prj;
	cerr << "--------------------------\n";
	cerr << "Cause Early Termination:\n";
	cerr << (fet_frame->exttype_early_term == EXT_X ? "X" : "R") << *fet_frame->ec_early_term << endl;
	cerr << *fet_frame->prj_ext_early_term;
	cerr << "--------------------------\n";
	cerr << "RExtensions:\n";
	const EdgeCode& ec_ext = shared->dfsc_[fet_frame->depth];
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
    template<class P>
    void debug_print(const FrameState<P>* frame, const SharedData* shared)
    {
	std::cerr << "===================================================================\n";
	std::cerr << "Frame: " << shared->num_calls() << std::endl;
	std::cerr << "DFSC:\n"
		  << shared->dfsc_;
	std::cerr << "REDGES---------------------------\n" << frame->r_edges;
	std::cerr << "XEDGES---------------------------\n" << frame->x_edges;
    }

    // *****************************************************************************
    //                         functions
    // *****************************************************************************

    template<class P>
    void enumerate_v3(FrameState<P>* frame, SharedData* shared)
    {
	typedef Extension<P, EdgeCodeCmpDfs> REdges;
	typedef Extension<P, EdgeCodeCmpLex> XEdges;

	REdges& r_edges = frame->r_edges;
	XEdges& x_edges = frame->x_edges;
	const P& prj = *frame->prj;
	const RMPath& rmpath = frame->rmpath;
	const DfscVI vi_dfsc_rmost   = rmpath.rightmost_vertex();
	const DfscVI vi_dfsc_new = vi_dfsc_rmost + 1;

	assert(rmpath.rightmost_vertex() == shared->dfsc_[rmpath.rightmost_edgeindex()].vi_dst());

	BOOST_FOREACH(SBG* sbg, prj)
	{
	    const Graph* g = sbg->get_graph();

	    //
	    // array of the graph VI, indexed by the dfsc VI
	    // so, sbg_vertices[dfsc_vi] == graph vi
	    //
	    const GraphVI* dfsc_to_graph_v = sbg->get_dfsc_to_graph_v();

	    //
	    // array of the dfsc VI, indexed by the graph VI
	    // so, sbg_vertices[graph_vi] == dfsc vi
	    //
	    DfscVI* graph_to_dfsc_v = 0; // created by first usage

	    VertexRMPathStatus vertex_rmpath_status(rmpath, sbg);

	    const GraphVI vi_rmost = dfsc_to_graph_v[rmpath.rightmost_vertex()];
	    
	    for (RMPath::const_iterator dfsc_vi_iter = rmpath.begin();
		 dfsc_vi_iter != rmpath.end(); ++dfsc_vi_iter)
	    {
		const DfscVI dfsc_vi = *dfsc_vi_iter;
		const GraphVI graph_vi = dfsc_to_graph_v[dfsc_vi];
		const bool rmpath_vertex = dfsc_vi_iter < rmpath.rmp_end();

		assert(graph_to_dfsc_v[graph_vi] == dfsc_vi);
		assert(graph_vi >= 0 && graph_vi < sbg->get_graph()->num_vertices());
		assert(graph_vi < SBG::ev_array_size(&g));

		if (sbg->has_no_extension(graph_vi))
		    continue;
		sbg->set_no_has_extension(graph_vi);
				
		const Graph::IncidentEdges& incid_edges = g->incident(graph_vi);
		Graph::IncidentEdgesIterator it = incid_edges.begin();
		const Graph::IncidentEdgesIterator it_end = incid_edges.end();
		for (; it != it_end; ++it)
		{
		    const Graph::Edge* e = *it;
		    if (sbg->has_edge(e->eid()))
			continue;
		    
		    sbg->set_has_extension(graph_vi);
		    
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

			if (! sbg->has_vertex(e->vi_dst()))
			{
			    // vi_dst is new vertex
			    // R forward
			    EdgeCode ec(dfsc_vi, vi_dfsc_new, e->vl_src(), e->el(), e->vl_dst(), true);
			    r_edges.insert(ec, *e, sbg);
			}
			else if (e->vi_src() == vi_rmost && vertex_rmpath_status[e->vi_dst()])
			{
			    // vi_src is rmost vertex AND vi_dst is any rmpath vertex
			    // R backward
			    if (!graph_to_dfsc_v)
				graph_to_dfsc_v =
				    shared->sbg_allocator.create_graph_to_dfsc_v(sbg, vi_dfsc_new);

			    DfscVI dfsc_vi_dst = graph_to_dfsc_v[e->vi_dst()];
			    assert(dfsc_vi_dst != vi_dfsc_new);
			    EdgeCode ec(dfsc_vi, dfsc_vi_dst, e->vl_src(), e->el(), e->vl_dst(), false);
			    r_edges.insert(ec, *e, sbg);
			}
			else
			{
			    if (! (e->vi_dst() == vi_rmost && vertex_rmpath_status[e->vi_src()]))
			    {
				// X backward
				if (!graph_to_dfsc_v)
				    graph_to_dfsc_v =
					shared->sbg_allocator.create_graph_to_dfsc_v(sbg, vi_dfsc_new);

				DfscVI dfsc_vi_dst = graph_to_dfsc_v[e->vi_dst()];
				assert(dfsc_vi_dst != vi_dfsc_new);
				EdgeCode ec(dfsc_vi, dfsc_vi_dst, e->vl_src(), e->el(), e->vl_dst(), false);
				x_edges.insert(ec, *e, sbg);
			    }
			}
		    }
		    else
		    {
			// from vertex not on rmpath
			if (!graph_to_dfsc_v)
			    graph_to_dfsc_v = shared->sbg_allocator.create_graph_to_dfsc_v(sbg, vi_dfsc_new);
			    
			DfscVI dfsc_vi_dst = graph_to_dfsc_v[e->vi_dst()];
			EdgeCode ec(dfsc_vi, dfsc_vi_dst,
				    e->vl_src(), e->el(), e->vl_dst(), ! sbg->has_vertex(e->vi_dst()));
			x_edges.insert(ec, *e, sbg);
		    }
		}
	    }
	    
	    if (graph_to_dfsc_v)
		shared->sbg_allocator.destroy_graph_to_dfsc_v(graph_to_dfsc_v, sbg);
	}

    }

    template<class ExtEdges>
    void remove_not_frequents(ExtEdges& edges, int minsup)
    {
	if (minsup > 1)
	    for (typename ExtEdges::iterator it = edges.begin(); it != edges.end();)
		if (it->second.support() < minsup)
		    edges.erase(it++);
		else
		    ++it;
    }
	

    template<class P>
    void trace_frame(FrameState<P>* frame, SharedData* shared)
    {
#ifdef GSPAN_TRACE
	if (frame->depth < shared->max_trace_depth_)
	    print_frame_trace(frame, shared);
#endif
    }
	
    template<class P>
    void result(FrameState<P>* frame, SharedData* shared)
    {
	if (frame->closed)
	{
	    (*shared->result_)(shared->dfsc_, frame->prj->container());
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


    bool fail_et(const SBG* sbg, DFSCode& dfsc, const RMPath& rmpath, ExtType exttype)
    {
	const Graph& g = *sbg->get_graph();
	std::vector<char> visited(g.num_vertices(), false);

	VertexRMPathStatus vertex_rmpath_status;
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


    bool fail_et(const ProjectedSimple<SBG>& prj, DFSCode& dfsc, const RMPath& rmpath, ExtType exttype)
    {
	BOOST_FOREACH(const SBG* sbg, prj)
	{
	    if (fail_et(sbg, dfsc, rmpath, exttype))
		return true;
	}
	return false;
    }


    bool find_f(EdgeCode& f, const DFSCode& dfsc, const RMPath& rmpath, DfscVI vi)
    {
	int n = dfsc.size();
	for (int i = 0; i < n; ++i)
	{
	    if (dfsc[i].vi_src() == vi && std::find(rmpath.begin(), rmpath.end(), i) != rmpath.end())
	    {
		f = dfsc[i];
		return true;
	    }
	}
	return false;
    }

    bool is_min_fast(SharedData* shared,
		     const DFSCode& dfsc, EdgeCode ec, const RMPath& rmpath, bool& last_result, EdgeCode& ec_last)
    {
	EdgeCode f;
	if (last_result && ec.is_forward() && ec_last.is_forward() && ec.vi_src() == ec_last.vi_src())
	    return true;
	else if (ec.is_forward() && find_f(f, dfsc, rmpath, ec.vi_dst()))
	{
	    EdgeCodeCmpDfs cmp;
	    return last_result = ! cmp(ec, f);
	}
	else
	{
	    last_result = is_min_iterative(dfsc, shared);
	    shared->on_ismin(last_result);
	    return last_result;
	}
    }

    template<class P>
    bool is_min(SharedData* shared, FrameState<P>* frame, const EdgeCode& ec, bool& last_result, EdgeCode& ec_last)
    {
	shared->dfsc_.push_back(ec);
	bool r1 = is_min_fast(shared, shared->dfsc_, ec, frame->rmpath, last_result, ec_last);
/*
  bool r2 = is_min_recursive(shared->dfsc_);
  if (r1 != r2)
  {
  std::cerr << "------ Testing is_min_fast ------------\n";
  std::cerr << " is_min_fast: " << r1 << std::endl;
  std::cerr << "      is_min: " << r2 << std::endl;
  std::cerr << "DFSCode:\n" << shared->dfsc_ << std::endl;
  std::cerr << "RMPath:\n" << frame->rmpath << std::endl;
  std::cerr << "Last Result: " << ec_last << std::endl;
  is_min_recursive(shared->dfsc_, true);
  exit(1);
  }
*/
	if (r1)
	    ec_last = ec;
	shared->dfsc_.pop_back();
	return r1;
    }

    template<class P, class ExtIter>
    void detect_nc(FrameState<P>* frame, SharedData* shared, const ExtIter& it, bool equiv, ExtType exttype)
    {
	if (frame->closed && (equiv || it->second.support() == frame->prj->support()))
	{
	    frame->closed = false;
	    assert(frame->exttype_notclose == EXT_NONE);
	    frame->exttype_notclose = exttype;
	    frame->ec_notclose = &it->first;
	    frame->prj_ext_notclose = &it->second;
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
		frame->early_term = ! fail_et(it->second, shared->dfsc_, frame->rmpath, exttype);
		if (frame->early_term)
		    shared->on_early_term_f(exttype);
	    }
	    else
	    {
		frame->early_term = it->second.size() == frame->prj->size();
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
		frame->prj_ext_early_term = &it->second;
	    }
	}
    }

    template<class P>
    void project(SharedData* shared, FrameState<P>* parent_frame, const P* prj)
    {
	typedef Extension<P, EdgeCodeCmpDfs> REdges;
	typedef Extension<P, EdgeCodeCmpLex> XEdges;

	typedef typename REdges::const_iterator REcIter;
	typedef typename XEdges::const_iterator XEcIter;

	shared->on_enter_frame();
	FrameState<P> frame(shared, parent_frame, prj);

	enumerate_v3(&frame, shared);
	
	remove_not_frequents(frame.r_edges, shared->minsup_);

	for (XEcIter it = frame.x_edges.begin(); it != frame.x_edges.end(); ++it)
	{
	    bool equiv = is_equal_occurrence(it->second, *frame.prj);
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
	    bool equiv = is_equal_occurrence(it->second, *frame.prj);
	    detect_nc(&frame, shared, it, equiv, EXT_R);
	    if (is_min(shared, &frame, it->first, is_min_last, ec_last))
	    {
		detect_et(&frame, shared, it, equiv, EXT_R);
		frame.children.push_back(std::pair<const EdgeCode*, const P*>(&it->first, &it->second));
	    }
	}
#else
	for (REcIter it = frame.r_edges.begin(); it != frame.r_edges.end(); ++it)
	{
	    bool equiv = is_equal_occurrence(it->second, *frame.prj);
	    detect_nc(&frame, shared, it, equiv, EXT_R);
	    if (!frame.early_term)
	    {
		if (is_min(shared, &frame, it->first, is_min_last, ec_last))
		{
		    detect_et(&frame, shared, it, equiv, EXT_R);
		    frame.children.push_back(std::pair<const EdgeCode*, const P*>(&it->first, &it->second));
		}
	    }
	    if (frame.early_term && frame.closed)
		break;	    
	}
#endif

	trace_frame(&frame, shared);
	result(&frame, shared);
	    
	typedef typename std::vector<std::pair<const EdgeCode*, const P*> >::const_iterator ChIter;
	for (ChIter it = frame.children.begin(); it != frame.children.end(); ++it)
	{
	    shared->dfsc_.push_back(*it->first);
	    project(shared, &frame, it->second);
	    shared->dfsc_.pop_back();

#ifdef CHECK_MODE
	    if (! frame.TERMINATED_BRANCH && it->first == frame.ec_early_term)
		frame.early_term_frame = &frame;
#endif
	}
    }
	
	
    template<class P>
    void run(SharedData* shared, const Extension<P, EdgeCodeCmpDfs>& r_edges)
    {
	for (typename Extension<P, EdgeCodeCmpDfs>::const_iterator i = r_edges.begin(); i != r_edges.end(); ++i)
	{
	    const P& prj = i->second;
	    if (prj.support() >= shared->minsup_)
	    {
		shared->dfsc_.push_back(i->first);
		project<P>(shared, 0, &prj);
		shared->dfsc_.pop_back();
	    }
	}
    }


    void closegraph(const Graph& graph, int minsup, GspanResult* result, int max_trace_depth)
    {
/*
	std::cerr << "SBG SIZE= " << sizeof(SBG) << std::endl;
	std::cerr << "= " << offsetof(SBG, automorph_next_);
	std::cerr << std::endl;
	return;
*/

	SharedData shared(minsup, result, max_trace_depth);
	typedef Projected<SubgraphsOfOneGraph> P;
	Extension<P, EdgeCodeCmpDfs> r_edges(&shared.sbg_allocator);
	enum_one_edges(r_edges, graph);
	run(&shared, r_edges);
    }

    void closegraph(const std::vector<const Graph*>& graphs, int minsup, GspanResult* result, int max_trace_depth)
    {
	SharedData shared(minsup, result, max_trace_depth);
	typedef Projected<SubgraphsOfManyGraph> P;
	Extension<P, EdgeCodeCmpDfs> r_edges(&shared.sbg_allocator);
	for (std::vector<const Graph*>::const_iterator i = graphs.begin(); i != graphs.end(); ++i)
	    enum_one_edges(r_edges, **i);
	run(&shared, r_edges);
    }
}
