#include "gspan2.hpp"

#include <boost/foreach.hpp>

#include <cassert>
#include <list>
#include <deque>	// used in RMPath
#include <memory>
#include <map>
#include <queue>
#include <algorithm>

int nnncc = 0;

namespace gSpan2
{
    
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
    VI max_vertex(const DFSCode& dfsc)
    {
	VI m = 0;
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
    class RMPathSimple
    {
    protected:
	std::deque<int> rmp_;
    public:
	RMPathSimple() {}
	explicit RMPathSimple(const DFSCode& dfsc);
	
	int operator[] (int i) const { return rmp_[i]; }
	typedef std::deque<int>::const_iterator const_iterator;
	const_iterator begin() const { return rmp_.begin(); }
	const_iterator end() const { return rmp_.end(); }
	int size() const { return rmp_.size(); }
	int rightmost() const { return rmp_.back(); }
	friend std::ostream& operator<<(std::ostream& out, const RMPathSimple& r)
	    {
		std::copy(r.rmp_.begin(), r.rmp_.end(),
			  std::ostream_iterator<int>(out, " "));
		return out;
	    }
    };

    RMPathSimple::RMPathSimple(const DFSCode& dfsc)
    {
	VI old_src = 0;
	for (int i = dfsc.size()-1; i >= 0; --i)
	    if (dfsc[i].is_forward() && (rmp_.empty() || old_src == dfsc[i].vi_dst()))
	    {
		rmp_.push_front(i);
		old_src = dfsc[i].vi_src();
	    }
    }


    typedef std::vector<char> RMPathVertices;

    class RMPath : public RMPathSimple
    {
	RMPathVertices vv_;
    public:
	RMPath() {}
	explicit RMPath(const DFSCode& dfsc);
	const RMPathVertices& dfs_rmpath_v() const { return vv_; }
    };
    

    RMPath::RMPath(const DFSCode& dfsc)
    {
	vv_.reserve(dfsc.size() * 2U);
	VI old_src = 0;
	
	for (int i = dfsc.size()-1; i >= 0; --i)
	{
	    if (dfsc[i].is_forward() && (rmp_.empty() || old_src == dfsc[i].vi_dst()))
	    {
		rmp_.push_front(i);

		if (vv_.size() <= dfsc[i].vi_dst())
		    vv_.resize(dfsc[i].vi_dst() + 1, false);
		vv_[dfsc[i].vi_src()] = true;
		vv_[dfsc[i].vi_dst()] = true;

		old_src = dfsc[i].vi_src();
	    }
	}
    }

    
    // *****************************************************************************
    //                          SBG
    // *****************************************************************************

    void SBG::init_e_() const
    {
	if (prev_)
	{
	    e_ = new std::vector<const Graph::Edge*>(prev_->get_e());
	    e_->push_back(&edge_);
	}
	else
	    e_ = new std::vector<const Graph::Edge*>(1, &edge_);
    }

    
    void SBG::init_vv_() const
    {
	if (prev_)
	    vv_ = new std::vector<short>(prev_->get_vv());
	else
	    vv_ = new std::vector<short>(graph_->num_vertices(), 0);
	
	++(*vv_)[edge_.vi_src()];
	++(*vv_)[edge_.vi_dst()];
    }

    void SBG::init_ee_() const
    {
	if (prev_)
	{
	    ee_ = new VecEE(prev_->get_ee());
	    sum_ = prev_->sum_;
	}
	else
	{
	    ee_ = new VecEE(graph_->num_edges(), false);
	    sum_ = 0;
	}
	(*ee_)[edge_.eid()] = true;
	sum_ += edge_.eid(); 
    }


    void SBG::init_s2g_vv_() const
    {
	if (prev_)
	{
	    s2g_vv_ = new std::vector<VI>(prev_->get_s2g_vv());
	    if (vi_src_dfsc_ >= s2g_vv_->size()) s2g_vv_->push_back(edge_.vi_src());	
	    if (vi_dst_dfsc_ >= s2g_vv_->size()) s2g_vv_->push_back(edge_.vi_dst());
	}
	else
	{
	    s2g_vv_ = new std::vector<VI>(2, VI_NULL);
	    (*s2g_vv_)[0] = edge_.vi_src();
	    (*s2g_vv_)[1] = edge_.vi_dst();
	}
    }

    std::ostream& operator<<(std::ostream& out, const SBG& sbg)
    {
	sbg.init_e();
	out << "sbg:";
	for (int i = 0; i < sbg.size(); ++i)
	    out << " " << *sbg[i];
	//out << " at address=" << &sbg << " parent=" << sbg.parent();
	return out;
    }


    VI SBG::dfsc_vindex(VI sbg_vi, VI vi_default) const
    {
	const SBG* s = this;
	do
	{
	    if (sbg_vi == s->edge_.vi_src())
		return s->vi_src_dfsc_;
	    if (sbg_vi == s->edge_.vi_dst())
		return s->vi_dst_dfsc_;
	    s = s->parent();
	} while (s);
	return vi_default;
    }

    void SBG::graph_to_dfsc_v(std::vector<VI>& vv, VI vi_default) const
    {
	vv.resize(graph_->num_vertices(), vi_default);
	const SBG* s = this;
	do
	{
	    vv[s->edge_.vi_src()] = s->vi_src_dfsc_;
	    vv[s->edge_.vi_dst()] = s->vi_dst_dfsc_;
	    s = s->parent();
	} while (s);
    }


    const SBG* SBG::find_sbg_by_depth(int depth) const
    {
	assert(depth <= depth_);
	const SBG* s = this;
	do
	{
	    if (depth == s->depth_)
		return s;
	    s = s->parent();
	} while (s);
	return 0;
    }
    

    void make_sbg_rmpath_vertices(RMPathVertices& sbg_rmpath_v,
				  const RMPathVertices& dfs_rmpath_v, const SBG* sbg)
    {
	sbg_rmpath_v.clear();
	const std::vector<VI>& s2g_vv = sbg->get_s2g_vv();
	sbg_rmpath_v.resize(sbg->get_graph()->num_vertices(), false);
	for (unsigned int i = 0; i < dfs_rmpath_v.size(); ++i)
	    sbg_rmpath_v[s2g_vv[i]] = dfs_rmpath_v[i];
    }

    // *****************************************************************************
    //                          Projected
    // *****************************************************************************
    //
    // Projected implementation part
    //

    class ProjectedSimple
    {
	typedef std::vector<SBG*> SBGS;
	SBGS sbgs_;
	typedef typename SBGS::iterator iterator;
	void release() { for (iterator i = sbgs_.begin(); i != sbgs_.end(); ++i) delete *i; }
    public:
	typedef typename SBGS::const_iterator const_iterator;
	~ProjectedSimple() { release(); }

	const_iterator begin() const	{ return sbgs_.begin(); }
	const_iterator end()   const	{ return sbgs_.end(); }
	SBG* push(SBG* s)		{ sbgs_.push_back(s); return s; }
	SBG* back()			{ return sbgs_.back(); }
	const SBG* back() const		{ return sbgs_.back(); }
	int size() const		{ return sbgs_.size(); }
	void clear()			{ release(); sbgs_.clear(); }
	bool empty() const		{ return sbgs_.empty(); }
	void swap(ProjectedSimple& r)	{ sbgs_.swap(r.sbgs_); }
    };


    std::ostream& operator<<(std::ostream& out, const ProjectedSimple& p);

    template<class S>
    class Projected : public ProjectedSimple
    {
	S sg_;
	typedef ProjectedSimple Base;
	Base& base()			{ return *this; }
	const Base& base() const	{ return *this; }
    public:
	int support() const		{ return sg_.support(); }
	int size() const		{ return base().size(); }
	int num_sbgs_uniq() const	{ return sg_.num_sbgs_uniq(); }
	void push(SBG* s)		{ sg_.push(base().push(s)); }
	const S& sg() const		{ return sg_; }
    };


    int calc_supp(const std::vector<const SBG*>& ss)
    {
	int support = std::numeric_limits<int>::max();
	const int n_ver_g = ss.front()->get_graph()->num_vertices();
	const int n_ver_s = ss.front()->num_vertices();
	const int n_ss = ss.size();
	for (int vi = 0; vi < n_ver_s; ++vi)
	{
	    int n = 0;
	    std::vector<bool> vvg(n_ver_g, false);
	    for (int i = 0; i < n_ss; ++i)
	    {
		if (! vvg[ss[i]->get_s2g_vv()[vi]])
		{
		    ++n;
		    vvg[ss[i]->get_s2g_vv()[vi]] = true;
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
		s = s->automorph_next;
	    } while (s != sbg);
	}

	support_ = calc_supp(ss);
	support_valid_ = true;
    }

    
    bool less_cmp_sbg(const SBG* sbg1, const SBG* sbg2)
    {
	return sbg1->get_s2g_vv() < sbg2->get_s2g_vv();
    }

    void insert_list_next(SBG* pos, SBG* s)
    {
	pos->automorph_next->automorph_prev = s;
	s->automorph_next = pos->automorph_next;
	s->automorph_prev = pos;
	pos->automorph_next = s;
    }

    int list_size(const SBG* sbg)
    {
	int n = 0;
	const SBG* s = sbg;
	do 
	{
	    ++n;
	    s = s->automorph_next;
	} while (s != sbg);
	return n;
    }

    void SubgraphsOfOneGraph::insert_sbgs_uniq_ptr(SBGS_PTR& cont, SBG* s)
    {
	s->init_ee();
	BOOST_FOREACH(SBG* sbg, cont)
	{
	    if (*sbg == *s)
	    {
		insert_list_next(sbg, s);
		return;
	    }
	}
	cont.push_back(s);
    }

    void SubgraphsOfOneGraph::push(SBG* s)
    {
	assert(s->automorph_next == s);
	assert(s->automorph_prev == s);
	
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
		s = s->automorph_next;
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


    std::ostream& operator<<(std::ostream& out, const ProjectedSimple& p)
    {
	for (ProjectedSimple::const_iterator i = p.begin(); i != p.end(); ++i)
	    out << '\t' << *i << std::endl;
	return out;
    }

    template<class S>
    std::ostream& operator<<(std::ostream& out, const Projected<S>& p)
    {
	const ProjectedSimple& prj_simple = p;	
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
	return new_prj.sg().is_equal_occurrence(prj.sg());
    }


    // *****************************************************************************
    //                         ExtEdges
    // *****************************************************************************
    template<class P>
    class R_ExtEdges : public std::map<EdgeCode, P, EdgeCodeCmpDfs>
    {
    public:
	void insert(const EdgeCode& ec, SBG* sbg) { (*this)[ec].push(sbg); }
    };

    template<class P>
    class X_ExtEdges : public std::map<EdgeCode, P, EdgeCodeCmpLex>
    {
    public:
	void insert(const EdgeCode& ec, SBG* sbg) { (*this)[ec].push(sbg); }
    };
    
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
	ProjectedSimple prj_;
    public:
	void insert(const EdgeCode& ec, SBG* sbg);

	const EdgeCode& get_ec() const		{ return ec_; }
	const ProjectedSimple& get_prj() const	{ return prj_; }
	bool empty() const			{ return prj_.empty(); }
    };


    void MostMinEdgeCode::insert(const EdgeCode& ec, SBG* sbg)
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
		    delete sbg;
	    }
	}
    }


    // *****************************************************************************
    //                         functions
    // *****************************************************************************
    template<class ExtEdges>
    void enum_one_edges(ExtEdges& m, const Graph& g)
    {
	const Graph::Edges& g_edges = g.edges();
	for (Graph::EdgesIterator it = g_edges.begin(); it != g_edges.end(); ++it)
	{
	    m.insert(EdgeCode(0, 1, it->vl_src(), it->el(), it->vl_dst(), true),
		     new SBG(&g, *it));
	}
    }

    
    void enumerate_min_bck(MostMinEdgeCode& most_min_b, const DFSCode& dfsc_min,
			   const ProjectedSimple& prj, const RMPathSimple& rmpath, const Graph& g)
    {
	VI vi_dfsc_rmost = dfsc_min[rmpath.rightmost()].vi_dst();
	VL vl_rmost = dfsc_min[rmpath.rightmost()].vl_dst();
	bool flg = false;

	const int N = rmpath.size() - 1;
	for (int i = 0; !flg && i < N; ++i)
	{
	    const EdgeCode& ec_rmpath = dfsc_min[rmpath[i]];
	    EL el_rmpath = ec_rmpath.el();
	    bool vl_less_eq = ec_rmpath.vl_dst() <= vl_rmost;

	    BOOST_FOREACH(const SBG* psbg, prj)
	    {
		const SBG& sbg = *psbg;
		sbg.init_e();
		sbg.init_ee();
		sbg.init_vv();

		const Graph::Edge* e_rmost = sbg[rmpath.rightmost()];
		
		const Graph::IncidentEdges& incid_edges = g.incident(sbg[rmpath[i]]->vi_src());
		Graph::IncidentEdgesIterator it = incid_edges.begin();
		Graph::IncidentEdgesIterator it_end = incid_edges.end();
		for (; it != it_end; ++it)
		{
		    const Graph::Edge* pe = *it;
		    if (sbg.has_edge(pe->eid()))
			continue;

		    if (pe->vi_dst() == e_rmost->vi_dst() &&
			((vl_less_eq && el_rmpath == pe->el()) || el_rmpath < pe->el()))
		    {
			Graph::Edge e = *pe;
			e.chgdir();
			EdgeCode ec(vi_dfsc_rmost, ec_rmpath.vi_src(), vl_rmost, e.el(), ec_rmpath.vl_src(), false);
			most_min_b.insert(ec, new SBG(&sbg, e, vi_dfsc_rmost, ec_rmpath.vi_src()));
			flg = true;
			break;
		    }
		}
	    }
	}
    }


    void enumerate_min_fwd(MostMinEdgeCode& most_min_f, const DFSCode& dfsc_min,
			   const ProjectedSimple& prj, const RMPathSimple& rmpath, const Graph& g)
    {
	VI vi_dfsc_rmost = dfsc_min[rmpath.rightmost()].vi_dst();
	VL vl_rmost = dfsc_min[rmpath.rightmost()].vl_dst();
	VL vl_minimum = dfsc_min[0].vl_src();
	bool flg = false;

	// forward pure
	BOOST_FOREACH(const SBG* psbg, prj)
	{
	    const SBG& sbg = *psbg;
	    sbg.init_e();
	    sbg.init_vv();
	    sbg.init_ee();

	    const Graph::Edge* e_rmost = sbg[rmpath.rightmost()];

	    const Graph::IncidentEdges& incid_edges = g.incident(e_rmost->vi_dst());
	    Graph::IncidentEdgesIterator it = incid_edges.begin();
	    Graph::IncidentEdgesIterator it_end = incid_edges.end();
	    for (; it != it_end; ++it)
	    {
		const Graph::Edge* e = *it;
		if (sbg.has_edge(e->eid()))
		    continue;
		if (! sbg.has_vertex(e->vi_dst()) && vl_minimum <= e->vl_dst())
		{
		    EdgeCode ec(vi_dfsc_rmost, vi_dfsc_rmost+1, vl_rmost, e->el(), e->vl_dst(), true);
		    most_min_f.insert(ec, new SBG(&sbg, *e, vi_dfsc_rmost, vi_dfsc_rmost+1));
		    flg = true;
		}
	    }
	}

	// forward rmpath
	for (int i = rmpath.size()-1; !flg && i >= 0; --i)
	{
	    const EdgeCode& ec_rmpath = dfsc_min[rmpath[i]];
	    BOOST_FOREACH(const SBG* psbg, prj)
	    {
		const SBG& sbg = *psbg;

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
			most_min_f.insert(ec, new SBG(&sbg, *e, ec_rmpath.vi_src(), vi_dfsc_rmost+1));
			flg = true;
		    }
		}
	    }
	}

    }


    bool is_min_iterative(const DFSCode& dfsc_tested, bool debug = false)
    {
	std::auto_ptr<Graph> graph(new Graph(dfsc_tested.begin(), dfsc_tested.end(), max_vertex(dfsc_tested)+1));
	std::list<MostMinEdgeCode> exts;
	
	exts.push_back(MostMinEdgeCode());
	MostMinEdgeCode* exts1 = &exts.back();

	enum_one_edges(*exts1, *graph);
	DFSCode dfsc_min;

	while (true)
	{
	    dfsc_min.push_back(exts1->get_ec());
	    if (dfsc_min[dfsc_min.size()-1] != dfsc_tested[dfsc_min.size()-1])
		return false;

	    RMPathSimple rmpath(dfsc_min);

	    exts.push_back(MostMinEdgeCode());
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


    namespace closegraph_detail
    {
	enum ExtType { EXT_X, EXT_R, EXT_NONE };

#ifdef GSPAN_WITH_STATISTICS
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

	    typedef R_ExtEdges<P> REdges;
	    typedef X_ExtEdges<P> XEdges;
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
	void enumerate_v2(FrameState<P>* frame, const SharedData* shared)
	{
	    typedef R_ExtEdges<P> REdges;
	    typedef X_ExtEdges<P> XEdges;

	    REdges& r_edges = frame->r_edges;
	    XEdges& x_edges = frame->x_edges;
	    const P& prj = *frame->prj;
	    const RMPath& rmpath = frame->rmpath;
	    const RMPathVertices& dfs_rmpath_v = rmpath.dfs_rmpath_v();
	    const DFSCode& dfsc = shared->dfsc_;
	    const VI vi_dfsc_rmost   = dfsc[rmpath.rightmost()].vi_dst();
	    const VI vi_dfsc_new = vi_dfsc_rmost + 1;

	    BOOST_FOREACH(const SBG* psbg, prj)
	    {
		const SBG& sbg = *psbg;
		//sbg.init_e();
		sbg.init_vv();
		sbg.init_ee();

		const Graph& g = *sbg.get_graph();

		//
		// array of the graph VI, indexed by the dfsc VI
		// so, sbg_vertices[dfsc_vi] == graph vi
		//
		const std::vector<VI>& sbg_vertices = sbg.get_s2g_vv();
		const VI NUM_SBG_VERTICES = sbg_vertices.size();

		//
		// array of the dfsc VI, indexed by the graph VI
		// so, sbg_vertices[graph_vi] == dfsc vi
		//
		std::vector<VI> dfsc_vertices;
		sbg.graph_to_dfsc_v(dfsc_vertices, vi_dfsc_new);
		
	        RMPathVertices sbg_rmpath_v;
		make_sbg_rmpath_vertices(sbg_rmpath_v, dfs_rmpath_v, &sbg);
		
		const Graph::Edge* e_rmost = &sbg.find_sbg_by_depth(rmpath.rightmost() + 1)->edge();
		const VI vi_rmost = e_rmost->vi_dst();

		for (VI dfsc_vi = 0; dfsc_vi < NUM_SBG_VERTICES; ++dfsc_vi)
		{
		    VI graph_vi = sbg_vertices[dfsc_vi];
		    const bool from_rmpath = sbg_rmpath_v[graph_vi];

		    assert(sbg_vertices[dfsc_vi] == graph_vi);
		    assert(dfsc_vertices[graph_vi] == dfsc_vi);

		    const Graph::IncidentEdges& incid_edges = g.incident(graph_vi);
		    Graph::IncidentEdgesIterator it = incid_edges.begin();
		    Graph::IncidentEdgesIterator it_end = incid_edges.end();
		    for (; it != it_end; ++it)
		    {
			const Graph::Edge* e = *it;
			if (sbg.has_edge(e->eid()))
			    continue;

			// -----------------------------------
			// R edges
			// -----------------------------------
			// 1) vi_src is rmpath vertex AND
			// 2) vi_dst is new vertex (forward) OR 
			//     (vi_src is rmost vertex AND vi_dst is any rmpath vertex)
			if (from_rmpath)
			{
			    // vi_src is rmpath vertex

			    if (! sbg.has_vertex(e->vi_dst()))
			    {
				// vi_dst is new vertex
				// R forward
				r_edges[EdgeCode(dfsc_vi, vi_dfsc_new,
						 e->vl_src(), e->el(), e->vl_dst(), true)].push(
						     new SBG(&sbg, *e, dfsc_vi, vi_dfsc_new));
			    }
			    else if (e->vi_src() == vi_rmost && sbg_rmpath_v[e->vi_dst()])
			    {
				// vi_src is rmost vertex AND vi_dst is any rmpath vertex
				// R backward
				VI dfsc_vi_dst = dfsc_vertices[e->vi_dst()];
				assert(dfsc_vi_dst != vi_dfsc_new);
				EdgeCode ec(dfsc_vi, dfsc_vi_dst, e->vl_src(), e->el(), e->vl_dst(), false);
				r_edges[ec].push(new SBG(&sbg, *e, dfsc_vi, dfsc_vi_dst));
			    }
			    else
			    {
				if (! (e->vi_dst() == vi_rmost && sbg_rmpath_v[e->vi_src()]))
				{
				    // X backward
				    VI dfsc_vi_dst = dfsc_vertices[e->vi_dst()];
				    assert(dfsc_vi_dst != vi_dfsc_new);
				    EdgeCode ec(dfsc_vi, dfsc_vi_dst, e->vl_src(), e->el(), e->vl_dst(), false);
				    x_edges[ec].push(new SBG(&sbg, *e, dfsc_vi, dfsc_vi_dst));
				}
			    }
			}
			else
			{
			    // from vertex not on rmpath
			    VI dfsc_vi_dst = dfsc_vertices[e->vi_dst()];
			    EdgeCode ec(dfsc_vi, dfsc_vi_dst,
					e->vl_src(), e->el(), e->vl_dst(), ! sbg.has_vertex(e->vi_dst()));
			    x_edges[ec].push(new SBG(&sbg, *e, dfsc_vi, dfsc_vi_dst));
			}
		    }
		} // for (VI dfsc_vi = 0; dfsc_vi < NUM_SBG_VERTICES; ++dfsc_vi)
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
		(*shared->result_)(shared->dfsc_, frame->prj->sg());
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


	bool fail_et(const SBG* sbg, const RMPathVertices& dfs_rmpath_v,
		     DFSCode& dfsc, int i_rmost, ExtType exttype)
	{
	    sbg->init_vv();
	    sbg->init_ee();

	    const Graph& g = *sbg->get_graph();
	    std::vector<char> visited(g.num_vertices(), false);

	    RMPathVertices sbg_rmpath_v;
	    if (exttype == EXT_R)
		make_sbg_rmpath_vertices(sbg_rmpath_v, dfs_rmpath_v, sbg);
	    
	    std::queue<VI> q;

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
		    if (exttype == EXT_X || sbg_rmpath_v[e->vi_dst()])
			return true;
		}
		else
		    q.push(e->vi_dst());
	    }
	    visited[sbg->edge().vi_dst()] = true;

	    bool failure = false;
	    while (!failure && !q.empty())
	    {
		VI vi = q.front();
		q.pop();
		visited[vi] = true;


		const Graph::IncidentEdges& incid_edges = g.incident(vi);
		Graph::IncidentEdgesIterator it = incid_edges.begin();
		Graph::IncidentEdgesIterator it_end = incid_edges.end();
		for (; it != it_end; ++it)
		{
		    const Graph::Edge* e = *it;
		    VI vi2 = e->vi_dst();
		    if (sbg->has_vertex(vi2))
		    {
			if (exttype == EXT_X || sbg_rmpath_v[vi2])
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


	bool fail_et(const ProjectedSimple& prj,
		     const RMPathVertices& rmpath_dfs_vertices,
		     DFSCode& dfsc,
		     int i_rmost,
		     ExtType exttype)
	{
	    BOOST_FOREACH(const SBG* sbg, prj)
	    {
		if (fail_et(sbg, rmpath_dfs_vertices, dfsc, i_rmost, exttype))
		    return true;
	    }
	    return false;
	}


	bool find_f(EdgeCode& f, const DFSCode& dfsc, const RMPath& rmpath, VI vi)
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
		last_result = is_min_iterative(dfsc);
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
		    frame->early_term = ! fail_et(it->second, frame->rmpath.dfs_rmpath_v(),
						  shared->dfsc_, frame->rmpath.rightmost(), exttype);
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
	    typedef R_ExtEdges<P> REdges;
	    typedef X_ExtEdges<P> XEdges;
	    typedef typename REdges::const_iterator REcIter;
	    typedef typename XEdges::const_iterator XEcIter;

	    shared->on_enter_frame();
	    FrameState<P> frame(shared, parent_frame, prj);

	    if (shared->num_calls() == 9)
		;//BR;

	    enumerate_v2(&frame, shared);
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
	void run(SharedData* shared, const R_ExtEdges<P>& r_edges)
	{
	    for (typename R_ExtEdges<P>::const_iterator i = r_edges.begin(); i != r_edges.end(); ++i)
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

    } // end: namespace closegraph_detail


    void closegraph(const Graph& graph, int minsup, GspanResult* result, int max_trace_depth)
    {
	closegraph_detail::SharedData shared(minsup, result, max_trace_depth);
	typedef Projected<SubgraphsOfOneGraph> P;
	R_ExtEdges<P> r_edges;
	enum_one_edges(r_edges, graph);
	run(&shared, r_edges);
    }

    void closegraph(const std::vector<const Graph*>& graphs, int minsup, GspanResult* result, int max_trace_depth)
    {
	closegraph_detail::SharedData shared(minsup, result, max_trace_depth);
	typedef Projected<SubgraphsOfManyGraph> P;
	R_ExtEdges<P> r_edges;
	for (std::vector<const Graph*>::const_iterator i = graphs.begin(); i != graphs.end(); ++i)
	    enum_one_edges(r_edges, **i);
	run(&shared, r_edges);
    }
}
