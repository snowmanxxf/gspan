#ifndef GSPAN_H_
#define GSPAN_H_

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

    template<class Policy>
    bool operator== (const EdgeCode<Policy>& ec1, const EdgeCode<Policy>& ec2)
    {
	return
	    ec1.vi_from == ec2.vi_from && ec1.vi_to == ec2.vi_to &&
	    ec1.vl_from == ec2.vl_from && ec1.vl_to == ec2.vl_to &&
	    ec1.el == ec2.el;
    }

    template<class Policy>
    bool operator!= (const EdgeCode<Policy>& ec1, const EdgeCode<Policy>& ec2) { return ! (ec1 == ec2); }

    template<class Policy>
    std::ostream& operator<<(std::ostream& out, const EdgeCode<Policy>& ec)
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
    //                          SBG
    // *****************************************************************************
    template<class Policy>
    class SBG
    {
    public:
	typedef typename Policy::Edge		Edge;
	typedef typename Policy::graph_t        G;
	typedef typename Policy::vertex_index_t VI;
	typedef typename Policy::edge_index_t   EI;

	SBG(const G* g, const Edge& e)
	    :prev_(0), edge_(e), rec_(0), graph_(g), depth_(1) { automorph_list = 0; }

	SBG(const SBG* s, const Edge& e)
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
	void init_rec() const				{ if (!rec_) init_rec_(); }
	const Edge* operator[] (int i) const		{ init_rec(); return rec_->e_[i]; }
	bool has_vertex(VI vi) const		        { init_rec(); return rec_->has_vertex(vi); }
	bool has_edge(EI ei) const       	        { init_rec(); return rec_->has_edge(ei); }
	const G* get_graph() const { return graph_; }
	const Edge& edge() const { return edge_; }
	const SBG* prev_sbg() const { return prev_; }
	SBG* automorph_list;

	static int num_common_edges(const SBG<Policy>& sbg1, const SBG<Policy>& sbg2);
    private:
	const SBG* prev_;
	Edge edge_;

	class R
	{
	public:
	    std::vector<const Edge*> e_;
	    std::vector<short> vv_;
	    std::vector<char> ee_;
	    R(const Edge& e, const G& g);
	    R(const Edge& e, const R* prev);
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
    SBG<Policy>::R::R(const Edge& e, const G& g)
	: e_(1, &e), vv_(Policy::nvertices(g), 0), ee_(Policy::nedges(g), false)
    {
	++vv_[e.vi_from];
	++vv_[e.vi_to];
	ee_[e.ei] = true;
    }

    template<class Policy>    
    SBG<Policy>::R::R(const Edge& e, const R* prev)
	: e_(prev->e_), vv_(prev->vv_), ee_(prev->ee_)
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
	    out << " " << *sbg[i];
	out << " at address: " << &sbg << " of the graph: " << sbg.get_graph();
	return out;
    }
 
    template<class Policy>
    int SBG<Policy>::num_common_edges(const SBG<Policy>& sbg1, const SBG<Policy>& sbg2)
    {
	sbg1.init_rec();
	sbg2.init_rec();
	const std::vector<char>& ee1 = sbg1.rec_->ee_;
	const std::vector<char>& ee2 = sbg2.rec_->ee_;
	int nce = 0;
	int num = ee1.size();
	for (int i = 0; i < num; ++i)
	{
	    if (ee1[i] && ee2[i])
		++nce;
	}
	return nce;
    }

    // *****************************************************************************
    //                          Projected
    // *****************************************************************************

    //
    // Projected part visible for user
    //

    template<class Policy>
    class SubgraphsOfOneGraph
    {	
	typedef std::vector<SBG<Policy>* > SBGS_PTR;
	SBGS_PTR sbgs_ptr_;
	int count_non_overlapped_;
	int num_sbgs_;
    public:
	SubgraphsOfOneGraph() :count_non_overlapped_(0), num_sbgs_(0) {}
	typedef typename SBGS_PTR::const_iterator const_iterator;
	const_iterator begin() const { return sbgs_ptr_.begin(); }
	const_iterator end()   const { return sbgs_ptr_.end(); }
	int support() const	{ return count_non_overlapped_; }
	int num_sbgs() const	{ return num_sbgs_; }
	void push(SBG<Policy>* s);
	void merge(SubgraphsOfOneGraph& sog);
    };

    template<class Policy>
    void SubgraphsOfOneGraph<Policy>::push(SBG<Policy>* s)
    {

	++num_sbgs_;
	const int ne = s->size();

	int count = 1;
	BOOST_FOREACH(SBG<Policy>* sbg, sbgs_ptr_)
	{
	    int nce = SBG<Policy>::num_common_edges(*s, *sbg);
	    if (nce == ne)
	    {
		assert(s->automorph_list == 0);
		s->automorph_list = sbg->automorph_list;
		sbg->automorph_list = s;
		return;
	    }

	    if (nce == 0)
		++count;
	}
	
	sbgs_ptr_.push_back(s);
	count_non_overlapped_ = count;

	//count_non_overlapped_ = 1;
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

    namespace detail
    {
	//
	// Projected implementation part
	//

	template<class Policy>
	class ProjectedSimple
	{
	    typedef std::list<SBG<Policy> > SBGS;
	    SBGS sbgs_;
	    int size_;
	public:
	    ProjectedSimple() :size_(0) {}
	    typedef typename SBGS::const_iterator const_iterator;
	    const_iterator begin() const	{ return sbgs_.begin(); }
	    const_iterator end()   const	{ return sbgs_.end(); }
	    SBG<Policy>& push(const SBG<Policy>& s) { sbgs_.push_back(s); ++size_; return sbgs_.back(); }
	    SBG<Policy>& back() { return sbgs_.back(); }
	    const SBG<Policy>& back() const { return sbgs_.back(); }
	    int size() const { return size_; }
	    void merge(ProjectedSimple& p) { size_ += p.size_; sbgs_.splice(sbgs_.begin(), p.sbgs_); p.size_ = 0; }
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
	//                          Maps
	// *****************************************************************************
	template<class Policy, class P>
	class MapTraits
	{
	    typedef typename Policy::vertex_index_t VI;
	    typedef typename Policy::vertex_label_t VL;
	    typedef typename Policy::edge_label_t EL;
	public:
	    typedef std::map<EL, P>			Map_EL_P;
	    typedef std::map<VL, P>			Map_VL_P;
	    typedef std::map<EL, Map_VL_P>		Map_EL_VL_P;
	    typedef std::map<VI, Map_EL_P>		Map_VI_EL_P;
	    typedef std::map<VI, Map_EL_VL_P>		Map_VI_EL_VL_P;
	    typedef std::map<VL, Map_EL_VL_P>		Map_VL_EL_VL_P;
	    typedef std::map<VI, Map_VL_EL_VL_P>	Map_VI_VL_EL_VL_P;
	    typedef std::map<VI, Map_VI_VL_EL_VL_P>	Map_VI_VI_VL_EL_VL_P;
	};
    
	// *****************************************************************************
	//                         functions
	// *****************************************************************************

	template<class Policy, class VLC>
	void vlab_count(VLC& vlc, const typename Policy::graph_t& g)
	{
	    for (typename Policy::AllVertexIt it(g); it.valid(); it.increment())
		++vlc[Policy::vlabel(it.vertex(), g)];
	}

	template<class Policy, class Map_VL_EL_VL_P>
	void enum_one_edges(Map_VL_EL_VL_P& m, const typename Policy::graph_t& g)
	{
	    for (typename Policy::AllEdgeIt it(g); it.valid(); it.increment())
	    {
		typename Policy::Edge e(it.edge());
		typename Policy::vertex_label_ref_t vl_from = Policy::vlabel(e.vi_from, g);
		typename Policy::vertex_label_ref_t vl_to   = Policy::vlabel(e.vi_to, g);
		typename Policy::edge_label_ref_t el        = Policy::elabel(e.ei, g);
		m[vl_from][el][vl_to].push(SBG<Policy>(&g, e));    
		e.chgdir();
		m[vl_to][el][vl_from].push(SBG<Policy>(&g, e));
	    }
	}

	template<class Policy>
	typename Policy::vertex_index_t dfsc_vindex(typename Policy::vertex_index_t sbg_vi,
						    const DFSCode<Policy>& dfsc,
						    const SBG<Policy>& sbg)
	{
	    const SBG<Policy>* s = &sbg;
	    for (int i = sbg.size()-1; i >= 0; --i, s = s->prev_sbg())
	    {
		if (sbg_vi == s->edge().vi_from)
		    return dfsc[i].vi_from;
		if (sbg_vi == s->edge().vi_to)
		    return dfsc[i].vi_to;
	    }
	    return Policy::void_vindex();
	}

	template<class Policy, class P>
	bool project_is_min(const P& projected, DFSCode<Policy>& dfsc_min, const DFSCode<Policy>& dfsc_tested)
	{
            typedef typename Policy::vertex_index_t VI;
            typedef typename Policy::vertex_label_ref_t VLR;
            typedef typename Policy::edge_label_ref_t ELR;
            typedef SBG<Policy> SBG;
            typedef typename Policy::Edge Edge;
            typedef EdgeCode<Policy> EdgeCode;

            if (dfsc_min[dfsc_min.size()-1] != dfsc_tested[dfsc_min.size()-1])
                return false;
                
            // --------------------------------------------------------------
            // enumerate
            typedef typename Policy::IncidEdgeIt EdgeIter;
            typedef typename MapTraits<Policy, P>::Map_EL_P    BEdges;
            typedef typename MapTraits<Policy, P>::Map_EL_VL_P FEdges;

            const VLR VL_NULL = Policy::vlabel_null();
            RMPath rmpath(dfsc_min);
            VI vi_dfsc_rmost   = dfsc_min[rmpath.rightmost()].vi_to;
            VLR vl_rmost = dfsc_min[rmpath.rightmost()].vl_to;

            const typename Policy::graph_t& g = *projected.back().get_graph();

            // backward
            {
                BEdges b_edges;
                VI newto = Policy::void_vindex();
                bool flg = false;
                for (int i = 0; !flg && i < rmpath.size() - 1; ++i)
                {
                    const EdgeCode& ec_rmpath = dfsc_min[rmpath[i]];
                    ELR el_rmpath = ec_rmpath.el;
                    bool vl_less_eq = ec_rmpath.vl_to <= vl_rmost;

                    BOOST_FOREACH(const SBG& sbg, projected)
                    {
                        const Edge* e_rmost = sbg[rmpath.rightmost()];
                        for (EdgeIter it(sbg[rmpath[i]]->vi_from, g); it.valid(); it.increment())
                        {
                            Edge e = it.edge();
                            if (sbg.has_edge(e.ei))
                                continue;

                            ELR el = Policy::elabel(e.ei, g);
                            if (e.vi_to == e_rmost->vi_to &&
                                (el_rmpath < el || (el_rmpath == el && vl_less_eq)))
                            {
                                b_edges[el].push(SBG(&sbg, e));
                                newto = ec_rmpath.vi_from;
                                flg = true;
                                break;
                            }
                        }
                    }
                }

                if (flg)
                {
                    typename BEdges::const_iterator i1 = b_edges.begin();
                    dfsc_min.push_back(EdgeCode(vi_dfsc_rmost, newto, VL_NULL, i1->first, VL_NULL));
                    return project_is_min(i1->second, dfsc_min, dfsc_tested);
                }
            }

            // forward
            {
                FEdges f_edges;
                VLR vl_minimum = dfsc_min[0].vl_from;
                VI newfrom = Policy::void_vindex();
                bool flg = false;
            
                // forward pure
                BOOST_FOREACH(const SBG& sbg, projected)
                {
                    const Edge* e_rmost = sbg[rmpath.rightmost()];
                    for (EdgeIter it(e_rmost->vi_to, g); it.valid(); it.increment())
                    {
                        Edge e = it.edge();
                        if (sbg.has_edge(e.ei))
                            continue;
                        VLR vl = Policy::vlabel(e.vi_to, g);
                        if (!sbg.has_vertex(e.vi_to) && vl_minimum <= vl)
                        {
                            ELR el = Policy::elabel(e.ei, g);
                            f_edges[el][vl].push(SBG(&sbg, e));
                            newfrom = vi_dfsc_rmost;
                            flg = true;
                        }
                    }
                }
            
                // forward rmpath
                for (int i = rmpath.size()-1; !flg && i >= 0; --i)
                {
                    const EdgeCode& ec_rmpath = dfsc_min[rmpath[i]];
                    BOOST_FOREACH(const SBG& sbg, projected)
                    {
                        for (EdgeIter it(sbg[rmpath[i]]->vi_from, g); it.valid(); it.increment())
                        {
                            Edge e = it.edge();
                            if (sbg.has_edge(e.ei))
                                continue;

                            VLR vl = Policy::vlabel(e.vi_to, g);
                            ELR el = Policy::elabel(e.ei, g);
                            if (!sbg.has_vertex(e.vi_to) && vl_minimum <= vl &&
                                (ec_rmpath.el < el || (ec_rmpath.el == el && ec_rmpath.vl_to <= vl)))
                            {
                                f_edges[el][vl].push(SBG(&sbg, e));
                                newfrom = dfsc_min[rmpath[i]].vi_from;
                                flg = true;
                            }
                        }
                    }
                }

                if (flg)
                {
                    typename FEdges::const_iterator i1 = f_edges.begin();
                    typename FEdges::mapped_type::const_iterator i2 = i1->second.begin();
                    dfsc_min.push_back(EdgeCode(newfrom, vi_dfsc_rmost+1, VL_NULL, i1->first, i2->first));
                    return project_is_min(i2->second, dfsc_min, dfsc_tested);
                }
            }

            return true;
	}

	template<class Policy>
	bool is_min(const DFSCode<Policy>& dfsc_tested)
	{
	    typedef ProjectedSimple<Policy> P;
	    typedef typename MapTraits<Policy, P>::Map_VL_EL_VL_P M3;
	    std::auto_ptr<typename Policy::graph_t> graph(Policy::create_graph(dfsc_tested));

	    M3 root;
	    enum_one_edges<Policy>(root, *graph);

	    DFSCode<Policy> dfsc_min;
	    typename M3::const_iterator i1 = root.begin();
	    typename M3::mapped_type::const_iterator i2 = i1->second.begin();
	    typename M3::mapped_type::mapped_type::const_iterator i3 = i2->second.begin();
	    dfsc_min.push_back(EdgeCode<Policy>(0, 1, i1->first, i2->first, i3->first));
	    return project_is_min(i3->second, dfsc_min, dfsc_tested);
	}


	// *****************************************************************************
	//			class Closegraph_alg
	// *****************************************************************************
	template<class Policy, class Output, class P>
	class Closegraph_alg
	{
	    typedef typename Policy::graph_t Graph;
	    typedef typename Policy::Edge Edge;
	    typedef typename Policy::vertex_index_t VI;
	    typedef typename Policy::edge_index_t EI;
	    typedef typename Policy::vertex_label_t VL;
	    typedef typename Policy::edge_label_t EL;
	    typedef typename Policy::vertex_label_ref_t VLR;
	    typedef typename Policy::edge_label_ref_t ELR;

	    typedef typename MapTraits<Policy, P>::Map_VL_EL_VL_P RootMap;
	    typedef typename MapTraits<Policy, P>::Map_VI_EL_P    BEdges;
	    typedef typename MapTraits<Policy, P>::Map_VI_EL_VL_P FEdges;
	    typedef typename MapTraits<Policy, P>::Map_VI_VI_VL_EL_VL_P XEdges;

	    const int minsup_;
	    Output& result_;
	    DFSCode<Policy> dfsc_;
	    std::vector<bool> early_termin_;
	    int call_depth_;

	    void run(RootMap& root, std::map<VL, int>& vlc);
	    void find_xedges_vito(const XEdges& x_edges, VI vi, VLR vl, ELR el, std::set<VI>& vii) const;
	    void fail_et(const EdgeCode<Policy>& ec, bool rightmost,
			 const BEdges& b_edges, const XEdges& x_edges, std::set<int>& vidx) const;
	    void project(const P* projected, const P* prev_projected);
	public:
	    template<class TGraphIterator>
	    Closegraph_alg(TGraphIterator tg_begin,
			   TGraphIterator tg_end,
			   int minsup,
			   Output& result);
	
	    Closegraph_alg(const Graph& graph, int minsup, Output& result);
	};



	template<class Policy, class Output, class P>
	void Closegraph_alg<Policy,Output,P>::run(RootMap& root, std::map<VL, int>& vlc)
	{
	    std::map<VL, bool> et1v;

	    typedef typename RootMap::const_iterator I1;
	    typedef typename RootMap::mapped_type::const_iterator I2;
	    typedef typename RootMap::mapped_type::mapped_type::const_iterator I3;
	    for (I1 i1 = root.begin(); i1 != root.end(); ++i1)
		for (I2 i2 = i1->second.begin(); i2 != i1->second.end() && !et1v[i1->first]; ++i2)
		    for (I3 i3 = i2->second.begin(); i3 != i2->second.end() && !et1v[i1->first]; ++i3)
		    {
			const P& projected = i3->second;
			if (projected.support() >=  minsup_)
			{
			    EdgeCode<Policy> ec(0, 1, i1->first, i2->first, i3->first);
#ifdef DEBUG_PRINT
			    std::cerr << "FIRST EC:" << ec << " support=" << projected.size() << " "
				      << "vlc[ec.vl_from]=" << vlc[ec.vl_from] << " "
				      << "vlc[ec.vl_to]=" << vlc[ec.vl_to] << std::endl;
#endif
			    if (projected.size() == vlc[ec.vl_from])
				et1v[ec.vl_from] = true;
			    if (projected.size() == vlc[ec.vl_to])
				et1v[ec.vl_to] = true;

			    dfsc_.push_back(ec);
			    project(&projected, 0);
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
	    
	    const VI VI_VOID = Policy::void_vindex();

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

	    std::map<VI, std::set<int> > vii_2;
	    for (typename std::set<VI>::const_iterator i = vii_1.begin(); i != vii_1.end(); ++i)
		for (unsigned int idx = 0; idx != dfsc_.size(); ++idx)
		{
		    if (dfsc_[idx].vi_from == *i)
			vii_2[dfsc_[idx].vi_to].insert(idx);
		    if (dfsc_[idx].vi_to == *i)
			vii_2[dfsc_[idx].vi_from].insert(idx);
		}

	    for (typename std::map<VI, std::set<int> >::const_iterator i = vii_2.begin(); i != vii_2.end(); ++i)
		for (unsigned int idx = 0; idx != dfsc_.size(); ++idx)
		{
		    if (dfsc_[idx].vi_from == i->first && i->second.count(idx) == 0)
			vidx.insert(idx);
		    if (dfsc_[idx].vi_to == i->first && i->second.count(idx) == 0)
			vidx.insert(idx);
		}
	}

	template<class Policy, class Output, class P>
	void Closegraph_alg<Policy,Output,P>::project(const P* projected, const P* prev_projected)
	{
	    ++call_depth_;

//#ifdef DEBUG_PRINT
	    if (call_depth_ < 9)
	    {
		std::cerr << call_depth_ << ":";
		for (int i = 0; i < call_depth_; ++i)
		    std::cerr << " ";
		std::cerr << dfsc_.back() << std::endl;
	    }
//#endif

	    early_termin_.push_back(false);

	    typedef SBG<Policy> SBG;
	    typedef EdgeCode<Policy> EdgeCode;
	    // --------------------------------------------------------------
	    //   discover all edge extensions
	    // --------------------------------------------------------------
	    typedef typename Policy::IncidEdgeIt EdgeIter;
	    BEdges b_edges;
	    FEdges f_edges;
	    XEdges x_edges;

	    const int NUM_EDGES = dfsc_.size(); 
	    RMPath rmpath(dfsc_);

	    VI vi_dfsc_rmost   = dfsc_[rmpath.rightmost()].vi_to;
	    VLR vl_rmost = dfsc_[rmpath.rightmost()].vl_to;
	    VLR vl_minimum = dfsc_[0].vl_from;

	    BOOST_FOREACH(const SBG& sbg, *projected)
	    {
		const Graph& g = *sbg.get_graph();

		std::vector<bool> r_extension(Policy::nedges(g), false);
		std::vector<bool> x_extension(Policy::nedges(g), false);

		// -----------------------------------
		// discover RMPath extension edges
		// -----------------------------------	
		const Edge* e_rmost = sbg[rmpath.rightmost()];

		// backward
		for (int i = 0; i < rmpath.size() - 1; ++i)
		{
		    const EdgeCode& ec_rmpath = dfsc_[rmpath[i]];
		    ELR el_rmpath = ec_rmpath.el;
		    bool vl_less_eq = ec_rmpath.vl_to <= vl_rmost;

		    for (EdgeIter it(sbg[rmpath[i]]->vi_from, g); it.valid(); it.increment())
		    {
			Edge e = it.edge();
			if (sbg.has_edge(e.ei))
			    continue;

			ELR el = Policy::elabel(e.ei, g);
			if (e.vi_to == e_rmost->vi_to &&
			    (el_rmpath < el || (el_rmpath == el && vl_less_eq)))
			{
			    // discover BACKWARD edge with DFSC order (B edge)
			    // from right most vertex (e_rmost->vi_to)
			    // to each vertices on right most path (e_sbg->vi_from)
			    assert(!r_extension[e.ei]);
			    b_edges[ec_rmpath.vi_from][el].push(SBG(&sbg, e));
			    r_extension[e.ei] = true;
			    break;
			}
		    }
		}

		// forward pure
		for (EdgeIter it(e_rmost->vi_to, g); it.valid(); it.increment())
		{
		    Edge e = it.edge();
		    if (sbg.has_edge(e.ei))
			continue;

		    VLR vl = Policy::vlabel(e.vi_to, g);
		    if (!sbg.has_vertex(e.vi_to) && vl_minimum <= vl)
		    {
			// discover FORWARD edge with DFSC order (F edge)
			// from right most vertex (e_rmmost->vi_to)
			assert(!r_extension[e.ei]);
			f_edges[vi_dfsc_rmost][Policy::elabel(e.ei, g)][vl].push(SBG(&sbg, e));
			r_extension[e.ei] = true;
		    }
		}

		
		// forward rmpath
		for (int i = rmpath.size()-1; i >= 0; --i)
		{
		    const EdgeCode& ec_rmpath = dfsc_[rmpath[i]];
		    ELR el_rmpath = ec_rmpath.el;

		    for (EdgeIter it(sbg[rmpath[i]]->vi_from, g); it.valid(); it.increment())
		    {
			Edge e = it.edge();
			if (sbg.has_edge(e.ei))
			    continue;

			VLR vl = Policy::vlabel(e.vi_to, g);
			ELR el = Policy::elabel(e.ei, g);
			if (!sbg.has_vertex(e.vi_to) && vl_minimum <= vl &&
			    (el_rmpath < el || (el_rmpath == el && ec_rmpath.vl_to <= vl)))
			{
			    // discover FORWARD edge with DFSC order (F edge)
			    // from right most path (e_rmpath->vi_from)
			    assert(!r_extension[e.ei]);
			    f_edges[ec_rmpath.vi_from][el][vl].push(SBG(&sbg, e));
			    r_extension[e.ei] = true;
			}
		    }
		}

		// -----------------------------------
		// discover NOT RMPath extension edges
		// -----------------------------------
		std::vector<bool> vv(Policy::nvertices(g), false);
		
		for (int i = 0; i < NUM_EDGES; ++i)
		{
		    VI vi = sbg[i]->vi_from;
		    if (!vv[vi])
		    {
			vv[vi] = true;
			for (EdgeIter it(vi, g); it.valid(); it.increment())
			{
			    Edge e = it.edge();
			    if (sbg.has_edge(e.ei))
				continue;
			
			    if (!r_extension[e.ei] && !x_extension[e.ei])
			    {
				VI vi_from = dfsc_vindex<Policy>(e.vi_from, dfsc_, sbg);
				VI vi_to = dfsc_vindex<Policy>(e.vi_to, dfsc_, sbg);
				VLR vl_from = Policy::vlabel(e.vi_from, g);
				VLR vl_to = Policy::vlabel(e.vi_to, g);
				ELR el = Policy::elabel(e.ei, g);

				x_edges[vi_from][vi_to][vl_from][el][vl_to].push(SBG(&sbg, e));
				x_extension[e.ei] = true;
			    }
			}
		    }

		    vi = sbg[i]->vi_to;
		    if (!vv[vi])
		    {
			vv[vi] = true;
			for (EdgeIter it(vi, g); it.valid(); it.increment())
			{
			    Edge e = it.edge();
			    if (sbg.has_edge(e.ei))
				continue;
			
			    if (!r_extension[e.ei] && !x_extension[e.ei])
			    {
				VI vi_from = dfsc_vindex<Policy>(e.vi_from, dfsc_, sbg);
				VI vi_to = dfsc_vindex<Policy>(e.vi_to, dfsc_, sbg);
				VLR vl_from = Policy::vlabel(e.vi_from, g);
				VLR vl_to = Policy::vlabel(e.vi_to, g);
				ELR el = Policy::elabel(e.ei, g);

				x_edges[vi_from][vi_to][vl_from][el][vl_to].push(SBG(&sbg, e));
				x_extension[e.ei] = true;
			    }
			}
		    }		    
		}

	    } // end: for each sbg


	    // --------------------------------------------------------------
	    //		recursive process SBG children
	    // --------------------------------------------------------------
	    bool r_closed = true;
	    typedef typename BEdges::iterator BckI1;
	    typedef typename BEdges::mapped_type::iterator BckI2;
	    for (BckI1 i1 = b_edges.begin(); i1 != b_edges.end() && !early_termin_.back(); ++i1)
		for (BckI2 i2 = i1->second.begin(); i2 != i1->second.end() && !early_termin_.back(); ++i2)
		{
		    EdgeCode ec(vi_dfsc_rmost, i1->first, Policy::vlabel_null(), i2->first, Policy::vlabel_null());
		    dfsc_.push_back(ec);
		    P& new_prj = i2->second;
		    
		    if (is_min(dfsc_))
		    {
			if (new_prj.support() >= minsup_)
			{
			    if (new_prj.support() == projected->support())
				r_closed = false; // projected is not closed
			    
			    // detect early termination
			    if (new_prj.size() == projected->size())
				early_termin_.back() = true;

			    project(&new_prj, projected);
			}
		    }
		    else
			x_edges[ec.vi_from][ec.vi_to][ec.vl_from][ec.el][ec.vl_to].merge(new_prj);

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
			EdgeCode ec(i1->first, vi_dfsc_rmost+1, Policy::vlabel_null(), i2->first, i3->first);
			dfsc_.push_back(ec);
			P& new_prj = i3->second;

			if (is_min(dfsc_))
			{
			    int new_prj_supp = new_prj.support();
			    if (new_prj_supp >= minsup_)
			    {
				if (new_prj_supp == projected->support())
				    r_closed = false; // projected is not closed

				// detect early termination
				if (new_prj.size() == projected->size())
				    early_termin_.back() = true;
				
				project(&new_prj, projected);
			    }
			    else
			    {

				std::set<int> idxs;
				fail_et(ec, ec.vi_from == vi_dfsc_rmost, b_edges, x_edges, idxs);

				for (std::set<int>::const_iterator ii = idxs.begin(); ii != idxs.end(); ++ii)
				    early_termin_[*ii] = false;
			    }
			}
			else
			    x_edges[ec.vi_from][ec.vi_to][ec.vl_from][ec.el][ec.vl_to].merge(new_prj);
			
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

	    --call_depth_;
	}

	template<class Policy, class Output, class P>
	template<class TGraphIterator>
	Closegraph_alg<Policy,Output,P>::Closegraph_alg(TGraphIterator tg_begin,
							TGraphIterator tg_end,
							int minsup,
							Output& result)
	    :minsup_(minsup), result_(result), call_depth_(0)
	{
	    std::map<VL, int> vlc;
	    RootMap root;
	    
	    for (; tg_begin != tg_end; ++tg_begin)
	    {
		vlab_count<Policy>(vlc, *tg_begin);
		enum_one_edges<Policy>(root, *tg_begin);
	    }
	    run(root, vlc);
	}

	template<class Policy, class Output, class P>
	Closegraph_alg<Policy,Output,P>::Closegraph_alg(const Graph& graph,
							int minsup,
							Output& result)
	    :minsup_(minsup), result_(result), call_depth_(0)
	{
	    std::map<VL, int> vlc;
	    RootMap root;

	    vlab_count<Policy>(vlc, graph);
	    enum_one_edges<Policy>(root, graph);

	    run(root, vlc);
	}

    } // end: namespace detail


    template<class Policy, class Graph, class Output>
    void closegraph(const Graph& graph, int minsup, Output& result)
    {
	detail::Closegraph_alg<Policy, Output, detail::Projected<Policy, SubgraphsOfOneGraph> >
	    (graph, minsup, result);
    }


    template<class Policy, class TGraphIterator, class Output>
    void closegraph(TGraphIterator tg_begin,
		    TGraphIterator tg_end,
		    int minsup,
		    Output& result)
    {
	detail::Closegraph_alg<Policy, Output, detail::Projected<Policy, SubgraphsOfManyGraph> >
	    (tg_begin, tg_end, minsup, result);
    }

} // end: namespace gSpan

#endif
