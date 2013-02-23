#ifndef CLOSEGRAPH_ST_H_
#define CLOSEGRAPH_ST_H_

#include "gspan.hpp"

extern unsigned long proj_calls;

namespace gSpan
{
    namespace st		// single threaded
    {
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
	   
	    typedef typename detail::MapTraits<Policy, P>::Map_VL_EL_VL_P RootMap;
	    typedef typename detail::MapTraits<Policy, P>::Map_VI_EL_P    BEdges;
	    typedef typename detail::MapTraits<Policy, P>::Map_VI_EL_VL_P FEdges;
	    typedef typename detail::MapTraits<Policy, P>::Map_VI_VI_VL_EL_VL_P XEdges;

	    const int minsup_;
	    Output& result_;
	    DFSCode<Policy> dfsc_;
	    std::vector<bool> early_termin_;
	    int call_depth_;

#ifdef GSPAN_WITH_STATISTICS
	    unsigned long num_project_calls_;
	    float num_ismin_calls_;
	    float num_ismin_true_ret_;
	    void statistics_init()
		{
		    num_project_calls_ = 0;
		    num_ismin_calls_ = num_ismin_true_ret_ = 0.0f;
		}

	    bool is_min(const DFSCode<Policy>& dfsc)
		{
		    bool ret = detail::is_min(dfsc);
		    num_ismin_calls_ += 1;
		    if (ret)
			num_ismin_true_ret_ += 1;
		    return ret;
		}

	    void statistics_report()
		{
		    float ismin_true_ret_prc = num_ismin_true_ret_;
		    ismin_true_ret_prc /= num_ismin_calls_;
		    std::cerr << "------- statistics -----------------------";
		    std::cerr << "\nnum_project_calls=" << num_project_calls_
			      << "\nnum_ismin_calls=" << num_ismin_calls_
			      << "\nismin_true_ret%=" << ismin_true_ret_prc*100.0f
			      << std::endl;
		    std::cerr << "------------------------------------------" << std::endl;
		}
#endif

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

	    ~Closegraph_alg()
		{
#ifdef GSPAN_WITH_STATISTICS
		    statistics_report();
#endif
		}
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
			    std::cerr << "FIRST EC:" << ec << " support=" << projected.support() << " "
				      << " size=" << projected.size() << " "
				      << "vlc[ec.vl_from]=" << vlc[ec.vl_from] << " "
				      << "vlc[ec.vl_to]=" << vlc[ec.vl_to] << std::endl;
#endif
			    if (projected.size() == vlc[ec.vl_from])
				et1v[ec.vl_from] = true;

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
	    using namespace detail;
	    ++call_depth_;

#ifdef GSPAN_WITH_STATISTICS
	    ++num_project_calls_;
#endif

#ifdef DEBUG_PRINT
	    static int ncall;
	    ++ncall;

	    std::cerr << ">--- project: NCALL=" << ncall << " depth=" << call_depth_ << "   "<< dfsc_ << std::endl
		      << "     projected.support=" << projected->support() << " size=" << projected->size()
		      << std::endl << *projected << std::endl;;
/*
	    if (call_depth_ < 3)
	    {
		std::cerr << call_depth_ << ":";
		for (int i = 0; i < call_depth_; ++i)
		    std::cerr << " ";
		std::cerr << dfsc_.back() << std::endl;
	    }
*/
#endif

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
	    detail::RMPath rmpath(dfsc_);

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
			    if (new_prj.size() == projected->size() && new_prj.support() == projected->support())
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
				if (new_prj.size() == projected->size() && new_prj_supp == projected->support())
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

#ifdef DEBUG_PRINT
	    std::cerr << "RET to " << call_depth_ << std::endl;;
#endif
	}

	template<class Policy, class Output, class P>
	template<class TGraphIterator>
	Closegraph_alg<Policy,Output,P>::Closegraph_alg(TGraphIterator tg_begin,
							TGraphIterator tg_end,
							int minsup,
							Output& result)
	    :minsup_(minsup), result_(result), call_depth_(0)
	{
#ifdef GSPAN_WITH_STATISTICS
		    statistics_init();
#endif
	    using namespace detail;
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
#ifdef GSPAN_WITH_STATISTICS
		    statistics_init();
#endif
	    using namespace detail;
	    std::map<VL, int> vlc;
	    RootMap root;

	    vlab_count<Policy>(vlc, graph);
	    enum_one_edges<Policy>(root, graph);

	    run(root, vlc);
	}

    } // end: namespace st


    template<class Policy, class Graph, class Output>
    void closegraph_st(const Graph& graph, int minsup, Output& result)
    {
	st::Closegraph_alg<Policy, Output, detail::Projected<Policy, SubgraphsOfOneGraph> >
	    (graph, minsup, result);
    }


    template<class Policy, class TGraphIterator, class Output>
    void closegraph_st(TGraphIterator tg_begin,
		       TGraphIterator tg_end,
		       int minsup,
		       Output& result)
    {
	st::Closegraph_alg<Policy, Output, detail::Projected<Policy, SubgraphsOfManyGraph> >
	    (tg_begin, tg_end, minsup, result);
    }

}

#endif
