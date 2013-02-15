#ifndef CLOSEGRAPH_MT_H_
#define CLOSEGRAPH_MT_H_

#include "gspan.hpp"
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/noncopyable.hpp>
#include <ctime>

namespace gSpan
{
    namespace mt		// multi threaded
    {
	const int MAX_THREAD_COUNT = 0;

	template<class T>
	struct ThreadTime
	{
	    ThreadTime()
		{
		    mtx_.lock();
		    std::cerr << "Thread " << boost::this_thread::get_id() << " started" << std::endl;
		    mtx_.unlock();
		    t_start_ = std::time(0);
		}
	    ~ThreadTime()
		{
		    unsigned long d = std::time(0) - t_start_;
		    mtx_.lock();
		    std::cerr << "Thread " << boost::this_thread::get_id() << "    exit"
			      << " " << d << " secs" << std::endl;
		    mtx_.unlock();
		}
	    std::time_t t_start_;
	    static boost::mutex mtx_;
	};

	template<class T>
	boost::mutex ThreadTime<T>::mtx_;

	// *****************************************************************************
	//			class CommonData
	// *****************************************************************************
	template<class Policy, class Output>
	class CommonData : private boost::noncopyable
	{
	    int minsup_;

	    Output& output_;
	    boost::mutex mtx_output_;
	    
	    const int max_thread_count_;
	    int thread_count_;
	    boost::mutex mtx_thread_count_;
	public:
	    const int avg_edges;
	    CommonData(int minsup, Output& output, int avg_edges,
		       int max_thread_count = MAX_THREAD_COUNT);
	    int minsup() const { return minsup_; }
	    template<template<class> class S>
	    void output(const DFSCode<Policy>& dfsc, const S<Policy>& sg);

	    template<typename F, class A>
	    boost::thread* create_thread(F threadfunc, A arg);
	    void on_exit_thread();
	};

	template<class Policy, class Output>
	CommonData<Policy,Output>::CommonData(int minsup, Output& output, int avg_edges,
					      int max_thread_count)
	    :minsup_(minsup),
	     output_(output),
	     max_thread_count_(max_thread_count),
	     thread_count_(0),
	     avg_edges(avg_edges)
	{
	}
	

	template<class Policy, class Output>
	template<typename F, class A>
	boost::thread* CommonData<Policy,Output>::create_thread(F threadfunc, A arg)
	{
	    boost::thread* thrd = 0;
	    mtx_thread_count_.lock();
	    assert(thread_count_ >= 0);
	    if (thread_count_ < max_thread_count_)
	    {
		thrd = new boost::thread(threadfunc, arg);
		++thread_count_;
	    }
	    mtx_thread_count_.unlock();
	    return thrd;
	}


	template<class Policy, class Output>
	void CommonData<Policy,Output>::on_exit_thread()
	{
	    mtx_thread_count_.lock();
	    assert(thread_count_ >= 1);
	    --thread_count_;
	    mtx_thread_count_.unlock();
	}
       
	
	template<class Policy, class Output>
	template<template<class> class S>
	void CommonData<Policy,Output>::output(const DFSCode<Policy>& dfsc, const S<Policy>& sg)
	{
	    mtx_output_.lock();
	    output_(dfsc, sg);
	    mtx_output_.unlock();
	}

	// *****************************************************************************
	//			class OnExitThread
	// *****************************************************************************
	template<class CmnData>
	class OnExitThread
	{
	    CmnData* dt_;
	public:
	    explicit OnExitThread(CmnData* dt) :dt_(dt) {}
	    ~OnExitThread() { dt_->on_exit_thread(); }
	};


	// *****************************************************************************
	//			class ET List
	//                early termination markers
	// *****************************************************************************
	class ET : private boost::noncopyable
	{
	    ET* prev_;
	    volatile bool et_;
	    boost::mutex mtx_;
	public:
	    ET() :prev_(0), et_(false) {}
	    explicit ET(ET* prev) :prev_(prev), et_(false) {}
	    volatile bool test() const	{ return et_; }
	    void set()			{ mtx_.lock(); et_ = true; mtx_.unlock(); }
	    void set_fast()		{ et_ = true; }
	};

	// *****************************************************************************
	//			class ThreadData
	// *****************************************************************************
	template<class Policy, class P, class CmnData>
	class ThreadData : private boost::noncopyable
	{
	    typedef typename Policy::graph_t Graph;
	    typedef typename Policy::Edge Edge;
	    typedef typename Policy::vertex_index_t VI;
	    typedef typename Policy::edge_index_t EI;
	    typedef typename Policy::vertex_label_t VL;
	    typedef typename Policy::edge_label_t EL;
	    typedef typename Policy::vertex_label_ref_t VLR;
	    typedef typename Policy::edge_label_ref_t ELR;
	    typedef DFSCode<Policy> DFSC;

	    typedef typename detail::MapTraits<Policy, P>::Map_VL_EL_VL_P RootMap;
	    typedef typename detail::MapTraits<Policy, P>::Map_VI_EL_P    BEdges;
	    typedef typename detail::MapTraits<Policy, P>::Map_VI_EL_VL_P FEdges;
	    typedef typename detail::MapTraits<Policy, P>::Map_VI_VI_VL_EL_VL_P XEdges;

	    CmnData* d_;
	    DFSC* dfsc_;
	    bool need_thread() const;
	public:
	    struct thread_param
	    {
		DFSC* dfsc;
		const P* projected;
		const P* prev_projected;
		ET* et;
		ThreadData* thrd_data_parent;
	    };
	    static void thread_func(thread_param* param);

	    explicit ThreadData(CmnData* d) :d_(d), dfsc_(new DFSC) {}
	    ThreadData(thread_param* param) :d_(param->thrd_data_parent->d_), dfsc_(param->dfsc) {}
	    ~ThreadData() { delete dfsc_; }

	    void run(RootMap& root, std::map<VL, int>& vlc);
	    void project(const P* projected, const P* prev_projected, ET* et);
	    void process(const P* projected, const P* prev_projected, ET* et, boost::thread_group* thrd_grp);
	};


	template<class Policy, class P, class CmnData>
	bool ThreadData<Policy,P,CmnData>::need_thread() const
	{
	    float num_edges = dfsc_->size();
	    //float rate = num_edges / d_->avg_edges;
	    //std::cerr << "num_edges=" << num_edges << " avg_edges=" << d_->avg_edges << " rate="<< rate << std::endl;
	    return num_edges < 10;
	}


	template<class Policy, class P, class CmnData>
	void ThreadData<Policy,P,CmnData>::thread_func(thread_param* param)
	{
	    ThreadTime<int> thrtime;
	    OnExitThread<CmnData> onexit(param->thrd_data_parent->d_);

	    ThreadData thrd_data(param);
	    const P* projected = param->projected;
	    const P* prev_projected = param->prev_projected;
	    ET* et = param->et;
	    delete param;
	    thrd_data.project(projected, prev_projected, et);
	}

	
	template<class Policy, class P, class CmnData>
	void ThreadData<Policy,P,CmnData>::run(RootMap& root, std::map<VL, int>& vlc)
	{
	    std::map<VL, bool> et1v;
	    const int minsup = d_->minsup();

	    boost::thread_group thrd_grp;
	    //int i = 0;
	    typedef typename RootMap::const_iterator I1;
	    typedef typename RootMap::mapped_type::const_iterator I2;
	    typedef typename RootMap::mapped_type::mapped_type::const_iterator I3;
	    for (I1 i1 = root.begin(); i1 != root.end(); ++i1)
	    {
		for (I2 i2 = i1->second.begin(); i2 != i1->second.end() && !et1v[i1->first]; ++i2)
		{
		    for (I3 i3 = i2->second.begin(); i3 != i2->second.end() && !et1v[i1->first]; ++i3)
		    {
			const P* projected = &i3->second;
			if (projected->support() >= minsup)
			{
			    EdgeCode<Policy> ec(0, 1, i1->first, i2->first, i3->first);
			    if (projected->size() == vlc[ec.vl_from])
				et1v[ec.vl_from] = true;
			    if (projected->size() == vlc[ec.vl_to])
				et1v[ec.vl_to] = true;
			    
			    dfsc_->push_back(ec);
			    process(projected, 0, 0, &thrd_grp);
			    dfsc_->pop_back();

			    //if (++i == 2)
			    //goto b;
			}
		    }
		}
	    }
	    //b:
	    thrd_grp.join_all();
	}

	template<class Policy, class P, class CmnData>
	void ThreadData<Policy,P,CmnData>::project(const P* projected, const P* prev_projected, ET* et_prev)
	{
	    using namespace detail;

	    typedef SBG<Policy> SBG;
	    typedef EdgeCode<Policy> EdgeCode;

	    ET et(et_prev);

	    // --------------------------------------------------------------
	    //   discover all edge extensions
	    // --------------------------------------------------------------
	    typedef typename Policy::IncidEdgeIt EdgeIter;
	    BEdges b_edges;
	    FEdges f_edges;
	    XEdges x_edges;
	    
	    DFSC& dfsc = *dfsc_;

	    const int NUM_EDGES = dfsc.size(); 
	    detail::RMPath rmpath(dfsc);

	    VI vi_dfsc_rmost   = dfsc[rmpath.rightmost()].vi_to;
	    VLR vl_rmost = dfsc[rmpath.rightmost()].vl_to;
	    VLR vl_minimum = dfsc[0].vl_from;

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
		    const EdgeCode& ec_rmpath = dfsc[rmpath[i]];
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
		    const EdgeCode& ec_rmpath = dfsc[rmpath[i]];
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
				VI vi_from = dfsc_vindex<Policy>(e.vi_from, dfsc, sbg);
				VI vi_to = dfsc_vindex<Policy>(e.vi_to, dfsc, sbg);
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
				VI vi_from = dfsc_vindex<Policy>(e.vi_from, dfsc, sbg);
				VI vi_to = dfsc_vindex<Policy>(e.vi_to, dfsc, sbg);
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
	    //		recursive process SBG extentions
	    // --------------------------------------------------------------
	    const int minsup = d_->minsup();
	    bool r_closed = true;
	    boost::thread_group thrd_grp;
	    // backward
	    typedef typename BEdges::iterator BckI1;
	    typedef typename BEdges::mapped_type::iterator BckI2;
	    for (BckI1 i1 = b_edges.begin(); i1 != b_edges.end() && !et.test(); ++i1)
		for (BckI2 i2 = i1->second.begin(); i2 != i1->second.end() && !et.test(); ++i2)
		{
		    EdgeCode ec(vi_dfsc_rmost, i1->first, Policy::vlabel_null(), i2->first, Policy::vlabel_null());
		    dfsc_->push_back(ec);
		    P& new_prj = i2->second;
		    
		    if (is_min(*dfsc_))
		    {
			if (new_prj.support() >= minsup)
			{
			    if (new_prj.support() == projected->support())
				r_closed = false; // projected is not closed
			    
			    // detect early termination
			    if (new_prj.size() == projected->size())
				et.set_fast();
			    
			    process(&new_prj, projected, &et, &thrd_grp);
			}
		    }
		    else
			x_edges[ec.vi_from][ec.vi_to][ec.vl_from][ec.el][ec.vl_to].merge(new_prj);

		    dfsc_->pop_back();
		}
	    thrd_grp.join_all();


	    // forward
	    typedef typename FEdges::reverse_iterator FwdRevI1;
	    typedef typename FEdges::mapped_type::iterator FwdI2;
	    typedef typename FEdges::mapped_type::mapped_type::iterator FwdI3;
	    for (FwdRevI1 i1 = f_edges.rbegin(); i1 != f_edges.rend() && !et.test(); ++i1)
		for (FwdI2 i2 = i1->second.begin(); i2 != i1->second.end() && !et.test(); ++i2)
		    for (FwdI3 i3 = i2->second.begin(); i3 != i2->second.end() && !et.test(); ++i3)
		    {
			EdgeCode ec(i1->first, vi_dfsc_rmost+1, Policy::vlabel_null(), i2->first, i3->first);
			dfsc_->push_back(ec);
			P& new_prj = i3->second;

			if (is_min(*dfsc_))
			{
			    int new_prj_supp = new_prj.support();
			    if (new_prj_supp >= minsup)
			    {
				if (new_prj_supp == projected->support())
				    r_closed = false; // projected is not closed

				// detect early termination
				if (new_prj.size() == projected->size())
				    et.set_fast();
				
				process(&new_prj, projected, &et, &thrd_grp);
			    }
			    else
			    {
/*
				std::set<int> idxs;
				fail_et(ec, ec.vi_from == vi_dfsc_rmost, b_edges, x_edges, idxs);

				for (std::set<int>::const_iterator ii = idxs.begin(); ii != idxs.end(); ++ii)
				    early_termin_[*ii] = false;
*/
			    }
			}
			else
			    x_edges[ec.vi_from][ec.vi_to][ec.vl_from][ec.el][ec.vl_to].merge(new_prj);
			
			dfsc_->pop_back();
		    }
	    thrd_grp.join_all();
	}


	template<class Policy, class P, class CmnData>
	void ThreadData<Policy,P,CmnData>::process(const P* projected, const P* prev_projected, ET* et,
						   boost::thread_group* thrd_grp)
	{
	    boost::thread* thrd = 0;
	    if (need_thread())
	    {
		thread_param* param = new thread_param;
		param->dfsc = new DFSC(*dfsc_); // make copy for thread
		param->projected = projected;
		param->prev_projected = prev_projected;
		param->et = et;
		param->thrd_data_parent = this;
		if ( (thrd = d_->create_thread(&thread_func, param)))
		    thrd_grp->add_thread(thrd);
		else
		{
		    delete param->dfsc;
		    delete param;
		}
	    }
	    
	    if (!thrd)
	    {
		project(projected, prev_projected, et);
	    }
	}



    } // end: namespace mt


    struct NumEdges
    {
	int count;
	NumEdges() :count(0) {}
	void operator() () { ++count; }
    };


    template<class Policy, class Graph, class Output>
    void closegraph_mt(const Graph& graph, int minsup, Output& result)
    {
	typedef detail::Projected<Policy, SubgraphsOfOneGraph> P;
	typedef typename detail::MapTraits<Policy, P>::Map_VL_EL_VL_P RootMap;
	typedef mt::CommonData<Policy,Output> CommonData;
	typedef mt::ThreadData<Policy,P,CommonData> ThreadData;

	std::map<typename Policy::vertex_label_t, int> vlc;
	RootMap root;

	detail::vlab_count<Policy>(vlc, graph);

	NumEdges n_edges;
	detail::enum_one_edges<Policy>(root, graph, n_edges);
	
	CommonData dt(minsup, result, n_edges.count);
	ThreadData thrd_data(&dt);
	thrd_data.run(root, vlc);
    }

    template<class Policy, class TGraphIterator, class Output>
    void closegraph_mt(TGraphIterator tg_begin,
		       TGraphIterator tg_end,
		       int minsup,
		       Output& result)
    {
    }
}

#endif
