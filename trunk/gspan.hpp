#ifndef GSPAN_H_
#define GSPAN_H_

#include <vector>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <limits>
#include <utility>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include <boost/foreach.hpp>

#define BR asm("int $3;")

extern std::map<int, char> mlv;
extern std::map<int, char> mle;

namespace graph_alg
{
    typedef unsigned int Vertex; // is the same as vertex_descriptor
    typedef unsigned int EdgeIndex;

    typedef int VertexLabel;		// vertex label
    typedef int EdgeLabel;		// edge label

    inline Vertex null_vertex() { return std::numeric_limits<Vertex>::max(); }

    // ****************************************************************************
    //				Edge
    // ****************************************************************************
    struct Edge
    {
	Vertex from, to;
	EdgeIndex eidx;
	Edge(Vertex v = 0, Vertex u = 0, EdgeIndex ei = 0) :from(v), to(u), eidx(ei) {}
    };

    template<class G>
    inline EdgeLabel fromlabel(const Edge& e, const G& g) { return g[e.from]; }

    template<class G>
    inline EdgeLabel tolabel(const Edge& e, const G& g) { return g[e.to]; }

    template<class G>
    inline EdgeLabel edgelabel(const Edge& e, const G& g) { return g[edge_from_index(e.eidx,g)]; }


    std::ostream& operator<<(std::ostream& out, const Edge& e)
    { return out << "Edge#" << e.eidx << ":" << "("<<e.from<<", "<<e.to<<")"; }
    
    template<class G>
    std::ostream& print_edge(std::ostream& out, const Edge& e, const G& g)
    { return out<<e<<" ["<<g[e.from]<<", "<<g[edge_from_index(e.eidx,g)]<<", "<<g[e.to]<< "]"; }

    // ****************************************************************************
    //				EdgeCode
    // ****************************************************************************
    struct EdgeCode
    {
	Vertex from, to;
	VertexLabel fromlabel, tolabel;
	EdgeLabel elabel;
	EdgeCode(Vertex from, Vertex to, VertexLabel fromlabel, EdgeLabel elabel, VertexLabel tolabel)
	    :from(from), to(to), fromlabel(fromlabel), tolabel(tolabel), elabel(elabel) {}
	bool is_forward() const { return from < to; }
	operator std::pair<Vertex,Vertex> () const { return std::pair<Vertex,Vertex>(from,to); }
	bool operator== (const EdgeCode& r) const
	    { return
		    from == r.from && to == r.to &&
		    fromlabel == r.fromlabel && tolabel == r.tolabel &&
		    elabel == r.elabel; }
	bool operator!= (const EdgeCode& r) const { return ! (*this == r); }	
    };

    inline std::ostream& operator<<(std::ostream& out, const EdgeCode& ec)
    {
	return out<<"("<<ec.from<<","<<ec.to<<", "<<ec.fromlabel
		  <<","<<ec.elabel<<","<<ec.tolabel<<")";
	//return out<<"("<<ec.from<<","<<ec.to<<", "<<mlv[ec.fromlabel]
	//<<","<<mle[ec.elabel]<<","<<mlv[ec.tolabel]<<")";
    }

    // ****************************************************************************
    //				DFSCode
    // ****************************************************************************
    class DFSCode
    {
	std::vector<EdgeCode> dfsc_;
    public:
	const EdgeCode& operator[] (int i) const { return dfsc_[i]; }
	bool empty() const { return dfsc_.empty(); }
	unsigned int size() const { return dfsc_.size(); }	
	typedef std::vector<EdgeCode>::const_iterator const_iterator;
	const_iterator begin() const { return dfsc_.begin(); }
	const_iterator end() const { return dfsc_.end(); }
	Vertex max_vertex() const
	    {
		Vertex m = 0;
		for (const_iterator i = begin(); i != end(); ++i)
		    m = std::max(m, std::max(i->from, i->to));
		return m;
	    }
	void push(const EdgeCode& ec) { dfsc_.push_back(ec); }
	void pop() { dfsc_.pop_back(); }
    };

    inline std::ostream& operator<<(std::ostream& out, const DFSCode& dfsc)
    {
	std::copy(dfsc.begin(), dfsc.end(), std::ostream_iterator<EdgeCode>(out, " "));
	return out;
    }

    // ****************************************************************************
    //				RMPath
    // ****************************************************************************
    class RMPath
    {
	std::deque<unsigned int> rmp_;
    public:
	RMPath(const DFSCode& dfsc)
	    {
		Vertex old_from = 0; // non initialized
		for (int i = dfsc.size()-1; i >= 0; --i)
		    if (dfsc[i].is_forward() && (rmp_.empty() || old_from == dfsc[i].to))
		    {
			rmp_.push_front(i);
			old_from = dfsc[i].from;
		    }
	    }

	unsigned int operator[] (int i) const { return rmp_[i]; }
	unsigned int size() const { return rmp_.size(); }
	unsigned int rightmost() const { return rmp_.back(); }
	friend std::ostream& operator<<(std::ostream& out, const RMPath& rpm);
    };

    inline std::ostream& operator<<(std::ostream& out, const RMPath& r)
    {
	std::copy(r.rmp_.begin(), r.rmp_.end(), std::ostream_iterator<unsigned int>(out, " "));
	return out;
    }


    // ****************************************************************************
    //				WGraph
    // ****************************************************************************
    template<class Directed>
    class WGraph
	: public boost::compressed_sparse_row_graph<Directed,
						    VertexLabel, EdgeLabel, boost::no_property,
						    Vertex, EdgeIndex>
    {
	typedef boost::compressed_sparse_row_graph<Directed,
						   VertexLabel, EdgeLabel, boost::no_property,
						   Vertex, EdgeIndex> Base;
    public:
	WGraph() {}
	WGraph(const DFSCode& dfsc);
    };

    template<class Directed>
    WGraph<Directed>::WGraph(const DFSCode& dfsc)
	:Base(boost::edges_are_unsorted_multi_pass_t(),
	      dfsc.begin(), dfsc.end(), dfsc.max_vertex() + 1)
    {
	// set labels
	DFSCode::const_iterator it = dfsc.begin();
	DFSCode::const_iterator it_end = dfsc.end();
	if (it != it_end)
	{
	    if (it->fromlabel != -1)
		(*this)[it->from] = it->fromlabel;
	    while (it != it_end)
	    {
		(*this)[boost::edge(it->from, it->to, *this).first] = it->elabel;
		if (it->is_forward())
		{
		    assert(it->tolabel != -1);
		    (*this)[it->to] = it->tolabel;
		}
		++it;
	    }
	}
    }


    // ****************************************************************************
    //				SubgraphOfTheGraph
    // ****************************************************************************
    template<class G>
    struct SubgraphOfTheGraph
    {
	std::vector<Edge> sbg_edges_;
	std::vector<int> vv_;
	std::vector<bool> ee_;
	const G* graph_;
    public:
	SubgraphOfTheGraph(const G* g) :graph_(g) {}
	void push(const Edge& edge);
	void pop();
	int size() const { return sbg_edges_.size(); }
	const Edge& operator[] (int i) const { return sbg_edges_[i]; }
	const Edge& get_edge(int i) const    { return sbg_edges_[i]; }
	bool has_vertex(Vertex v) const { return vv_[v]; }
	bool has_edge(const Edge& edge) const { return ee_[edge.eidx]; }
	const G& get_graph() const { return *graph_; }
	const G* get_graph_p() const { return graph_; }

	template<class G_>
	friend std::ostream& operator<<(std::ostream& out, const SubgraphOfTheGraph<G_>& r);
    };

    template<class G>
    void SubgraphOfTheGraph<G>::push(const Edge& edge)
    {
	if (sbg_edges_.empty())
	{
	    sbg_edges_.reserve(num_edges(*graph_));
	    vv_.resize(num_vertices(*graph_), 0);
	    ee_.resize(num_edges(*graph_), false);
	}
	sbg_edges_.push_back(edge);
	++vv_[edge.from];
	++vv_[edge.to];
	ee_[edge.eidx] = true;
    }

    template<class G>
    void SubgraphOfTheGraph<G>::pop()
    {
	Edge edge = sbg_edges_.back();
	sbg_edges_.pop_back();
	--vv_[edge.from];
	--vv_[edge.to];
	ee_[edge.eidx] = false;

	assert(vv_[edge.from] >= 0);
	assert(vv_[edge.to] >= 0);
    }

    template<class G>
    std::ostream& operator<<(std::ostream& out, const SubgraphOfTheGraph<G>& r)
    {
	std::copy(r.sbg_edges_.begin(), r.sbg_edges_.end(), std::ostream_iterator<Edge>(out, " "));
	out << std::endl;
	std::copy(r.vv_.begin(), r.vv_.end(), std::ostream_iterator<bool>(out, " "));	
	return out;
    }

    // ****************************************************************************
    //				EdgeIterator
    // ****************************************************************************
    template<class G, class DirTag>
    class EdgeIterator;

    template<class G>
    struct EdgeIteratorSelect {
	typedef EdgeIterator<G, 
			     typename boost::graph_traits<G>::directed_category
			     //boost::bidirectionalS
			     > Type;
    };
    

    // -----------------------------------------------------
    //		EdgeIterator for boost::directedS
    // -----------------------------------------------------
    template<class G>
    class EdgeIterator<G, boost::directedS>
    {
    };

    // -----------------------------------------------------
    //		EdgeIterator for boost::bidirectionalS
    // -----------------------------------------------------
    template<class G>
    class EdgeIterator<G, boost::bidirectional_tag>
    {
    public:
	EdgeIterator(Vertex source, const G& g);
	bool is_end() const { return state_ == END; }
	const Edge& get_edge() const { return curr_edge_; }
	void increment();
    private:
	typedef typename boost::graph_traits<G>::edge_descriptor edge_descriptor;
	typedef typename boost::graph_traits<G>::out_edge_iterator out_edge_iterator;
	typedef typename boost::graph_traits<G>::in_edge_iterator in_edge_iterator;
	std::pair<out_edge_iterator,out_edge_iterator> ou_iters_;
	std::pair<in_edge_iterator,in_edge_iterator>   in_iters_;
	const G& graph_;
	enum State { OUT_ITER, IN_ITER, END };
	State state_;
	void set_edge_ou()
	    {
		curr_edge_.from = source(*ou_iters_.first, graph_);
		curr_edge_.to   = target(*ou_iters_.first, graph_);
		curr_edge_.eidx = get_edge_index(*ou_iters_.first);
	    }
	void set_edge_in()
	    {
		curr_edge_.from = target(*in_iters_.first, graph_);
		curr_edge_.to   = source(*in_iters_.first, graph_);
		curr_edge_.eidx = get_edge_index(*in_iters_.first);
	    }
	EdgeIndex get_edge_index(edge_descriptor e) const
	    { return get(boost::edge_index_t(), graph_, e); }
    protected:
	Edge curr_edge_;
    };

    template<class G>
    EdgeIterator<G, boost::bidirectional_tag>::EdgeIterator(Vertex v, const G& g)
	:ou_iters_(out_edges(v, g)),
	 in_iters_(in_edges(v, g)), graph_(g)
    {
	if (ou_iters_.first != ou_iters_.second)
	{
	    set_edge_ou();
	    state_ = OUT_ITER;
	}
	else if (in_iters_.first != in_iters_.second)
	{
	    set_edge_in();
	    state_ = IN_ITER;
	}
	else
	    state_ = END;
    }

    template<class G>
    void EdgeIterator<G, boost::bidirectional_tag>::increment()
    {
	switch (state_)
	{
	case OUT_ITER:
	    if (++ou_iters_.first != ou_iters_.second)
		set_edge_ou();
	    else
		if (in_iters_.first != in_iters_.second)
		{
		    set_edge_in();
		    state_ = IN_ITER;
		}
		else state_ = END;
	    break;

	case IN_ITER: 
	    if (++in_iters_.first != in_iters_.second)
		set_edge_in();
	    else
		state_ = END;
	    break;

	case END: assert(0); break;
	}
    }

    // -----------------------------------------------------
    //		EdgeIterator forward pure
    // -----------------------------------------------------
    template<class G>
    class ForwardPureEdgeIterator : public EdgeIteratorSelect<G>::Type
    {
	typedef typename EdgeIteratorSelect<G>::Type Base;
	VertexLabel minlabel_;
	const SubgraphOfTheGraph<G>* sbg_;
	void next_();
    public:
	ForwardPureEdgeIterator(const Edge& e, VertexLabel minlabel,
				const SubgraphOfTheGraph<G>* sbg)
	    :Base(e.to, sbg->get_graph()), minlabel_(minlabel), sbg_(sbg) { next_(); }
	void next() { this->increment(); next_(); }
    };

    template<class G>
    void ForwardPureEdgeIterator<G>::next_()
    {
	const G& g = sbg_->get_graph();
	for (; !this->is_end(); this->increment())
	{
	    const Edge& e = this->curr_edge_;
	    if (minlabel_ <= tolabel(e,g) && !sbg_->has_vertex(e.to))
		break;
	}
    }


    // -----------------------------------------------------
    //		EdgeIterator forward rmpath
    // -----------------------------------------------------
    template<class G>
    class ForwardRMPathEdgeIterator : public EdgeIteratorSelect<G>::Type
    {
	typedef typename EdgeIteratorSelect<G>::Type Base;
	VertexLabel minlabel_;
	const SubgraphOfTheGraph<G>* sbg_;
	Edge edge_;
	void next_();
    public:
	ForwardRMPathEdgeIterator(const Edge& e, VertexLabel minlabel,
				  const SubgraphOfTheGraph<G>* sbg)
	    :Base(e.from, sbg->get_graph()), minlabel_(minlabel), sbg_(sbg), edge_(e)
	    { next_(); }
	void next() { this->increment(); next_(); }
    };

    template<class G>
    void ForwardRMPathEdgeIterator<G>::next_()
    {
	const G& g = sbg_->get_graph();
	for (; !this->is_end(); this->increment())
	{
	    const Edge& e = this->curr_edge_;
	    if (e.eidx==4 and e.from==2)
		;//std::cerr << *sbg_ << std::endl;
	    VertexLabel tolabel2 = tolabel(e,g);
	    
	    if (sbg_->has_edge(e) || minlabel_ > tolabel2 || sbg_->has_vertex(e.to))
		continue;
	    
	    EdgeLabel elabel = edgelabel(edge_,g);
	    EdgeLabel eilabel = edgelabel(e,g);
	    if (elabel < eilabel || (elabel == eilabel && tolabel(edge_,g) <= tolabel2) )
		break;
	}
    }


    // -----------------------------------------------------
    //		EdgeIterator backward
    // -----------------------------------------------------
    template<class G>
    class BackwardEdgeIterator : public EdgeIteratorSelect<G>::Type
    {
	typedef typename EdgeIteratorSelect<G>::Type Base;
    public:
	BackwardEdgeIterator(const Edge& e1, const Edge& e2, const SubgraphOfTheGraph<G>* sbg);
    };

    template<class G>
    BackwardEdgeIterator<G>::BackwardEdgeIterator(const Edge& e1, const Edge& e2,
						  const SubgraphOfTheGraph<G>* sbg)
	:Base(e2.to, sbg->get_graph())
    {
	assert(e1.eidx != e2.eidx);
	const G& g = sbg->get_graph();
	EdgeLabel e1label = edgelabel(e1, g);
	bool b = tolabel(e1,g) <= tolabel(e2,g);
	for (; !this->is_end(); this->increment())
	{
	    const Edge& e = this->curr_edge_;
	    if (sbg->has_edge(e))
		continue;
	    EdgeLabel elabel = edgelabel(e,g);
	    if (e.to == e1.from && (e1label < elabel || (e1label == elabel && b) ) )
		break;
	}
    }


    // ****************************************************************************
    //				Projected
    // ****************************************************************************

    template<class G>
    class Projected : public std::vector<SubgraphOfTheGraph<G>* >
    {
    public:
	typedef std::vector<SubgraphOfTheGraph<G>*> Base;
	typedef SubgraphOfTheGraph<G> SBG;
	typedef std::pair<SBG*,Edge> ExtEdge;
	typedef std::list<ExtEdge> ExtEdgeList;

	Projected() {}

	Projected(const ExtEdgeList& exl)
	    {
		reserve(exl.size());
		for (typename ExtEdgeList::const_iterator i = exl.begin(); i != exl.end(); ++i)
		{
		    i->first->push(i->second);
		    push_back(i->first);
		}
	    }

	~Projected()
	    {
		for (typename Base::iterator i = this->begin(); i != this->end(); ++i)
		    (*i)->pop();
	    }

	int num_graphs() const
	    {
		std::set<const G*> n;
		BOOST_FOREACH(const SBG* p, *this) n.insert(p->get_graph_p());
		return n.size();
	    }
    };

    template<class G>
    std::ostream& operator<<(std::ostream& out, const Projected<G>& r)
    {
	typedef SubgraphOfTheGraph<G> SBG;
	BOOST_FOREACH(const SBG* p, r) out << *p << std::endl;
	return out;
    }


    // ****************************************************************************
    //				Maps
    // ****************************************************************************
    template<class K1, class T>
    struct Map1 : public std::map<K1, T>
    {
	typedef std::map<K1, T> Base;
	typedef typename Base::const_iterator M1_iterator;
    };

    template<class K1, class K2, class T>
    struct Map2 : public std::map<K1, std::map<K2, T> >
    {
	typedef std::map<K1, std::map<K2, T> > Base;
	typedef typename Base::const_iterator M1_iterator;
	typedef typename Base::mapped_type::const_iterator M2_iterator;
    };

    template<class K1, class K2, class K3, class T>
    struct Map3 : public std::map<K1, std::map<K2, std::map<K1, T> > >
    {
	typedef std::map<K1, std::map<K2, std::map<K1, T> > > Base;
	typedef typename Base::const_iterator				M1_iterator;
	typedef typename Base::const_reverse_iterator			M1_riterator;
	typedef typename Base::mapped_type::const_iterator		M2_iterator;
	typedef typename Base::mapped_type::mapped_type::const_iterator M3_iterator;
    };


    // ****************************************************************************
    //				gspan function
    // ****************************************************************************

    template<class G>
    void enumerate_one(Map3<VertexLabel,EdgeLabel,VertexLabel, Projected<G> >& m3,
		       std::list<SubgraphOfTheGraph<G> >& sbgs,
		       const G& g)
    {
	typedef SubgraphOfTheGraph<G> SBG;
	typename boost::graph_traits<G>::vertex_iterator vi, viend;
	for (boost::tie(vi,viend) = vertices(g); vi != viend; ++vi)
	{
	    typedef typename EdgeIteratorSelect<G>::Type EdgeIter;
	    for (EdgeIter iter(*vi, g); !iter.is_end(); iter.increment())
	    {
		Edge e = iter.get_edge();
		sbgs.push_back(SBG(&g));
		SBG* sbg = &sbgs.back();
		sbg->push(e);
		VertexLabel fromlab = fromlabel(e, g);
		VertexLabel tolab   = tolabel(e, g);
		EdgeLabel elab      = edgelabel(e, g);
		m3[fromlab][elab][tolab].push_back(sbg);
	    }
	}
    }


    template<class G>
    int support(const Projected<G>& projected)
    {
	return projected.num_graphs();
    }

    
    template<class G, class Output>
    void report(const Projected<G>& projected, const DFSCode& dfsc, Output& result)
    {
	std::cout << "report(): resulted " << dfsc << std::endl;
	result(projected, dfsc);
    }


    template<class G>
    bool project_is_min(const Projected<G>& projected, DFSCode& dfsc_min, const DFSCode& dfsc_tested)
    {
	typedef SubgraphOfTheGraph<G> SBG;
	typedef std::pair<SBG*,Edge> ExtEdge;
	typedef std::list<ExtEdge> ExtEdgeList;
	typedef Map1<EdgeLabel, ExtEdgeList> M1;
	typedef Map2<EdgeLabel, VertexLabel, ExtEdgeList> M2;

	RMPath rmpath(dfsc_min);
	Vertex maxtoc = dfsc_min[rmpath.rightmost()].to;

	// ---------------------------------------------------
	// backward
	{
	    M1 bck_edges;
	    Vertex newto = null_vertex();
	    bool flg = false;
	    for (unsigned int i = 0; !flg && i < rmpath.size() - 1; ++i)
	    {
		BOOST_FOREACH(SBG* sbg, projected)
		{
		    SBG& s = *sbg;
		    const G& g = s.get_graph();
		    BackwardEdgeIterator<G> iter(s[rmpath[i]], s[rmpath.rightmost()], sbg);
		    if (!iter.is_end())
		    {
			const Edge& e = iter.get_edge();
			ExtEdge ext(sbg, e);		    
			bck_edges[edgelabel(e,g)].push_back(ext);
			newto = dfsc_min[rmpath[i]].from;
			flg = true;
		    }
		}
	    }

	    if (flg)
	    {
		typename M1::M1_iterator i1 = bck_edges.begin();
		dfsc_min.push(EdgeCode(maxtoc, newto, -1, i1->first, -1));
		if (dfsc_min[dfsc_min.size()-1] != dfsc_tested[dfsc_min.size()-1])
		    return false;
		Projected<G> next_projected(i1->second);
		return project_is_min(next_projected, dfsc_min, dfsc_tested);
	    }
	}
	
	// ---------------------------------------------------
	// forward
	{
	    M2 fwd_edges;
	    VertexLabel minlabel = dfsc_min[0].fromlabel;
	    Vertex newfrom = null_vertex();
	    bool flg = false;
	    
	    // forward pure
	    BOOST_FOREACH(SBG* sbg, projected)
	    {
		SBG& s = *sbg;
		const G& g = s.get_graph();
		for (ForwardPureEdgeIterator<G> iter(s[rmpath.rightmost()], minlabel, sbg);
		     !iter.is_end(); iter.next())
		{
		    const Edge& e = iter.get_edge();
		    ExtEdge ext(sbg, e);
		    fwd_edges[edgelabel(e,g)][tolabel(e,g)].push_back(ext);
		    newfrom = maxtoc;
		    flg = true;
		}
	    }

	    // forward rmpath
	    for (int i = rmpath.size()-1; !flg && i >= 0; --i)
	    {
		BOOST_FOREACH(SBG* sbg, projected)
		{
		    SBG& s = *sbg;
		    const G& g = s.get_graph();
		    for (ForwardRMPathEdgeIterator<G> iter(s[rmpath[i]], minlabel, sbg);
			 !iter.is_end(); iter.next())
		    {
			const Edge& e = iter.get_edge();
			ExtEdge ext(sbg, e);
			fwd_edges[edgelabel(e,g)][tolabel(e,g)].push_back(ext);
			newfrom = dfsc_min[rmpath[i]].from;
			flg = true;
		    }
		}
	    }

	    if (flg)
	    {
		typename M2::M1_iterator i1 = fwd_edges.begin();
		typename M2::M2_iterator i2 = i1->second.begin();
		dfsc_min.push(EdgeCode(newfrom, maxtoc+1, -1, i1->first, i2->first));
		if (dfsc_min[dfsc_min.size()-1] != dfsc_tested[dfsc_min.size()-1])
		    return false;
		Projected<G> next_projected(i2->second);
		return project_is_min(next_projected, dfsc_min, dfsc_tested);
	    }
	}
	return true;
    }


    template<class G>
    bool is_min(const DFSCode& dfsc_tested)
    {
	typedef SubgraphOfTheGraph<G> SBG;
	typedef Map3<VertexLabel,EdgeLabel,VertexLabel, Projected<G> > M3;

	G graph(dfsc_tested);
	std::list<SBG> sbgs;
	M3 root;
	enumerate_one(root, sbgs, graph);

	DFSCode dfsc_min;
	typename M3::M1_iterator i1 = root.begin();
	typename M3::M2_iterator i2 = i1->second.begin();
	typename M3::M3_iterator i3 = i2->second.begin();
	dfsc_min.push(EdgeCode(0, 1, i1->first, i2->first, i3->first));
	return project_is_min(i3->second, dfsc_min, dfsc_tested);
    }


    template<class G, class Output>
    void project(const Projected<G>& projected, DFSCode& dfsc, int minsup, Output& result)
    {
	int sup = support(projected);
	if (sup < minsup)
	    return;

	if (! is_min<G>(dfsc))
	    return;

	report(projected, dfsc, result);

	// ---------------------------------------------------
	// enumerate
	typedef SubgraphOfTheGraph<G> SBG;
	typedef std::pair<SBG*,Edge> ExtEdge;
	typedef std::list<ExtEdge> ExtEdgeList;
	typedef Map2<Vertex, EdgeLabel, ExtEdgeList> M2;
	typedef Map3<Vertex, EdgeLabel, VertexLabel, ExtEdgeList> M3;

	RMPath rmpath(dfsc);

	Vertex maxtoc = dfsc[rmpath.rightmost()].to;
	VertexLabel minlabel = dfsc[0].fromlabel;

	M2 bck_edges;
	M3 fwd_edges;

	BOOST_FOREACH(SBG* sbg, projected)
	{
	    SBG& s = *sbg;
	    const G& g = s.get_graph();

	    // backward
	    for (unsigned int i = 0; i < rmpath.size() - 1; ++i)
	    {
		BackwardEdgeIterator<G> iter(s[rmpath[i]], s[rmpath.rightmost()], sbg);
		if (!iter.is_end())
		{
		    const Edge& e = iter.get_edge();
		    ExtEdge ext(sbg, e);
		    bck_edges[dfsc[rmpath[i]].from][edgelabel(e,g)].push_back(ext);
		}
	    }

	    // forward
	    for (ForwardPureEdgeIterator<G> iter(s[rmpath.rightmost()], minlabel, sbg);
		 !iter.is_end(); iter.next())
	    {
		const Edge& e = iter.get_edge();
		ExtEdge ext(sbg, e);
		fwd_edges[maxtoc][edgelabel(e,g)][tolabel(e,g)].push_back(ext);
	    }

	    for (int i = rmpath.size()-1; i >= 0; --i)
		for (ForwardRMPathEdgeIterator<G> iter(s[rmpath[i]], minlabel, sbg);
		     !iter.is_end(); iter.next())
		{
		    const Edge& e = iter.get_edge();
		    if (e.eidx==4 && e.from==2)
			;//std::cerr << "TTT: " << *sbg << std::endl;

		    ExtEdge ext(sbg, e);
		    fwd_edges[dfsc[rmpath[i]].from][edgelabel(e,g)][tolabel(e,g)].push_back(ext);
		}
	}

	// ---------------------------------------------------
	// test SBG + extended edge children

	// backward
	for (typename M2::M1_iterator it1 = bck_edges.begin(); it1 != bck_edges.end(); ++it1)
	    for (typename M2::M2_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
	    {
		EdgeCode ec(maxtoc, it1->first, -1, it2->first, -1);
		dfsc.push(ec);
		Projected<G> next_projected(it2->second);
		project(next_projected, dfsc, minsup, result);
		dfsc.pop();
	    }

	// forward
	for (typename M3::M1_riterator it1 = fwd_edges.rbegin(); it1 != fwd_edges.rend(); ++it1)
	    for (typename M3::M2_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
		for (typename M3::M3_iterator it3 = it2->second.begin(); it3 != it2->second.end(); ++it3)
		{
		    EdgeCode ec(it1->first, maxtoc+1, -1, it2->first, it3->first);
		    dfsc.push(ec);
		    Projected<G> next_projected(it3->second);
		    project(next_projected, dfsc, minsup, result);
		    dfsc.pop();
		}
	return;
    }


    template<class TGraphIterator, class Output>
    void gspan(TGraphIterator tg_begin, TGraphIterator tg_end, int minsup, Output& result)
    {
	typedef typename std::iterator_traits<TGraphIterator>::value_type G;
	typedef SubgraphOfTheGraph<G> SBG;
	typedef Map3<VertexLabel,EdgeLabel,VertexLabel, Projected<G> > M3;

	std::list<SBG> sbgs;
	M3 root;
	for (; tg_begin != tg_end; ++tg_begin)
	    enumerate_one(root, sbgs, *tg_begin);

	DFSCode dfsc;
	for (typename M3::M1_iterator i1 = root.begin(); i1 != root.end(); ++i1)
	    for (typename M3::M2_iterator i2 = i1->second.begin(); i2 != i1->second.end(); ++i2)
		for (typename M3::M3_iterator i3 = i2->second.begin(); i3 != i2->second.end(); ++i3)
		{
		    EdgeCode ec(0, 1, i1->first, i2->first, i3->first);
		    dfsc.push(ec);
		    project(i3->second, dfsc, minsup, result);
		    dfsc.pop();
		}
    }
}
#endif
