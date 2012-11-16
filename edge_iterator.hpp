#ifndef EDGE_ITERATOR_H_
#define EDGE_ITERATOR_H_

#include <boost/graph/graph_traits.hpp>

template<class T>
struct Deref
{
    template<class I>
    T operator() (const I& i) const { return *i; }
};

struct RetTrue
{
    template<class T>
    bool operator() (const T&) { return true; }
};

template<class G,
	 class E = typename boost::graph_traits<G>::edge_descriptor,
	 class GetE = Deref<E>,
	 class DirTag = typename boost::graph_traits<G>::directed_category>
class EdgeIterator;

template<class G, class E, class GetE>
class EdgeIterator<G,E,GetE,boost::bidirectional_tag>
{
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G>::edge_descriptor   edge_descriptor;
    typedef typename boost::graph_traits<G>::out_edge_iterator out_edge_iterator;
    typedef typename boost::graph_traits<G>::in_edge_iterator  in_edge_iterator;
    std::pair<out_edge_iterator,out_edge_iterator> ou_iters_;
    std::pair<in_edge_iterator,in_edge_iterator>   in_iters_;
    enum State { OUT_ITER, IN_ITER, END } state_;
    const GetE get_e_;
    E edge_;
    void init();
public:
    typedef E value_type;

    EdgeIterator(vertex_descriptor v, const G& g, const GetE& get_e = Deref<E>())
	:ou_iters_(out_edges(v,g)), in_iters_(in_edges(v,g)),
	 get_e_(get_e) { init(); }

    bool is_end() const
	{
	    return state_ == END;
	}

    const value_type& edge() const { return edge_; }
    void increment();
};

template<class G, class E, class GetE>
void EdgeIterator<G, E, GetE, boost::bidirectional_tag>::init()
{
    if (ou_iters_.first != ou_iters_.second)
    {
	edge_ = get_e_(ou_iters_.first);
	state_ = OUT_ITER;
    }
    else if (in_iters_.first != in_iters_.second)
    {
	edge_ = get_e_(in_iters_.first);
	state_ = IN_ITER;
    }
    else
	state_ = END;
}

template<class G, class E, class GetE>
void EdgeIterator<G, E, GetE, boost::bidirectional_tag>::increment()
{
    switch (state_)
    {
    case OUT_ITER:
	if (++ou_iters_.first != ou_iters_.second)
	    edge_ = get_e_(ou_iters_.first);
	else if (in_iters_.first != in_iters_.second)
	{
	    edge_ = get_e_(in_iters_.first);
	    state_ = IN_ITER;
	}
	else
	    state_ = END;
	break;
    case IN_ITER: 
	if (++in_iters_.first != in_iters_.second)
	    edge_ = get_e_(in_iters_.first);
	else
	    state_ = END;
	break;
    case END: assert(0); break;
    }
}



template<class Iter, class Pred, class OnFalse>
class EdgeConditionalIterator : private Iter
{
    Pred pred_;
    OnFalse onfalse_;
    void increment_();
public:
    EdgeConditionalIterator(const Iter& iter,
			    const Pred& pred,
			    OnFalse onfalse)
	:Iter(iter), pred_(pred), onfalse_(onfalse)
	{ increment_(); }

    bool is_end() const { return Iter::is_end(); }
    const typename Iter::value_type& edge() const { return Iter::edge(); }
    void increment() { Iter::increment(); increment_(); }
};

template<class Iter, class Pred, class OnFalse>
void EdgeConditionalIterator<Iter,Pred,OnFalse>::increment_()
{
    for (; !Iter::is_end(); Iter::increment())
    {
	if (pred_(Iter::edge()))
	    break;
	onfalse_(Iter::edge());
    }
}

#endif
