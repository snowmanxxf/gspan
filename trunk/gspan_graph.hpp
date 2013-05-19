#ifndef GSPAN_GRAPH_H_
#define GSPAN_GRAPH_H_

#include "gspan_allocator.hpp"

#include <vector>
#include <list>
#include <limits>
#include <iterator>             // std::distance()
#include <algorithm>            // std::max()
#include <iostream>
#include <functional>           // std::less()

#include <boost/intrusive/set.hpp>
#include <boost/intrusive/list.hpp>
#include <boost/noncopyable.hpp>

#ifndef BR
#define BR asm volatile ("int3;")
#endif

namespace gSpan
{

    namespace detail
    {
        typedef unsigned short VI;
        typedef unsigned short EI;
    }

    typedef short VL;
    typedef short EL;

    const detail::VI VI_NULL = std::numeric_limits<detail::VI>::max();
    const detail::EI EI_NULL = std::numeric_limits<detail::EI>::max();
    const VL VL_NULL = -1;
    const EL EL_NULL = -1;

#ifdef TYPE_CHECK
    template<class T, int TypeID>
    class Value
    {
        T x_;
    public:
        Value() :x_() {}
        Value(const T x) :x_(x) {}
        operator T() const { return x_; }

        const T operator= (const T& x) { x_ = x; return x_; }
        const T operator++ () { return ++x_; }
        const T operator++ (int) { return x_++; }
        const T operator-- () { return --x_; }
        const T operator-- (int) { return x_--; }
    };

    typedef Value<detail::VI, 1>        GraphVI;
    typedef Value<detail::EI, 2>        GraphEI;
    typedef Value<detail::VI, 3>        DfscVI;
    typedef Value<detail::EI, 4>        DfscEI;
#else
    typedef detail::VI  GraphVI;
    typedef detail::EI  GraphEI;
    typedef detail::VI  DfscVI;
    typedef detail::EI  DfscEI;
#endif
    class Graph;

    namespace detail
    {
        // *****************************************************************************
        //                          Edge
        // *****************************************************************************

        class Edge : private boost::noncopyable
        {
            friend class gSpan::Graph;
        public:
            GraphVI vi_src() const              { return vi_src_; }
            GraphVI vi_dst() const              { return vi_dst_; }
            VL vl_src() const                   { return vl_src_; }
            VL vl_dst() const                   { return vl_dst_; }
            GraphEI eid() const                 { return eid_; }
            EL el() const                       { return el_; }
            const Edge* reverse() const         { return reverse_edge_; }
            const Edge& operator- () const      { return *reverse_edge_; }

            const Edge* next_adjacent() const   { return next_adjacent_; }
            bool operator< (const Edge& e) const
                {
                    if (vl_src_ < e.vl_src_)
                        return true;
                    if (vl_src_ == e.vl_src_ && el_ < e.el_)
                        return true;
                    if (vl_src_ == e.vl_src_ && el_ == e.el_ && vl_dst_ < e.vl_dst_)
                        return true;
                    return false;
                }
        private:
            GraphEI eid_;
            GraphVI vi_src_;
            GraphVI vi_dst_;
            VL vl_src_;
            VL vl_dst_;
            EL el_;

            Edge()
                :eid_(EI_NULL),
                 vi_src_(VI_NULL),
                 vi_dst_(VI_NULL),
                 vl_src_(VL_NULL),
                 vl_dst_(VL_NULL),
                 el_(EL_NULL),
                 next_adjacent_(0)
                {}

            Edge(GraphEI eid, GraphVI src, GraphVI dst, VL srclab, VL dstlab, EL elab)
                :eid_(eid),
                 vi_src_(src),
                 vi_dst_(dst),
                 vl_src_(srclab),
                 vl_dst_(dstlab),
                 el_(elab),
                 next_adjacent_(0)
                {}

        public:           
            // intrusive set hooks
            typedef 
            boost::intrusive::set_member_hook<
            //boost::intrusive::link_mode<boost::intrusive::normal_link>,
            boost::intrusive::optimize_size<true>
            > SetHook;
            SetHook edge_set_hook_;
            SetHook incid_set_hook_;
        private:
            Edge* reverse_edge_;
            Edge* next_adjacent_;
        };
    }

    // *****************************************************************************
    //                          Graph
    // *****************************************************************************

    template<class EdgeCodeIterator>
    std::size_t calc_num_vertices(EdgeCodeIterator it, const EdgeCodeIterator it_end)
    {
        std::size_t max_vi = 0U;
        while (it != it_end)
        {
            std::size_t m = std::max(it->vi_src(), it->vi_dst());
            if (max_vi < m)
                max_vi = m;
            ++it;
        }
        return max_vi + 1;
    }

    class Graph : private boost::noncopyable
    {
    public:

        // ----------------------------------------------
        // preferable for transactional graph creation
        // ----------------------------------------------
        template<class EdgeCodeIterator>
        Graph(EdgeCodeIterator it, const EdgeCodeIterator it_end)
            :num_vertices_(calc_num_vertices(it, it_end)),
             num_edges_(0),
             vertices_(new IncidEdges[num_vertices_]),
             edge_allocator_(sizeof(detail::Edge))
            {
                vertices_ = new IncidEdges[num_vertices_];                
                while (it != it_end) push_edge(*it++);
            }

        // preferable for internal use
        Graph(std::size_t max_num_vertices)
            :num_vertices_(0),
             num_edges_(0),
             vertices_(new IncidEdges[max_num_vertices]),
             edge_allocator_(sizeof(detail::Edge))
            { 
                float f = (max_num_vertices * max_num_vertices) * 0.6;
                history_.reserve(std::size_t(f));
            }

        ~Graph();

        std::size_t num_vertices() const                { return num_vertices_; }
        std::size_t num_edges() const                   { return num_edges_; }

        typedef detail::Edge            Edge;

        // Intrusive Sets
        typedef boost::intrusive::member_hook<Edge, Edge::SetHook, &Edge::edge_set_hook_> EdgeSetHookMemberOption;
        typedef boost::intrusive::member_hook<Edge, Edge::SetHook, &Edge::incid_set_hook_> IncidSetHookMemberOption;
        typedef boost::intrusive::multiset<Edge, EdgeSetHookMemberOption> EdgesSet;
        typedef boost::intrusive::multiset<Edge, IncidSetHookMemberOption> IncidentEdgesSet;        
        typedef EdgesSet::const_iterator EdgesIterator;

        const EdgesSet& edges() const         { return edges_; }

        // return pointer to first edge in the list
        const Edge* incident(GraphVI vi) const   { return vertices_[vi].adj_list(); }


        // preferable for internal use
        template<class ECode>
        void push_edge(const ECode& ec) { push_edge(ec.vi_src(), ec.vi_dst(), ec.vl_src(), ec.vl_dst(), ec.el()); }
        void push_edge(DfscVI src, DfscVI dst, VL srclab, VL dstlab, EL elab);
        void pop_edge();
    private:
        std::size_t num_vertices_;
        std::size_t num_edges_;

        class IncidEdges
        {
            IncidentEdgesSet*   set_;
            Edge*               adjacent_list_;
        public:
            IncidEdges() :set_(new IncidentEdgesSet), adjacent_list_(0) {}
            ~IncidEdges() { delete set_; }
            void add_edge(Edge* e);
            void del_edge(Edge* e);
            const Edge* adj_list() const { return adjacent_list_; }
        };

        IncidEdges* vertices_;  // array
        EdgesSet edges_;

        struct Record
        {
            Edge* edge;
            std::size_t num_vertices;
            Record(Edge* e, std::size_t nver) : edge(e), num_vertices(nver) {}
        };
        std::vector<Record> history_;

        FixedAllocator edge_allocator_;
    };


    inline std::ostream& operator<<(std::ostream& out, const Graph::Edge& e)
    {
        return out << e.eid() << "(" << e.vi_src() << "," << e.vi_dst()<<")";
    }
}
#endif
