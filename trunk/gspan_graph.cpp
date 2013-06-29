
#include "gspan_graph.hpp"


namespace gSpan
{
    void Graph::IncidEdges::add_edge(Edge* e)
    {
        assert(0 == e->next_adjacent_);

        const IncidentEdgesSet::iterator edge_iter = set_.insert(*e);
        IncidentEdgesSet::iterator right_iter = edge_iter; ++right_iter;
        IncidentEdgesSet::iterator left_iter = edge_iter;
        if (left_iter != set_.begin())
        {
            // insert after left_iter

            --left_iter;
            
            assert(right_iter != set_.end() ?
                   left_iter->next_adjacent_ == &*right_iter :
                   left_iter->next_adjacent_ == 0
                );

            e->next_adjacent_ = left_iter->next_adjacent_;
            left_iter->next_adjacent_ = e;
        }
        else
        {
            // insert to the beginning
            e->next_adjacent_ = adjacent_list_;
            adjacent_list_ = e;
        }
    }


    void Graph::IncidEdges::del_edge(Edge* e)
    {
        const IncidentEdgesSet::iterator edge_iter = set_.iterator_to(*e);
        IncidentEdgesSet::iterator left_iter = edge_iter;
        if (left_iter != set_.begin())
        {
            --left_iter;
            left_iter->next_adjacent_ = e->next_adjacent_;
        }
        else
        {
            adjacent_list_ = e->next_adjacent_;
        }

        set_.erase(edge_iter);
    }


    Graph::~Graph()
    {
        while (num_edges_ != 0)
            pop_edge();        
        delete[] history_;
        delete[] vertices_;
    }

    GraphEI Graph::push_edge(DfscVI src, DfscVI dst, VL srclab, VL dstlab, EL elab)
    {
        assert(src != VI_NULL);
        assert(dst != VI_NULL);

	const GraphEI eid = num_edges_;

	Edge* edge1 = new (edge_allocator_.allocate()) Edge(eid, GraphVI(src), GraphVI(dst), srclab, dstlab, elab);
	Edge* edge2 = new (edge_allocator_.allocate()) Edge(eid, GraphVI(dst), GraphVI(src), dstlab, srclab, elab);
        edge1->reverse_edge_ = edge2;
        edge2->reverse_edge_ = edge1;

        // add to EdgesSet
	edges_.insert(*edge1);
	edges_.insert(*edge2);

        // add to IncidEdges
        vertices_[src].add_edge(edge1);
        vertices_[dst].add_edge(edge2);

        // push to history stack and update counters
        history_[num_edges_++] = Record(edge1, num_vertices_);

        std::size_t max_vi = std::max(src, dst);
        if (num_vertices_ < max_vi + 1)
            num_vertices_ = max_vi + 1;

        return eid;
    }


    void Graph::pop_edge()
    {
        using namespace boost::intrusive;

        Edge* edge1 = history_[num_edges_ - 1].edge;
        Edge* edge2 = edge1->reverse_edge_;

        // remove from IncidEdges
        vertices_[edge1->vi_src()].del_edge(edge1);
        vertices_[edge2->vi_src()].del_edge(edge2);

        // remove from EdgesSet
        edges_.erase(edges_.iterator_to(*edge1));
        edges_.erase(edges_.iterator_to(*edge2));        

        edge1->~Edge();
        edge2->~Edge();
        edge_allocator_.deallocate(edge1);
        edge_allocator_.deallocate(edge2);
        
        assert(num_edges_ > 0);
        assert(num_vertices_ > 0);
        
        //  pop history, update counters
        num_vertices_ = history_[--num_edges_].num_vertices;
    }

}
