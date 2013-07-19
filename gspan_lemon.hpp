#ifndef GSPAN_LEMON_H_
#define GSPAN_LEMON_H_

#include "gspan.hpp"
#include <lemon/concepts/graph.h>
#include <lemon/maps.h>
#include <vector>

template<typename GR, typename NM, typename EM>
void create_gspan_graph(gSpan::Graph& gspan_graph,
                        std::vector<typename GR::Node>& node_ref,
                        std::vector<typename GR::Edge>& edge_ref,
                        const GR& lemon_graph,
                        const NM& nm,
                        const EM& em)
{
    int node_count = countNodes(lemon_graph);
    int edge_count = countEdges(lemon_graph);
    gspan_graph.init(node_count, edge_count);
    node_ref.resize(node_count);
    edge_ref.resize(edge_count);

    lemon::RangeIdMap<GR,typename GR::Node> node_id(lemon_graph);
    for (int i = 0; i < node_count; ++i)
        node_ref[i] = node_id(i);

    lemon::RangeIdMap<GR,typename GR::Edge> edge_id(lemon_graph);
    for (int i = 0; i < edge_count; ++i)
    {
        typename GR::Edge e = edge_id(i);
        edge_ref[i] = e;
        typename GR::template Node u = lemon_graph.u(e);
        typename GR::template Node v = lemon_graph.v(e);
        gspan_graph.push_edge(node_id[u], node_id[v], nm[u], nm[v], em[e]);
    }
}


template<typename GR, typename NM, typename EM>
void create_lemon_graph(const gSpan::Graph& gspan_graph,
                        std::vector<typename GR::Node>& node_ref,
                        std::vector<typename GR::Edge>& edge_ref,
                        GR& lemon_graph,
                        NM& nm,
                        EM& em)
{
    int num_vertices    = gspan_graph.num_vertices();
    int num_edges       = gspan_graph.num_edges();
    node_ref.resize(num_vertices);
    edge_ref.resize(num_edges);
    for (int i = 0; i < num_vertices; ++i)
        node_ref[i] = lemon_graph.addNode();

    const gSpan::Graph::EdgesSet& eset = gspan_graph.edges();
    const gSpan::Graph::Edge* opp = eset.begin()->reverse();
    for (gSpan::Graph::EdgesSet::const_iterator i = eset.begin(); i != eset.end(); ++i)
    {
        if (&*i == opp && i != eset.begin())
            break;
        typename GR::Edge e = lemon_graph.addEdge(node_ref[i->vi_src()],
                                                  node_ref[i->vi_dst()]);
        edge_ref[i->eid()] = e;

        nm.set(node_ref[i->vi_src()], i->vl_src());
        nm.set(node_ref[i->vi_dst()], i->vl_dst());
        em.set(e, i->el());
    }
}


template<typename M>
void make_embd_node_map(M& m,
                        const gSpan::SBG* sbg,
                        const std::vector<typename M::Key>& node_ref1,
                        const std::vector<typename M::Value>& node_ref2)
{
    const gSpan::GraphVI* vi_dfsc_to_graph = sbg->get_dfsc_to_graph_v();
    std::size_t num_vertices = sbg->num_vertices();
    
    for (gSpan::DfscVI i = 0; i < num_vertices; ++i)
        m.set(node_ref1[vi_dfsc_to_graph[i]], node_ref2[i]);
}


template<typename M>
void make_embd_edge_map(M& m,
                        const gSpan::SBG* sbg,
                        const std::vector<typename M::Key>& edge_ref1,
                        const std::vector<typename M::Value>& edge_ref2)
{
    std::vector<const gSpan::SBG*> ch;
    gSpan::get_chain(ch, sbg);

    int n = ch.size();
    for (int i = 0; i < n; ++i)
        m.set(edge_ref1[ch[i]], edge_ref2[i]);
}


#endif
