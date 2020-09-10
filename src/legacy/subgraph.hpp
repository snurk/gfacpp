#pragma once

#include "graph.hpp"

static void Go(const Graph &g, size_t v_id, uint32_t curr_depth, uint32_t max_depth, std::vector<uint32_t> &considered) {
    assert(v_id < g.node_cnt());
    //std::cout << "In node " << g.nodes[v_id] << " on depth " << curr_depth << std::endl;
    considered[v_id] = curr_depth;
    if (curr_depth == max_depth) {
        //std::cout << "Hit max depth " << max_depth << std::endl;
        return;
    }
    for (size_t n : g.adjacent(v_id)) {
        //not considered or considered at a higher depth
        if (considered[n] > curr_depth + 1) {
            Go(g, n, curr_depth + 1, max_depth, considered);
        }
    }
}

static std::set<std::string> CollectNeighborhood(const Graph &g, const std::vector<std::string> &nodes_of_interest, uint32_t max_depth) {
    std::vector<size_t> ids_of_interest;
    ids_of_interest.reserve(nodes_of_interest.size());
    for (const std::string &n : nodes_of_interest) {
        ids_of_interest.push_back(g.id(n));
    }

    std::cout << "Searching for neighbourhood" << std::endl;
    std::vector<uint32_t> considered(g.node_cnt(), std::numeric_limits<uint32_t>::max());
    for (size_t n : ids_of_interest) {
        considered[n] = 0;
    }
    for (size_t n : ids_of_interest) {
        assert(n < g.node_cnt());
        std::cout << "Starting from node " << g.str(n) << std::endl;
        Go(g, n, 0, max_depth, considered);
        //std::cout << "Node " << g.node(n) << " processed" << std::endl;
    }

    std::cout << "Extracting the subgraph" << std::endl;
    std::set<std::string> reached_segments;
    for (size_t i = 0, n = considered.size(); i < n; ++i)
        if (considered[i] < std::numeric_limits<uint32_t>::max())
            reached_segments.insert(g.str(i));

    return reached_segments;
}
