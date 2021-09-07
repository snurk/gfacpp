#pragma once
#include "wrapper.hpp"
#include "utils.hpp"

namespace gfa {

inline Path UnambiguousPathForward(const Graph &g, DirectedSegment v) {
    Path p(v);
    std::set<LinkInfo> used_links;
    while (g.unique_outgoing(v)) {
        auto l = *g.outgoing_begin(v);
        if (used_links.count(l)) {
            DEBUG("Loop detected");
            break;
        }
        used_links.insert(l);
        p.Extend(l);
        v = l.end;
    }
    return p;
}

}
