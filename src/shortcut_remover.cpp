#include "wrapper.hpp"
#include "utils.hpp"

#include <vector>
#include <set>
#include <cassert>
#include <iostream>

inline bool UnambiguousBackwardPath(const gfa::Graph &g, gfa::DirectedSegment w, gfa::DirectedSegment v,
                                    const std::unordered_map<std::string, uint32_t> &segment_cov,
                                    uint32_t min_path_coverage) {
    assert(w != v);
    std::set<gfa::LinkInfo> used_links;
    while (g.unique_incoming(w)
            && utils::get(segment_cov, g.segment_name(w)) >= min_path_coverage
            && w != v) {
        auto l = *g.incoming_begin(w);
        if (used_links.count(l)) {
            DEBUG("Loop detected");
            return false;
        }
        //w = (*g.incoming_begin(w)).start;
        used_links.insert(l);
        w = l.start;
    }
    return w == v;
}

inline bool UnambiguousBackwardAlternative(const gfa::Graph &g, gfa::DirectedSegment w, gfa::DirectedSegment v,
                                           const std::unordered_map<std::string, uint32_t> &segment_cov,
                                           uint32_t min_path_coverage) {
    for (auto l : g.incoming_links(w)) {
        assert(l.end == w);
        auto w1 = l.start;
        if (w1 == v || g.outgoing_link_cnt(w1) > 1)
            continue;
        if (UnambiguousBackwardPath(g, w1, v, segment_cov, min_path_coverage)) {
            return true;
        }
    }
    return false;
}

int main(int argc, char *argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <gfa input> <gfa output> <node coverage file> <max base coverage (int)> <min path coverage (int)>" << std::endl;
        std::cerr << "Removing 'shortcut' links if the connected segments have coverage less than max_base_coverage "
                     "and 'start' can be accessed by an unambiguous path back passing over the nodes of coverage no less than min_path_coverage "
                     "from the 'end' after exclusion of the 'shortcut'" << std::endl;
        exit(239);
    }

    const std::string in_fn(argv[1]);
    const std::string out_fn(argv[2]);
    const std::string coverage_fn(argv[3]);
    const uint32_t max_base_coverage = std::stoi(argv[4]);
    const uint32_t min_path_coverage = std::stoi(argv[5]);

    std::cout << "Max base segment coverage set to " << max_base_coverage << std::endl;

    std::cout << "Reading coverage from " << coverage_fn << std::endl;
    const auto segment_cov = utils::ReadCoverage(coverage_fn);

    gfa::Graph g;
    std::cout << "Loading graph from GFA file " << in_fn << std::endl;
    g.open(in_fn);
    std::cout << "Segment cnt: " << g.segment_cnt() << "; link cnt: " << g.link_cnt() << std::endl;

    //std::set<std::string> neighbourhood;

    size_t l_ndel = 0;
    for (gfa::DirectedSegment v : g.directed_segments()) {
        if (g.outgoing_link_cnt(v) < 2)
            continue;

        DEBUG("Looking at directed node v " << g.str(v));
        if (utils::get(segment_cov, g.segment_name(v)) >= max_base_coverage) {
            DEBUG("Coverage of node v is too high");
            continue;
        }

        for (auto l : g.outgoing_links(v)) {
            assert(l.start == v);
            auto w = l.end;

            DEBUG("Looking at link to w " << g.str(w) << " (overlap size " << l.overlap() << ")");

            if (utils::get(segment_cov, g.segment_name(v)) >= max_base_coverage) {
                DEBUG("Coverage of node w is too high");
                continue;
            }

            if (UnambiguousBackwardAlternative(g, w, v, segment_cov, min_path_coverage)) {
                DEBUG("Unambiguous backward alternative found");
                std::cout << "Removing link " << g.str(v) << "," << g.str(w) << std::endl;
                //std::cout << "Removing link " << g.str(v) << " -> " << g.str(w) << std::endl;
                //TODO some links are counted 'twice' along with its conjugate
                ++l_ndel;
                g.DeleteLink(l);
            }
        }
    }

    std::cout << "Total of " << l_ndel << " links removed" << std::endl;
    if (l_ndel > 0)
        g.Cleanup();

    std::cout << "Writing output to " << out_fn << std::endl;
    g.write(out_fn);
    std::cout << "END" << std::endl;
}

