#include "wrapper.hpp"
#include "utils.hpp"

#include <vector>
#include <set>
#include <cassert>
#include <iostream>

int main(int argc, char *argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <gfa input> <gfa output> <node coverage file> <max base coverage (int)>" << std::endl;
        std::cerr << "Loop link will be killed if coverage of the node is < max_base_coverage and other links are present" << std::endl;
        exit(239);
    }

    const std::string in_fn(argv[1]);
    const std::string out_fn(argv[2]);
    const std::string coverage_fn(argv[3]);
    const uint32_t max_base_coverage = std::stoi(argv[4]);
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
        //TODO add forward-only iterator
        if (v.direction == gfa::Direction::REVERSE)
            continue;
        if (g.outgoing_link_cnt(v) < 2 || g.incoming_link_cnt(v) < 2)
            continue;

        DEBUG("Looking at directed node " << g.str(v));
        for (auto l : g.outgoing_links(v)) {
            assert(l.start == v);
            if (l.end != v)
                continue;
            if (utils::get(segment_cov, g.segment_name(v)) < max_base_coverage) {
                std::cout << "Removing loop link from segment " << g.str(v) << '\n';
                g.DeleteLink(l);
                ++l_ndel;
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
