#include "wrapper.hpp"
#include "utils.hpp"

#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <cassert>

//TODO consider making iterative right here after I can compress and track reads here
//TODO put coverage into GFA (check support in parcer, etc)
int main(int argc, char *argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <gfa input> <gfa output> <node coverage file> <coverage ratio threshold>" << std::endl;
        exit(239);
    }

    const std::string in_fn(argv[1]);
    const std::string out_fn(argv[2]);
    const std::string coverage_fn(argv[3]);
    const double coverage_ratio = std::stod(argv[4]);

    std::cout << "Reading coverage from " << coverage_fn << std::endl;
    const auto segment_cov = utils::ReadCoverage(coverage_fn);

    gfa::Graph g;
    std::cout << "Loading graph from GFA file " << in_fn << std::endl;
    g.open(in_fn);
    std::cout << "Segment cnt: " << g.segment_cnt() << "; link cnt: " << g.link_cnt() << std::endl;

    size_t ndel = 0;
    for (gfa::DirectedSegment ds : g.directed_segments()) {
        //std::cout << "Considering directed segment " << g.segment(ds.segment_id).name << " " << gfa::PrintDirection(ds.direction) << std::endl;

        uint32_t max_onb_cov = 0;
        for (gfa::LinkInfo l : g.outgoing_links(ds)) {
            //std::cerr << "Considering link between " <<
            //    g.segment(l.start.segment_id).name << " (" << gfa::PrintDirection(l.start.direction) << ") and " <<
            //    g.segment(l.end.segment_id).name << " (" << gfa::PrintDirection(l.end.direction) << "). Overlaps " <<
            //    l.start_overlap << " and " << l.end_overlap << std::endl;
            //std::cerr << "Use overlap " << ovl << std::endl;
            //std::cout << "Checking name " << g.segment(l.end.segment_id).name << std::endl;
            assert(l.start == ds);
            assert(segment_cov.find(g.segment_name(l.end)) != segment_cov.end());

            auto nb_cov = segment_cov.find(g.segment_name(l.end))->second;
            if (nb_cov > max_onb_cov)
                max_onb_cov = nb_cov;
        }

        uint32_t baseline_cov = segment_cov.find(g.segment_name(ds))->second;

        for (gfa::LinkInfo l : g.outgoing_links(ds)) {
            auto nb_cov = segment_cov.find(g.segment_name(l.end))->second;
            if (nb_cov > uint32_t(std::floor(coverage_ratio * baseline_cov)))
                continue;

            if (nb_cov == max_onb_cov)
                continue;

            //std::cout << "Removing link between " <<
            //    g.segment(l.start.segment_id).name << " (" << gfa::PrintDirection(l.start.direction) << ") and " <<
            //    g.segment(l.end.segment_id).name << " (" << gfa::PrintDirection(l.end.direction) << "). Overlaps " <<
            //    l.start_overlap << " and " << l.end_overlap << std::endl;
            //TODO optimize?
            g.DeleteLink(l);
            ndel++;
        }
    }

    std::cout << "Total of " << ndel << " links removed" << std::endl;
    if (ndel > 0)
        g.Cleanup();

    std::cout << "Writing output to " << out_fn << std::endl;
    g.write(out_fn);
}
