#include "wrapper.hpp"

#include <iostream>
#include <vector>
#include <set>
#include <cassert>

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <gfa input> <gfa output> <min overlap threshold>" << std::endl;
        exit(239);
    }

    const std::string in_fn(argv[1]);
    const std::string out_fn(argv[2]);
    uint32_t min_overlap = std::stoi(argv[3]);

    gfa::Graph g;
    std::cout << "Loading graph from GFA file " << in_fn << std::endl;
    g.open(in_fn);
    std::cout << "Segment cnt: " << g.segment_cnt() << "; link cnt: " << g.link_cnt() << std::endl;

    size_t ndel = 0;
    for (gfa::DirectedSegment ds : g.directed_segments()) {
        if (g.outgoing_link_cnt(ds) > 1) {
            for (gfa::LinkInfo l : g.outgoing_links(ds)) {
                if (std::max(l.start_overlap, l.end_overlap) < min_overlap) {
                    std::cout << "Removing link between " <<
                        g.segment(l.start.segment_id).name << " (" << gfa::PrintDirection(l.start.direction) << ") and " <<
                        g.segment(l.end.segment_id).name << " (" << gfa::PrintDirection(l.end.direction) << "). Overlaps " <<
                        l.start_overlap << " and " << l.end_overlap << std::endl;
                    //TODO optimize?
                    g.DeleteLink(l);
                    ndel++;
                }
            }
        }
    }
    std::cout << "Total of " << ndel << " links removed" << std::endl;
    if (ndel > 0)
        g.Cleanup();

    std::cout << "Writing output to " << out_fn << std::endl;
    g.write(out_fn);
}
