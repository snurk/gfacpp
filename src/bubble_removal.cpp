#include "superbubbles.hpp"
#include "utils.hpp"

#include <vector>
#include <set>
#include <cassert>
#include <iostream>

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <gfa input> <gfa output> [max_length (20000)] [max_diff (2000)]" << std::endl;
        //std::cerr << "Usage: " << argv[0] << " <gfa input>" << std::endl;
        exit(239);
    }

    const std::string in_fn(argv[1]);
    const std::string out_fn(argv[2]);

    size_t max_length = 20000;
    if (argc > 3) {
        max_length = std::stoi(argv[3]);
    }
    std::cout << "Max length set to " << max_length << std::endl;

    size_t max_diff = 2000;
    if (argc > 4) {
        max_diff = std::stoi(argv[4]);
    }
    std::cout << "Max length diff set to " << max_diff << std::endl;

    gfa::Graph g;
    std::cout << "Loading graph from GFA file " << in_fn << std::endl;
    g.open(in_fn);
    std::cout << "Segment cnt: " << g.segment_cnt() << "; link cnt: " << g.link_cnt() << std::endl;

    //std::set<std::string> neighbourhood;

    std::cout << "Searching for bubbles" << std::endl;
    std::set<gfa::DirectedSegment> v_in_bubble;
    std::set<std::pair<gfa::DirectedSegment, gfa::DirectedSegment>> l_in_bubble;
    //consist of heaviest paths of outermost bubbles
    std::set<uint32_t> segments_to_keep;
    std::set<std::pair<gfa::DirectedSegment, gfa::DirectedSegment>> links_to_keep;
    //FIXME replace with iterator
    for (gfa::DirectedSegment v : g.directed_segments()) {
        DEBUG("Looking at directed node " << g.str(v));
        if (v_in_bubble.count(v) != 0) {
            DEBUG("Not considering. Was part of bubble.");
            continue;
        }
        bubbles::SuperbubbleFinder finder(g, v, max_length, max_diff);
        if (finder.FindSuperbubble()) {
            std::cout << "Found superbubble between " << g.str(finder.start_vertex()) << " and " << g.str(finder.end_vertex()) << std::endl;
            for (gfa::DirectedSegment v : finder.segments()) {
                std::cout << g.str(v) << '\n';
                //Updating sets of segments and links belonging to all bubbles
                //And resetting the 'keep' marks within within the bubble

                //end vertex can be start of a different bubble
                if (v != finder.end_vertex()) {
                    v_in_bubble.insert(v);
                    for (auto l : g.outgoing_links(v)) {
                        //TODO move to canonical?!
                        l_in_bubble.insert(std::make_pair(l.start, l.end));
                        l_in_bubble.insert(std::make_pair(l.end.Complement(), l.start.Complement()));

                        links_to_keep.erase(std::make_pair(l.start, l.end));
                        links_to_keep.erase(std::make_pair(l.end.Complement(), l.start.Complement()));
                    }
                }

                //complement of the start vertex vertex can be start of a different bubble
                if (v != finder.start_vertex()) {
                    v_in_bubble.insert(v.Complement());
                }

                segments_to_keep.erase(v.segment_id);
            }

            if (finder.segments().size() == finder.HeaviestPath().segment_cnt()) {
                std::cout << "New processing only" << std::endl;
            }

            //Putting new 'keep' marks
            gfa::Path heaviest_path = finder.HeaviestPath();
            for (gfa::DirectedSegment v : heaviest_path.segments) {
                std::cout << "Keeping node " << g.str(v) << std::endl;
                segments_to_keep.insert(v.segment_id);
            }
            for (gfa::LinkInfo l : heaviest_path.links) {
                std::cout << "Keeping link " << g.str(l) << std::endl;
                //TODO move to canonical?!
                links_to_keep.insert(std::make_pair(l.start, l.end));
                links_to_keep.insert(std::make_pair(l.end.Complement(), l.start.Complement()));
            }
        }
    }

    //for (size_t id = 0, n = g.node_cnt(); id < n; ++id) {
    //    neighbourhood.insert(g.str(id));
    //}

    size_t l_ndel = 0;
    size_t v_ndel = 0;

    //TODO currently same link will be marked for deletion twice
    for (auto l_desc : l_in_bubble) {
        if (links_to_keep.count(l_desc) == 0) {
            ++l_ndel;
            g.DeleteLink(l_desc.first, l_desc.second);
        }
    }

    //TODO currently same segment can be marked for deletion twice
    for (auto v : v_in_bubble) {
        if (segments_to_keep.count(v.segment_id) == 0) {
            ++v_ndel;
            g.DeleteSegment(v);
        }
    }

    std::cout << "Total of " << l_ndel << " links and " << v_ndel << " segments removed" << std::endl;
    if (l_ndel > 0 || v_ndel > 0)
        g.Cleanup();

    std::cout << "Writing output to " << out_fn << std::endl;
    g.write(out_fn);
    std::cout << "END" << std::endl;
}
