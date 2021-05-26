#include "superbubbles.hpp"
#include "tooling.hpp"

#include <vector>
#include <set>
#include <functional>
#include <algorithm>
#include <cassert>
#include <iostream>

struct cmd_cfg: public tooling::cmd_cfg_base {
    size_t max_length = 0;
    size_t max_diff = 0;
    bool use_coverage = false;
};

static void process_cmdline(int argc, char **argv, cmd_cfg &cfg) {
    using namespace clipp;

    auto cli = (tooling::BaseCfg(cfg), (
                (option("--max-length") & integer("value", cfg.max_length)) % "max (additional) bubble path length (default 20000)",
                (option("--max-diff") & integer("value", cfg.max_diff)) % "max bubble path length difference (default: 2000)",
                option("--use-coverage").set(cfg.use_coverage) % "use coverage instead of overlap sizes (default: false)"
                //option("--use-cov-ratios").set(cfg.use_cov_ratios) % "enable procedures based on unitig coverage ratios (default: false)",
    ) % "algorithm settings");

    auto result = parse(argc, argv, cli);
    if (cfg.use_coverage && cfg.coverage.empty()) {
        std::cerr << "Option to use coverage values was enabled, but coverage file wasn't provided" << std::endl;
        exit(2);
    }

    if (!result) {
        std::cerr << "Super-bubble removal" << std::endl;
        std::cerr << make_man_page(cli, argv[0]);
        exit(1);
    }
}

int main(int argc, char *argv[]) {
    cmd_cfg cfg;
    process_cmdline(argc, argv, cfg);
    std::cout << "Max length set to " << cfg.max_length << std::endl;
    std::cout << "Max length diff set to " << cfg.max_diff << std::endl;

    std::unique_ptr<utils::SegmentCoverageMap> segment_cov_ptr;
    if (!cfg.coverage.empty()) {
        std::cout << "Reading coverage from " << cfg.coverage << std::endl;
        segment_cov_ptr = std::make_unique<utils::SegmentCoverageMap>(utils::ReadCoverage(cfg.coverage));
    }

    gfa::Graph g;
    std::cout << "Loading graph from GFA file " << cfg.graph_in << std::endl;
    g.open(cfg.graph_in);
    std::cout << "Segment cnt: " << g.segment_cnt() << "; link cnt: " << g.link_cnt() << std::endl;

    bubbles::SuperbubbleFinder::SegmentCoverageF segment_cov_f;
    if (cfg.use_coverage) {
        assert(segment_cov_ptr);
        segment_cov_f = [&](gfa::DirectedSegment v) {
            return double(utils::get(*segment_cov_ptr, g.segment_name(v)));
        };
    }

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
        bubbles::SuperbubbleFinder finder(g, v, segment_cov_f, cfg.max_length, cfg.max_diff);
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

    //std::cout << "Total of " << l_ndel << " links and " << v_ndel << " segments removed" << std::endl;
    tooling::OutputGraph(g, cfg, (l_ndel + v_ndel) == 0 ? 0 : size_t(-1), segment_cov_ptr.get());
    std::cout << "END" << std::endl;
}
