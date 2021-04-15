#include "tooling.hpp"

#include <iostream>
#include <cassert>
#include <vector>
#include <set>

struct cmd_cfg: public tooling::cmd_cfg_base {
    //coverage ratio threshold
    int32_t min_overlap = 0;
    bool prevent_deadends = false;
};

static void process_cmdline(int argc, char **argv, cmd_cfg &cfg) {
    using namespace clipp;

    auto cli = (std::move(tooling::BaseCfg(cfg)), (
            (required("--min-overlap") & integer("value", cfg.min_overlap)) % "overlap size threshold (default: 0)",
            option("--prevent-deadends").set(cfg.prevent_deadends) % "check that no new dead-ends are formed (default: false)"
    ) % "algorithm settings");

    auto result = parse(argc, argv, cli);
    if (!result) {
        std::cerr << "Removing overlaps weaker than the provided threshold" << std::endl;
        std::cerr << make_man_page(cli, argv[0]);
        exit(1);
    }
}

//todo replace with std::any_of?
inline bool CheckHasStrongIncoming(const gfa::Graph &g, gfa::DirectedSegment v, int32_t weak_thr) {
    for (const auto &l : g.incoming_links(v)) {
        auto ovl = std::max(l.start_overlap, l.end_overlap);
        if (ovl >= weak_thr) {
            return true;
        }
    }
    return false;
}

//TODO consider making iterative right here after I can compress and track reads here
//TODO put coverage into GFA (check support in parser, etc)
int main(int argc, char *argv[]) {
    cmd_cfg cfg;
    process_cmdline(argc, argv, cfg);

    std::unique_ptr<utils::SegmentCoverageMap> segment_cov_ptr;
    if (!cfg.coverage.empty()) {
        std::cout << "Reading coverage from " << cfg.coverage << std::endl;
        segment_cov_ptr = std::make_unique<utils::SegmentCoverageMap>(utils::ReadCoverage(cfg.coverage));
    }

    gfa::Graph g;
    std::cout << "Loading graph from GFA file " << cfg.graph_in << std::endl;
    g.open(cfg.graph_in);
    std::cout << "Segment cnt: " << g.segment_cnt() << "; link cnt: " << g.link_cnt() << std::endl;

    size_t ndel = 0;
    for (gfa::DirectedSegment ds : g.directed_segments()) {
        //std::cerr << "Considering directed segment " << g.segment(ds.segment_id).name << " " << gfa::PrintDirection(ds.direction) << std::endl;

        int32_t max_ovl = 0;
        for (gfa::LinkInfo l : g.outgoing_links(ds)) {
            //std::cerr << "Considering link between " <<
            //    g.segment(l.start.segment_id).name << " (" << gfa::PrintDirection(l.start.direction) << ") and " <<
            //    g.segment(l.end.segment_id).name << " (" << gfa::PrintDirection(l.end.direction) << "). Overlaps " <<
            //    l.start_overlap << " and " << l.end_overlap << std::endl;
            //std::cerr << "Use overlap " << ovl << std::endl;
            auto ovl = std::max(l.start_overlap, l.end_overlap);
            if (ovl > max_ovl)
                max_ovl = ovl;
        }

        for (gfa::LinkInfo l : g.outgoing_links(ds)) {
            auto ovl = std::max(l.start_overlap, l.end_overlap);
            if (ovl >= cfg.min_overlap)
                continue;
            if (max_ovl < cfg.min_overlap && ovl == max_ovl)
                continue;

            if (cfg.prevent_deadends) {
                if (!CheckHasStrongIncoming(g, l.end, cfg.min_overlap)) {
                    DEBUG("Not removing link " << g.str(l) << " because end has no strong alternatives");
                    continue;
                }
            }

            INFO("Removing link " << g.str(l) << ". Overlaps " <<
                l.start_overlap << " and " << l.end_overlap);
            //TODO optimize?
            g.DeleteLink(l);
            ndel++;
        }
    }

    tooling::OutputGraph(g, cfg, ndel, segment_cov_ptr.get());
    std::cout << "END" << std::endl;
}
