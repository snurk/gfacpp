#include "tooling.hpp"

#include <iostream>
#include <cassert>
#include <set>

struct cmd_cfg: public tooling::cmd_cfg_base {
    double max_base_coverage = 0.;
    double min_path_coverage = 0.;
};

static void process_cmdline(int argc, char **argv, cmd_cfg &cfg) {
    using namespace clipp;

    auto cli = (std::move(tooling::BaseCfg(cfg)), (
            (required("--max-base-cov") & number("value", cfg.max_base_coverage)) % "maximal coverage of the surrounding nodes (required)",
            (required("--min-path-cov") & number("value", cfg.min_path_coverage)) % "minimal coverage along 'alternative' path (required)"
    ) % "algorithm settings");


    auto result = parse(argc, argv, cli);
    assert(!cfg.coverage.empty());
    if (!result) {
        std::cerr << "Removing 'shortcut' links if the connected segments have coverage less than max_base_coverage "
                     "and 'start' can be accessed by an unambiguous path back passing over the nodes of coverage no less than min_path_coverage "
                     "from the 'end' after exclusion of the 'shortcut'" << std::endl;
        std::cerr << make_man_page(cli, argv[0]);
        exit(1);
    }
}

inline bool UnambiguousBackwardPath(const gfa::Graph &g, gfa::DirectedSegment w, gfa::DirectedSegment v,
                                    const utils::SegmentCoverageMap &segment_cov,
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
                                           const utils::SegmentCoverageMap &segment_cov,
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
    cmd_cfg cfg;
    process_cmdline(argc, argv, cfg);

    std::cout << "Max base segment coverage set to " << cfg.max_base_coverage << std::endl;

    std::cout << "Reading coverage from " << cfg.coverage << std::endl;
    const auto segment_cov = utils::ReadCoverage(cfg.coverage);

    gfa::Graph g;
    std::cout << "Loading graph from GFA file " << cfg.graph_in << std::endl;
    g.open(cfg.graph_in);
    std::cout << "Segment cnt: " << g.segment_cnt() << "; link cnt: " << g.link_cnt() << std::endl;

    //std::set<std::string> neighbourhood;

    size_t l_ndel = 0;
    for (gfa::DirectedSegment v : g.directed_segments()) {
        if (g.outgoing_link_cnt(v) < 2)
            continue;

        DEBUG("Looking at directed node v " << g.str(v));
        if (utils::get(segment_cov, g.segment_name(v)) >= cfg.max_base_coverage) {
            DEBUG("Coverage of node v is too high");
            continue;
        }

        for (auto l : g.outgoing_links(v)) {
            assert(l.start == v);
            auto w = l.end;

            DEBUG("Looking at link to w " << g.str(w) << " (overlap size " << l.overlap() << ")");

            if (utils::get(segment_cov, g.segment_name(v)) >= cfg.max_base_coverage) {
                DEBUG("Coverage of node w is too high");
                continue;
            }

            if (UnambiguousBackwardAlternative(g, w, v, segment_cov, cfg.min_path_coverage)) {
                DEBUG("Unambiguous backward alternative found");
                std::cout << "Removing link " << g.str(v) << "," << g.str(w) << std::endl;
                //std::cout << "Removing link " << g.str(v) << " -> " << g.str(w) << std::endl;
                //TODO some links are counted 'twice' along with its conjugate
                ++l_ndel;
                g.DeleteLink(l);
            }
        }
    }

    tooling::OutputGraph(g, cfg, l_ndel, &segment_cov);
    std::cout << "END" << std::endl;
}

