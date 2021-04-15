#include "tooling.hpp"

#include <iostream>
#include <cassert>

struct cmd_cfg: public tooling::cmd_cfg_base {
    //coverage ratio threshold
    double coverage_ratio = 0.;
};

static void process_cmdline(int argc, char **argv, cmd_cfg &cfg) {
    using namespace clipp;

    auto cli = (std::move(tooling::BaseCfg(cfg)), (
        (required("--cov-ratio") & number("value", cfg.coverage_ratio)) % "minimal coverage ratio to keep "
                                                                          "; should be < 1., (default: 0. -- effectively disabled)"
    ) % "algorithm settings");

    auto result = parse(argc, argv, cli);
    assert(!cfg.coverage.empty());
    if (!result) {
        std::cerr << make_man_page(cli, argv[0]);
        exit(1);
    }
}

//TODO consider making iterative right here after I can compress and track reads here
//TODO put coverage into GFA (check support in parcer, etc)
int main(int argc, char *argv[]) {
    cmd_cfg cfg;
    process_cmdline(argc, argv, cfg);

    assert(cfg.coverage_ratio <= 1.);
    std::cout << "Removing links with node coverage ratio less than " << cfg.coverage_ratio << std::endl;

    std::cout << "Reading coverage from " << cfg.coverage << std::endl;
    const auto segment_cov = utils::ReadCoverage(cfg.coverage);

    gfa::Graph g;
    std::cout << "Loading graph from GFA file " << cfg.graph_in << std::endl;
    g.open(cfg.graph_in);
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

        uint32_t baseline_cov = utils::get(segment_cov, g.segment_name(ds));

        for (gfa::LinkInfo l : g.outgoing_links(ds)) {
            auto nb_cov = utils::get(segment_cov, g.segment_name(l.end));
            if (double(nb_cov) / baseline_cov > cfg.coverage_ratio - 1e-5)
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

    tooling::OutputGraph(g, cfg, ndel, &segment_cov);
    std::cout << "END" << std::endl;
}
