#include "tooling.hpp"

#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <cassert>

struct cmd_cfg: public tooling::cmd_cfg_base {
    //coverage ratio threshold
    double coverage_ratio = 0.;
};

static void process_cmdline(int argc, char **argv, cmd_cfg &cfg) {
    using namespace clipp;

    auto cli = (
            cfg.graph_in << value("input file in GFA (ending with .gfa)"),
            cfg.graph_out << value("output file"),
            (required("--coverage") & value("file", cfg.coverage)) % "file with coverage information",
            (required("--cov-ratio") & number("value", cfg.coverage_ratio)) % "minimal coverage ratio to keep "
                                                                              "; should be < 1., (default: 0. -- effectively disabled)",
            option("--compact").set(cfg.compact) % "compact the graph after cleaning (default: false)",
            (option("--id-mapping") & value("file", cfg.id_mapping)) % "file with compacted segment id mapping",
            (option("--prefix") & value("vale", cfg.compacted_prefix)) % "prefix used to form compacted segment names",
            option("--drop-sequence").set(cfg.drop_sequence) % "flag to drop sequences even if present in original file (default: false)"
            //option("--use-cov-ratios").set(cfg.use_cov_ratios) % "enable procedures based on unitig coverage ratios (default: false)",
            //(required("-k") & integer("value", cfg.k)) % "k-mer length to use",
    );


    auto result = parse(argc, argv, cli);
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

        uint32_t baseline_cov = segment_cov.find(g.segment_name(ds))->second;

        for (gfa::LinkInfo l : g.outgoing_links(ds)) {
            auto nb_cov = segment_cov.find(g.segment_name(l.end))->second;
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
