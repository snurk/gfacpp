#include "tooling.hpp"

#include <iostream>
#include <cassert>
#include <compact.hpp>

static void process_cmdline(int argc, char **argv,
                            tooling::cmd_cfg_base &cfg) {
    using namespace clipp;

    auto cli = tooling::BaseCfg(cfg);

    auto result = parse(argc, argv, cli);
    cfg.compact = true;

    if (!result) {
        std::cerr << "Graph normalizer & compactifier (contrary to doc '--compact' always enabled)" << std::endl;
        std::cerr << make_man_page(cli, argv[0]);
        exit(1);
    }
}

int main(int argc, char *argv[]) {
    tooling::cmd_cfg_base cfg;
    process_cmdline(argc, argv, cfg);

    const std::string in_fn(cfg.graph_in);
    const std::string out_fn(cfg.graph_out);

    std::unique_ptr<utils::SegmentCoverageMap> segment_cov_ptr;
    if (!cfg.coverage.empty()) {
        std::cout << "Reading coverage from " << cfg.coverage << std::endl;
        segment_cov_ptr = std::make_unique<utils::SegmentCoverageMap>(utils::ReadCoverage(cfg.coverage));
    }

    gfa::Graph g;
    std::cout << "Loading graph from GFA file " << in_fn << std::endl;
    g.open(in_fn);

    std::cout << "Segment cnt: " << g.segment_cnt() << "; link cnt: " << g.link_cnt() << std::endl;
    //gfa::CompactAndWrite(g, out_fn);
    gfa::Compactifier compactifier(g, cfg.compacted_prefix, segment_cov_ptr.get(), cfg.dbg_k, /*normalize overlaps*/true);
    std::cout << "Writing compacted graph to " << out_fn << std::endl;
    compactifier.Compact(out_fn, cfg.id_mapping, cfg.drop_sequence, cfg.rename_all);
    std::cout << "Writing complete" << std::endl;
    std::cout << "Done" << std::endl;
}
