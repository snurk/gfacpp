#include "tooling.hpp"

#include <iostream>
#include <cassert>
#include <compact.hpp>

static void process_cmdline(int argc, char **argv,
                            tooling::cmd_cfg_base &cfg) {
    using namespace clipp;

    auto cli = (
            cfg.graph_in << value("input file in GFA (ending with .gfa)"),
            cfg.graph_out << value("output file"),
            (option("--coverage") & value("file", cfg.coverage)) % "file with coverage information",
            (option("--id-mapping") & value("file", cfg.id_mapping)) % "file with compacted segment id mapping",
            (option("--prefix") & value("vale", cfg.compacted_prefix)) % "prefix used to form compacted segment names",
            option("--drop-sequence").set(cfg.drop_sequence) % "flag to drop sequences even if present in original file (default: false)"
            //(required("-k") & integer("value", cfg.k)) % "k-mer length to use",
    );

    cfg.compact = true;

    auto result = parse(argc, argv, cli);
    if (!result) {
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
    gfa::Compactifier compactifier(g, cfg.compacted_prefix, segment_cov_ptr.get(), /*normalize overlaps*/true);
    std::cout << "Writing compacted graph to " << out_fn << std::endl;
    compactifier.Compact(out_fn, cfg.id_mapping, cfg.drop_sequence);
    std::cout << "Writing complete" << std::endl;
    std::cout << "Done" << std::endl;
}
