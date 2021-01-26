#include "tooling.hpp"

#include <vector>
#include <set>
#include <cassert>

struct cmd_cfg: public tooling::cmd_cfg_base {
    //coverage ratio threshold
    size_t max_length = 10000;
    double cov_thr = -1.;
};

static void process_cmdline(int argc, char **argv, cmd_cfg &cfg) {
    using namespace clipp;

    auto cli = (
            cfg.graph_in << value("input file in GFA (ending with .gfa)"),
            cfg.graph_out << value("output file"),
            (option("--max-length") & integer("value", cfg.max_length)) % "length threshold (default: 10Kb)",
            (option("--cov-thr") & number("value", cfg.cov_thr)) % "coverage upper bound (exclusive, default: -1. -- disabled)",
            (option("--coverage") & value("file", cfg.coverage)) % "file with coverage information",
            option("--compact").set(cfg.compact) % "compact the graph after cleaning (default: false)",
            (option("--id-mapping") & value("file", cfg.id_mapping)) % "file with compacted segment id mapping",
            (option("--prefix") & value("vale", cfg.compacted_prefix)) % "prefix used to form compacted segment names",
            option("--drop-sequence").set(cfg.drop_sequence) % "flag to drop sequences even if present in original file (default: false)"
            //option("--use-cov-ratios").set(cfg.use_cov_ratios) % "enable procedures based on unitig coverage ratios (default: false)",
            //(required("-k") & integer("value", cfg.k)) % "k-mer length to use",
    );

    auto result = parse(argc, argv, cli);
    if (!result) {
        std::cerr << "Removing isolated nodes shorter than max-length" << std::endl;
        std::cerr << make_man_page(cli, argv[0]);
        exit(1);
    }

    if (cfg.cov_thr >= 0.) {
        if (cfg.coverage.empty()) {
            std::cerr << "Provide --coverage file\n";
            exit(2);
        }
    }
}

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
    std::cout << "Isolated segments shorter than " << cfg.max_length << "bp will be removed" << std::endl;

    if (cfg.cov_thr >= 0.)
        std::cout << "Only segments with coverage below " << cfg.cov_thr << " will be considered" << std::endl;

    for (gfa::DirectedSegment ds : g.directed_segments()) {
        if (ds.direction == gfa::Direction::REVERSE)
            continue;

        if (g.incoming_link_cnt(ds) + g.outgoing_link_cnt(ds) == 0 &&
                g.segment_length(ds) < cfg.max_length) {
            if (cfg.cov_thr >= 0. && utils::get(*segment_cov_ptr, g.segment_name(ds)) >= cfg.cov_thr) {
                DEBUG("Coverage of segment " << g.str(ds) << " exceeded upper bound");
            } else {
                INFO("Will remove node " << g.str(ds));
                g.DeleteSegment(ds);
                ndel++;
            }
        }
    }

    tooling::OutputGraph(g, cfg, ndel, segment_cov_ptr.get());
    std::cout << "END" << std::endl;
}
