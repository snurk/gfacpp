#include "tooling.hpp"

#include <vector>
#include <set>
#include <cassert>

struct cmd_cfg: public tooling::cmd_cfg_base {
    //coverage ratio threshold
    size_t max_length = 10000;
    double cov_thr = 0.;
};

static void process_cmdline(int argc, char **argv, cmd_cfg &cfg) {
    using namespace clipp;

    auto cli = (
            cfg.graph_in << value("input file in GFA (ending with .gfa)"),
            cfg.graph_out << value("output file"),
            (option("--max-length") & integer("value", cfg.max_length)) % "tip length threshold (default: 10Kb)",
            (option("--cov-thr") & number("value", cfg.cov_thr)) % "coverage threshold (default: 0.)",
            (required("--coverage") & value("file", cfg.coverage)) % "file with coverage information",
            option("--compact").set(cfg.compact) % "compact the graph after cleaning (default: false)",
            (option("--id-mapping") & value("file", cfg.id_mapping)) % "file with compacted segment id mapping",
            (option("--prefix") & value("vale", cfg.compacted_prefix)) % "prefix used to form compacted segment names",
            option("--drop-sequence").set(cfg.drop_sequence) % "flag to drop sequences even if present in original file (default: false)"
            //option("--use-cov-ratios").set(cfg.use_cov_ratios) % "enable procedures based on unitig coverage ratios (default: false)",
            //(required("-k") & integer("value", cfg.k)) % "k-mer length to use",
    );


    auto result = parse(argc, argv, cli);
    if (!result) {
        std::cerr << "Removing nodes shorter than max-length with coverage below cov-thr" << std::endl;
        std::cerr << make_man_page(cli, argv[0]);
        exit(1);
    }
}

int main(int argc, char *argv[]) {
    cmd_cfg cfg;
    process_cmdline(argc, argv, cfg);

    std::cout << "Reading coverage from " << cfg.coverage << std::endl;
    const auto segment_cov = utils::ReadCoverage(cfg.coverage);

    gfa::Graph g;
    std::cout << "Loading graph from GFA file " << cfg.graph_in << std::endl;
    g.open(cfg.graph_in);
    std::cout << "Segment cnt: " << g.segment_cnt() << "; link cnt: " << g.link_cnt() << std::endl;

    size_t ndel = 0;

    std::cout << "Nodes with coverage below " << cfg.cov_thr <<
            " no longer than " << cfg.max_length << "bp will be removed" << std::endl;

    for (gfa::DirectedSegment ds : g.directed_segments()) {
        //DEBUG("Processing segment " << g.str(ds) << " of length " << g.segment_length(ds));
        if (ds.direction == gfa::Direction::REVERSE || g.segment_length(ds) > cfg.max_length)
            continue;

        double cov = utils::get(segment_cov, g.segment_name(ds));
        //DEBUG("Segment " << g.str(ds) << " with coverage " << cov);
        if (cov < cfg.cov_thr) {
            INFO("Will remove node " << g.str(ds) << " of length " << g.segment_length(ds) << " with coverage " << cov);
            g.DeleteSegment(ds);
            ndel++;
        }
    }

    tooling::OutputGraph(g, cfg, ndel, &segment_cov);
    std::cout << "END" << std::endl;
}
