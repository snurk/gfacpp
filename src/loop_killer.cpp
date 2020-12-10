#include "tooling.hpp"

#include <iostream>
#include <cassert>

struct cmd_cfg: public tooling::cmd_cfg_base {
    //coverage ratio threshold
    double max_base_coverage = 0.;
};

static void process_cmdline(int argc, char **argv, cmd_cfg &cfg) {
    using namespace clipp;

    auto cli = (
            cfg.graph_in << value("input file in GFA (ending with .gfa)"),
            cfg.graph_out << value("output file"),
            (required("--coverage") & value("file", cfg.coverage)) % "file with coverage information",
            (required("--max-base-cov") & number("value", cfg.max_base_coverage)) % "maximal coverage of the node with loop",
            option("--compact").set(cfg.compact) % "compact the graph after cleaning (default: false)",
            (option("--id-mapping") & value("file", cfg.id_mapping)) % "file with compacted segment id mapping",
            (option("--prefix") & value("vale", cfg.compacted_prefix)) % "prefix used to form compacted segment names",
            option("--drop-sequence").set(cfg.drop_sequence) % "flag to drop sequences even if present in original file (default: false)"
            //option("--use-cov-ratios").set(cfg.use_cov_ratios) % "enable procedures based on unitig coverage ratios (default: false)",
            //(required("-k") & integer("value", cfg.k)) % "k-mer length to use",
    );


    auto result = parse(argc, argv, cli);
    if (!result) {
        std::cerr << "Loop link will be killed if coverage of the node doesn't exceed 'max_base_coverage' and other links are present" << std::endl;
        std::cerr << make_man_page(cli, argv[0]);
        exit(1);
    }
}

//TODO consider making iterative right here after I can compress and track reads here
//TODO put coverage into GFA (check support in parser, etc)
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
        //TODO add forward-only iterator
        if (v.direction == gfa::Direction::REVERSE)
            continue;
        if (g.outgoing_link_cnt(v) < 2 || g.incoming_link_cnt(v) < 2)
            continue;

        DEBUG("Looking at directed node " << g.str(v));
        for (auto l : g.outgoing_links(v)) {
            assert(l.start == v);
            if (l.end != v)
                continue;
            if (utils::get(segment_cov, g.segment_name(v)) <= cfg.max_base_coverage) {
                std::cout << "Removing loop link from segment " << g.str(v) << '\n';
                g.DeleteLink(l);
                ++l_ndel;
            }
        }
    }

    tooling::OutputGraph(g, cfg, l_ndel, &segment_cov);
    std::cout << "END" << std::endl;
}
