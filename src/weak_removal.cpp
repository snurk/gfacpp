#include "tooling.hpp"

#include <iostream>
#include <cassert>
#include <vector>
#include <set>

struct cmd_cfg: public tooling::cmd_cfg_base {
    //coverage ratio threshold
    int32_t min_overlap = 0;
};

static void process_cmdline(int argc, char **argv, cmd_cfg &cfg) {
    using namespace clipp;

    auto cli = (
            cfg.graph_in << value("input file in GFA (ending with .gfa)"),
            cfg.graph_out << value("output file"),
            (option("--coverage") & value("file", cfg.coverage)) % "file with coverage information",
            (required("--min-overlap") & integer("value", cfg.min_overlap)) % "overlap size threshold (default: 0)",
            option("--compact").set(cfg.compact) % "compact the graph after cleaning (default: false)",
            (option("--id-mapping") & value("file", cfg.id_mapping)) % "file with compacted segment id mapping",
            (option("--prefix") & value("vale", cfg.compacted_prefix)) % "prefix used to form compacted segment names",
            option("--drop-sequence").set(cfg.drop_sequence) % "flag to drop sequences even if present in original file (default: false)"
            //option("--use-cov-ratios").set(cfg.use_cov_ratios) % "enable procedures based on unitig coverage ratios (default: false)",
            //(required("-k") & integer("value", cfg.k)) % "k-mer length to use",
    );

    auto result = parse(argc, argv, cli);
    if (!result) {
        std::cerr << "Removing overlaps weaker than the provided threshold" << std::endl;
        std::cerr << make_man_page(cli, argv[0]);
        exit(1);
    }
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

            std::cout << "Removing link between " <<
                g.segment(l.start.segment_id).name << " (" << gfa::PrintDirection(l.start.direction) << ") and " <<
                g.segment(l.end.segment_id).name << " (" << gfa::PrintDirection(l.end.direction) << "). Overlaps " <<
                l.start_overlap << " and " << l.end_overlap << std::endl;
            //TODO optimize?
            g.DeleteLink(l);
            ndel++;
        }
    }

    tooling::OutputGraph(g, cfg, ndel, segment_cov_ptr.get());
    std::cout << "END" << std::endl;
}
