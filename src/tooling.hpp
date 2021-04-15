#pragma once

#include "compact.hpp"
#include "utils.hpp"
#include "clipp.h"

namespace tooling {

struct cmd_cfg_base {
    //input file
    std::string graph_in;

    //output file
    std::string graph_out;

    //optional file with coverage
    std::string coverage;

    //compact the graph after cleaning
    bool compact = false;

    //optional file with id mapping
    std::string id_mapping;

    //prefix used to form compacted segment names
    std::string compacted_prefix = "m_";

    //flag to drop sequences even if present in original file
    bool drop_sequence = false;

    //flag to drop sequences even if present in original file
    bool rename_all = false;

    //DBG vertex size, enables DBG mode of coverage transformation
    int32_t dbg_k = 0;
};

inline
clipp::group BaseCfg(cmd_cfg_base &cfg) {
    using namespace clipp;

    return ( cfg.graph_in << value("input file in GFA (ending with .gfa)"),
            cfg.graph_out << value("output file"),
            (option("--coverage") & value("file", cfg.coverage)) % "file with coverage information",
            option("--compact").set(cfg.compact) % "compact the graph after cleaning (default: false)",
            (option("--id-mapping") & value("file", cfg.id_mapping)) % "file with compacted segment id mapping",
            (option("--prefix") & value("vale", cfg.compacted_prefix)) % "prefix used to form compacted segment names (default: m_, use _ for empty)",
            option("--drop-sequence").set(cfg.drop_sequence) % "flag to drop sequences even if present in original file (default: false)",
            option("--rename-all").set(cfg.rename_all) % "flag to rename all segments (default: false)",
            (option("--dbg-k") & integer("value", cfg.dbg_k)) % "DBG k-mer length to use in coverage transformation (default: 0 -- disabled)"
    ) % "common settings";
}

void OutputGraph(gfa::Graph &g,
                 const cmd_cfg_base &cfg,
                 size_t ndel = size_t(-1),
                 const utils::SegmentCoverageMap *segment_cov_ptr = nullptr) {
    if (ndel != size_t(-1))
        std::cout << "Triggered " << ndel << " times" << std::endl;

    if (ndel > 0) {
        std::cout << "Cleanup" << std::endl;
        g.Cleanup();
    }

    assert(g.CheckNoDeadLinks());

    std::cout << "Writing complete" << std::endl;
    if (ndel > 0 && cfg.compact) {
        gfa::Compactifier compactifier(g, cfg.compacted_prefix, segment_cov_ptr);
        std::cout << "Writing compacted graph to " << cfg.graph_out << std::endl;
        compactifier.Compact(cfg.graph_out, cfg.id_mapping, cfg.drop_sequence);
    } else {
        std::cout << "Writing output to " << cfg.graph_out << std::endl;
        g.write(cfg.graph_out, cfg.drop_sequence);
    }

    std::cout << "Writing complete" << std::endl;
}

}
