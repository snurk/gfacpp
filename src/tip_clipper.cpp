#include "tooling.hpp"

#include <vector>
#include <set>
#include <cassert>

struct cmd_cfg: public tooling::cmd_cfg_base {
    //coverage ratio threshold
    size_t max_length = 10000;
    std::string read_cnt_file;
    uint32_t max_read_cnt = 0;
};

static void process_cmdline(int argc, char **argv, cmd_cfg &cfg) {
    using namespace clipp;
    //<read counts file> <max read count>
    auto cli = (
            cfg.graph_in << value("input file in GFA (ending with .gfa)"),
            cfg.graph_out << value("output file"),
            (option("--max-length") & integer("value", cfg.max_length)) % "tip length threshold (default: 10Kb)",
            (option("--read-cnt-file") & value("value", cfg.read_cnt_file)) % "file with read counts",
            (option("--max-read-cnt") & integer("value", cfg.max_read_cnt)) % "max read count (default: 0)",
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
        std::cerr << "Removing tips based on length and optionally read counts" << std::endl;
        std::cerr << make_man_page(cli, argv[0]);
        exit(1);
    }

    if (cfg.max_read_cnt > 0 && cfg.read_cnt_file.empty()) {
        std::cerr << "Non-trivial threshold on read counts requires read-cnt-file to be provided" << std::endl;
        exit(2);
    }
}

int main(int argc, char *argv[]) {
    cmd_cfg cfg;
    process_cmdline(argc, argv, cfg);

    std::unique_ptr<utils::SegmentCoverageMap> read_cnt_ptr;
    if (!cfg.read_cnt_file.empty()) {
        std::cout << "Reading read counts from " << cfg.read_cnt_file << std::endl;
        read_cnt_ptr = std::make_unique<utils::SegmentCoverageMap>(utils::ReadCoverage(cfg.read_cnt_file));
    }

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

    std::cout << "Searching for tips" << std::endl;

//if max_length == 0 check returns false
//TODO improve to support multiple outgoing links
//TODO introduce coverage threshold
    auto is_tip = [&] (const gfa::DirectedSegment &v) {
        if (cfg.max_length == 0)
            return false;

        if (g.incoming_link_cnt(v) > 0 || g.outgoing_link_cnt(v) != 1)
            return false;

        auto l = *g.outgoing_begin(v);
        auto n = l.end;
        assert(g.incoming_link_cnt(n) > 0);
        if (g.incoming_link_cnt(n) == 1)
            return false;

        if (!cfg.read_cnt_file.empty()) {
            uint32_t read_cnt = utils::get(*read_cnt_ptr, g.segment_name(v));
            if (read_cnt > cfg.max_read_cnt) {
                DEBUG("Node " << g.str(v) << " consisting of too many reads: " << read_cnt);
                return false;
            }
        }

        return g.segment_length(v) < cfg.max_length + l.start_overlap;
    };

    for (gfa::DirectedSegment ds : g.directed_segments()) {
        DEBUG("Looking at vertex " << g.str(v));
        if (is_tip(ds)) {
            INFO("Found tip " << g.str(ds));
            g.DeleteSegment(ds);
            ndel++;
        }
    }

    tooling::OutputGraph(g, cfg, ndel, segment_cov_ptr.get());
    std::cout << "END" << std::endl;
}