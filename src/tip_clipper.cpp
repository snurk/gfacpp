#include "tooling.hpp"

#include <vector>
#include <set>
#include <cassert>

struct cmd_cfg: public tooling::cmd_cfg_base {
    //coverage ratio threshold
    size_t max_length = 10000;
    std::string read_cnt_file;
    uint32_t max_read_cnt = 0;
    double cov_thr = -1.;
};

static void process_cmdline(int argc, char **argv, cmd_cfg &cfg) {
    using namespace clipp;

    auto cli = (tooling::BaseCfg(cfg), (
            (option("--max-length") & integer("value", cfg.max_length)) % "tip length threshold (default: 10Kb)",
            (option("--read-cnt-file") & value("value", cfg.read_cnt_file)) % "file with read counts",
            (option("--max-read-cnt") & integer("value", cfg.max_read_cnt)) % "max read count (default: 0)",
            (option("--cov-thr") & number("value", cfg.cov_thr)) % "coverage upper bound (exclusive, default: -1. -- disabled)"
                //option("--use-cov-ratios").set(cfg.use_cov_ratios) % "enable procedures based on unitig coverage ratios (default: false)",
    ) % "algorithm settings");

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

    std::cout << "Searching for tips with (adjusted for overlap size) length below " << cfg.max_length << std::endl;

    if (cfg.cov_thr >= 0.)
        std::cout << "Only segments with coverage below " << cfg.cov_thr << " will be considered" << std::endl;

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

        if (g.segment_length(v) >= cfg.max_length + l.start_overlap) {
            DEBUG("Length of segment " << g.str(v) << " " <<
                    g.segment_length(v) << "bp (adjusted for overlap size " <<
                    l.start_overlap << "bp) exceeded tip length threshold of " <<
                    cfg.max_length << "bp");
            return false;
        }

        if (!cfg.read_cnt_file.empty()) {
            uint32_t read_cnt = utils::get(*read_cnt_ptr, g.segment_name(v));
            if (read_cnt > cfg.max_read_cnt) {
                DEBUG("Segment " << g.str(v) << " consisting of too many reads: " << read_cnt);
                return false;
            }
        }

        if (cfg.cov_thr >= 0. && utils::get(*segment_cov_ptr, g.segment_name(v)) >= cfg.cov_thr) {
            DEBUG("Coverage of segment " << g.str(v) << " exceeded upper bound");
            return false;
        }

        return true;
    };

    for (gfa::DirectedSegment ds : g.directed_segments()) {
        DEBUG("Looking at node " << g.str(ds));
        if (is_tip(ds)) {
            INFO("Found tip " << g.str(ds));
            INFO("Removing segment " << g.str(ds));
            g.DeleteSegment(ds);
            ndel++;
        }
    }

    tooling::OutputGraph(g, cfg, ndel, segment_cov_ptr.get());
    std::cout << "END" << std::endl;
}
