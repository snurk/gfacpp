#include "clipp.h"
#include "wrapper.hpp"
#include "utils.hpp"

#include <vector>
#include <set>
#include <map>

static void Go(const gfa::Graph &g, gfa::DirectedSegment v,
                uint32_t curr_depth, const uint32_t max_depth,
                std::map<gfa::DirectedSegment, uint32_t> &considered) {
    //std::cout << "In node " << g.nodes[v_id] << " on depth " << curr_depth << std::endl;
    considered[v] = curr_depth;
    considered[v.Complement()] = curr_depth;
    if (curr_depth == max_depth) {
        //std::cout << "Hit max depth " << max_depth << std::endl;
        return;
    }
    for (const auto& l : g.outgoing_links(v)) {
        //not considered or considered at a higher depth
        if (considered.count(l.end) == 0 || considered[l.end] > curr_depth + 1) {
            Go(g, l.end, curr_depth + 1, max_depth, considered);
        }
    }
    for (const auto& l : g.incoming_links(v)) {
        //not considered or considered at a higher depth
        if (considered.count(l.start) == 0 || considered[l.start] > curr_depth + 1) {
            Go(g, l.start, curr_depth + 1, max_depth, considered);
        }
    }
}

static std::map<gfa::DirectedSegment, uint32_t>
CollectNeighborhood(const gfa::Graph &g,
                    const std::set<std::string> &nodes_of_interest,
                    const uint32_t max_depth) {

    INFO("Searching for neighbourhood");
    std::map<gfa::DirectedSegment, uint32_t> considered;
    std::vector<gfa::DirectedSegment> ids_of_interest;
    ids_of_interest.reserve(nodes_of_interest.size() * 2);

    for (auto ds : g.directed_segments()) {
        if (ds.direction == gfa::Direction::REVERSE) {
            DEBUG("Skipping " << g.str(ds));
            continue;
        }

        if (nodes_of_interest.count(g.segment_name(ds))) {
            ids_of_interest.push_back(ds);
            considered[ds] = 0;
        }
    }

    for (auto ds : ids_of_interest) {
        INFO("Starting from vertex " << g.str(ds));
        Go(g, ds, 0, max_depth, considered);
    }

    //std::cout << "Extracting the subgraph" << std::endl;
    //std::set<gfa::DirectedSegment> reached;
    //for (const auto &k_v: considered)
    //    reached.insert(k_v.first);

    //return reached_segments;
    return considered;
}

struct cmd_cfg {
    //input file
    std::string graph_in;

    //output file
    std::string graph_out;

    //optional file with coverage
    std::string coverage;

    //file with nodes of interest
    std::string nodes;

    //flag to drop sequences even if present in original file
    bool drop_sequence = false;

    //neighborhood radius
    size_t radius = 10;

};

static void process_cmdline(int argc, char **argv, cmd_cfg &cfg) {
  using namespace clipp;

  auto cli = ( cfg.graph_in << value("input file in GFA (ending with .gfa)"),
         cfg.graph_out << value("output file"),
         (required("-n", "--nodes") & value("file", cfg.nodes)) % "file with nodes ids of interest",
         (option("-c", "--coverage") & value("file", cfg.coverage)) % "file with coverage information",
         option("--drop-sequence").set(cfg.drop_sequence) % "flag to drop sequences even if present in original file (default: false)",
         (option("-r", "--radius") & integer("value", cfg.radius)) % "neighborhood radius (default: 10)"
  ) % "algorithm settings";

  auto result = parse(argc, argv, cli);

  if (!result) {
      std::cerr << make_man_page(cli, argv[0]);
      exit(1);
  }
}

int main(int argc, char *argv[]) {
    cmd_cfg cfg;
    process_cmdline(argc, argv, cfg);

    std::unique_ptr<utils::SegmentCoverageMap> segment_cov_ptr;
    if (!cfg.coverage.empty()) {
        INFO("Reading coverage from " << cfg.coverage);
        segment_cov_ptr = std::make_unique<utils::SegmentCoverageMap>(utils::ReadCoverage(cfg.coverage));
    }

    gfa::Graph g;
    INFO("Loading graph from GFA file " << cfg.graph_in);
    g.open(cfg.graph_in);
    INFO("Segment cnt: " << g.segment_cnt() << "; link cnt: " << g.link_cnt());

    //std::set<std::string> neighbourhood;

    //min(unique_left_ovl, unique_right_ovl) & segment_id
    std::set<std::string> nodes_of_interest;
    utils::ReadSet(cfg.nodes, nodes_of_interest);

    auto neighborhood = CollectNeighborhood(g, nodes_of_interest, cfg.radius);

    std::ofstream out(cfg.graph_out);

    for (gfa::DirectedSegment v : g.directed_segments()) {
        DEBUG("Considering vertex " << g.str(v));
        auto seg = g.segment(v);

        //TODO add forward-only iterator
        if (v.direction == gfa::Direction::REVERSE) {
            DEBUG("Skipping " << g.str(v));
            continue;
        }

        if (neighborhood.count(v) == 0) {
            continue;
        }

        if (seg.removed()) {
            DEBUG("Removed segment " << seg.name);
            assert(false);
            continue;
        }

        std::string s = (cfg.drop_sequence || !seg.sequence) ? "" : std::string(seg.sequence);
        out << "S\t" << seg.name <<
            "\t" << (s.empty() ? "*" : s) <<
            "\tLN:i:" << std::to_string(seg.length);

        if (segment_cov_ptr) {
            double cov = double(utils::get(*segment_cov_ptr, g.segment_name(v)));
            //adding Mikko-style output to simplify scripting
            out << "\tRC:i:" << uint64_t(std::round(cov * seg.length));
            out << "\tll:f:" << std::round(cov * 1000) / 1000;
        }
        out << "\n";
    }

    for (auto v : g.directed_segments()) {
        for (auto l : g.outgoing_links(v)) {
            if (!l.IsCanonical())
                continue;

            if (neighborhood.count(l.start) == 0
                    || neighborhood.count(l.end) == 0) {
                continue;
            }

            //TODO support CIGAR?
            out <<"L\t" << g.str(l.start, "\t")
                << "\t" << g.str(l.end, "\t")
                << "\t" << l.overlap() << "M" << "\n";
        }
    }

    INFO("Finished");
}
