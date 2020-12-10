#include "tooling.hpp"

#include <vector>
#include <functional>
#include <memory>
#include <algorithm>
#include <set>
#include <cassert>
#include <iostream>

//NB. By default nothing is unique and everything is reliable!
struct cmd_cfg: public tooling::cmd_cfg_base {
    size_t unique_len = std::numeric_limits<size_t>::max();

    //FIXME initialization!
    double max_unique_cov = -1.;

    size_t reliable_len = 0;

    double reliable_cov = -1.;

    int32_t reliable_ovl = 0;

    bool require_both_sides = false;
};

static void process_cmdline(int argc, char **argv, cmd_cfg &cfg) {
  using namespace clipp;

  auto cli = (
      cfg.graph_in << value("input file in GFA (ending with .gfa)"),
      cfg.graph_out << value("output file"),
      (option("--coverage") & value("file", cfg.coverage)) % "file with coverage information",
      option("--compact").set(cfg.compact) % "compact the graph after cleaning (default: false)",
      (option("--id-mapping") & value("file", cfg.id_mapping)) % "file with compacted segment id mapping",
      (option("--prefix") & value("vale", cfg.compacted_prefix)) % "prefix used to form compacted segment names",
      option("--drop-sequence").set(cfg.drop_sequence) % "flag to drop sequences even if present in original file (default: false)",
      (option("--unique-len") & integer("length", cfg.unique_len)) % "longer nodes are considered unique",
      (option("--max-unique-cov") & number("value", cfg.max_unique_cov)) % "node below this coverage is likely unique (coverage must be available)",
      (option("--reliable-cov") & number("value", cfg.reliable_cov)) % "node above this coverage is more likely to be correct (coverage must be available)",
      (option("--reliable-len") & integer("value", cfg.reliable_len)) % "node above this length is more likely to be correct",
      (option("--reliable-ovl") & integer("value", cfg.reliable_ovl)) % "only overlaps exceeding this size are considered reliable",
      option("--both-sides").set(cfg.require_both_sides) % "only remove links that don't look genomic from both sides (default: false)"
      //(required("-k") & integer("value", cfg.k)) % "k-mer length to use",
  );

  auto result = parse(argc, argv, cli);
  if (!result) {
      std::cerr << "Removing links that are unlikely to be part of 'genomic' traversal" << std::endl;
      std::cout << make_man_page(cli, argv[0]);
      exit(1);
  }

  if (cfg.max_unique_cov > -1. || cfg.reliable_cov > -1.) {
      if (cfg.coverage.empty()) {
          std::cerr << "Provide --coverage file" << std::endl;
          exit(2);
      }
  }
}

typedef std::function<bool (gfa::SegmentId)> UniquenessF;

//If a unique node has two unambiguously incoming then those are marked suspicious (the ones shorter than nonsuspicious_length)
//TODO maybe check for reliable coverage
inline std::set<gfa::SegmentId> FindSuspicious(const gfa::Graph &g,
                                          UniquenessF uniqueness_f,
                                          std::set<gfa::SegmentId> &suspected_false) {
   std::set<gfa::SegmentId> suspected_repeats;
    for (auto w: g.directed_segments()) {
        if (!uniqueness_f(w.segment_id))
            continue;

        std::vector<gfa::DirectedSegment> unambiguously_incoming;
        for (const auto &l: g.incoming_links(w))
            if (g.unique_outgoing(l.start))
                unambiguously_incoming.push_back(l.start);

        if (unambiguously_incoming.size() > 1) {
            for (auto v: unambiguously_incoming) {
                DEBUG("Segment " << g.str(v.segment_id) << " is suspected to be false");
                suspected_false.insert(v.segment_id);
            }
            DEBUG("Segment " << g.str(w.segment_id) << " is suspected to be repeat");
            suspected_repeats.insert(w.segment_id);
        }
    }
    return suspected_repeats;
}

inline std::set<gfa::SegmentId> FindDeadends(const gfa::Graph &g) {
    std::set<gfa::SegmentId> answer;
    for (auto v: g.directed_segments()) {
        if (g.no_outgoing(v)) {
            answer.insert(v.segment_id);
        }
    }
    return answer;
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

    auto initial_deadends = FindDeadends(g);

    auto cov_f = [&](gfa::SegmentId s) {
        assert(segment_cov_ptr);
        return double(utils::get(*segment_cov_ptr, g.segment_name(s)));
    };

    //Segments for which coverage is not enough to be considered unique
    std::set<gfa::SegmentId> suspected_repeats;
    //Segments for which coverage is not enough to be considered reliable
    std::set<gfa::SegmentId> suspected_false;
    auto uniqueness_f = [&](gfa::SegmentId s) {
        if (g.segment_length(s) > cfg.unique_len)
            return true;
        if (segment_cov_ptr && !suspected_repeats.count(s))
            if (cov_f(s) < cfg.max_unique_cov + 1e-5)
                return true;
        return false;
    };

    suspected_repeats = FindSuspicious(g, uniqueness_f, suspected_false);
    //NB. From here on uniqueness_f starts to check suspected repeats
    //NB. longer repeat nodes with 'single' coverage can be result of heterozygosity loss

    auto check_reliable_ext = [&](gfa::LinkInfo l) {
        auto w = l.end;
        if (l.overlap() < cfg.reliable_ovl)
            return false;
        if (g.segment_length(w) >= cfg.reliable_len) {
            return true;
        }
        if (segment_cov_ptr && !suspected_false.count(w.segment_id))
            if (cov_f(w.segment_id) > cfg.reliable_cov - 1e-5)
                return true;
        return false;
    };

    auto has_nongenomic_start = [&](gfa::LinkInfo l) {
        //DEBUG("Checking start of link " << g.str(l) << " for not being genomic");
        auto v = l.start;
        auto w = l.end;

        if (g.unique_incoming(w)) {
            //DEBUG("Link end has single incoming");
            return false;
        }
        DEBUG("Checking start of link " << g.str(l) << " for not being genomic");

        if (!uniqueness_f(v.segment_id)) {
            DEBUG("Link start doesn't seem unique");
            return false;
        }

        //TODO extra point if w is also single copy?

        for (const auto &l1: g.outgoing_links(v)) {
            if (l1.end == w) {
                assert(l == l1);
                continue;
            }
            if (g.unique_incoming(l1.end)) {
                if (check_reliable_ext(l1)) {
                    //unambiguous for some reliable extension
                    DEBUG("Reliable unambiguous extension " << g.str(l1) << " found");
                    return true;
                } else {
                    DEBUG("Unambiguous extension " << g.str(l1) << " was unreliable");
                }
            }
        }
        return false;
    };

    size_t l_ndel = 0;
    for (auto v: g.directed_segments()) {
        DEBUG("Looking at directed node " << g.str(v));
        for (const auto &l: g.outgoing_links(v)) {
            if (has_nongenomic_start(l)) {
                DEBUG("Start of the link " << g.str(l) << " doesn't look genomic");
                if (!cfg.require_both_sides || has_nongenomic_start(l.Complement())) {
                    INFO("Removing link " << g.str(l));
                    g.DeleteLink(l);
                    ++l_ndel;
                } else {
                    DEBUG("End of the link " << g.str(l) << " looked genomic");
                }
            }
        }
    }

    //    if (protected_segments.count(seg_id)) {
    //        DEBUG("Segment " << g.str(seg_id) << " is protected");
    //        continue;
    //    }
    //    if (FormsSimpleBulge(g, gfa::DirectedSegment::Forward(seg_id),
    //                max_length, bulge_check_f, protected_segments)
    //        || FormsSimpleBulge(g, gfa::DirectedSegment::Reverse(seg_id),
    //            max_length, bulge_check_f, protected_segments)) {
    //        std::cout << "Removing simple bulge " << g.str(seg_id) << std::endl;
    //        g.DeleteSegment(seg_id);
    //        ++ndel;
    //    }
    //}

    tooling::OutputGraph(g, cfg, l_ndel, segment_cov_ptr.get());

    for (auto s_id: FindDeadends(g)) {
        if (initial_deadends.count(s_id) == 0) {
            WARN("New deadend was formed! Node: " << g.str(s_id));
        }
    }
    std::cout << "END" << std::endl;
}
