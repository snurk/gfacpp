#include "tooling.hpp"
#include "wrapper.hpp"
#include "utils.hpp"

#include <vector>
#include <functional>
#include <memory>
#include <algorithm>
#include <set>
#include <cassert>
#include <iostream>

struct cmd_cfg: public tooling::cmd_cfg_base {
    //check that length contributed by the node is below
    size_t max_length = std::numeric_limits<size_t>::max();

    //check that length difference between the paths is below
    size_t max_diff = std::numeric_limits<size_t>::max();

    ////check that relative length difference between the paths is below frac
    //double max_rel_diff = std::numeric_limits<double>::max();

    //check that alt is no more than value shorter (set to 0 to force only longer alternatives)
    size_t max_shortening = std::numeric_limits<size_t>::max();

    ////check that alternative is no more than frac longer (set to 0. to force only shorter alternatives)
    //double max_rel_alt_increase = std::numeric_limits<double>::max();

    //check that minimal overlap on alternative path is above or equal to value
    int32_t min_alt_overlap = 0;

    //check that first & last node are below this coverage (coverage must be provided)
    double max_unique_cov = std::numeric_limits<double>::max();

    //check that [node coverage / min coverage on alternative path] is below (coverage must be provided)
    //if set to << 1. will only remove low coverage variants
    //if set to > 1. ensures that the coverage of the removed node is no more than X times higher than of the retained one
    double max_coverage_ratio = std::numeric_limits<double>::max();

    //flag to use coverage values instead of overlap sizes to prioritize segments
    bool use_coverage = false;
};

static void process_cmdline(int argc, char **argv, cmd_cfg &cfg) {
  using namespace clipp;

  auto cli = (tooling::BaseCfg(cfg), (
      (option("-l", "--max-length") & integer("length", cfg.max_length)) % "check that length contributed by the node is below <length>",
      (option("-d", "--max-diff") & integer("value", cfg.max_diff)) % "check that length difference between paths is below <value>",
      //(option("--max-rel-diff") & number("frac", cfg.max_rel_diff)) % "check that relative length difference between paths is below <frac>",
      (option("--max-shortening") & integer("value", cfg.max_shortening)) % "check that base is no more than <value>bp longer (setting to 0 forces only longer alternatives)",
      //(option("--max-rel-shortening") & number("frac", cfg.max_rel_alt_increase)) % "check that base is no more than <frac> longer (setting to 0. forces only longer alternatives)",
      (option("--min-alt-ovl") & integer("value", cfg.min_alt_overlap)) % "check that minimal overlap on alternative path is stronger than on base path or is >= <value>",
      (option("--max-unique-cov") & number("value", cfg.max_unique_cov)) % "check that first & last node are below this coverage (coverage must be available)",
      (option("--max-cov-ratio") & number("value", cfg.max_coverage_ratio)) % "check that [node coverage / min coverage on alternative path] is below <value> (coverage must be provided)."
                                                                          "If set to << 1. (e.g. to 0.2) will only remove low coverage (5x less covered) variants,"
                                                                          "if set to > 1. ensures that the coverage of the removed node is no more than <value> times higher than of the alternative",
      option("--use-coverage").set(cfg.use_coverage) % "use coverage instead of overlap sizes (default: false)"
  ) % "algorithm settings");

  auto result = parse(argc, argv, cli);
  if (cfg.use_coverage && cfg.coverage.empty()) {
      std::cerr << "Option to use coverage values was enabled, but coverage file wasn't provided" << std::endl;
      exit(2);
  }

  if (!result) {
      std::cerr << make_man_page(cli, argv[0]);
      exit(1);
  }

  if (cfg.max_unique_cov != std::numeric_limits<double>::max() ||
      cfg.max_coverage_ratio != std::numeric_limits<double>::max()) {
      if (cfg.coverage.empty()) {
          std::cerr << "Provide --coverage file\n";
          exit(2);
      }
  }
}

//TODO maybe consider all paths rather than unambiguous
//Outputting unambiguous path from w to v or empty path
inline gfa::Path UnambiguousBackwardPath(const gfa::Graph &g, gfa::DirectedSegment w, gfa::DirectedSegment v) {
    assert(w != v);
    std::vector<gfa::LinkInfo> rev_links;
    std::set<gfa::LinkInfo> used_links;
    while (g.unique_incoming(w) && w != v) {
        auto l = *g.incoming_begin(w);
        if (used_links.count(l)) {
            DEBUG("Loop detected");
            return gfa::Path();
        }
        rev_links.push_back(l);
        used_links.insert(l);
        w = l.start;
    }
    if (w == v) {
        return gfa::Path(std::vector<gfa::LinkInfo>(rev_links.rbegin(), rev_links.rend()));
    }
    //Allow unambiguous path to stop one step away from v
    //for (const auto &l : g.incoming_links(w)) {
    //    if (l.start == v) {
    //        rev_links.push_back(l);
    //        return gfa::Path(std::vector<gfa::LinkInfo>(rev_links.rbegin(), rev_links.rend()));
    //    }
    //}
    return gfa::Path();
}

inline bool CheckNotInPath(const gfa::Path &p, gfa::DirectedSegment n) {
    for (const auto &v : p.segments) {
        if (v == n || v == n.Complement()) {
            return false;
        }
    }
    return true;
}

typedef std::function<bool (const gfa::Path &base, const gfa::Path &alt)> BulgeCheckF;

inline bool FormsSimpleBulge(const gfa::Graph &g, gfa::DirectedSegment n,
                             size_t max_length,
                             const BulgeCheckF &check_f,
                             std::set<gfa::SegmentId> &protected_segments) {
    assert(g.unique_incoming(n) && g.unique_outgoing(n));
    DEBUG("Considering node " << g.str(n));
    gfa::Path p({*g.incoming_begin(n), *g.outgoing_begin(n)});
    const size_t total_len = g.total_length(p);
    gfa::DirectedSegment v = p.segments.front();
    gfa::DirectedSegment w = p.segments.back();

    if (v == n || v == n.Complement() || w == n || w == n.Complement()) {
        DEBUG("Unambiguously links to self or complement");
        return false;
    }

    //todo simplify? maybe use int64_t instead of size_t?
    if (total_len > g.segment_length(v) + g.segment_length(w) &&
        total_len - g.segment_length(v) - g.segment_length(w) > max_length) {
        DEBUG("Base path too long (" << (total_len - g.segment_length(v) - g.segment_length(w)) << "bp of 'internal' bases)");
        return false;
    }

    for (auto l : g.incoming_links(w)) {
        auto w1 = l.start;
        if (w1 == n || w1 == v || w1 == n.Complement())
            continue;
        auto alt_p = UnambiguousBackwardPath(g, w1, v);
        if (!CheckNotInPath(alt_p, n)) {
            DEBUG("Alternative path hit 'base' node");
            continue;
        }
        if (!alt_p.empty()) {
            alt_p.Extend(l);
            DEBUG("Found suitable 'backward' path to check " << g.str(alt_p));

            if (check_f(p, alt_p)) {
                DEBUG("Check successful, removing node " << g.str(n));
                for (const auto &a : alt_p.segments) {
                    DEBUG("Marking alternative path node " << g.str(a) << " as protected");
                    protected_segments.insert(a.segment_id);
                }
                return true;
            }
        }
    }
    return false;
}

//NB apply only to the graph after weak link removal rounds!
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

    //std::set<std::string> neighbourhood;
    std::set<gfa::SegmentId> protected_segments;

    auto cov_f = [&](gfa::DirectedSegment v) {
        assert(segment_cov_ptr);
        return double(utils::get(*segment_cov_ptr, g.segment_name(v)));
    };

    //min(unique_left_ovl, unique_right_ovl) & segment_id
    std::vector<std::pair<double, gfa::SegmentId>> segments_of_interest;
    segments_of_interest.reserve(g.segment_cnt());
    for (gfa::SegmentId s = 0; s < g.segment_cnt(); ++s) {
        gfa::DirectedSegment v(s, gfa::Direction::FORWARD);
        if (g.unique_outgoing(v) && g.unique_incoming(v)) {
            double weight;
            if (cfg.use_coverage) {
                weight = cov_f(v);
            } else {
                int32_t min_ovl = std::min((*g.outgoing_begin(v)).start_overlap, (*g.incoming_begin(v)).end_overlap);
                weight = double(min_ovl);
            }
            segments_of_interest.push_back(std::make_pair(weight, s));
        }
    }

    //sort in order of increased min overlaps or coverage
    std::sort(segments_of_interest.begin(), segments_of_interest.end());

    auto inner_cov_f = [&](const gfa::Path &p) {
        assert(p.segment_cnt() >= 3);
        double min_cov = std::numeric_limits<double>::max();
        for (size_t i = 1; i < p.segment_cnt() - 1; ++i) {
            min_cov = std::min(min_cov, cov_f(p.segments[i]));
        }
        return min_cov;
    };

    auto bulge_check_f = [&](const gfa::Path &base, const gfa::Path &alt) {
        assert(base.segment_cnt() == 3 && alt.segment_cnt() >= 3);
        auto diff = utils::abs_diff(g.total_length(alt), g.total_length(base));

        if (diff > cfg.max_diff) {
            DEBUG(diff << "bp diff in length between 'alt' and 'base' paths exceeded max_diff=" << cfg.max_diff);
            return false;
        }

        if (g.total_length(base) > g.total_length(alt) && diff > cfg.max_shortening) {
            DEBUG("'Alt' length was " << diff << "bp shorter than 'base', which exceeded max_shortening threshold=" << cfg.max_shortening);
            return false;
        }

        //assert(alt.min_overlap() > 0 && base.min_overlap() > 0);
        if (alt.min_overlap() < base.min_overlap() && alt.min_overlap() < cfg.min_alt_overlap) {
            DEBUG("Minimal overlap along the 'alt' path " << g.str(alt)
                << " was shorter than for the 'base' path " << g.str(base)
                << " and shorter than " << cfg.min_alt_overlap << " threshold");
            return false;
        }

        if (segment_cov_ptr) {
            //todo can be optimized if the check is actually disabled
            auto v = base.segments.front();
            auto w = base.segments.back();
            if (cov_f(v) > cfg.max_unique_cov + 1e-5 || cov_f(w) > cfg.max_unique_cov + 1e-5) {
                DEBUG("Coverage on one of the sides exceeded 'uniqueness' threshold=" << cfg.max_unique_cov);
                return false;
            }

            DEBUG("Cov base " << inner_cov_f(base) << " cov alt " << inner_cov_f(alt));
            if (inner_cov_f(alt) < 1e-5 || (inner_cov_f(base) / inner_cov_f(alt)) > cfg.max_coverage_ratio) {
                DEBUG("Ratio between estimated coverage of the node and alternative path exceeded specified ratio threshold=" << cfg.max_coverage_ratio);
                return false;
            }
        }
        return true;
    };

    size_t ndel = 0;
    for (auto ovl_s : segments_of_interest) {
        gfa::SegmentId seg_id = ovl_s.second;
        std::cout << "Considering segment " << g.str(seg_id) << ". Min overlap " << ovl_s.first << std::endl;

        if (protected_segments.count(seg_id)) {
            DEBUG("Segment " << g.str(seg_id) << " is protected");
            continue;
        }
        if (FormsSimpleBulge(g, gfa::DirectedSegment::Forward(seg_id),
                    cfg.max_length, bulge_check_f, protected_segments)
            || FormsSimpleBulge(g, gfa::DirectedSegment::Reverse(seg_id),
                    cfg.max_length, bulge_check_f, protected_segments)) {
            std::cout << "Removing simple bulge " << g.str(seg_id) << std::endl;
            g.DeleteSegment(seg_id);
            ++ndel;
        }
    }

    tooling::OutputGraph(g, cfg, ndel, segment_cov_ptr.get());

    std::cout << "END" << std::endl;
}
