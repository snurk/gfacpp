#include "wrapper.hpp"
#include "utils.hpp"

#include <vector>
#include <functional>
#include <algorithm>
#include <set>
#include <cassert>
#include <iostream>

//static void process_cmdline(int argc, char **argv, gcfg &cfg) {
//  using namespace clipp;
//
//  auto cli = (
//      cfg.graph << value("graph. In GFA (ending with .gfa) or prefix to SPAdes graph pack"),
//      cfg.outfile << value("output prefix"),
//      option("--gfa").set(cfg.save_gfa) % "produce GFA output (default: true)",
//      option("--spades-gp").set(cfg.save_gp) % "produce output graph pack in SPAdes internal format (default: false). "
//                                                      "Recommended if bulges are removed to improve further read mapping. "
//                                                      "In case GFA output is required with graph pack specify '--gfa'",
//      option("--use-cov-ratios").set(cfg.use_cov_ratios) % "enable procedures based on unitig coverage ratios (default: false)",
//      (required("-k") & integer("value", cfg.k)) % "k-mer length to use",
//      (required("--read-length") & integer("value", cfg.RL)) % "read length",
//      (option("-c", "--coverage") & value("coverage", cfg.bin_cov_str)) % "estimated average (k+1-mer) bin coverage (default: 0.) "
//                                                                          "or 'auto' (works only with '-d/--dead-ends' provided)",
//      (option("-t", "--threads") & integer("value", cfg.nthreads)) % "# of threads to use (default: max_threads / 2)",
//      (option("-p", "--profile") & value("file", cfg.edge_profile_fn)) % "file with edge coverage profiles across multiple samples",
//      (option("-s", "--stop-codons") & value("file", cfg.stop_codons_fn)) % "file stop codon positions",
//      (option("-d", "--dead-ends") & value("file", cfg.deadends_fn)) % "while processing a subgraph -- file listing edges which are dead-ends in the original graph",
//      (option("--tmpdir") & value("dir", cfg.tmpdir)) % "scratch directory to use (default: <output prefix>.tmp)"
//  );
//
//  auto result = parse(argc, argv, cli);
//  if (!result) {
//      std::cout << make_man_page(cli, argv[0]);
//      exit(1);
//  }
//}

//TODO maybe consider all paths rather than unambiguous
inline gfa::Path UnambiguousBackwardPath(const gfa::Graph &g, gfa::DirectedSegment w, gfa::DirectedSegment v) {
    assert(w != v);
    std::vector<gfa::LinkInfo> rev_links;
    while (g.unique_incoming(w) && w != v) {
        auto l = *g.incoming_begin(w);
        rev_links.push_back(l);
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

template<typename T>
T abs_diff(T a, T b) {
    return (a > b) ? a-b : b-a;
}

inline bool CheckNotInPath(const gfa::Path &p, gfa::DirectedSegment n) {
    for (const auto &v : p.segments) {
        if (v == n || v == n.Complement()) {
            return false;
        }
    }
    return true;
}

inline bool FormsSimpleBubble(const gfa::Graph &g, gfa::DirectedSegment n,
                              uint32_t max_length, uint32_t max_diff,
                              uint32_t min_alt_overlap,
                              std::set<uint32_t> &protected_segments) {
    assert(g.unique_incoming(n) && g.unique_outgoing(n));
    DEBUG("Considering node " << g.str(n));
    gfa::Path p({*g.incoming_begin(n), *g.outgoing_begin(n)});
    const int32_t min_ovl = p.min_overlap();
    const size_t total_len = g.total_length(p);
    gfa::DirectedSegment v = p.segments.front();
    gfa::DirectedSegment w = p.segments.back();

    if (v == n || v == n.Complement() || w == n || w == n.Complement()) {
        DEBUG("Unambiguously links to self or complement");
        return false;
    }

    if (total_len > max_length
            + g.segment_length(v) + g.segment_length(w)) {
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
            DEBUG("Found suitable 'backward' path " << g.str(alt_p));
            if (abs_diff(g.total_length(alt_p), total_len) > max_diff) {
                DEBUG("Difference in length between 'alt' and 'base' paths exceeded max_diff=" << max_diff);
                continue;
            }
            if (alt_p.min_overlap() < min_ovl && alt_p.min_overlap() < min_alt_overlap) {
                DEBUG("Minimal overlap along the 'alt' path " << g.str(alt_p)
                    << " was shorter than for the 'base' path via " << g.str(n)
                    << " and shorted than " << min_alt_overlap << " threshold");
                continue;
            }
            for (const auto &a : alt_p.segments) {
                DEBUG("Marking alternative path node " << g.str(a) << " as protected");
                protected_segments.insert(a.segment_id);
            }
            return true;
        }
    }
    return false;
}

//NB apply only to the graph after weak link removal rounds!
int main(int argc, char *argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <gfa input> <gfa output> <max length> <max diff> <min alt overlap>" << std::endl;
        std::cerr << "Removing simple bubble of simple node against 'unambiguous backward' alternative path" << std::endl;
        std::cerr << "Removal only occurs if minimal overlap on the alternative path is bigger than min_alt_overlap" << std::endl;
        exit(239);
    }

    const std::string in_fn(argv[1]);
    const std::string out_fn(argv[2]);
    const uint32_t max_length = std::stoi(argv[3]);
    const uint32_t max_diff = std::stoi(argv[4]);
    const uint32_t min_alt_overlap = std::stoi(argv[5]);

    gfa::Graph g;
    std::cout << "Loading graph from GFA file " << in_fn << '\n';
    g.open(in_fn);
    std::cout << "Segment cnt: " << g.segment_cnt() << "; link cnt: " << g.link_cnt() << '\n';

    //std::set<std::string> neighbourhood;
    std::set<uint32_t> protected_segments;

    //min(unique_left_ovl, unique_right_ovl) & segment_id
    std::vector<std::pair<int32_t, uint32_t>> segments_of_interest;
    segments_of_interest.reserve(g.segment_cnt());
    for (uint32_t s = 0; s < g.segment_cnt(); ++s) {
        gfa::DirectedSegment v(s, gfa::Direction::FORWARD);
        if (g.unique_outgoing(v) && g.unique_incoming(v)) {
            int32_t min_ovl = std::min((*g.outgoing_begin(v)).start_overlap, (*g.incoming_begin(v)).end_overlap);
            segments_of_interest.push_back(std::make_pair(min_ovl, s));
        }
    }

    //sort in descending order
    std::sort(segments_of_interest.begin(), segments_of_interest.end());

    size_t ndel = 0;
    for (auto ovl_s : segments_of_interest) {
        uint32_t seg_id = ovl_s.second;
        std::cout << "Considering segment " << g.str(seg_id) << ". Min overlap " << ovl_s.first << '\n';

        if (protected_segments.count(seg_id)) {
            DEBUG("Segment " << g.str(seg_id) << " is protected");
            continue;
        }
        if (FormsSimpleBubble(g, gfa::DirectedSegment::Forward(seg_id),
                    max_length, max_diff, min_alt_overlap, protected_segments)
            || FormsSimpleBubble(g, gfa::DirectedSegment::Reverse(seg_id),
                max_length, max_diff, min_alt_overlap, protected_segments)) {
            std::cout << "Removing simple bulge " << g.str(seg_id) << '\n';
            g.DeleteSegment(seg_id);
            ++ndel;
        }
    }

    std::cout << "Writing output to " << out_fn << std::endl;
    g.write(out_fn);
    std::cout << "END" << std::endl;
}


