#pragma once

#include "wrapper.hpp"
#include "utils.hpp"

#include <iostream>
#include <string>
#include <map>
#include <set>
#include <functional>
#include <cmath>

namespace gfa {

//TODO move to cpp
class Compactifier {
    const Graph &g_;
    const std::string name_prefix_;
    std::function<double (const std::string&)> coverage_f_;

    LinkInfo NonbranchingExtension(DirectedSegment v) const {
        if (g_.unique_outgoing(v)) {
            auto l = *g_.outgoing_begin(v);
            if (g_.unique_incoming(l.end)) {
                DEBUG("Non-branching extension link " << g_.str(l)
                    << " found for segment " << g_.str(v));
                return l;
            }
        }
        return LinkInfo();
    }

    Path NonbranchingForward(DirectedSegment v_init) const {
        DEBUG("Searching for non-branching path forward from " << g_.str(v_init));
        std::set<SegmentId> used;
        used.insert(v_init.segment_id);
        Path nb_path(v_init);
        LinkInfo l;
        auto v = v_init;
        while ((l = NonbranchingExtension(v)) != LinkInfo()) {
            v = l.end;
            if (used.count(v.segment_id) == 0) {
                DEBUG("Extending path by " << g_.str(l));
                nb_path.Extend(l);
            } else {
                DEBUG("Cycle detected, couldn't extend beyond " << g_.str(l.start));
                break;
            }
        }
        return nb_path;
    }

    Path NonbranchingPath(DirectedSegment v_init) const {
        //TODO optimize by picking appropriate starts?
        DEBUG("Searching around " << g_.str(v_init));
        DEBUG("First searching backward");
        auto nb_path = NonbranchingForward(v_init.Complement()).Complement();
        assert(nb_path.segments.back() == v_init);
        if (g_.unique_outgoing(v_init) &&
                (*g_.outgoing_begin(v_init)).end == nb_path.segments.front()) {
            DEBUG("Loop detected, not traversing twice");
            return nb_path;
        }

        DEBUG("Then searching forward");
        for (auto l : NonbranchingForward(v_init).links) {
            nb_path.Extend(l);
        }
        assert(!nb_path.empty());
        return nb_path;
    }

    std::tuple<std::string, std::size_t, double> CompactedSequence(const Path &p,
                                                                   bool drop_sequence) const {
        assert(!p.segments.empty());
        const auto first_seg = g_.segment(p.segments.front());

        std::string result = (drop_sequence || !first_seg.sequence) ? "" : std::string(first_seg.sequence);
        std::size_t total_len = first_seg.length;

        //TODO think more about coverage averaging computation
        std::size_t len_sum = total_len;
        double coverage = 0.;
        if (coverage_f_) {
            coverage += coverage_f_(first_seg.name) * first_seg.length;
        }

        for (auto l : p.links) {
            const auto seg_info = g_.segment(l.end);
            assert(l.end_overlap >= 0);
            if (uint32_t(l.end_overlap) >= seg_info.length) {
                WARN("Overlap is longer than (or equal to) segment");
            }
            auto trim = std::min(seg_info.length, uint32_t(l.end_overlap));
            total_len += seg_info.length - trim;
            len_sum += seg_info.length;
            if (coverage_f_) {
                coverage += coverage_f_(seg_info.name) * seg_info.length;
            }
            if (!result.empty()) {
                std::string trimmed(seg_info.sequence + trim);
                assert(trimmed.size() == seg_info.length - trim);
                result += trimmed;
            }
        }
        return std::make_tuple(result, total_len, coverage / len_sum);
    }

    std::string FormName(std::size_t compact_cnt) const {
        return name_prefix_ + std::to_string(compact_cnt);
    }

public:
    Compactifier(const Graph &g,
                 std::string name_prefix = "m_",
                 const utils::SegmentCoverageMap *segment_cov_ptr = nullptr):
        g_(g), name_prefix_(std::move(name_prefix)) {
        if (segment_cov_ptr) {
            coverage_f_ = [&](const std::string &s_name) {
                assert(segment_cov_ptr);
                return double(utils::get(*segment_cov_ptr, s_name));
            };
        }
    }

    void Compact(const std::string &out_fn,
                 const std::string &mapping_fn = "",
                 bool drop_sequence = false) const {

        //Do not forget to check both link and complement
        std::set<LinkInfo> inner_links;
        //original segment name to compacted count and orientation 'match'
        std::map<SegmentId, std::pair<std::string, bool>> orig2new;

        std::ofstream out(out_fn);
        std::ofstream mapping_out;
        if (!mapping_fn.empty()) {
            mapping_out.open(mapping_fn, std::ios_base::app);
        }

        out << "H\tVN:Z:1.0" << "\n";
        //S       utg000026l      *       LN:i:18541      RC:i:166869
        //L       m113_3  +       utg511904l      -       9240M

        {
            std::set<SegmentId> used_segments;
            size_t compact_cnt = 0;
            for (gfa::DirectedSegment v : g_.directed_segments()) {
                //TODO add forward-only iterator
                if (v.direction == gfa::Direction::REVERSE || used_segments.count(v.segment_id))
                    continue;

                if (g_.segment(v).removed()) {
                    WARN("Graph had removed segments");
                    continue;
                }

                auto nb_path = NonbranchingPath(v);
                for (auto ds: nb_path.segments) {
                    used_segments.insert(ds.segment_id);
                    assert(!g_.segment(ds).removed());
                }

                for (auto l: nb_path.links) {
                    inner_links.insert(l);
                }

                auto start = nb_path.segments.front();
                auto end = nb_path.segments.back();

                assert(start != end || nb_path.links.empty());

                std::string name;
                if (nb_path.links.empty()) {
                    //keeping the name if trivial path
                    name = g_.segment_name(start);
                } else {
                    name = FormName(compact_cnt++);
                    INFO("Compacting path " << g_.str(nb_path) << " into " << name);
                    if (!mapping_fn.empty()) {
                        mapping_out << name << " " << g_.str(nb_path, ",") << "\n";
                    }
                }

                orig2new[start.segment_id] =
                        std::make_pair(name, start.direction == Direction::FORWARD);

                if (!nb_path.links.empty()) {
                    orig2new[end.segment_id] =
                            std::make_pair(name, end.direction == Direction::FORWARD);
                }

                //compacted seq/len/cov
                std::string cs;
                std::size_t cl;
                double cc;
                std::tie(cs, cl, cc) = CompactedSequence(nb_path, drop_sequence);

                out << "S\t" << name <<
                    "\t" << (cs.empty() ? "*" : cs) <<
                    "\tLN:i:" << std::to_string(cl);

                if (coverage_f_) {
                    out << "\tRC:i:" << uint64_t(std::round(cc * cl));
                }
                out << "\n";
            }
        }

        auto get_compacted = [&] (DirectedSegment v) {
            auto new_id_o = utils::get(orig2new, v.segment_id);
            auto d = v.direction;
            if (!new_id_o.second)
                d = Swap(d);

            return new_id_o.first + "\t" + PrintDirection(d);
        };

        for (DirectedSegment v : g_.directed_segments()) {
            if (g_.segment(v).removed()) {
                //WARN("Graph had removed segments");
                continue;
            }
            for (auto l : g_.outgoing_links(v)) {
                if (!l.IsCanonical())
                    continue;

                if (inner_links.count(l) || inner_links.count(l.Complement()))
                    continue;
                DEBUG("Processing link " << g_.str(l));

                //TODO support CIGAR?
                out <<"L\t" << get_compacted(l.start)
                    << "\t" << get_compacted(l.end)
                    << "\t" << l.overlap() << "M" << "\n";
            }
        }
    }
};

//inline void CompactAndWrite(const Graph &g, const std::string &fn) {
//    std::cout << "Writing compacted graph to " << fn << std::endl;
//    Compactifier compactifier(g);
//    compactifier.Compact(fn);
//    std::cout << "Writing complete" << std::endl;
//}

}