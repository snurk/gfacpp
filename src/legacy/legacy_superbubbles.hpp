#pragma once

#include "graph.hpp"
#include "range.hpp"
#include "math.hpp"

#include <utility>
#include <map>
#include <set>

//TODO update to pseudo-code from miniasm paper
class SuperbubbleFinder {
    typedef Graph::DirectedNode DirectedNode;

    const Graph& g_;
    DirectedNode start_vertex_;
    size_t max_length_;
    size_t max_diff_;
    size_t max_count_;

    size_t cnt_;
    //TODO think of alternative definitions of weight (currently: total k-mer multiplicity)
    //vertex to heaviest path weight / path length range
    std::map<DirectedNode, std::pair<int64_t, Range>> superbubble_vertices_;
    std::map<DirectedNode, DirectedNode> heaviest_backtrace_;
    DirectedNode end_vertex_;

    bool CheckCanBeProcessed(DirectedNode v) const {
        DEBUG("Check if vertex " << g_.str(v) << " is dominated close neighbour");
        for (const Graph::Link &l : g_.incoming(v)) {
            DirectedNode neighbour_v = l.dn;
            if (superbubble_vertices_.count(neighbour_v) == 0) {
                DEBUG("Blocked by external vertex " << g_.str(neighbour_v));
                DEBUG("Check fail");
                return false;
            }
        }
        DEBUG("Check ok");
        return true;
    }

    void UpdateCanBeProcessed(DirectedNode v,
                              std::set<DirectedNode>& can_be_processed,
                              std::set<DirectedNode>& border) const {
        DEBUG("Updating can be processed");
        for (const Graph::Link &l : g_.outgoing(v)) {
            DirectedNode neighbour_v = l.dn;
            DEBUG("Considering neighbor " << g_.str(neighbour_v));
            if (neighbour_v == start_vertex_) {
                assert(v == start_vertex_);
                continue;
            }
            assert(superbubble_vertices_.count(neighbour_v) == 0);
            DEBUG("Adding vertex " << g_.str(neighbour_v) << " to border");
            border.insert(neighbour_v);
            if (CheckCanBeProcessed(neighbour_v)) {
                DEBUG("Adding vertex " << g_.str(neighbour_v) << " to 'can be processed' set");
                can_be_processed.insert(neighbour_v);
            }
        }
    }

    bool CheckNoEdgeToStart(DirectedNode v) {
        for (const Graph::Link &l : g_.outgoing(v)) {
            DirectedNode neighbour_v = l.dn;
            if (neighbour_v == start_vertex_) {
                return false;
            }
        }
        return true;
    }

public:
    SuperbubbleFinder(const Graph& g, DirectedNode v,
                      size_t max_length = -1ull, size_t max_diff = -1ull, size_t max_count = -1ull)
            : g_(g),
              start_vertex_(v),
              max_length_(max_length),
              max_diff_(max_diff),
              max_count_(max_count),
              cnt_(0) {

    }

    //todo handle case when first/last vertex have other outgoing/incoming edges
    //true if no thresholds exceeded
    bool FindSuperbubble() {
        if (g_.outgoing_cnt(start_vertex_) < 2) {
            return false;
        }
        DEBUG("Adding starting vertex " << g_.str(start_vertex_) << " to dominated set");
        superbubble_vertices_[start_vertex_] = std::make_pair(std::numeric_limits<int64_t>::max(), Range(0, 0));
        heaviest_backtrace_[start_vertex_] = DirectedNode();
        cnt_++;
        std::set<DirectedNode> can_be_processed;
        std::set<DirectedNode> border;
        UpdateCanBeProcessed(start_vertex_, can_be_processed, border);
        //prevents the corner case of 'bubble' of single link and loop on the start node
        bool nontrivial = false;
        while (true) {
            //finish after checks and adding the vertex
            //bool is_end = (border.size() == 1 && can_be_processed.size() == 1);
            const bool is_end = (border.size() == 1);
            DEBUG("is_end: " << is_end);
            if (++cnt_ > max_count_) {
                break;
            }
            DirectedNode v;
            if (!is_end) {
                if (can_be_processed.empty()) {
                    DEBUG("No more nodes could be added");
                    break;
                }
                v = *can_be_processed.begin();
            } else {
                v= *border.begin();
                DEBUG("End node mode activated for vertex " << g_.str(v));
            }

            can_be_processed.erase(v);

            DEBUG("Counting distance range for vertex " << g_.str(v));
            size_t min_d = std::numeric_limits<size_t>::max();
            size_t max_d = 0;
            int64_t max_w = 0;
            DirectedNode entry;

            assert(g_.incoming_cnt(v) > 0);
            assert(is_end || CheckCanBeProcessed(v));

            uint32_t used_incoming_cnt = 0;
            for (const Graph::Link &l : g_.incoming(v)) {
                DirectedNode neighbour_v = l.dn;
                //in case of dominated_only == false
                assert(is_end || superbubble_vertices_.count(neighbour_v));
                if (is_end && superbubble_vertices_.count(neighbour_v) == 0) {
                    DEBUG("Incoming link into end node from the node " << g_.str(v) << " not part of the bubble");
                    continue;
                }
                ++used_incoming_cnt;

                int64_t weight;
                Range range;
                std::tie(weight, range) = utils::get(superbubble_vertices_, neighbour_v);
                range.shift(l.cigar.size_on_second() < g_.node_length(v) ? (int64_t) g_.node_length(v) - l.cigar.size_on_second() : 1);
                DEBUG("Link to " << g_.str(l.dn) << " cigar " << l.cigar.str() << " provide distance range " << range);
                if (range.start_pos < min_d)
                    min_d = range.start_pos;
                if (range.end_pos > max_d)
                    max_d = range.end_pos;

                //path weight is the minimal overlap size along the path
                weight = std::min(weight, int64_t(l.cigar.size_on_second()));

                if (weight > max_w) {
                    max_w = weight;
                    entry = neighbour_v;
                }
            }

            nontrivial |= (used_incoming_cnt > 1);

            assert(entry != DirectedNode());
            assert((max_d > 0) && (min_d < std::numeric_limits<size_t>::max()) && (min_d <= max_d));

            DEBUG("Range " << Range(min_d, max_d));
            Range r(min_d, max_d);
            //FIXME re-introduce early check somehow
            //if (r.start_pos > max_length_) {
            //    return false;
            //}
            //Inner vertices cannot have edge to start vertex
            //TODO Also all added nodes have to have an outgoing edge
            if (!is_end && (!CheckNoEdgeToStart(v) || g_.outgoing_cnt(v) == 0)) {
                break;
            }

            DEBUG("Adding vertex " << g_.str(v) << " to dominated set");
            superbubble_vertices_[v] = std::make_pair(max_w, r);
            heaviest_backtrace_[v] = entry;
            DEBUG("Backtrace " << g_.str(entry));
            border.erase(v);
            if (is_end) {
                //FIXME it seems like only start_pos is ever checked
                //can not simplify check since max_length_ default is close to overflow
                if (r.start_pos > g_.node_length(v) && (r.start_pos - g_.node_length(v)) > max_length_) {
                    DEBUG("Length of minimal additional sequence " << (r.start_pos - g_.node_length(v)) << " exceeded limit " << max_length_);
                    break;
                }
                if (r.size() > max_diff_) {
                    DEBUG("Minimal and maximal lengths differed by " << r.size() << " exceeded limit " << max_diff_);
                    break;
                }
                if (!nontrivial) {
                    DEBUG("Trivial bubble component");
                    break;
                }
                end_vertex_ = v;
                return true;
            } else {
                UpdateCanBeProcessed(v, can_be_processed, border);
            }
        }
        DEBUG("Finished search for starting vertex " << g_.str(start_vertex_));
        return false;
    }

    //const std::map<DirectedNode, std::pair<int64_t, Range>>& nodes_info() const {
    //    return superbubble_vertices_;
    //}

    std::set<DirectedNode> nodes() const {
        std::set<DirectedNode> answer;
        for (const auto n_i : superbubble_vertices_) {
            answer.insert(n_i.first);
        }
        return answer;
    }

    Range PathLengthRange() const {
        return end_vertex_ == DirectedNode() ? Range() :
               utils::get(superbubble_vertices_, end_vertex_).second;
    }

    DirectedNode start_vertex() const {
        return start_vertex_;
    }

    DirectedNode end_vertex() const {
        return end_vertex_;
    }

    const std::vector<DirectedNode> HeaviestPath() const {
        assert(end_vertex_ != DirectedNode());
        std::vector<DirectedNode> rev_nodes;
        DirectedNode v = end_vertex_;
        while (v != DirectedNode()) {
            DEBUG("Added to heaviest path node " << g_.str(v));
            rev_nodes.push_back(v);
            v = utils::get(heaviest_backtrace_, v);
        }
        return std::vector<DirectedNode>(rev_nodes.rbegin(), rev_nodes.rend());
    }

};
