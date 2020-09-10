#pragma once

#include "gfakluge.hpp"
#include "utils.hpp"

#include <set>
#include <vector>

class Cigar {
    std::string str_;
    bool defined_;
    std::vector<std::pair<size_t, char>> events_;

    static char Complement(char c) {
        switch (c) {
            case 'I': return 'D';
            case 'D': return 'I';
            case 'M':
            case 'X':
            case '=':
                return c;
            default:
                assert(false);
                return 'M';
        }
    }

public:

    explicit Cigar(const std::string& cigar = "*") {
        str_ = cigar;
        if (cigar == "*") {
            defined_ = false;
        } else {
            defined_ = true;
            size_t prev_pos = 0;
            size_t pos;
            while ((pos = cigar.find_first_of("MIDX=", prev_pos)) != std::string::npos) {
                if (pos > prev_pos)
                    events_.emplace_back(std::stoul(cigar.substr(prev_pos, pos - prev_pos)), cigar[pos]);
                prev_pos = pos + 1;
            }
            assert(prev_pos == cigar.length());
        }
    }

    size_t size_on_first() const {
        size_t sum = 0;
        for (const auto &l_t : events_) {
            if (l_t.second != 'I') {
                sum += l_t.first;
            }
        }
        return sum;
    }

    size_t size_on_second() const {
        size_t sum = 0;
        for (const auto &l_t : events_) {
            if (l_t.second != 'D') {
                sum += l_t.first;
            }
        }
        return sum;
    }

    //Cigar for reverse-complement pair
    //Switches I/D and reverses the order
    Cigar Complement() const {
        std::stringstream ss;
        for (auto it = events_.rbegin(), end = events_.rend(); it != end; ++it) {
            ss << it->first << Complement(it->second);
        }
        return Cigar(ss.str());
    }

    const std::string &str() const {
        return str_;
    }

};

class Graph {
public:
    typedef size_t VertexId;

    enum class NodeOrientation {
        FORWARD,
        REVERSE
    };

    static NodeOrientation OppositeOrientation(NodeOrientation o) {
        return o == NodeOrientation::FORWARD ? NodeOrientation::REVERSE : NodeOrientation::FORWARD;
    }

    static std::string OrientationStr(NodeOrientation o) {
        return o == NodeOrientation::FORWARD ? "+" : "-";
    }

    static std::vector<NodeOrientation> PossibleOrientations() {
        return {NodeOrientation::FORWARD, NodeOrientation::REVERSE};
    }

    struct NodeInfo {
        std::string name;
        size_t length;
        std::string sequence;

        NodeInfo(const std::string &name_, size_t length_) :
            name(name_), length(length_), sequence("*") {}

        NodeInfo(const std::string &name_, const std::string &sequence_) :
            name(name_), length(sequence_.size()), sequence(sequence_) {}

        NodeInfo() : NodeInfo("", 0) {}

    };

    struct DirectedNode {
        VertexId n;
        NodeOrientation o;

        explicit DirectedNode(VertexId n_, NodeOrientation o_ = NodeOrientation::FORWARD): n(n_), o(o_) {}

        DirectedNode(): DirectedNode(-1) {}

        DirectedNode Complement() const {
            return DirectedNode(n, OppositeOrientation(o));
        }

        bool operator==(const DirectedNode &other) const {
            return n == other.n && o == other.o;
        }

        bool operator!=(const DirectedNode &other) const {
            return !(*this == other);
        }

        bool operator<(const DirectedNode &other) const {
            return n == other.n ? (o == NodeOrientation::FORWARD && other.o == NodeOrientation::REVERSE) : n < other.n;
        }
    };

    //static std::vector<DirectedNode> SwitchOrientation(const std::vector<DirectedNode> &ns) {
    //    std::vector<DirectedNode> answer;
    //    answer.reserve(ns.size());
    //    for (const auto &n : ns) {
    //        answer.emplace_back(n.Opposite());
    //    }
    //    return answer;
    //}

    struct Link {
        DirectedNode dn;
        Cigar cigar;

        explicit Link(DirectedNode dn_, const Cigar &cigar_ = Cigar()) : dn(dn_), cigar(cigar_) {}

        Link() {}

        explicit Link(VertexId n_, NodeOrientation o_, const Cigar &cigar_ = Cigar()): dn(n_, o_), cigar(cigar_) {}

        Link Complement() const {
            return Link(dn.Complement(), cigar.Complement());
        }
    };

    //static Graph FromGFAKluge(gfak::GFAKluge &gg) {
    //}

    void FillFromGFAKluge(gfak::GFAKluge &gg) {
        const auto name_to_seq = gg.get_name_to_seq();
        nodes_.reserve(name_to_seq.size());
        for (const auto &n2s : name_to_seq) {
            const std::string &name = n2s.first;
            const gfak::sequence_elem &seq_info = n2s.second;
            //const gfak::sequence_elem &seq = n2s.second;
            //std::cout << "Name " << name << " name " << seq.name << " id " << seq.id << " seq " << seq.sequence << " len " << seq.length << std::endl;
            node_2_id_[name] = nodes_.size();
            nodes_.push_back((seq_info.sequence == "*") ? NodeInfo(name, seq_info.length) : NodeInfo(name, seq_info.sequence));
        }

        outgoing_lists_.resize(nodes_.size());
        incoming_lists_.resize(nodes_.size());
        const auto seq_to_edges = gg.get_seq_to_edges();
        for (const auto &s : gg.get_name_to_seq()) {
            //std::cout << "Seq " << s.first << std::endl;
            for (const auto &link: seq_to_edges.find(s.first)->second){
                //std::cout << " source " << link.source_name << " sink " << link.sink_name << " sof " << link.source_orientation_forward << " sif " << link.sink_orientation_forward << " cigar " << link.alignment << std::endl;
                VertexId sink_id = utils::get(node_2_id_, link.sink_name);
                VertexId source_id = utils::get(node_2_id_, link.source_name);
                Cigar cigar(link.alignment);

                if (link.source_orientation_forward) {
                    if (link.sink_orientation_forward) {
                        outgoing_lists_[source_id].emplace_back(sink_id, NodeOrientation::FORWARD, cigar);
                        incoming_lists_[sink_id].emplace_back(source_id, NodeOrientation::FORWARD, cigar);
                    } else {
                        outgoing_lists_[source_id].emplace_back(sink_id, NodeOrientation::REVERSE, cigar);
                        outgoing_lists_[sink_id].emplace_back(source_id, NodeOrientation::REVERSE, cigar.Complement());
                    }
                } else {
                    if (link.sink_orientation_forward) {
                        incoming_lists_[source_id].emplace_back(sink_id, NodeOrientation::REVERSE, cigar.Complement());
                        incoming_lists_[sink_id].emplace_back(source_id, NodeOrientation::REVERSE, cigar);
                    } else {
                        incoming_lists_[source_id].emplace_back(sink_id, NodeOrientation::FORWARD, cigar.Complement());
                        outgoing_lists_[sink_id].emplace_back(source_id, NodeOrientation::FORWARD, cigar.Complement());
                    }
                }
            }
        }
    }

    bool contains(const std::string &node) const {
        return node_2_id_.count(node);
    }

    const NodeInfo &node(VertexId id) const {
        return nodes_.at(id);
    }

    std::string str(VertexId id) const {
        return node(id).name;
    }

    std::string str(DirectedNode dn) const {
        return str(dn.n) + OrientationStr(dn.o);
    }

    VertexId id(const std::string &node) const {
        return utils::get(node_2_id_, node);
    }

    size_t node_cnt() const {
        return nodes_.size();
    }

    size_t node_length(VertexId id) const {
        return node(id).length;
    }

    size_t node_length(DirectedNode dn) const {
        return node_length(dn.n);
    }

    std::vector<Link> incoming(DirectedNode dn) const {
        return dn.o == NodeOrientation::FORWARD ? incoming_lists_.at(dn.n) : ComplementLinks(outgoing_lists_.at(dn.n));
    }

    size_t incoming_cnt(DirectedNode dn) const {
        return dn.o == NodeOrientation::FORWARD ? incoming_lists_.at(dn.n).size() : outgoing_lists_.at(dn.n).size();
    }

    std::vector<Link> outgoing(DirectedNode dn) const {
        return dn.o == NodeOrientation::FORWARD ? outgoing_lists_.at(dn.n) : ComplementLinks(incoming_lists_.at(dn.n));
    }

    size_t outgoing_cnt(DirectedNode dn) const {
        return dn.o == NodeOrientation::FORWARD ? outgoing_lists_.at(dn.n).size() : incoming_lists_.at(dn.n).size();
    }

    std::vector<VertexId> adjacent(VertexId id) const {
        assert(id < nodes_.size());
        std::vector<VertexId> answer;
        for (const auto &l : incoming_lists_[id]) {
            answer.push_back(l.dn.n);
        }
        for (const auto &l : outgoing_lists_[id]) {
            answer.push_back(l.dn.n);
        }
        return answer;
    }

private:

    //Transforms list of incoming/outgoing links into list of outgoing/incoming links of the reverse node
    static std::vector<Link> ComplementLinks(const std::vector<Link> &links) {
        std::vector<Link> answer;
        answer.reserve(links.size());
        for (const auto &l : links) {
            answer.emplace_back(l.Complement());
        }
        return answer;
    }

    std::vector<NodeInfo> nodes_;
    std::map<std::string, VertexId> node_2_id_;
    std::vector<std::vector<Link>> outgoing_lists_;
    std::vector<std::vector<Link>> incoming_lists_;
};

static void ReadFromFile(const std::string &fn,
                         Graph &g) {
    auto gg = gfak::GFAKluge();
    gg.parse_gfa_file(fn);
    g.FillFromGFAKluge(gg);
}

static void OutputGraph(gfak::GFAKluge &g, const std::string &fn) {
    std::cout << "Outputting the graph to " << fn << std::endl;
    std::ofstream os(fn);
    const std::string header_string = g.header_string();
    if (header_string.size() > 0){
        os << header_string << "\n";
    }

    auto name_to_seq = g.get_name_to_seq();
    for (auto &kv : name_to_seq) {
        os << kv.second.to_string_1() << "\n";
    }

    auto seq_to_edges = g.get_seq_to_edges();
    for (auto &kv : seq_to_edges){
        for (auto &l : kv.second){
            os << l.to_string_1() << "\n";
        }
    }
}

static void OldGetSubgraph(const std::string &i_fn,
                        const std::string &o_fn,
                        const std::set<std::string> &subgraph_segments) {
    auto gg = gfak::GFAKluge();
    gg.parse_gfa_file(i_fn);

    auto subgraph = gfak::GFAKluge();

    std::cout << "Adding segments" << std::endl;
    const auto name_to_seq = gg.get_name_to_seq();
    for (const auto &n2s : name_to_seq) {
        const std::string name = n2s.first;
        if (subgraph_segments.count(name))
            subgraph.add_sequence(n2s.second);
    }

    std::cout << "Adding links" << std::endl;
    const auto seq_to_edges = gg.get_seq_to_edges();
    for (const auto &s : gg.get_name_to_seq()) {
        for (const auto &link: seq_to_edges.find(s.first)->second){
            //std::cout << " source " << link.source_name << " sink " << link.sink_name << " sof " << link.source_orientation_forward << " sif " << link.sink_orientation_forward << " cigar " << link.alignment << std::endl;
            if (subgraph_segments.count(link.source_name) && subgraph_segments.count(link.sink_name))
                subgraph.add_edge(link.source_name, link);
        }
    }

    OutputGraph(subgraph, o_fn);
    std::cout << "Done" << std::endl;
}

static void GetSubgraph(const std::string &i_fn,
                        const std::string &o_fn,
                        const std::set<std::string> &subgraph_segments) {
    //std::cout << "subgraph segments cnt " << subgraph_segments.size() << std::endl;
    auto gg = gfak::GFAKluge();
    std::ofstream os(o_fn);
    os << "H\tVN:Z:1.0" << std::endl;
    gg.for_each_sequence_line_in_file(i_fn.c_str(), [&](gfak::sequence_elem s) {
            if (subgraph_segments.count(s.name)) {
                os << s.to_string_1() << std::endl;
            }
        });

    gg.for_each_edge_line_in_file(i_fn.c_str(), [&](gfak::edge_elem e) {
            //std::cout << "source name " << e.source_name << std::endl;
            //std::cout << "sink name " << e.sink_name << std::endl;
            if (subgraph_segments.count(e.source_name) && subgraph_segments.count(e.sink_name)) {
                os << e.to_string_1() << std::endl;
            } else {
                //std::cout << "Ignored" << std::endl;
            }
        });
}
