#pragma once

#include <memory>
#include <string>
#include "gfa.h"
#include "gfa-priv.h"

namespace gfa {

enum class Direction {
    FORWARD, REVERSE
};

inline Direction Swap(Direction d) {
    return d == Direction::FORWARD ? Direction::REVERSE : Direction::FORWARD;
}

std::string PrintDirection(Direction d) {
    return d == Direction::FORWARD ? "forward" : "reverse";
}

struct DirectedSegment {
    uint32_t segment_id;
    Direction direction;

    DirectedSegment(uint32_t segment_id_, Direction direction_) : segment_id(segment_id_), direction(direction_) { }

    DirectedSegment() : segment_id(-1u), direction(Direction::FORWARD) { }

    //DirectedSegment& operator=(const DirectedSegment &ds) {segment_id = ds.segment_id;direction = ds.direction; return *this;}// = default;

    static DirectedSegment FromInnerVertexT(uint32_t v) {
        return DirectedSegment(v >> 1, (v & 1) == 0 ? Direction::FORWARD : Direction::REVERSE);
    }

    uint32_t AsInnerVertexT() const {
        //TODO simplify?
        return segment_id << 1 | (direction == Direction::FORWARD ? 0 : 1);
    }

    DirectedSegment Complement() const {
        return DirectedSegment(segment_id, Swap(direction));
    }

};

//typedef struct {
//	uint32_t m_aux, l_aux;
//	uint8_t *aux;
//} gfa_aux_t;
//
//typedef struct {
//	uint32_t start, end; // start: starting vertex in the string graph; end: ending vertex
//	uint32_t len2, dummy; // len_r: the other length of the unitig
//	uint32_t m, n; // number of reads
//	uint64_t *a; // list of reads
//	char **name;
//} gfa_utg_t;
//
//typedef struct {
//	int32_t len;
//	uint32_t del:16, circ:16;
//	int32_t snid; // stable name ID
//	int32_t soff; // stable start position
//	int32_t rank; // stable rank
//	char *name, *seq;
//	gfa_utg_t *utg;
//	gfa_aux_t aux;
//} gfa_seg_t;
//
//typedef struct {
//	int32_t len, snid, soff, rank;
//	uint64_t end[2];
//	char *seq;
//} gfa_sfa_t;
//
//typedef struct {
//	char *name;
//	int32_t min, max, rank;
//} gfa_sseq_t;
//
//#define gfa_n_vtx(g) ((g)->n_seg << 1)
//
//typedef struct {
//	// segments
//	uint32_t m_seg, n_seg, max_rank;
//	gfa_seg_t *seg;
//	void *h_names;
//	// persistent names
//	uint32_t m_sseq, n_sseq;
//	gfa_sseq_t *sseq;
//	void *h_snames;
//	// links
//	uint64_t m_arc, n_arc:62, is_srt:1, is_symm:1;
//	gfa_arc_t *arc;
//	gfa_aux_t *link_aux;
//	uint64_t *idx;
//} gfa_t;

//FIXME consider making a wrapper over original type rather than copying fields
struct SegmentInfo {
    //std::string sequence;
    char* sequence;
    //FIXME what about long segments?
    uint32_t length;
    char* name;
    bool removed;

    static SegmentInfo FromInnerSegT(gfa_seg_t &s) {
        return SegmentInfo(s);
    }

//        gfa_seg_t *seg = gfa_->seg + i;
//        uint8_t *kc = gfa_aux_get(seg->aux.l_aux, seg->aux.aux, "KC");
//        unsigned cov = 0;
//        if (kc && kc[0] == 'i')
//            cov = *(int32_t*)(kc+1);

private:

    explicit SegmentInfo(gfa_seg_t &s) :
        sequence(s.seq),
        length(s.len),
        name(s.name),
        removed(s.del) {}
};

//TODO upgrade iterators
class SegmentIterator {
    gfa_seg_t *seg_ptr_;

public:
    explicit SegmentIterator(gfa_seg_t *seg_ptr) :
        seg_ptr_(seg_ptr) {}

    SegmentIterator& operator++() {
        ++seg_ptr_;
        return *this;
    }

    SegmentIterator operator++(int) {
        auto retval = *this;
        ++(*this);
        return retval;
    }

    bool operator==(SegmentIterator other) const {
        return seg_ptr_ == other.seg_ptr_;
    }

    bool operator!=(SegmentIterator other) const {
        return !(*this == other);
    }

    SegmentInfo operator*() const {
        return SegmentInfo::FromInnerSegT(*seg_ptr_);
    }

    // iterator traits
    using difference_type = ptrdiff_t;
    using value_type = SegmentInfo;
    using pointer = const SegmentInfo*;
    using reference = const SegmentInfo&;
    using iterator_category = std::forward_iterator_tag;
};

class DirectedSegmentIterator {
    uint32_t inner_v_;

public:
    explicit DirectedSegmentIterator(uint32_t inner_v) :
        inner_v_(inner_v) {}

    DirectedSegmentIterator& operator++() {
        ++inner_v_;
        return *this;
    }

    DirectedSegmentIterator operator++(int) {
        auto retval = *this;
        ++inner_v_;
        return retval;
    }

    bool operator==(DirectedSegmentIterator other) const {
        return inner_v_ == other.inner_v_;
    }

    bool operator!=(DirectedSegmentIterator other) const {
        return !(*this == other);
    }

    DirectedSegment operator*() const {
        return DirectedSegment::FromInnerVertexT(inner_v_);
    }

    // iterator traits
    using difference_type = uint32_t;
    using value_type = DirectedSegment;
    using pointer = const DirectedSegment*;
    using reference = const DirectedSegment&;
    using iterator_category = std::forward_iterator_tag;
};

struct LinkInfo {
    DirectedSegment start;
    DirectedSegment end;
    int32_t start_overlap;
    int32_t end_overlap;
    uint64_t id; // link_id: a pair of dual arcs are supposed to have the same link_id
    bool removed;
    bool complement;

    LinkInfo():
        start_overlap(0), end_overlap(0),
        id(-1ull), removed(false), complement(false) {}

    static LinkInfo FromInnerArcT(const gfa_arc_t &a) {
        return LinkInfo(a);
    }

    LinkInfo Complement() const {
        LinkInfo answer(*this);
        //std::swap(answer.start, answer.end);
        DirectedSegment tmp = start;
        //start = end;
        //end = tmp;
        std::swap(answer.start_overlap, answer.end_overlap);
        answer.complement = !answer.complement;
        return answer;
    }

private:

    LinkInfo(const gfa_arc_t &a) :
        start(DirectedSegment::FromInnerVertexT(a.v_lv >> 32)),
        end(DirectedSegment::FromInnerVertexT(a.w)),
        start_overlap(a.ov), end_overlap(a.ow),
        id(a.link_id), removed(a.del), complement(a.comp) {}
};

class LinkIterator {
    gfa_arc_t *arc_ptr_;
    bool complement_;

public:
    explicit LinkIterator(gfa_arc_t *arc_ptr, bool complement = false) :
        arc_ptr_(arc_ptr), complement_(complement) {}

    LinkIterator& operator++() {
        ++arc_ptr_;
        return *this;
    }

    LinkIterator operator++(int) {
        auto retval = *this;
        ++(*this);
        return retval;
    }

    bool operator==(LinkIterator other) const {
        return arc_ptr_ == other.arc_ptr_;
    }

    bool operator!=(LinkIterator other) const {
        return !(*this == other);
    }

    LinkInfo operator*() const {
        return complement_ ? LinkInfo::FromInnerArcT(*arc_ptr_).Complement() : LinkInfo::FromInnerArcT(*arc_ptr_);
    }

    // iterator traits
    using difference_type = ptrdiff_t;
    using value_type = LinkInfo;
    using pointer = const LinkInfo*;
    using reference = const LinkInfo&;
    using iterator_category = std::forward_iterator_tag;
};

template<class It>
class ProxyContainer {
    It begin_;
    It end_;
public:
    ProxyContainer(It begin, It end): begin_(begin), end_(end) {}

    It begin() const {
        return begin_;
    }

    It end() const {
        return end_;
    }
};

class Graph {
    std::unique_ptr<gfa_t, void(*)(gfa_t*)> g_ptr_;

public:
    Graph(): g_ptr_(nullptr, gfa_destroy) {}

    Graph(const std::string &filename)
        : g_ptr_(gfa_read(filename.c_str()), gfa_destroy) {}

    gfa_t *get() const { return g_ptr_.get(); }

    uint32_t segment_cnt() const { return g_ptr_->n_seg; }
    uint64_t link_cnt() const { return g_ptr_->n_arc; }

    int32_t id(const std::string &name) const {
        return gfa_name2id(get(), name.c_str());
    }

    bool open(const std::string &filename) {
        g_ptr_.reset(gfa_read(filename.c_str()));
        return (bool)g_ptr_;
    }

    void write(const std::string &filename) const {
        FILE *f = fopen(filename.c_str(), "w");
        gfa_print(get(), f, 0);
        fclose(f);
    }

    bool valid() const { return bool(g_ptr_); }

    void DeleteSegment(uint32_t segment_id) { gfa_seg_del(get(), segment_id); }

    void Cleanup() {
        gfa_cleanup(get());
        gfa_fix_symm_del(get());
    }

    //FIXME make faster by including a pointer to an arc into LinkInfo instead of copying
    void DeleteLink(LinkInfo l) {
        DeleteLink(l.start, l.end);
    }

    void DeleteLink(DirectedSegment v, DirectedSegment w) {
        gfa_arc_del(get(), v.AsInnerVertexT(), w.AsInnerVertexT(), 1);
    }

    uint32_t outgoing_link_cnt(DirectedSegment v) const {
        return gfa_arc_n(get(), v.AsInnerVertexT());
    }

    ProxyContainer<LinkIterator> outgoing_links(DirectedSegment v) const {
        return ProxyContainer<LinkIterator>(outgoing_begin(v), outgoing_end(v));
    }

    LinkIterator outgoing_begin(DirectedSegment v) const {
        return LinkIterator(gfa_arc_a(get(), v.AsInnerVertexT()));
    }

    LinkIterator outgoing_end(DirectedSegment v) const {
        return LinkIterator(gfa_arc_a(get(), v.AsInnerVertexT()) + gfa_arc_n(get(), v.AsInnerVertexT()));
    }

    ProxyContainer<LinkIterator> incoming_links(DirectedSegment v) const {
        return ProxyContainer<LinkIterator>(incoming_begin(v), incoming_end(v));
    }

    LinkIterator incoming_begin(DirectedSegment v) const {
        return LinkIterator(gfa_arc_a(get(), v.Complement().AsInnerVertexT()), /*complement*/true);
    }

    LinkIterator incoming_end(DirectedSegment v) const {
        auto inner_v = v.Complement().AsInnerVertexT();
        return LinkIterator(gfa_arc_a(get(), inner_v) + gfa_arc_n(get(), inner_v), /*complement*/true);
    }

    SegmentInfo segment(uint32_t segment_id) const {
        assert(segment_id < segment_cnt());
        return SegmentInfo::FromInnerSegT(get()->seg[segment_id]);
    }

    SegmentIterator segment_begin() const {
        return SegmentIterator(get()->seg);
    }

    SegmentIterator segment_end() const {
        return SegmentIterator(get()->seg + segment_cnt());
    }

    ProxyContainer<SegmentIterator> segments() const {
        return ProxyContainer<SegmentIterator>(segment_begin(), segment_end());
    }

    DirectedSegmentIterator directed_segment_begin() const {
        return DirectedSegmentIterator(0);
    }

    DirectedSegmentIterator directed_segment_end() const {
        return DirectedSegmentIterator(gfa_n_vtx(get()));
    }

    ProxyContainer<DirectedSegmentIterator> directed_segments() const {
        return ProxyContainer<DirectedSegmentIterator>(directed_segment_begin(), directed_segment_end());
    }

};
}
