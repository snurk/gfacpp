#pragma once

#include <cassert>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <fstream>

//#define DEBUG_LOGGING 1
//#define TRACE_LOGGING 1

#define LOG_MSG(msg)                                                    \
  do {                                                                  \
     std::cout << msg << std::endl;                                     \
  } while(0);

#ifdef DEBUG_LOGGING
# define DEBUG(message)                      LOG_MSG(message)
#else
# define DEBUG(message)                      {}/* No trace */
#endif

#if defined DEBUG_LOGGING && defined TRACE_LOGGING
# define TRACE(message)                      LOG_MSG(message)
#else
# define TRACE(message)                      {}/* No trace */
#endif

# define INFO(message)                      LOG_MSG(message)
# define WARN(message)                      LOG_MSG("WARNING: " << message)

namespace utils {

template<class MapT>
const typename MapT::mapped_type &get(const MapT &from, const typename MapT::key_type &key) {
    auto it = from.find(key);
    assert(it != from.end());
    return it->second;
}

template<class MapT>
typename MapT::mapped_type &get(MapT &from, const typename MapT::key_type &key) {
    auto it = from.find(key);
    assert(it != from.end());
    return it->second;
}

template<class MMapT>
const std::vector<typename MMapT::mapped_type> get_all(const MMapT &from, const typename MMapT::key_type &key) {
    std::vector<typename MMapT::mapped_type> answer;
    for (auto it = from.lower_bound(key); it != from.upper_bound(key); ++it) {
        answer.push_back(it->second);
    }
    return answer;
}

template<typename T>
T abs_diff(T a, T b) {
    return (a > b) ? a-b : b-a;
}

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

template<class Vec>
void ReadVec(const std::string &fn, Vec &vec) {
    std::string seg_name;
    typename Vec::value_type val;

    std::ifstream is(fn);

    while (is >> val) {
        vec.push_back(val);
    }
}

template<class Set>
void ReadSet(const std::string &fn, Set &set) {
    std::string seg_name;
    typename Set::value_type val;

    std::ifstream is(fn);

    while (is >> val) {
        set.insert(val);
    }
}

template<class Map>
void ReadMap(const std::string &fn, Map &map) {
    std::string seg_name;
    typename Map::mapped_type val;

    std::ifstream is(fn);

    while (is >> seg_name >> val) {
        map[seg_name] = val;

        TRACE("Populating map with '" << seg_name << "' and " << cov);
    }
}

typedef std::unordered_map<std::string, double> SegmentCoverageMap;

inline
SegmentCoverageMap ReadCoverage(std::string fn) {
    SegmentCoverageMap segment_coverage;
    ReadMap(fn, segment_coverage);
    //std::string seg_name;
    //double cov;

    //std::ifstream is(fn);

    //while (is >> seg_name >> cov) {
    //    segment_coverage[seg_name] = cov;

    //    TRACE("Populating coverage with '" << seg_name << "' and " << cov);
    //}

    return segment_coverage;
}

}
