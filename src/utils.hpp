#pragma once

#include <cassert>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <fstream>

//#define DEBUG_LOGGING 1

#define LOG_MSG(msg)                                                    \
  do {                                                                  \
     std::cout << msg << std::endl;                                     \
  } while(0);

#ifdef DEBUG_LOGGING
# define DEBUG(message)                      LOG_MSG(message)
#else
# define DEBUG(message)                      {}/* No trace */
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

typedef std::unordered_map<std::string, uint32_t> SegmentCoverageMap;
inline auto ReadCoverage(std::string fn) {
    SegmentCoverageMap segment_coverage;
    std::string seg_name;
    uint32_t cov;

    std::ifstream is(fn);

    while (is >> seg_name >> cov) {
        segment_coverage[seg_name] = cov;

        //std::cout << "Populating coverage with " << seg_name << " and " << cov << std::endl;
    }

    return segment_coverage;
}

}