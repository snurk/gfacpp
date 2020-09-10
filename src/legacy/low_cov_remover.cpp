#include "utils.hpp"
#include "subgraph.hpp"

#include <vector>
#include <set>
#include <cassert>

std::unordered_map<std::string, uint32_t> ReadCoverage(const std::string &fn) {
    std::unordered_map<std::string, uint32_t> segment_coverage;
    std::string seg_name;
    uint32_t cov;

    std::ifstream is(fn);

    assert(is);
    while (is >> seg_name >> cov) {
        segment_coverage[seg_name] = cov;

        //std::cout << "Populating coverage with " << seg_name << " and " << cov << std::endl;
    }

    return segment_coverage;
}

int32_t Coverage(const Graph &g, const std::unordered_map<std::string, uint32_t> &cov_map, Graph::VertexId n) {
    //std::cout << "Checking for coverage of " << g.str(n) << std::endl;
    return cov_map.find(g.str(n))->second;
}

int main(int argc, char *argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <gfa input> <gfa output> <node coverage file> <coverage threshold (int)> [max length (default: 10Kb)]" << std::endl;
        std::cerr << "Nodes with coverage below <coverage threshold> and shorter than [max length] will be removed" << std::endl;
        exit(239);
    }

    const std::string in_fn(argv[1]);
    const std::string out_fn(argv[2]);
    const std::string cov_fn(argv[3]);
    const uint32_t cov_thr = std::stoi(argv[4]);

    size_t max_length = 10000;
    if (argc > 5) {
        max_length = std::stol(argv[5]);
    }

    std::cout << "Reading coverage info from " << cov_fn << std::endl;

    auto node_cov = ReadCoverage(cov_fn);

    //TODO optimize
    std::set<std::string> neighbourhood;

    {
        Graph g;
        std::cout << "Reading graph in GFA from " << in_fn << std::endl;
        ReadFromFile(in_fn, g);

        std::cout << "Searching for short low coverage nodes" << std::endl;
        std::cout << "Nodes with coverage below " << cov_thr << " and shorter than " << max_length << "bp will be removed" << std::endl;
        std::set<Graph::VertexId> to_remove;
        //FIXME replace with iterator
        for (size_t id = 0, n = g.node_cnt(); id < n; ++id) {
            if (Coverage(g, node_cov, id) < cov_thr && g.node_length(id) < max_length) {
                std::cout << "Will remove node " << g.str(id) << std::endl;
                to_remove.insert(id);
            }
        }

        for (size_t id = 0, n = g.node_cnt(); id < n; ++id) {
            neighbourhood.insert(g.str(id));
        }

        for (Graph::VertexId n : to_remove) {
            neighbourhood.erase(g.str(n));
        }
    }

    std::cout << "Outputting resulting subgraph" << std::endl;
    GetSubgraph(in_fn, out_fn, neighbourhood);

    std::cout << "END" << std::endl;
}
