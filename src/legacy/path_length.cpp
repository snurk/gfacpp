#include "utils.hpp"
#include "subgraph.hpp"

#include <vector>
#include <set>
#include <cassert>

//if max_length == 0 check returns false
static size_t PathLength(const Graph &g, const std::vector<Graph::DirectedNode> &path) {
    if (path.empty())
        return 0;

    size_t total = g.node_length(path.front());

    for (size_t i = 0; i < path.size() - 1; ++i) {
        auto v = path[i];
        auto next = path[i + 1];
        bool found = false;
        for (auto l : g.outgoing(v)) {
            if (l.dn == next) {
                total += g.node_length(next) - l.cigar.size_on_second();
                found = true;
                break;
            }
        }
        if (!found) {
            std::cerr << "Couldn't find necessary link" << std::endl;
            assert(false);
        }
    }
    return total;
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <gfa input> <node+-, node+-, ...>" << std::endl;
        //std::cerr << "Usage: " << argv[0] << " <gfa input>" << std::endl;
        exit(239);
    }

    const std::string in_fn(argv[1]);
    Graph g;
    std::cout << "Reading graph in GFA from " << in_fn << std::endl;
    ReadFromFile(in_fn, g);

    std::vector<Graph::DirectedNode> path;

    for (size_t i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        std::string name = arg.substr(0, arg.size() - 1);
        std::cout << name << std::endl;
        path.emplace_back(g.id(name),
                arg[arg.size() -1] == '+' ? Graph::NodeOrientation::FORWARD : Graph::NodeOrientation::REVERSE);
    }

    std::cout << PathLength(g, path) << std::endl;

}
