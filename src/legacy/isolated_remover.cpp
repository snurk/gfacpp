#include "utils.hpp"
#include "subgraph.hpp"

#include <vector>
#include <set>
#include <cassert>

//if max_length == 0 check returns false
static bool IsIsolated(const Graph &g, const Graph::DirectedNode &v, size_t max_length) {
    if (max_length == 0)
        return false;

    if (g.incoming_cnt(v) + g.outgoing_cnt(v) > 0)
        return false;

    return g.node_length(v) < max_length;
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <gfa input> <gfa output> [max length (default: 10Kb)]" << std::endl;
        //std::cerr << "Usage: " << argv[0] << " <gfa input>" << std::endl;
        exit(239);
    }

    const std::string in_fn(argv[1]);
    const std::string out_fn(argv[2]);

    size_t max_length = 10000;
    if (argc > 3) {
        max_length = std::stol(argv[3]);
    }

    //TODO optimize
    std::set<std::string> neighbourhood;

    {
        Graph g;
        std::cout << "Reading graph in GFA from " << in_fn << std::endl;
        ReadFromFile(in_fn, g);

        std::cout << "Searching for isolated nodes" << std::endl;
        std::set<Graph::DirectedNode> isolated;
        //consist of heaviest paths of outermost bubbles
        std::set<Graph::VertexId> keep;
        //FIXME replace with iterator
        for (size_t id = 0, n = g.node_cnt(); id < n; ++id) {
            for (auto o : Graph::PossibleOrientations()) {
                Graph::DirectedNode v(id, o);
                DEBUG("Looking at directed node " << g.str(v));
                if (IsIsolated(g, v, max_length)) {
                    std::cout << "Found isolated node " << g.str(v) << std::endl;
                    isolated.insert(v);
                }
            }
        }

        for (size_t id = 0, n = g.node_cnt(); id < n; ++id) {
            neighbourhood.insert(g.str(id));
        }

        for (auto dn : isolated) {
            neighbourhood.erase(g.str(dn.n));
        }
    }

    std::cout << "Outputting isolated node free subgraph" << std::endl;
    GetSubgraph(in_fn, out_fn, neighbourhood);

    std::cout << "END" << std::endl;
}
