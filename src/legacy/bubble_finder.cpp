#include "legacy_superbubbles.hpp"
#include "utils.hpp"
#include "subgraph.hpp"

#include <vector>
#include <set>
#include <cassert>

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <gfa input> <gfa output> [max_length (20000)] [max_diff (2000)]" << std::endl;
        //std::cerr << "Usage: " << argv[0] << " <gfa input>" << std::endl;
        exit(239);
    }

    const std::string in_fn(argv[1]);
    const std::string out_fn(argv[2]);

    size_t max_length = 20000;
    if (argc > 3) {
        max_length = std::stoi(argv[3]);
    }
    std::cout << "Max length set to " << max_length << std::endl;

    size_t max_diff = 2000;
    if (argc > 4) {
        max_diff = std::stoi(argv[4]);
    }
    std::cout << "Max length diff set to " << max_diff << std::endl;

    std::set<std::string> neighbourhood;
    {
        Graph g;
        std::cout << "Reading graph in GFA from " << in_fn << std::endl;
        ReadFromFile(in_fn, g);

        std::cout << "Searching for bubbles" << std::endl;
        std::set<Graph::DirectedNode> in_bubble;
        //consist of heaviest paths of outermost bubbles
        std::set<Graph::VertexId> keep;
        //FIXME replace with iterator
        for (size_t id = 0, n = g.node_cnt(); id < n; ++id) {
            for (auto o : Graph::PossibleOrientations()) {
                Graph::DirectedNode v(id, o);
                DEBUG("Looking at directed node " << g.str(v));
                if (in_bubble.count(v) != 0) {
                    DEBUG("Not considering. Was part of bubble.");
                    continue;
                }
                SuperbubbleFinder finder(g, v, max_length, max_diff);
                if (finder.FindSuperbubble()) {
                    std::cout << "Found superbubble between " << g.str(finder.start_vertex()) << " and " << g.str(finder.end_vertex()) << std::endl;
                    for (Graph::DirectedNode dn : finder.nodes()) {
                        std::cout << g.str(dn) << '\n';
                        //end vertex can be start of a different bubble
                        if (dn != finder.end_vertex()) {
                            in_bubble.insert(dn);
                        }
                        //complement of the start vertex vertex can be start of a different bubble
                        if (dn != finder.start_vertex()) {
                            in_bubble.insert(dn.Complement());
                        }

                        keep.erase(dn.n);
                    }
                    for (Graph::DirectedNode dn : finder.HeaviestPath()) {
                        std::cout << "Keeping node " << g.str(dn) << std::endl;
                        keep.insert(dn.n);
                    }
                }
            }
        }

        for (size_t id = 0, n = g.node_cnt(); id < n; ++id) {
            neighbourhood.insert(g.str(id));
        }

        for (auto dn : in_bubble) {
            size_t id = dn.n;
            if (keep.count(id) == 0) {
                neighbourhood.erase(g.str(id));
            }
        }
    }

    std::cout << "Outputting bubble-free subgraph" << std::endl;
    GetSubgraph(in_fn, out_fn, neighbourhood);

    std::cout << "END" << std::endl;
}

