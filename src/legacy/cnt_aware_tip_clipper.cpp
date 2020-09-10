#include "utils.hpp"
#include "subgraph.hpp"

#include <unordered_map>
#include <vector>
#include <set>
#include <cassert>

static std::unordered_map<std::string, uint32_t> ReadCnts(std::string fn) {
    std::unordered_map<std::string, uint32_t> read_cnts;
    std::string seg_name;
    uint32_t cnt;

    std::ifstream is(fn);

    while (is >> seg_name >> cnt) {
        read_cnts[seg_name] = cnt;

        //std::cout << "Populating coverage with " << seg_name << " and " << cov << std::endl;
    }

    return read_cnts;
}

//if max_length == 0 check returns false
static bool IsTip(const Graph &g, const Graph::DirectedNode &v, size_t max_length) {
    if (max_length == 0)
        return false;

    if (g.incoming_cnt(v) > 0 || g.outgoing_cnt(v) != 1)
        return false;

    Graph::Link l = g.outgoing(v).front();
    Graph::DirectedNode n = l.dn;
    assert(g.incoming_cnt(n) > 0);
    if (g.incoming_cnt(n) == 1)
        return false;

    return g.node_length(v) < max_length + l.cigar.size_on_first();
}

int main(int argc, char *argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <gfa input> <gfa output> <read counts file> <max read count> [max length (default: 10Kb)]" << std::endl;
        //std::cerr << "Usage: " << argv[0] << " <gfa input>" << std::endl;
        exit(239);
    }

    const std::string in_fn(argv[1]);
    const std::string out_fn(argv[2]);

    const std::string read_cnts_fn(argv[3]);

    size_t max_read_cnt = std::stol(argv[4]);

    size_t max_length = 10000;
    if (argc > 5) {
        max_length = std::stol(argv[5]);
    }

    const auto read_cnts = ReadCnts(read_cnts_fn);

    //TODO optimize
    std::set<std::string> neighbourhood;

    {
        Graph g;
        std::cout << "Reading graph in GFA from " << in_fn << std::endl;
        ReadFromFile(in_fn, g);

        std::cout << "Searching for tips" << std::endl;
        std::set<Graph::DirectedNode> tips;
        //FIXME replace with iterator
        for (size_t id = 0, n = g.node_cnt(); id < n; ++id) {
            auto it = read_cnts.find(g.str(id));
            assert(it != read_cnts.end());
            if (it->second > max_read_cnt) {
                DEBUG("Node " << g.str(id) << " consisting of too many reads: " << it->second);
                continue;
            }
            for (auto o : Graph::PossibleOrientations()) {
                Graph::DirectedNode v(id, o);
                DEBUG("Looking at directed node " << g.str(v));

                if (IsTip(g, v, max_length)) {
                    std::cout << "Found tip " << g.str(v) << std::endl;
                    tips.insert(v);
                }
            }
        }

        for (size_t id = 0, n = g.node_cnt(); id < n; ++id) {
            neighbourhood.insert(g.str(id));
        }

        for (auto dn : tips) {
            neighbourhood.erase(g.str(dn.n));
        }
    }

    std::cout << "Outputting tip-free subgraph" << std::endl;
    GetSubgraph(in_fn, out_fn, neighbourhood);

    std::cout << "END" << std::endl;
}

