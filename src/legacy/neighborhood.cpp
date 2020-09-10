#include "subgraph.hpp"
#include "utils.hpp"

#include <vector>
#include <set>
#include <cassert>

std::vector<std::string> ReadNodes(const std::string &fn) {
    std::vector<std::string> answer;
    std::ifstream is(fn);
    std::string node;
    while (is >> node) {
        answer.push_back(node);
    }
    return answer;
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <gfa input> <gfa output> <nodes of interest file> [neighbourhood depth (default 10)]" << std::endl;
        exit(239);
    }

    const std::string in_fn(argv[1]); // = "/gpfs/gsfs10/users/Phillippy/projects/gfa_works/defensine/processed.gfa";
    const std::string out_fn(argv[2]); // = "/gpfs/gsfs10/users/Phillippy/projects/gfa_works/defensine/subgraph.gfa";
    const std::string nodes_fn(argv[3]); // = "/gpfs/gsfs10/users/Phillippy/projects/gfa_works/defensine/reads.txt";

    uint32_t max_depth = 10;
    if (argc > 4) {
        max_depth = std::stoi(argv[4]);
    }
    std::cout << "Neighbourhood depth set to " << max_depth << std::endl;

    Graph g;
    std::cout << "Reading graph in GFA from " << in_fn << std::endl;
    ReadFromFile(in_fn, g);

    std::cout << "Reading list of nodes of interest from " << nodes_fn << std::endl;
    std::vector<std::string> nodes_of_interest;
    for (const auto &node_str : ReadNodes(nodes_fn)) {
        if (!g.contains(node_str)) {
            std::cerr << "Couldn't find the node " << node_str << std::endl;
            continue;
        }
        nodes_of_interest.push_back(node_str);
    }

    std::set<std::string> neighbourhood = CollectNeighborhood(g, nodes_of_interest, max_depth);
    GetSubgraph(in_fn, out_fn, neighbourhood);

    //std::cout << "HERE " << std::endl;
    //auto links = gg.get_seq_to_link();
    //std::cout << "HERE " << std::endl;
    //for (const auto &s2l : links) {
    //    const std::string seq = s2l.first;
    //    std::cout << "Seq " << seq << std::endl;
    //    for (const gfak::link_elem link : s2l.second) {
    //        std::cout << " source " << link.source_name << " sink " << link.sink_name << " sof " << link.source_orientation_forward << " sif " << link.sink_orientation_forward << " cigar " << link.cigar << std::endl;
    //    }
    ////        std::map<std::string, std::vector<link_elem> > seq_to_link;
    ////struct link_elem{
    ////    std::string source_name;
    ////    std::string sink_name;
    ////    bool source_orientation_forward;
    ////    bool sink_orientation_forward;
    ////    std::string cigar;
    //}

    std::cout << "END" << std::endl;
}
