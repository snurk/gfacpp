#include "tooling.hpp"

#include <iostream>
#include <cassert>
#include <compact.hpp>

int main(int argc, char *argv[]) {
    using namespace clipp;
    std::string graph_in;
    std::string fasta_out;

    auto cli = (
            graph_in << value("input file in GFA (ending with .gfa)"),
            fasta_out << value("output file")
    );

    auto result = parse(argc, argv, cli);

    if (!result) {
        std::cerr << "Outputing unambiguously extended node sequences" << std::endl;
        std::cerr << make_man_page(cli, argv[0]);
        exit(1);
    }

    gfa::Graph g;
    std::cout << "Loading graph from GFA file " << graph_in << std::endl;
    g.open(graph_in);

    std::cout << "Segment cnt: " << g.segment_cnt() << "; link cnt: " << g.link_cnt() << std::endl;

    gfa::UnambiguousFinder ua_finder(g);
    ua_finder.OutputUnambiguous(fasta_out);

    std::cout << "Done" << std::endl;
}

