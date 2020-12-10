#include "wrapper.hpp"
#include "gfa-priv.h"

namespace gfa {

void Graph::Cleanup() {
    gfa_cleanup(get());
    gfa_fix_symm_del(get());
}

bool Graph::CheckNoDeadLinks() const {
    for (uint64_t k = 0; k < get()->n_arc; ++k) {
        const gfa_arc_t *a = &get()->arc[k];
        if (a->del)
            return false;
    }
    return true;
}

}

