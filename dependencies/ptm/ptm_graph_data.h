#ifndef PTM_GRAPH_DATA_H
#define PTM_GRAPH_DATA_H

#include <stdint.h>
#include "ptm_constants.h"

namespace ptm {

typedef struct
{
    int id;
    uint64_t hash;
    int automorphism_index;
    int num_automorphisms;
    int8_t canonical_labelling[PTM_MAX_POINTS];
    int8_t facets[PTM_MAX_FACETS][3];
} graph_t;

#define NUM_SC_GRAPHS 1
#define NUM_ICO_GRAPHS 1
#define NUM_FCC_GRAPHS 8
#define NUM_HCP_GRAPHS 16
#define NUM_BCC_GRAPHS 218
#define NUM_DCUB_GRAPHS 12
#define NUM_DHEX_GRAPHS 24

extern int8_t automorphisms[][PTM_MAX_POINTS];

extern graph_t graphs_sc[NUM_SC_GRAPHS];
extern graph_t graphs_fcc[NUM_FCC_GRAPHS];
extern graph_t graphs_hcp[NUM_HCP_GRAPHS];
extern graph_t graphs_ico[NUM_ICO_GRAPHS];
extern graph_t graphs_bcc[NUM_BCC_GRAPHS];
extern graph_t graphs_dcub[NUM_DCUB_GRAPHS];
extern graph_t graphs_dhex[NUM_DHEX_GRAPHS];

//------------------------------------------------------------------------------
// Sorted hash index over each structure's graph table.
//
// check_graphs() needs every graph whose canonical-form hash equals the hash of
// the candidate's convex hull. The obvious implementation is a linear scan over
// refdata_t::graphs, but graph_t is ~128 bytes (it carries canonical_labelling
// and a 28x3 facet table), so scanning BCC's 218 entries touches roughly 28 KB /
// 436 cache lines per candidate atom just to read 218 8-byte hashes.
//
// These arrays are (hash, graph_index) pairs sorted by hash, then by graph_index
// so that graphs sharing a hash are still visited in their original order and
// results stay bit-identical. 16 bytes per entry, so the whole BCC index is
// 3.5 KB and a lookup touches a handful of lines.
//
// Built by ptm_initialize_global(), indexed by refdata_t::type (1..8).
//------------------------------------------------------------------------------
typedef struct
{
    uint64_t hash;
    int32_t graph_index;
} graph_hash_entry_t;

extern const graph_hash_entry_t* graph_hash_index[9];
extern int graph_hash_index_size[9];

}

#endif

