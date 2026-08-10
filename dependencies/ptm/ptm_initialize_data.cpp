#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <algorithm>
#include "ptm_initialize_data.h"


static void make_facets_clockwise(int num_facets, int8_t (*facets)[3], const double (*points)[3])
{
    double plane_normal[3];
    double origin[3] = {0, 0, 0};

    for (int i = 0;i<num_facets;i++)
        ptm::add_facet(points, facets[i][0], facets[i][1], facets[i][2], facets[i], plane_normal, origin, 0, NULL);
}

static int initialize_graphs(const ptm::refdata_t* s, int8_t* colours)
{
    for (int i = 0;i<s->num_graphs;i++)
    {
        int8_t code[2 * PTM_MAX_EDGES];
        int8_t degree[PTM_MAX_NBRS];
        int _max_degree = ptm::graph_degree(s->num_facets, s->graphs[i].facets, s->num_nbrs, degree);
        assert(_max_degree <= s->max_degree);

        make_facets_clockwise(s->num_facets, s->graphs[i].facets, &s->points[0][1]);
        int ret = ptm::canonical_form_coloured(s->num_facets, s->graphs[i].facets, s->num_nbrs, degree, colours, s->graphs[i].canonical_labelling, (int8_t*)&code[0], &s->graphs[i].hash);
        if (ret != 0)
            return ret;
    }

    return PTM_NO_ERROR;
}

namespace ptm {

// Storage for the sorted hash indices declared in ptm_graph_data.h. Sized to the
// largest graph table so one static buffer per structure suffices; no allocation.
static graph_hash_entry_t hash_index_sc[NUM_SC_GRAPHS];
static graph_hash_entry_t hash_index_fcc[NUM_FCC_GRAPHS];
static graph_hash_entry_t hash_index_hcp[NUM_HCP_GRAPHS];
static graph_hash_entry_t hash_index_ico[NUM_ICO_GRAPHS];
static graph_hash_entry_t hash_index_bcc[NUM_BCC_GRAPHS];
static graph_hash_entry_t hash_index_dcub[NUM_DCUB_GRAPHS];
static graph_hash_entry_t hash_index_dhex[NUM_DHEX_GRAPHS];

const graph_hash_entry_t* graph_hash_index[9] = {NULL};
int graph_hash_index_size[9] = {0};

// std::sort is not stable, so sort on (hash, graph_index) to make the order of
// graphs sharing a hash deterministic and identical to the old linear scan.
static bool hash_entry_less(const graph_hash_entry_t& a, const graph_hash_entry_t& b)
{
    if (a.hash != b.hash)
        return a.hash < b.hash;
    return a.graph_index < b.graph_index;
}

static void build_hash_index(const refdata_t* s, graph_hash_entry_t* storage)
{
    for (int i = 0; i < s->num_graphs; i++)
    {
        storage[i].hash = s->graphs[i].hash;
        storage[i].graph_index = i;
    }
    std::sort(storage, storage + s->num_graphs, hash_entry_less);
    graph_hash_index[s->type] = storage;
    graph_hash_index_size[s->type] = s->num_graphs;
}

}

bool ptm_initialized = false;
int ptm_initialize_global()
{
    if (ptm_initialized)
        return PTM_NO_ERROR;

    int8_t colours[PTM_MAX_POINTS] = {0};
    int8_t dcolours[PTM_MAX_POINTS] = {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    int ret = initialize_graphs(&ptm::structure_sc, colours);
    ret |= initialize_graphs(&ptm::structure_fcc, colours);
    ret |= initialize_graphs(&ptm::structure_hcp, colours);
    ret |= initialize_graphs(&ptm::structure_ico, colours);
    ret |= initialize_graphs(&ptm::structure_bcc, colours);
    ret |= initialize_graphs(&ptm::structure_dcub, dcolours);
    ret |= initialize_graphs(&ptm::structure_dhex, dcolours);

    if (ret != PTM_NO_ERROR)
        return ret;

    // Hashes only exist once initialize_graphs() has run.
    ptm::build_hash_index(&ptm::structure_sc,   ptm::hash_index_sc);
    ptm::build_hash_index(&ptm::structure_fcc,  ptm::hash_index_fcc);
    ptm::build_hash_index(&ptm::structure_hcp,  ptm::hash_index_hcp);
    ptm::build_hash_index(&ptm::structure_ico,  ptm::hash_index_ico);
    ptm::build_hash_index(&ptm::structure_bcc,  ptm::hash_index_bcc);
    ptm::build_hash_index(&ptm::structure_dcub, ptm::hash_index_dcub);
    ptm::build_hash_index(&ptm::structure_dhex, ptm::hash_index_dhex);

    ptm_initialized = true;
    return PTM_NO_ERROR;
}

ptm_local_handle_t ptm_initialize_local()
{
    assert(ptm_initialized);
    return (ptm_local_handle_t)ptm::voronoi_initialize_local();
}

void ptm_uninitialize_local(ptm_local_handle_t ptr)
{
    ptm::voronoi_uninitialize_local(ptr);
}

