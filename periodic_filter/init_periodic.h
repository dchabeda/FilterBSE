#include "fd.h"

void init_periodic(
    lattice_st*    lattice,
    vector*        G_vecs,
    vector*        k_vecs,
    grid_st*       grid, 
    index_st*      ist, 
    par_st*        par, 
    flag_st*       flag, 
    parallel_st*   parallel
);

void gen_recip_lat_vecs(
    lattice_st*    lattice, 
    index_st*      ist, 
    par_st*        par, 
    flag_st*       flag, 
    parallel_st*   parallel
);

void gen_G_vecs(
    vector*        G_vecs, 
    grid_st*       grid, 
    index_st*      ist, 
    par_st*        par, 
    flag_st*       flag, 
    parallel_st*   parallel
);

void gen_k_vecs(
    vector*        k_vecs, 
    lattice_st*    lattice, 
    index_st*      ist, 
    par_st*        par, 
    flag_st*       flag, 
    parallel_st*   parallel
);

void read_k_path(
    vector**       k_vecs, 
    lattice_st*    lattice, 
    index_st*      ist, 
    par_st*        par, 
    flag_st*       flag, 
    parallel_st*   parallel
);
