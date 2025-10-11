#include "fd.h"
#include "aux.h"
#include "vector.h"

void init_grid_params(
    grid_st*      grid,
    xyz_st*       R,
    index_st*     ist,
    par_st*       par,
    flag_st*      flag,
    parallel_st*  parallel);

void build_grid_ksqr(
    double*       ksqr, 
    xyz_st*       R, 
    grid_st*      grid, 
    index_st*     ist, 
    par_st*       par, 
    flag_st*      flag, 
    parallel_st*  parallel);

void set_ene_targets(
    double*       ene_targets, 
    index_st*     ist, 
    par_st*       par, 
    flag_st*      flag, 
    parallel_st*  parallel);

void build_local_pot(
    double*       pot_local,
    pot_st*       pot, 
    xyz_st*       R,
    atom_info*    atom,
    grid_st*      grid,
    index_st*     ist,
    par_st*       par,
    flag_st*      flag,
    parallel_st*  parallel
);

void init_SO_projectors(
    double*       SO_projectors,
    grid_st*      grid, 
    xyz_st*       R, 
    atom_info*    atm, 
    index_st*     ist, 
    par_st*       par, 
    flag_st*      flag, 
    parallel_st*  parallel
);

void init_NL_projectors(
    nlc_st*       nlc, 
    long*         nl, 
    double*       SO_projectors, 
    grid_st*      grid, 
    xyz_st*       R, 
    atom_info*    atm, 
    index_st*     ist, 
    par_st*       par, 
    flag_st*      flag, 
    parallel_st*  parallel
);

void init_filter_states(
    double*       psi_rank, 
    zomplex*      psi, 
    grid_st*      grid, 
    long*          rand_seed, 
    index_st*     ist, 
    par_st*       par, 
    flag_st*      flag, 
    parallel_st*  parallel
);
