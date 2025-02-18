#include "fd.h"
#include "hamiltonian.h"
#include "energy.h"
#include "aux.h"

void run_filter_cycle(
    double*       psi_rank, 
    double*       pot_local,
    zomplex*      LS,
    nlc_st*       nlc, 
    long*         nl, 
    double*       ksqr, 
    zomplex*      an, 
    double*       zn, 
    double*       ene_targets, 
    grid_st*      grid, 
    index_st*     ist, 
    par_st*       par, 
    flag_st*      flag, 
    parallel_st*  parallel);

void time_hamiltonian(
    zomplex*      phi, 
    zomplex*      psi, 
    double*       pot_local, 
    zomplex*      LS,
    nlc_st*       nlc, 
    long*         nl, 
    double*       ksqr,
    index_st*     ist, 
    par_st*       par, 
    flag_st*      flag, 
    parallel_st*  parallel
);
