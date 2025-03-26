#include "fd.h"
#include "init.h"
#include "energy.h"
#include "hamiltonian.h"
#include "filter.h"

void mod_filter(
    double*       psi_rank,
    zomplex*      psi,
    zomplex*      phi,
    double*       pot_local,
    grid_st*      grid,
    zomplex*      LS,
    nlc_st*       nlc,
    long*         nl,
    zomplex*      an,
    double*       zn,
    double*       ene_targets,
    double*       ksqr,
    index_st*     ist,
    par_st*       par,
    flag_st*      flag,
    parallel_st*  parallel
);

void gen_newton_coeff(
    zomplex*      an, 
    double*       samp, 
    double*       ene_targets, 
    index_st*     ist, 
    par_st*       par, 
    parallel_st*  parallel
);



