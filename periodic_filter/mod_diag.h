#include "fd.h"
#include "aux.h"
#include "energy.h"

void mod_diag(
    double*       psitot,
    double*       pot_local,
    double*       eig_vals,
    double*       sigma_E,
    vector*       G_vecs,
    vector*       k_vecs,
    grid_st*      grid,
    zomplex*      LS,
    nlc_st*       nlc,
    long*         nl,
    zomplex*      an,
    double*       zn,
    double*       ene_targets,
    index_st*     ist,
    par_st*       par,
    flag_st*      flag,
    parallel_st*  parallel
);

void diag_H(
    double*        psitot, 
    double*        pot_local, 
    zomplex*       LS, 
    nlc_st*        nlc, 
    long*          nl, 
    double*        ksqr, 
    double*        eval, 
    index_st*      ist, 
    par_st*        par, 
    flag_st*       flag, 
    parallel_st*   parallel
);
