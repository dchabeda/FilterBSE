#include "fd.h"
#include "Hmat.h"
#include "energy.h"
#include "aux.h"

void mod_diag(
    double*       psitot,
    double*       pot_local,
    double*       eig_vals,
    double*       sigma_E,
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

