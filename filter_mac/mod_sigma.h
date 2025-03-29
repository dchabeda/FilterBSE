#include "fd.h"
#include "energy.h"
#include "aux.h"

void mod_sigma(
    double*       psitot,
    double*       pot_local,
    double*       eig_vals,
    double*       sigma_E,
    grid_st*      grid,
    zomplex*      LS,
    nlc_st*       nlc,
    long*         nl,
    double*       ksqr,
    index_st*     ist,
    par_st*       par,
    flag_st*      flag,
    parallel_st*  parallel
);


void restart_from_sigma(
    double**      psitot,
    double*       pot_local,
    double*       eig_vals,
    grid_st*      grid,
    zomplex*      LS,
    nlc_st*       nlc,
    long*         nl,
    double*       ksqr,
    index_st*     ist,
    par_st*       par,
    flag_st*      flag,
    parallel_st*  parallel
);

