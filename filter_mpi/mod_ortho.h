#include "fd.h"
#include "ortho.h"
#include "aux.h"

void mod_ortho(
    double*       psitot,
    double*       pot_local,
    grid_st*      grid,
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