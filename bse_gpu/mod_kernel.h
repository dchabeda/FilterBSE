#include "fd.h"
#include "coulomb.h"

void mod_kernel(
    double*         psi_qp,
    zomplex**       direct,
    zomplex**       exchange,
    double*         pot_bare,
    double*         pot_screened,
    index_st*       ist,
    par_st*         par,
    flag_st*        flag,
    parallel_st*    parallel
);
