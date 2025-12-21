#include "fd.h"
#include "coulomb.h"

void mod_kernel(
    double *psi_qp,
    double **direct,
    double **exchange,
    double complex *pot_bare,
    double complex *pot_screened,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel);
