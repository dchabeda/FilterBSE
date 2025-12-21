#include "fd.h"

void bethe_salpeter(
    double *direct,
    double *exchange,
    double *bsmat,
    double *bs_coeff,
    double *h0mat,
    double *xton_ene,
    grid_st *grid,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel);

void diag(
    const int n,
    int nthreads,
    double *mat,
    double *eval);
