#include "fd.h"
#include "hamiltonian.h"


long ortho_real(
    double *psitot,
    double dv,
    long n_cols,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel
);

long ortho_cplx(
    MKL_Complex16 *psitot,
    double dv,
    long n_cols,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel
);

void restart_from_ortho(
    double**      psitot,
    index_st*     ist,
    par_st*       par,
    flag_st*      flag,
    parallel_st*  parallel
);
