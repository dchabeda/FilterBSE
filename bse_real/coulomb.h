#include "fd.h"
#include "aux.h"
#include "hartree.h"
#include "write.h"

void calc_eh_kernel_cplx(
    double *psi_qp,
    double complex *pot_bare,
    double complex *pot_screened,
    double *direct,
    double *exchange,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel);

int load_coulomb_mat(
    double *mat,
    char *fileName,
    long *a,
    long *b,
    long *i,
    long *j,
    index_st *ist);
