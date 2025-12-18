#include "fd.h"
#include "bethe-salpeter.h"

void mod_bse(
    double *psi_qp,
    double *direct,
    double *exchange,
    double **bsmat,
    double **bs_coeff,
    double **h0mat,
    double **xton_ene,
    double *eig_vals,
    grid_st *grid,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel);

void build_h0_mat(
    double *h0mat,
    double *eval,
    index_st *ist);

void build_BSE_mat(
    double *bsmat,
    double *direct,
    double *exchange,
    index_st *ist);
