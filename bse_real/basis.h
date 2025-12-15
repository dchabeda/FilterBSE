#include "fd.h"

/************************************************************/

void get_qp_basis_indices(
    double *eig_vals,
    double *sigma_E,
    long **eval_hole_idxs,
    long **eval_elec_idxs,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel);

/************************************************************/

void get_qp_basis(
    double *psi_qp,
    double *eig_vals,
    double *sigma_E,
    index_st *ist,
    par_st *par,
    flag_st *flag);

/************************************************************/
