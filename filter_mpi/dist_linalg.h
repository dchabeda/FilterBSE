#pragma once
#include "fd.h"
#include "aux.h"
#include "hamiltonian.h"

/*****************************************************************************/
/* Fully distributed (MPI) post-filter pipeline for the NON-PERIODIC path.   */
/*                                                                           */
/* The filtered states stay column-distributed exactly as the filter step    */
/* leaves them: MPI rank p owns a contiguous block of whole states (full grid */
/* each). Nothing is ever gathered to a single rank, so the tall-skinny       */
/* matrix A (Ngrid x Nstates, with Ngrid ~ 1e8) is never materialised in one  */
/* place.                                                                     */
/*                                                                           */
/* Two ring-communication primitives do all the heavy lifting:               */
/*   - dist_gram_*      : S_ij = dv <a_i|b_j>, the N x N small matrix         */
/*                        (overlap for ortho, <u_i|H u_j> for diag).          */
/*   - dist_backtransform_*: out = A W, applying a small N x r basis-change   */
/*                        matrix W to the big distributed vectors.            */
/* Both circulate the (large) column blocks around the rank ring while doing  */
/* local BLAS3 (z/dgemm) products; only tiny N x N matrices are reduced.      */
/*                                                                           */
/* Ortho is an SVD via the Gram matrix (S = A^H A = V Lambda V^H, singular    */
/* values sigma = sqrt(lambda)); the SVDEPS cutoff is applied on sigma and    */
/* the orthonormal basis is U = A V Lambda^{-1/2}. A second pass restores      */
/* orthonormality to machine precision (CholeskyQR2-style). Diag builds the    */
/* Hamiltonian in the orthonormal basis, diagonalises the small matrix, and    */
/* back-transforms the eigenvectors to the grid basis. The eigenvectors are    */
/* written to psi.dat with collective MPI-IO at each state's global offset.    */
/*****************************************************************************/

void run_dist_postfilter(
    double **psi_rank,     /* in/out: this rank's block of filtered states    */
    double *pot_local,
    xyz_st *R,
    grid_st *grid,
    zomplex *LS,
    nlc_st *nlc,
    long *nl,
    double *ksqr,
    double *eig_vals,
    double *sigma_E,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel);
