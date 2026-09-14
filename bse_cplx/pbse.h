/*****************************************************************************/
#pragma once
#include "fd.h"

/*****************************************************************************
 * Distributed (ScaLAPACK) build + solve of the dense Bethe-Salpeter          *
 * eigenproblem, block-cyclic across ALL MPI ranks.                            *
 *                                                                            *
 * The serial path (bethe-salpeter.c) diagonalizes the N x N (N = n_xton)      *
 * Hermitian BSE matrix H = h0 - (direct + exchange) with LAPACKE_zheev on a   *
 * single rank, and the Coulomb kernel used to build the full N x N direct /   *
 * exchange matrices on EVERY rank (16.8 GB each at n_xton = 32400) -> host     *
 * OOM at production sizes.                                                    *
 *                                                                            *
 * The distributed path (mpi_size > 1) never materializes the full matrix on   *
 * any rank: the kernel writes matrix elements straight into a 2D block-cyclic *
 * layout (one ScaLAPACK descriptor, set up once via bse_setup_blockcyclic and *
 * stored on parallel_st), so each rank holds only ~N^2 / P, and pzheevd       *
 * diagonalizes it in place. Only the lower triangle (ibs >= jbs) is stored;    *
 * pzheevd/pzhemm are called with uplo = 'L' so no symmetrization is needed.    *
 *                                                                            *
 * Numerics: the distributed reductions (pzheevd, pzhemm) sum in a different    *
 * order than the serial LAPACK path, so eigenvalues / <H_*> match to solver    *
 * tolerance (~1e-8), not bit-identically. Eigenvectors additionally carry an   *
 * arbitrary per-vector phase; physical observables (eigenvalues, |coeff|^2,    *
 * oscillator strengths, <H_*>) are reproduced.                                *
 *****************************************************************************/

/*---------------------------------------------------------------------------*
 * Block-cyclic grid / descriptor setup (fills the bc_* fields of parallel).   *
 * Idempotent: a second call with the same N is a no-op. Must be called by ALL *
 * ranks before any router / redistribution / solve call.                      *
 *---------------------------------------------------------------------------*/
void bse_setup_blockcyclic(parallel_st *parallel, long N);
void bse_teardown_blockcyclic(parallel_st *parallel);

/* Local element count of this rank's block-cyclic tile (mloc*nloc, min 1). */
long bc_local_size(const parallel_st *p);

/* Owner MPI rank and column-major local index of global element (Ig, Jg).
 * (Global indices are named Ig/Jg, not I/J: <complex.h> makes I a macro.) */
void bc_map(const parallel_st *p, long Ig, long Jg, int *owner, long *locidx);

/* Inverse of the local->global row map: global row index of local row lr. */
long bc_global_row(const parallel_st *p, long lr);
long bc_global_col(const parallel_st *p, long lc);

/*---------------------------------------------------------------------------*
 * Router: bins locally-computed matrix elements by their block-cyclic owner   *
 * (payload = owner-local index + value), then a single MPI_Alltoallv places   *
 * them into every rank's local tile. Each global element is produced by       *
 * exactly one rank, so placement is a store (no summation).                   *
 *---------------------------------------------------------------------------*/
typedef struct
{
  long loc;       /* column-major local index in the owner's tile */
  double re, im;  /* matrix element value                          */
} bc_elem;

typedef struct
{
  bc_elem **buf; /* [P] per-destination growable arrays */
  long *cnt;     /* [P] element counts                  */
  long *cap;     /* [P] capacities                      */
  int P;
} bc_router;

void router_init(bc_router *R, int P);
void router_push(bc_router *R, const parallel_st *p, long ibs, long jbs, double complex v);
void router_flush(bc_router *R, double complex *A_loc, MPI_Comm comm); /* frees per-dest buffers */
void router_free(bc_router *R);

/* Restart helpers (distribute-on-read). router_load_file re-pushes a rank's own
 * checkpoint dump into the router (owners re-resolved); bc_load_owned scans a
 * merged .dat once and keeps only the lower-triangle elements this rank owns. */
long router_load_file(bc_router *R, const parallel_st *p, const char *fileName,
                      long *a_max, long *b_max, long *i_max, long *j_max, index_st *ist);
int bc_load_owned(double complex *A_loc, const char *fileName, const parallel_st *p, index_st *ist);

/*---------------------------------------------------------------------------*
 * Distributed solve. All inputs are this rank's block-cyclic tile (descA in   *
 * parallel->bc_desc), storing only the lower triangle (ibs >= jbs):           *
 *   H_loc   in/out : BSE matrix H = h0 - direct - exchange; overwritten.       *
 *   dir_loc in     : direct kernel   (for <H_dir>).                            *
 *   exc_loc in     : exchange kernel (for <H_exc>).                            *
 *   eig_vals in    : QP eigenvalues  (for the diagonal <H_0>).                 *
 * Outputs: xton_ene[0..N) ascending eigenvalues on every rank; rank 0's        *
 * bs_coeff filled row-major bs_coeff[ibs*N + j] = component ibs of exciton j.  *
 *---------------------------------------------------------------------------*/
void bethe_salpeter_dist(
    double complex *H_loc,
    double complex *dir_loc,
    double complex *exc_loc,
    double complex *bs_coeff,
    double *eig_vals,
    double *xton_ene,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel);
