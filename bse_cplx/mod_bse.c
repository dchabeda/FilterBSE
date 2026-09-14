#include "mod_bse.h"
#ifdef USE_SCALAPACK
#include "pbse.h"
#endif

/***************************************************************************************/

void mod_bse(
    double complex *psi_qp,
    double complex *direct,
    double complex *exchange,
    double complex **bsmat,
    double complex **bs_coeff,
    double **h0mat,
    double **xton_ene,
    double *eig_vals,
    grid_st *grid,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel)
{

  int mpir = parallel->mpi_rank;

  if (mpir == 0)
  {
    write_separation(stdout, "T");
    printf("\n4.\tSOLVING BETHE-SALPETER EQUATION | %s\n", get_time());
    write_separation(stdout, "B");
    fflush(stdout);
  }

  /* xton_ene (the eigenvalues) is small and needed on every rank. */
  ALLOCATE(xton_ene, ist->n_xton, "xton_ene");

#ifdef USE_SCALAPACK
  const int dist = (parallel->mpi_size > 1);
#else
  const int dist = 0;
#endif

#ifdef USE_SCALAPACK
  if (dist)
  {
    /* Distributed path: `direct` / `exchange` are this rank's block-cyclic tile
     * (descA), storing only the lower triangle. Assemble H = h0 - direct -
     * exchange locally (no full N x N anywhere): off-diagonal H_loc = -(dir+exc);
     * on owned diagonal entries add h0 = eval[a]-eval[i] and force imag=0
     * (matching build_BSE_mat). bsmat/h0mat are unused; bs_coeff is gathered
     * full onto rank 0 for the downstream optical/angular properties. */
    bse_setup_blockcyclic(parallel, ist->n_xton); /* idempotent (kernel set it up) */
    const long ls = bc_local_size(parallel);
    double complex *H_loc = NULL;
    ALLOCATE(&H_loc, ls, "H_loc (block-cyclic tile)");
    for (long k = 0; k < ls; k++)
      H_loc[k] = -(direct[k] + exchange[k]);

    const long N = ist->n_xton, n_ho = ist->n_holes, lidx = ist->lumo_idx;
    for (long Ig = 0; Ig < N; Ig++)
    {
      int owner;
      long loc;
      bc_map(parallel, Ig, Ig, &owner, &loc);
      if (owner == mpir)
      {
        double h0 = eig_vals[lidx + Ig / n_ho] - eig_vals[Ig % n_ho];
        H_loc[loc] = creal(h0 - direct[loc] - exchange[loc]) + 0.0 * I;
      }
    }

    *bsmat = NULL;
    *h0mat = NULL;
    if (mpir == 0)
      ALLOCATE(bs_coeff, ist->n_xton * ist->n_xton, "bs_coeff (rank 0 gather)");
    else
      *bs_coeff = NULL;

    bethe_salpeter_dist(H_loc, direct, exchange, *bs_coeff, eig_vals, *xton_ene,
                        ist, par, flag, parallel);
    free(H_loc);
  }
  else
#endif
  {
    ALLOCATE(bsmat, ist->n_xton * ist->n_xton, "bsmat");
    ALLOCATE(h0mat, ist->n_xton * ist->n_xton, "h0mat");
    ALLOCATE(bs_coeff, ist->n_xton * ist->n_xton, "bs_coeff");

    build_BSE_mat(*bsmat, direct, exchange, ist);
    build_h0_mat(*h0mat, eig_vals, ist);

    bethe_salpeter(
        direct, exchange, *bsmat, *bs_coeff, *h0mat, *xton_ene,
        grid, ist, par, flag, parallel);
  }

  if (mpir == 0)
  {
    printf("\nDone solving BSE | %s\n", get_time());
    fflush(0);
  }

  return;
}

/***************************************************************************************/

void build_h0_mat(
    double *h0mat,
    double *eval,
    index_st *ist)
{

  long a, i, j, ibs;

  FILE *ppsi;
  ibs = 0UL;

  for (a = ist->lumo_idx; a < ist->lumo_idx + ist->n_elecs; a++)
  {
    for (i = 0; i < ist->n_holes; i++, ibs++)
    {
      h0mat[ibs * ist->n_xton + ibs] = eval[a] - eval[i];
    }
  }

  ppsi = fopen("h0.dat", "w");
  for (i = 0; i < ist->n_xton; i++, fprintf(ppsi, "\n"))
  {
    for (j = 0; j < ist->n_xton; j++)
    {
      // fprintf(ppsi,"%.*g ", PR_LEN, h0mat[i*ist->n_xton+j]);
      fprintf(ppsi, "%.6g ", h0mat[i * ist->n_xton + j]);
    }
  }
  fclose(ppsi);

  return;
}

/***************************************************************************************/

void build_BSE_mat(
    double complex *bsmat,
    double complex *direct,
    double complex *exchange,
    index_st *ist)
{

  FILE *ppsi;
  long ibs, jbs;
  long i, j;
  long ut; // upper triangle
  long lt; // lower triangle

  // Construct the BSE matrix from the exchange and direct kernels
  for (ibs = 0; ibs < ist->n_xton; ibs++)
  {
    for (jbs = 0; jbs <= ibs; jbs++)
    {
      ut = ibs * ist->n_xton + jbs;
      lt = jbs * ist->n_xton + ibs;
      // Symmetrize the matrices
      direct[lt] = conj(direct[ut]);

      exchange[lt] = conj(exchange[ut]);

      // Collect values for bsmat
      bsmat[ut] = direct[ut] + exchange[ut];
      bsmat[lt] = direct[lt] + exchange[lt];

      // Enforce Hermitivity by setting imag part of diag elems to 0.0
      if (ibs == jbs)
      {
        bsmat[ut] = creal(bsmat[ut]) + 0.0 * I; //
      }
    }
  }

  ppsi = fopen("bsRE.dat", "w");
  for (i = 0; i < ist->n_xton; i++, fprintf(ppsi, "\n"))
  {
    for (j = 0; j < ist->n_xton; j++)
    {
      fprintf(ppsi, "%.*g ", PR_LEN, creal(bsmat[i * ist->n_xton + j]));
    }
  }
  fclose(ppsi);

  ppsi = fopen("bsIM.dat", "w");
  for (i = 0; i < ist->n_xton; i++, fprintf(ppsi, "\n"))
  {
    for (j = 0; j < ist->n_xton; j++)
    {
      fprintf(ppsi, "%.*g ", PR_LEN, cimag(bsmat[i * ist->n_xton + j]));
    }
  }
  fclose(ppsi);

  return;
}

/*****************************************************************************/
