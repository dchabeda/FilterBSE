#include "dist_linalg.h"

/* Gram-matrix (normal-equations) SVD resolves singular values only down to
 * ~sqrt(machine_eps) ~ 1.5e-8 relative to sigma_max: below that the squared
 * spectrum sinks under the eigensolver noise floor. Directions retained below
 * this level would be numerical noise (non-orthonormal). We therefore clamp the
 * SVDEPS cutoff to this resolvable floor. Effective cutoff = max(SVDEPS,
 * DIST_SVD_FLOOR) on the singular values. A finer cutoff at scale needs a QR-
 * based (TSQR) orthogonalization instead; see the header/notes. */
#define DIST_SVD_FLOOR 1.0e-7

/*****************************************************************************/
/* Small utilities                                                           */
/*****************************************************************************/

/* Block-partition Ntot items over P ranks: rank i gets cnt[i] contiguous
 * items starting at global offset off[i]. Remainder spread over the low ranks. */
static void block_partition(long Ntot, int P, long *cnt, long *off)
{
  long base = Ntot / P;
  long rem = Ntot % P;
  long o = 0;
  for (int i = 0; i < P; i++)
  {
    cnt[i] = base + ((long)i < rem ? 1 : 0);
    off[i] = o;
    o += cnt[i];
  }
}

static long lmax_arr(const long *a, int P)
{
  long mx = 0;
  for (int i = 0; i < P; i++)
    if (a[i] > mx)
      mx = a[i];
  return mx;
}

/* complex zero literal usable inside expressions */
static inline MKL_Complex16 zero_c(void)
{
  MKL_Complex16 z = {0.0, 0.0};
  return z;
}

/*****************************************************************************/
/* Ring Gram matrix (complex):  S[i,j] = dv * <Aleft_i | Bright_j>            */
/*                                                                           */
/* Aleft and Bright are column-distributed identically (cnt/off, N cols       */
/* total, m complex rows each). Aleft blocks circulate the ring; each rank     */
/* accumulates the column panel S[:, its own cols] and Allgathers the panels   */
/* into the full N x N column-major S (Hermitian) held on every rank.          */
/*****************************************************************************/
static void dist_gram_cplx(
    const MKL_Complex16 *Aleft, const MKL_Complex16 *Bright,
    long m, const long *cnt, const long *off, long N, int P, int rank,
    double dv, MPI_Comm comm, MPI_Datatype state_t, MKL_Complex16 *S)
{
  const long cmax = lmax_arr(cnt, P);
  const long nloc = cnt[rank];

  MKL_Complex16 *buf = (MKL_Complex16 *)malloc((size_t)cmax * m * sizeof(MKL_Complex16));
  MKL_Complex16 *Cblk = (MKL_Complex16 *)malloc((size_t)cmax * (nloc > 0 ? nloc : 1) * sizeof(MKL_Complex16));
  MKL_Complex16 *panel = (MKL_Complex16 *)calloc((size_t)N * (nloc > 0 ? nloc : 1), sizeof(MKL_Complex16));
  if (!buf || !Cblk || !panel)
  {
    fprintf(stderr, "OUT OF MEMORY in dist_gram_cplx (rank %d)\n", rank);
    MPI_Abort(comm, EXIT_FAILURE);
  }

  memcpy(buf, Aleft, (size_t)nloc * m * sizeof(MKL_Complex16));

  const MKL_Complex16 alpha = {dv, 0.0};
  const MKL_Complex16 zero = {0.0, 0.0};
  int cur_src = rank;

  for (int s = 0; s < P; s++)
  {
    const long cs = cnt[cur_src];
    if (nloc > 0 && cs > 0)
    {
      /* Cblk(cs x nloc) = dv * buf^H * Bright  -> S[off[cur_src].., off[rank]..] */
      cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans,
                  cs, nloc, m, &alpha, buf, m, Bright, m, &zero, Cblk, cs);
      for (long j = 0; j < nloc; j++)
        memcpy(&panel[off[cur_src] + j * N], &Cblk[j * cs], (size_t)cs * sizeof(MKL_Complex16));
    }
    if (rank == 0 && ((s == 0) || (0 == (s % (P / 4 + 1))) || (s == P - 1)))
      print_progress_bar(s + 1, P, -1);
    if (s < P - 1)
    {
      const int next = (rank + 1) % P;
      const int prev = (rank - 1 + P) % P;
      MPI_Sendrecv_replace(buf, (int)cmax, state_t, next, 77, prev, 77, comm, MPI_STATUS_IGNORE);
      cur_src = (cur_src - 1 + P) % P;
    }
  }

  /* Gather column panels into the full N x N (column-major) matrix. */
  MPI_Datatype col_t;
  MPI_Type_contiguous((int)(2 * N), MPI_DOUBLE, &col_t);
  MPI_Type_commit(&col_t);
  int *rc = (int *)malloc(P * sizeof(int));
  int *dp = (int *)malloc(P * sizeof(int));
  for (int i = 0; i < P; i++)
  {
    rc[i] = (int)cnt[i];
    dp[i] = (int)off[i];
  }
  MPI_Allgatherv(panel, (int)nloc, col_t, S, rc, dp, col_t, comm);
  MPI_Type_free(&col_t);

  free(rc);
  free(dp);
  free(buf);
  free(Cblk);
  free(panel);
}

/* Ring back-transform (complex): out_local (m x cnt_out[rank]) = A * W[:, out-block]
 * A is column-distributed (cnt_in/off_in, Nin cols); W is Nin x Nout column-major,
 * replicated on all ranks. Output is column-distributed (cnt_out/off_out, Nout cols). */
static void dist_backtransform_cplx(
    const MKL_Complex16 *A_local, long m,
    const long *cnt_in, const long *off_in, long Nin, int P, int rank,
    const MKL_Complex16 *W, const long *cnt_out, const long *off_out, long Nout,
    MPI_Comm comm, MPI_Datatype state_t, MKL_Complex16 *out_local)
{
  (void)Nout;
  const long cmax = lmax_arr(cnt_in, P);
  const long nout = cnt_out[rank];

  MKL_Complex16 *buf = (MKL_Complex16 *)malloc((size_t)cmax * m * sizeof(MKL_Complex16));
  if (!buf)
  {
    fprintf(stderr, "OUT OF MEMORY in dist_backtransform_cplx (rank %d)\n", rank);
    MPI_Abort(comm, EXIT_FAILURE);
  }
  memcpy(buf, A_local, (size_t)cnt_in[rank] * m * sizeof(MKL_Complex16));

  const MKL_Complex16 one = {1.0, 0.0};
  int cur_src = rank;

  for (int s = 0; s < P; s++)
  {
    const long cs = cnt_in[cur_src];
    if (nout > 0 && cs > 0)
    {
      const MKL_Complex16 *Wsub = &W[off_in[cur_src] + off_out[rank] * Nin]; /* ld = Nin */
      const MKL_Complex16 beta = (s == 0) ? zero_c() : one;
      cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  m, nout, cs, &one, buf, m, Wsub, Nin, &beta, out_local, m);
    }
    else if (s == 0 && nout > 0)
    {
      memset(out_local, 0, (size_t)m * nout * sizeof(MKL_Complex16));
    }
    if (rank == 0 && ((s == 0) || (0 == (s % (P / 4 + 1))) || (s == P - 1)))
      print_progress_bar(s + 1, P, -1);
    if (s < P - 1)
    {
      const int next = (rank + 1) % P;
      const int prev = (rank - 1 + P) % P;
      MPI_Sendrecv_replace(buf, (int)cmax, state_t, next, 78, prev, 78, comm, MPI_STATUS_IGNORE);
      cur_src = (cur_src - 1 + P) % P;
    }
  }
  free(buf);
}

/*****************************************************************************/
/* Ring Gram matrix / back-transform (real, double)                          */
/*****************************************************************************/
static void dist_gram_real(
    const double *Aleft, const double *Bright,
    long m, const long *cnt, const long *off, long N, int P, int rank,
    double dv, MPI_Comm comm, MPI_Datatype state_t, double *S)
{
  const long cmax = lmax_arr(cnt, P);
  const long nloc = cnt[rank];

  double *buf = (double *)malloc((size_t)cmax * m * sizeof(double));
  double *Cblk = (double *)malloc((size_t)cmax * (nloc > 0 ? nloc : 1) * sizeof(double));
  double *panel = (double *)calloc((size_t)N * (nloc > 0 ? nloc : 1), sizeof(double));
  if (!buf || !Cblk || !panel)
  {
    fprintf(stderr, "OUT OF MEMORY in dist_gram_real (rank %d)\n", rank);
    MPI_Abort(comm, EXIT_FAILURE);
  }
  memcpy(buf, Aleft, (size_t)nloc * m * sizeof(double));

  int cur_src = rank;
  for (int s = 0; s < P; s++)
  {
    const long cs = cnt[cur_src];
    if (nloc > 0 && cs > 0)
    {
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                  cs, nloc, m, dv, buf, m, Bright, m, 0.0, Cblk, cs);
      for (long j = 0; j < nloc; j++)
        memcpy(&panel[off[cur_src] + j * N], &Cblk[j * cs], (size_t)cs * sizeof(double));
    }
    if (rank == 0 && ((s == 0) || (0 == (s % (P / 4 + 1))) || (s == P - 1)))
      print_progress_bar(s + 1, P, -1);
    if (s < P - 1)
    {
      const int next = (rank + 1) % P;
      const int prev = (rank - 1 + P) % P;
      MPI_Sendrecv_replace(buf, (int)cmax, state_t, next, 77, prev, 77, comm, MPI_STATUS_IGNORE);
      cur_src = (cur_src - 1 + P) % P;
    }
  }

  MPI_Datatype col_t;
  MPI_Type_contiguous((int)N, MPI_DOUBLE, &col_t);
  MPI_Type_commit(&col_t);
  int *rc = (int *)malloc(P * sizeof(int));
  int *dp = (int *)malloc(P * sizeof(int));
  for (int i = 0; i < P; i++)
  {
    rc[i] = (int)cnt[i];
    dp[i] = (int)off[i];
  }
  MPI_Allgatherv(panel, (int)nloc, col_t, S, rc, dp, col_t, comm);
  MPI_Type_free(&col_t);
  free(rc);
  free(dp);
  free(buf);
  free(Cblk);
  free(panel);
}

static void dist_backtransform_real(
    const double *A_local, long m,
    const long *cnt_in, const long *off_in, long Nin, int P, int rank,
    const double *W, const long *cnt_out, const long *off_out, long Nout,
    MPI_Comm comm, MPI_Datatype state_t, double *out_local)
{
  (void)Nout;
  const long cmax = lmax_arr(cnt_in, P);
  const long nout = cnt_out[rank];
  double *buf = (double *)malloc((size_t)cmax * m * sizeof(double));
  if (!buf)
  {
    fprintf(stderr, "OUT OF MEMORY in dist_backtransform_real (rank %d)\n", rank);
    MPI_Abort(comm, EXIT_FAILURE);
  }
  memcpy(buf, A_local, (size_t)cnt_in[rank] * m * sizeof(double));

  int cur_src = rank;
  for (int s = 0; s < P; s++)
  {
    const long cs = cnt_in[cur_src];
    if (nout > 0 && cs > 0)
    {
      const double *Wsub = &W[off_in[cur_src] + off_out[rank] * Nin];
      const double beta = (s == 0) ? 0.0 : 1.0;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  m, nout, cs, 1.0, buf, m, Wsub, Nin, beta, out_local, m);
    }
    else if (s == 0 && nout > 0)
    {
      memset(out_local, 0, (size_t)m * nout * sizeof(double));
    }
    if (rank == 0 && ((s == 0) || (0 == (s % (P / 4 + 1))) || (s == P - 1)))
      print_progress_bar(s + 1, P, -1);
    if (s < P - 1)
    {
      const int next = (rank + 1) % P;
      const int prev = (rank - 1 + P) % P;
      MPI_Sendrecv_replace(buf, (int)cmax, state_t, next, 78, prev, 78, comm, MPI_STATUS_IGNORE);
      cur_src = (cur_src - 1 + P) % P;
    }
  }
  free(buf);
}

/*****************************************************************************/
/* Small dense eigensolves (root only), and SVD-cutoff W builders            */
/*****************************************************************************/

/* Hermitian eigensolve of an N x N column-major matrix (overwrites S with
 * eigenvectors in columns, ascending eigenvalues in w). LAPACK integers are
 * 64-bit under -DMKL_ILP64, matching the long long usage elsewhere. */
static void heev_cplx(MKL_Complex16 *S, double *w, long N)
{
  long long n = (long long)N, lwork = 3 * (long long)N, info = 0;
  MKL_Complex16 *work = (MKL_Complex16 *)malloc((size_t)lwork * sizeof(MKL_Complex16));
  double *rwork = (double *)malloc((size_t)(3 * N) * sizeof(double));
  zheev_("V", "U", &n, S, &n, w, work, &lwork, rwork, &info);
  if (info)
  {
    fprintf(stderr, "error in zheev_ (dist) info=%lld\n", info);
    exit(EXIT_FAILURE);
  }
  free(work);
  free(rwork);
}

static void syev_real(double *S, double *w, long N)
{
  long long n = (long long)N, lwork = 3 * (long long)N, info = 0;
  double *work = (double *)malloc((size_t)lwork * sizeof(double));
  dsyev_("V", "U", &n, S, &n, w, work, &lwork, &info);
  if (info)
  {
    fprintf(stderr, "error in dsyev_ (dist) info=%lld\n", info);
    exit(EXIT_FAILURE);
  }
  free(work);
}

/* Given the eigendecomposition (eigvecs S, ascending eigvals w) of the Gram
 * matrix, build the SVD back-transform W = V Lambda^{-1/2} for the retained
 * directions (sigma_i/sigma_max >= SVDEPS, i.e. lambda_i/lambda_max >= SVDEPS^2),
 * ordered by decreasing singular value. Returns the retained count r. */
static long build_svd_W_cplx(const MKL_Complex16 *S, const double *w, long N, MKL_Complex16 **Wout)
{
  const double eff = (SVDEPS > DIST_SVD_FLOOR) ? SVDEPS : DIST_SVD_FLOOR;
  const double lmax = w[N - 1];
  const double thresh = lmax * (eff * eff);
  long r = 0;
  for (long t = 0; t < N; t++)
  {
    long e = N - 1 - t;
    if (w[e] > 0.0 && w[e] >= thresh)
      r++;
    else
      break;
  }
  MKL_Complex16 *W = (MKL_Complex16 *)malloc((size_t)N * r * sizeof(MKL_Complex16));
  for (long t = 0; t < r; t++)
  {
    long e = N - 1 - t;
    double inv = 1.0 / sqrt(w[e]);
    for (long i = 0; i < N; i++)
    {
      W[i + t * N].real = S[i + e * N].real * inv;
      W[i + t * N].imag = S[i + e * N].imag * inv;
    }
  }
  *Wout = W;
  return r;
}

static long build_svd_W_real(const double *S, const double *w, long N, double **Wout)
{
  const double eff = (SVDEPS > DIST_SVD_FLOOR) ? SVDEPS : DIST_SVD_FLOOR;
  const double lmax = w[N - 1];
  const double thresh = lmax * (eff * eff);
  long r = 0;
  for (long t = 0; t < N; t++)
  {
    long e = N - 1 - t;
    if (w[e] > 0.0 && w[e] >= thresh)
      r++;
    else
      break;
  }
  double *W = (double *)malloc((size_t)N * r * sizeof(double));
  for (long t = 0; t < r; t++)
  {
    long e = N - 1 - t;
    double inv = 1.0 / sqrt(w[e]);
    for (long i = 0; i < N; i++)
      W[i + t * N] = S[i + e * N] * inv;
  }
  *Wout = W;
  return r;
}

/*****************************************************************************/
/* Distributed SVD orthogonalization                                          */
/*                                                                           */
/* On entry *pA is this rank's block of cnt_in[rank] input states (uniform     */
/* Nin total). On return *pA is realloced to this rank's block of the r        */
/* orthonormal basis vectors, distributed by cnt_out/off_out (block layout).   */
/* Returns the global rank r of the basis.                                    */
/*****************************************************************************/
static long dist_ortho(
    double **pA, long m, int complex_idx, double dv,
    const long *cnt_in, const long *off_in, long Nin,
    long *cnt_out, long *off_out,
    int P, int rank, MPI_Comm comm, MPI_Datatype state_t)
{
  const int isC = (complex_idx == 2);
  long r = 0;

  /* Current distribution of *pA. Starts at the (uniform) input distribution and
   * is updated to the output partition after each pass. Kept in dedicated
   * buffers so a pass never aliases its own input/output distribution. */
  long *cin = (long *)malloc(P * sizeof(long));
  long *oin = (long *)malloc(P * sizeof(long));
  for (int i = 0; i < P; i++) { cin[i] = cnt_in[i]; oin[i] = off_in[i]; }
  long Ncur = Nin;

  for (int pass = 0; pass < 2; pass++)
  {
    if (rank == 0)
    {
      printf("  ortho pass %d: Gram + back-transform (%ld states, ring over %d blocks)\n",
             pass, Ncur, P);
      fflush(stdout);
    }
    /* --- Gram matrix S = dv A^H A (N x N, on every rank) --- */
    long rnew;
    if (isC)
    {
      MKL_Complex16 *S = (MKL_Complex16 *)malloc((size_t)Ncur * Ncur * sizeof(MKL_Complex16));
      dist_gram_cplx((MKL_Complex16 *)(*pA), (MKL_Complex16 *)(*pA), m, cin, oin, Ncur, P, rank, dv, comm, state_t, S);

      double *w = (double *)malloc((size_t)Ncur * sizeof(double));
      MKL_Complex16 *W = NULL;
      long rr = 0;
      if (rank == 0)
      {
        heev_cplx(S, w, Ncur);
        rr = build_svd_W_cplx(S, w, Ncur, &W);
      }
      MPI_Bcast(&rr, 1, MPI_LONG, 0, comm);
      if (rank != 0)
        W = (MKL_Complex16 *)malloc((size_t)Ncur * rr * sizeof(MKL_Complex16));
      MPI_Bcast(W, (int)(2 * Ncur * rr), MPI_DOUBLE, 0, comm);
      free(S);
      free(w);
      rnew = rr;

      block_partition(rnew, P, cnt_out, off_out);
      MKL_Complex16 *U = (MKL_Complex16 *)malloc((size_t)(cnt_out[rank] > 0 ? cnt_out[rank] : 1) * m * sizeof(MKL_Complex16));
      dist_backtransform_cplx((MKL_Complex16 *)(*pA), m, cin, oin, Ncur, P, rank, W, cnt_out, off_out, rnew, comm, state_t, U);
      free(W);
      free(*pA);
      *pA = (double *)U;
    }
    else
    {
      double *S = (double *)malloc((size_t)Ncur * Ncur * sizeof(double));
      dist_gram_real(*pA, *pA, m, cin, oin, Ncur, P, rank, dv, comm, state_t, S);

      double *w = (double *)malloc((size_t)Ncur * sizeof(double));
      double *W = NULL;
      long rr = 0;
      if (rank == 0)
      {
        syev_real(S, w, Ncur);
        rr = build_svd_W_real(S, w, Ncur, &W);
      }
      MPI_Bcast(&rr, 1, MPI_LONG, 0, comm);
      if (rank != 0)
        W = (double *)malloc((size_t)Ncur * rr * sizeof(double));
      MPI_Bcast(W, (int)(Ncur * rr), MPI_DOUBLE, 0, comm);
      free(S);
      free(w);
      rnew = rr;

      block_partition(rnew, P, cnt_out, off_out);
      double *U = (double *)malloc((size_t)(cnt_out[rank] > 0 ? cnt_out[rank] : 1) * m * sizeof(double));
      dist_backtransform_real(*pA, m, cin, oin, Ncur, P, rank, W, cnt_out, off_out, rnew, comm, state_t, U);
      free(W);
      free(*pA);
      *pA = U;
    }

    if (rank == 0)
      printf("  dist-ortho pass %d: %ld -> %ld orthonormal states\n", pass, Ncur, rnew);
    fflush(stdout);

    /* the output partition becomes the input distribution for the next pass */
    for (int i = 0; i < P; i++) { cin[i] = cnt_out[i]; oin[i] = off_out[i]; }
    Ncur = rnew;
    r = rnew;
  }

  free(cin);
  free(oin);
  return r;
}

/*****************************************************************************/
/* Per-rank Hamiltonian application over the rank's local states             */
/*   HU_local[:,j] = H U_local[:,j], using threaded p_hamiltonian + FFT.      */
/*****************************************************************************/
static void apply_H_local(
    const double *U_local, double *HU_local, long nloc, long m, int complex_idx,
    double *pot_local, zomplex *LS, nlc_st *nlc, long *nl, double *ksqr,
    index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel)
{
  const int isC = (complex_idx == 2);
  const long stlen = complex_idx * m;

  fftw_plan_loc planfw, planbw;
  fftw_complex *fftwpsi;
  fftw_init_threads();
  fftw_plan_with_nthreads(par->ham_threads);
  create_fft_plans(&planfw, &planbw, &fftwpsi, ist, FFTW_MEASURE);

  zomplex *psi = (zomplex *)calloc(ist->nspinngrid, sizeof(zomplex));
  zomplex *phi = (zomplex *)calloc(ist->nspinngrid, sizeof(zomplex));
  if (!psi || !phi)
  {
    fprintf(stderr, "OUT OF MEMORY: psi/phi in apply_H_local\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  omp_set_dynamic(0);
  omp_set_num_threads(parallel->nthreads);

  for (long j = 0; j < nloc; j++)
  {
    const long base = j * stlen;
    if (isC)
    {
      for (long g = 0; g < ist->nspinngrid; g++)
      {
        psi[g].re = U_local[base + 2 * g];
        psi[g].im = U_local[base + 2 * g + 1];
      }
    }
    else
    {
      for (long g = 0; g < ist->nspinngrid; g++)
      {
        psi[g].re = U_local[base + g];
        psi[g].im = 0.0;
      }
    }

    memcpy(phi, psi, ist->nspinngrid * sizeof(zomplex));
    p_hamiltonian(phi, psi, pot_local, LS, nlc, nl, ksqr, ist, par, flag,
                  planfw, planbw, fftwpsi, parallel->nthreads);

    if (isC)
    {
      for (long g = 0; g < ist->nspinngrid; g++)
      {
        HU_local[base + 2 * g] = phi[g].re;
        HU_local[base + 2 * g + 1] = phi[g].im;
      }
    }
    else
    {
      for (long g = 0; g < ist->nspinngrid; g++)
        HU_local[base + g] = phi[g].re;
    }

    // Progress bar (rank 0's local share of the H-applications, ~every 25%).
    // The loop is purely local (no collectives), so printing here is safe.
    if ((parallel->mpi_rank == 0) &&
        ((j == 0) || (0 == (j % (nloc / 4 + 1))) || (j == nloc - 1)))
      print_progress_bar((int)(j + 1), (int)nloc, -1);
  }

  free(psi);
  free(phi);
  destroy_fft_plans(planfw, planbw, fftwpsi);
}

/*****************************************************************************/
/* Per-rank eigenvalue variance sigma_E for the local eigenstates            */
/*   sigma_i^2 = <psi_i|H^2|psi_i> - <psi_i|H|psi_i>^2  (all local)           */
/*****************************************************************************/
static void dist_sigma_local(
    const double *Psi_local, long nloc, long m, int complex_idx,
    double *pot_local, zomplex *LS, nlc_st *nlc, long *nl, double *ksqr,
    double *sigma_local, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel)
{
  const int isC = (complex_idx == 2);
  const long stlen = complex_idx * m;
  const double dv = par->dv;

  fftw_plan_loc planfw, planbw;
  fftw_complex *fftwpsi;
  fftw_init_threads();
  fftw_plan_with_nthreads(par->ham_threads);
  create_fft_plans(&planfw, &planbw, &fftwpsi, ist, FFTW_MEASURE);

  zomplex *psi = (zomplex *)calloc(ist->nspinngrid, sizeof(zomplex));
  zomplex *phi = (zomplex *)calloc(ist->nspinngrid, sizeof(zomplex));

  omp_set_dynamic(0);
  omp_set_num_threads(parallel->nthreads);

  for (long j = 0; j < nloc; j++)
  {
    const long base = j * stlen;
    if (isC)
      for (long g = 0; g < ist->nspinngrid; g++)
      {
        psi[g].re = Psi_local[base + 2 * g];
        psi[g].im = Psi_local[base + 2 * g + 1];
      }
    else
      for (long g = 0; g < ist->nspinngrid; g++)
      {
        psi[g].re = Psi_local[base + g];
        psi[g].im = 0.0;
      }

    memcpy(phi, psi, ist->nspinngrid * sizeof(zomplex));
    p_hamiltonian(phi, psi, pot_local, LS, nlc, nl, ksqr, ist, par, flag,
                  planfw, planbw, fftwpsi, parallel->nthreads);
    /* phi = H|psi> now */

    /* eval  = <psi|H|psi> = <psi|phi> */
    double eval = 0.0;
    for (long g = 0; g < ist->nspinngrid; g++)
      eval += psi[g].re * phi[g].re + psi[g].im * phi[g].im;
    eval *= dv;

    /* eval2 = <psi|H^2|psi> = <H psi|H psi> = ||phi||^2  (H Hermitian) */
    double eval2 = 0.0;
    for (long g = 0; g < ist->nspinngrid; g++)
      eval2 += phi[g].re * phi[g].re + phi[g].im * phi[g].im;
    eval2 *= dv;

    sigma_local[j] = sqrt(fabs(eval2 - eval * eval));

    // Progress bar (rank 0's local share, ~every 25%). Local loop, no collectives.
    if ((parallel->mpi_rank == 0) &&
        ((j == 0) || (0 == (j % (nloc / 4 + 1))) || (j == nloc - 1)))
      print_progress_bar((int)(j + 1), (int)nloc, -1);
  }

  free(psi);
  free(phi);
  destroy_fft_plans(planfw, planbw, fftwpsi);
}

/*****************************************************************************/
/* Distributed diagonalization: build H in the orthonormal basis, diagonalise */
/* the small matrix, back-transform eigenvectors to the grid basis.           */
/* On return *pU holds this rank's block of eigenvectors (same cnt/off), and   */
/* eig_vals[0..r) is filled (ascending) on every rank.                        */
/*****************************************************************************/
static void dist_diag(
    double **pU, long m, int complex_idx, double dv,
    const long *cnt, const long *off, long r, int P, int rank,
    MPI_Comm comm, MPI_Datatype state_t,
    double *pot_local, zomplex *LS, nlc_st *nlc, long *nl, double *ksqr,
    double *eig_vals, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel)
{
  const int isC = (complex_idx == 2);
  const long nloc = cnt[rank];
  const long stlen = complex_idx * m;

  /* H U (per-rank, local FFT-based Hamiltonian) */
  double *HU = (double *)malloc((size_t)(nloc > 0 ? nloc : 1) * stlen * sizeof(double));
  if (!HU)
  {
    fprintf(stderr, "OUT OF MEMORY: HU in dist_diag (rank %d)\n", rank);
    MPI_Abort(comm, EXIT_FAILURE);
  }
  if (rank == 0)
  {
    printf("  dist-diag: applying H to %ld orthonormal states | %s\n", r, get_time());
    fflush(stdout);
  }
  apply_H_local(*pU, HU, nloc, m, complex_idx, pot_local, LS, nlc, nl, ksqr,
                ist, par, flag, parallel);

  if (rank == 0)
  {
    printf("  dist-diag: reducing H to %ld x %ld (Gram, ring over %d blocks)\n", r, r, P);
    fflush(stdout);
  }

  /* Hmat = dv U^H (H U)  (r x r, on every rank) */
  if (isC)
  {
    MKL_Complex16 *Hm = (MKL_Complex16 *)malloc((size_t)r * r * sizeof(MKL_Complex16));
    dist_gram_cplx((MKL_Complex16 *)(*pU), (MKL_Complex16 *)HU, m, cnt, off, r, P, rank, dv, comm, state_t, Hm);
    free(HU);

    MKL_Complex16 *C = NULL;
    if (rank == 0)
    {
      printf("  dist-diag: diagonalizing %ld x %ld Hamiltonian | %s\n", r, r, get_time());
      fflush(stdout);
      heev_cplx(Hm, eig_vals, r); /* Hm overwritten with eigenvectors, ascending eig_vals */
      C = Hm;
    }
    else
      C = (MKL_Complex16 *)malloc((size_t)r * r * sizeof(MKL_Complex16));
    MPI_Bcast(eig_vals, (int)r, MPI_DOUBLE, 0, comm);
    MPI_Bcast(C, (int)(2 * r * r), MPI_DOUBLE, 0, comm);

    if (rank == 0)
    {
      printf("  dist-diag: back-transforming eigenvectors to grid (ring over %d blocks)\n", P);
      fflush(stdout);
    }
    MKL_Complex16 *Psi = (MKL_Complex16 *)malloc((size_t)(nloc > 0 ? nloc : 1) * m * sizeof(MKL_Complex16));
    dist_backtransform_cplx((MKL_Complex16 *)(*pU), m, cnt, off, r, P, rank, C, cnt, off, r, comm, state_t, Psi);
    free(C);
    free(*pU);
    *pU = (double *)Psi;
  }
  else
  {
    double *Hm = (double *)malloc((size_t)r * r * sizeof(double));
    dist_gram_real(*pU, HU, m, cnt, off, r, P, rank, dv, comm, state_t, Hm);
    free(HU);

    double *C = NULL;
    if (rank == 0)
    {
      printf("  dist-diag: diagonalizing %ld x %ld Hamiltonian | %s\n", r, r, get_time());
      fflush(stdout);
      syev_real(Hm, eig_vals, r);
      C = Hm;
    }
    else
      C = (double *)malloc((size_t)r * r * sizeof(double));
    MPI_Bcast(eig_vals, (int)r, MPI_DOUBLE, 0, comm);
    MPI_Bcast(C, (int)(r * r), MPI_DOUBLE, 0, comm);

    if (rank == 0)
    {
      printf("  dist-diag: back-transforming eigenvectors to grid (ring over %d blocks)\n", P);
      fflush(stdout);
    }
    double *Psi = (double *)malloc((size_t)(nloc > 0 ? nloc : 1) * m * sizeof(double));
    dist_backtransform_real(*pU, m, cnt, off, r, P, rank, C, cnt, off, r, comm, state_t, Psi);
    free(C);
    free(*pU);
    *pU = Psi;
  }
}

/*****************************************************************************/
/* Collective MPI-IO write of the grid-basis eigenvectors to psi.dat.        */
/* State g lives at byte offset g * stlen * sizeof(double), matching the       */
/* serial writer and the get_n_states.x / BSE readers.                         */
/*****************************************************************************/
static void dist_write_psi_dat(
    const char *fname, const double *Psi_local, long stlen,
    const long *cnt, const long *off, int rank, MPI_Comm comm)
{
  MPI_Datatype state_t;
  MPI_Type_contiguous((int)stlen, MPI_DOUBLE, &state_t);
  MPI_Type_commit(&state_t);

  MPI_File fh;
  MPI_File_delete(fname, MPI_INFO_NULL); /* start fresh (ignore error if absent) */
  MPI_Barrier(comm);
  int ierr = MPI_File_open(comm, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  if (ierr != MPI_SUCCESS)
  {
    fprintf(stderr, "ERROR: MPI_File_open(%s) failed on rank %d\n", fname, rank);
    MPI_Abort(comm, EXIT_FAILURE);
  }

  const MPI_Offset base = (MPI_Offset)off[rank] * (MPI_Offset)stlen * (MPI_Offset)sizeof(double);
  if (cnt[rank] > 0)
    MPI_File_write_at(fh, base, Psi_local, (int)cnt[rank], state_t, MPI_STATUS_IGNORE);

  MPI_File_close(&fh);
  MPI_Type_free(&state_t);
}

/*****************************************************************************/
/* Driver                                                                     */
/*****************************************************************************/
void run_dist_postfilter(
    double **psi_rank, double *pot_local, xyz_st *R, grid_st *grid,
    zomplex *LS, nlc_st *nlc, long *nl, double *ksqr,
    double *eig_vals, double *sigma_E,
    index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel)
{
  const int P = parallel->mpi_size;
  const int rank = parallel->mpi_rank;
  const int complex_idx = ist->complex_idx;
  const long m = ist->nspinngrid;             /* complex vector length (rows) */
  const long stlen = complex_idx * m;         /* doubles per state            */
  const double dv = par->dv;
  MPI_Comm comm = MPI_COMM_WORLD;

  double t_clock = (double)clock();
  double t_wall = (double)time(NULL);

  if (rank == 0)
  {
    write_separation(stdout, "T");
    printf("\n5-7. DISTRIBUTED ORTHO / DIAG / SIGMA (MPI, %d ranks) | %s\n", P, get_time());
    write_separation(stdout, "B");
    fflush(stdout);
  }

  /* one MPI datatype for a whole state (stlen doubles) used by every ring */
  MPI_Datatype state_t;
  MPI_Type_contiguous((int)stlen, MPI_DOUBLE, &state_t);
  MPI_Type_commit(&state_t);

  /* --- Time-reverse spinors (local, doubles the states on each rank) --- */
  long nloc_in = ist->n_states_per_rank;
  if ((1 == flag->useSpinors) && (1 != flag->noTimeRev))
  {
    *psi_rank = (double *)realloc(*psi_rank, (size_t)par->t_rev_factor * nloc_in * stlen * sizeof(double));
    if (!*psi_rank)
    {
      fprintf(stderr, "OUT OF MEMORY: psi_rank realloc for time reversal (rank %d)\n", rank);
      MPI_Abort(comm, EXIT_FAILURE);
    }
    time_reverse_all(*psi_rank, &(*psi_rank)[nloc_in * stlen], nloc_in, ist, parallel);
    nloc_in *= par->t_rev_factor;
  }
  normalize_all(*psi_rank, nloc_in, ist, par, flag, parallel);

  /* Uniform input distribution: every rank owns nloc_in states. */
  long *cnt_in = (long *)malloc(P * sizeof(long));
  long *off_in = (long *)malloc(P * sizeof(long));
  for (int i = 0; i < P; i++)
  {
    cnt_in[i] = nloc_in;
    off_in[i] = (long)i * nloc_in;
  }
  const long Nin = (long)P * nloc_in;

  if (rank == 0)
  {
    printf("  states for ortho: %ld total (%ld per rank) | %s\n", Nin, nloc_in, get_time());
    fflush(stdout);
  }

  /* --- Distributed SVD orthogonalization --- */
  long *cnt = (long *)malloc(P * sizeof(long));
  long *off = (long *)malloc(P * sizeof(long));
  long r = dist_ortho(psi_rank, m, complex_idx, dv, cnt_in, off_in, Nin,
                      cnt, off, P, rank, comm, state_t);
  free(cnt_in);
  free(off_in);

  ist->mn_states_tot = r; /* global count of orthonormal/eigen states */
  if (rank == 0)
  {
    printf("\n  SVD cutoff (orthonormal states) r = %ld | %s\n", r, get_time());
    fflush(stdout);
  }

  /* --- Distributed diagonalization (eigenvectors returned in grid basis) --- */
  dist_diag(psi_rank, m, complex_idx, dv, cnt, off, r, P, rank, comm, state_t,
            pot_local, LS, nlc, nl, ksqr, eig_vals, ist, par, flag, parallel);

  normalize_all(*psi_rank, cnt[rank], ist, par, flag, parallel);

  /* --- Eigenvalue variance (ghost-state diagnostic), per-rank local --- */
  if (rank == 0)
  {
    printf("  dist-sigma: computing sigma_E for %ld states | %s\n", r, get_time());
    fflush(stdout);
  }
  double *sigma_local = (double *)calloc((size_t)(cnt[rank] > 0 ? cnt[rank] : 1), sizeof(double));
  dist_sigma_local(*psi_rank, cnt[rank], m, complex_idx, pot_local, LS, nlc, nl, ksqr,
                   sigma_local, ist, par, flag, parallel);

  /* gather sigma_E to root (eig_vals is already global) */
  int *rc = (int *)malloc(P * sizeof(int));
  int *dp = (int *)malloc(P * sizeof(int));
  for (int i = 0; i < P; i++)
  {
    rc[i] = (int)cnt[i];
    dp[i] = (int)off[i];
  }
  MPI_Gatherv(sigma_local, (int)cnt[rank], MPI_DOUBLE, sigma_E, rc, dp, MPI_DOUBLE, 0, comm);
  free(rc);
  free(dp);
  free(sigma_local);

  /* --- Write eigenvectors (MPI-IO) and evals (root) --- */
  if (rank == 0)
  {
    printf("  writing %ld eigenvectors to psi.dat (MPI-IO) | %s\n", r, get_time());
    fflush(stdout);
  }
  dist_write_psi_dat("psi.dat", *psi_rank, stlen, cnt, off, rank, comm);

  if (rank == 0)
  {
    write_eval_dat(eig_vals, sigma_E, r, "eval.dat");
    printf("  wrote eval.dat (%ld eigenvalues)\n", r);
    printf("\n  NOTE: distributed path writes psi.dat + eval.dat. The monolithic\n"
           "        output.dat is not built at this scale; run BSE from the\n"
           "        psi.par/eval.par/conf.par (unsafe_input) path.\n");
    fflush(stdout);
  }

  free(cnt);
  free(off);
  MPI_Type_free(&state_t);

  if (rank == 0)
  {
    printf("\ndone with distributed ortho/diag/sigma, CPU time (sec) %g, wall run time (sec) %g\n",
           ((double)clock() - t_clock) / (double)CLOCKS_PER_SEC, (double)time(NULL) - t_wall);
    fflush(stdout);
  }

  (void)R;
  (void)grid;
}

/*****************************************************************************/
/* Self-test (compile with -DDIST_LINALG_TEST as a standalone executable).   */
/* Validates the ring Gram matrix and the distributed SVD orthogonalization  */
/* against a dense single-rank reference on a small random complex matrix.   */
/*****************************************************************************/
#ifdef DIST_LINALG_TEST

/* Deterministic pseudo-random complex entry for global column g, row i, so
 * every rank agrees on the same global matrix A while storing only its cols. */
static MKL_Complex16 gen_entry(long g, long i)
{
  unsigned long h = (unsigned long)(g * 2654435761UL) ^ (unsigned long)(i * 40503UL + 12345UL);
  h ^= h >> 13;
  h *= 0x5bd1e995UL;
  h ^= h >> 15;
  double re = ((double)(h & 0xffff) / 65535.0) - 0.5;
  double im = ((double)((h >> 16) & 0xffff) / 65535.0) - 0.5;
  MKL_Complex16 z = {re, im};
  return z;
}

/* Rank-deficient generator: columns >= K are damped to 1e-13. */
static MKL_Complex16 gen_scaled(long g, long i, long K)
{
  MKL_Complex16 z = gen_entry(g, i);
  if (g >= K)
  {
    z.real *= 1e-13;
    z.imag *= 1e-13;
  }
  return z;
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  int P, rank;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &P);
  MPI_Comm_rank(comm, &rank);

  const long m = 2000; /* grid rows */
  const long N = 12 * P; /* columns, divisible by P */
  const long nloc = N / P;
  const double dv = 0.7;
  /* Rank-deficient test matrix: the first KRANK columns are O(1) and
   * independent; the rest are scaled to 1e-13 (singular values well below
   * SVDEPS=1e-10, with a clean gap) so the SVD cutoff must return r = KRANK.
   * KRANK=17 is coprime to typical rank counts, exercising the uneven
   * block_partition / uneven-ring paths. */
  const long KRANK = (N >= 17) ? 17 : N;

  MPI_Datatype state_t;
  MPI_Type_contiguous((int)(2 * m), MPI_DOUBLE, &state_t);
  MPI_Type_commit(&state_t);

  long *cnt = malloc(P * sizeof(long)), *off = malloc(P * sizeof(long));
  for (int i = 0; i < P; i++) { cnt[i] = nloc; off[i] = (long)i * nloc; }

  /* local columns of A */
  MKL_Complex16 *A = malloc((size_t)nloc * m * sizeof(MKL_Complex16));
  for (long j = 0; j < nloc; j++)
    for (long i = 0; i < m; i++)
      A[j * m + i] = gen_scaled(off[rank] + j, i, KRANK);

  /* ---- Test 1: ring Gram vs dense reference ---- */
  MKL_Complex16 *S = malloc((size_t)N * N * sizeof(MKL_Complex16));
  dist_gram_cplx(A, A, m, cnt, off, N, P, rank, dv, comm, state_t, S);

  double gerr = 0.0;
  if (rank == 0)
  {
    for (long a = 0; a < N; a++)
      for (long b = 0; b < N; b++)
      {
        double sre = 0, sim = 0;
        for (long i = 0; i < m; i++)
        {
          MKL_Complex16 ca = gen_scaled(a, i, KRANK), cb = gen_scaled(b, i, KRANK);
          /* conj(ca)*cb */
          sre += ca.real * cb.real + ca.imag * cb.imag;
          sim += ca.real * cb.imag - ca.imag * cb.real;
        }
        sre *= dv; sim *= dv;
        double dr = S[a + b * N].real - sre, di = S[a + b * N].imag - sim;
        double e = sqrt(dr * dr + di * di);
        if (e > gerr) gerr = e;
      }
    printf("[test] max |S_ring - S_ref| = %.3e (N=%ld, m=%ld, P=%d)\n", gerr, N, m, P);
  }

  /* ---- Test 2: distributed SVD orthogonalization ---- */
  MKL_Complex16 *Acopy = malloc((size_t)nloc * m * sizeof(MKL_Complex16));
  memcpy(Acopy, A, (size_t)nloc * m * sizeof(MKL_Complex16));
  double *Aptr = (double *)Acopy;
  long *co = malloc(P * sizeof(long)), *oo = malloc(P * sizeof(long));
  long r = dist_ortho(&Aptr, m, 2, dv, cnt, off, N, co, oo, P, rank, comm, state_t);

  /* gather U to rank 0 and check dv*U^H U ~ I_r */
  MPI_Datatype col_t;
  MPI_Type_contiguous((int)(2 * m), MPI_DOUBLE, &col_t);
  MPI_Type_commit(&col_t);
  MKL_Complex16 *U0 = (rank == 0) ? malloc((size_t)r * m * sizeof(MKL_Complex16)) : NULL;
  int *rc = malloc(P * sizeof(int)), *dp = malloc(P * sizeof(int));
  for (int i = 0; i < P; i++) { rc[i] = (int)co[i]; dp[i] = (int)oo[i]; }
  MPI_Gatherv(Aptr, (int)co[rank], col_t, U0, rc, dp, col_t, 0, comm);

  if (rank == 0)
  {
    double oerr = 0.0;
    for (long a = 0; a < r; a++)
      for (long b = 0; b < r; b++)
      {
        double sre = 0, sim = 0;
        for (long i = 0; i < m; i++)
        {
          MKL_Complex16 ca = U0[a * m + i], cb = U0[b * m + i];
          sre += ca.real * cb.real + ca.imag * cb.imag;
          sim += ca.real * cb.imag - ca.imag * cb.real;
        }
        sre *= dv; sim *= dv;
        double target = (a == b) ? 1.0 : 0.0;
        double e = sqrt((sre - target) * (sre - target) + sim * sim);
        if (e > oerr) oerr = e;
      }
    printf("[test] r = %ld (expected %ld = KRANK)\n", r, KRANK);
    printf("[test] max |dv U^H U - I| = %.3e\n", oerr);
    printf("[test] %s\n", (gerr < 1e-9 && oerr < 1e-8 && r == KRANK) ? "PASS" : "FAIL");
  }

  MPI_Finalize();
  return 0;
}
#endif
