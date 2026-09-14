/*****************************************************************************/
#include "pbse.h"
#include <limits.h>

/*****************************************************************************
 * ScaLAPACK / PBLAS / BLACS prototypes.                                      *
 * cray-libsci ships these in libsci_<compiler>_mpi[_mp] but provides no C     *
 * header, so we declare the Fortran symbols (trailing underscore, all args    *
 * by pointer) and the C BLACS wrappers we use. Integers are the default       *
 * 32-bit Fortran INTEGER (cray-libsci is NOT ILP64 in the NVHPC build), which *
 * is safe here: only 2D (i,j) indices and per-rank block sizes are ever       *
 * passed as ints -- the flat N*N element count is never handed to ScaLAPACK.  *
 * double complex (C99) and Fortran COMPLEX*16 share layout, so we pass         *
 * (double complex *) directly.                                               *
 *****************************************************************************/

/* BLACS C wrappers */
extern void Cblacs_pinfo(int *mypnum, int *nprocs);
extern void Cblacs_get(int ictxt, int what, int *val);
extern void Cblacs_gridinit(int *ictxt, const char *order, int nprow, int npcol);
extern void Cblacs_gridinfo(int ictxt, int *nprow, int *npcol, int *myrow, int *mycol);
extern void Cblacs_gridexit(int ictxt);

/* ScaLAPACK tools (Fortran) */
extern int numroc_(const int *n, const int *nb, const int *iproc,
                   const int *isrcproc, const int *nprocs);
extern int indxl2g_(const int *indxloc, const int *nb, const int *iproc,
                    const int *isrcproc, const int *nprocs);
extern void descinit_(int *desc, const int *m, const int *n, const int *mb,
                      const int *nb, const int *irsrc, const int *icsrc,
                      const int *ictxt, const int *lld, int *info);

/* PBLAS / ScaLAPACK routines */
extern void pzgemr2d_(const int *m, const int *n,
                      const double complex *a, const int *ia, const int *ja, const int *desca,
                      double complex *b, const int *ib, const int *jb, const int *descb,
                      const int *gcontext);
extern void pzhemm_(const char *side, const char *uplo,
                    const int *m, const int *n,
                    const double complex *alpha,
                    const double complex *a, const int *ia, const int *ja, const int *desca,
                    const double complex *b, const int *ib, const int *jb, const int *descb,
                    const double complex *beta,
                    double complex *c, const int *ic, const int *jc, const int *descc);
extern void pzheevd_(const char *jobz, const char *uplo, const int *n,
                     double complex *a, const int *ia, const int *ja, const int *desca,
                     double *w,
                     double complex *z, const int *iz, const int *jz, const int *descz,
                     double complex *work, const int *lwork,
                     double *rwork, const int *lrwork,
                     int *iwork, const int *liwork, int *info);

/* descriptor field indices (0-based, ScaLAPACK convention) */
#define DTYPE_ 0
#define CTXT_ 1
#define M_ 2
#define N_ 3
#define MB_ 4
#define NB_ 5
#define RSRC_ 6
#define CSRC_ 7
#define LLD_ 8

#define BSE_BLOCK 64 /* block-cyclic block size (MB = NB) */

/*****************************************************************************/
/* Factor P into a near-square process grid (nprow >= npcol, nprow*npcol=P).  */
/*****************************************************************************/
static void factor_grid(int P, int *nprow, int *npcol)
{
  int r = (int)(sqrt((double)P) + 1e-9);
  while (r > 1 && (P % r != 0))
    r--;
  if (r < 1)
    r = 1;
  *npcol = r;
  *nprow = P / r;
}

/* C reimplementation of ScaLAPACK numroc_ with isrcproc = 0: number of local
 * rows/cols of an n-element dimension, block nb, owned by process iproc of
 * nprocs. Lets the sender compute an owner's local leading dimension without a
 * BLACS call in the hot path. */
static long numroc0(long n, int nb, int iproc, int nprocs)
{
  long nblocks = n / nb;
  long loc = (nblocks / nprocs) * nb;
  long extra = nblocks % nprocs;
  if (iproc < extra)
    loc += nb;
  else if (iproc == extra)
    loc += n % nb;
  return loc;
}

/*****************************************************************************/
/* Grid / descriptor setup. Idempotent for a fixed N.                         */
/*****************************************************************************/
void bse_setup_blockcyclic(parallel_st *parallel, long N)
{
  if (parallel->bc_ready && parallel->bc_N == N)
    return;

  const int P = parallel->mpi_size;
  const int izero = 0;
  const int iN = (int)N;
  if ((long)iN != N)
  {
    if (parallel->mpi_rank == 0)
      fprintf(stderr, "ERROR: n_xton = %ld exceeds 32-bit ScaLAPACK index range\n", N);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  int nprow, npcol;
  factor_grid(P, &nprow, &npcol); /* nprow*npcol == P, so every rank is in-grid */
  int ctxt;
  Cblacs_get(0, 0, &ctxt);
  Cblacs_gridinit(&ctxt, "Row", nprow, npcol);
  int myrow, mycol;
  Cblacs_gridinfo(ctxt, &nprow, &npcol, &myrow, &mycol);

  const int nb = (N < BSE_BLOCK) ? iN : BSE_BLOCK;
  int mloc = numroc_(&iN, &nb, &myrow, &izero, &nprow);
  int nloc = numroc_(&iN, &nb, &mycol, &izero, &npcol);
  int lld = (mloc > 0) ? mloc : 1;
  int info = 0;
  descinit_(parallel->bc_desc, &iN, &iN, &nb, &nb, &izero, &izero, &ctxt, &lld, &info);
  if (info != 0 && parallel->mpi_rank == 0)
    fprintf(stderr, "WARNING: descinit returned info = %d\n", info);

  parallel->bc_ctxt = ctxt;
  parallel->bc_nprow = nprow;
  parallel->bc_npcol = npcol;
  parallel->bc_myrow = myrow;
  parallel->bc_mycol = mycol;
  parallel->bc_nb = nb;
  parallel->bc_mloc = mloc;
  parallel->bc_nloc = nloc;
  parallel->bc_lld = lld;
  parallel->bc_N = N;
  parallel->bc_ready = 1;

  if (parallel->mpi_rank == 0)
  {
    printf("  Block-cyclic BSE layout: N = %ld, grid %d x %d, block %d "
           "(<= %.2f GB/rank per N^2 matrix) | %s\n",
           N, nprow, npcol, nb,
           (double)mloc * nloc * sizeof(double complex) / 1e9, get_time());
    fflush(stdout);
  }
}

void bse_teardown_blockcyclic(parallel_st *parallel)
{
  if (!parallel->bc_ready)
    return;
  Cblacs_gridexit(parallel->bc_ctxt);
  parallel->bc_ready = 0;
}

long bc_local_size(const parallel_st *p)
{
  long m = (p->bc_mloc > 0) ? p->bc_mloc : 1;
  long n = (p->bc_nloc > 0) ? p->bc_nloc : 1;
  return m * n;
}

void bc_map(const parallel_st *p, long Ig, long Jg, int *owner, long *locidx)
{
  const int nb = p->bc_nb;
  const int nprow = p->bc_nprow, npcol = p->bc_npcol;
  const long pr = (Ig / nb) % nprow;
  const long pc = (Jg / nb) % npcol;
  *owner = (int)(pr * npcol + pc);
  const long lr = (Ig / ((long)nb * nprow)) * nb + (Ig % nb);
  const long lc = (Jg / ((long)nb * npcol)) * nb + (Jg % nb);
  long owner_mloc = numroc0(p->bc_N, nb, (int)pr, nprow);
  if (owner_mloc < 1)
    owner_mloc = 1;
  *locidx = lc * owner_mloc + lr;
}

long bc_global_row(const parallel_st *p, long lr)
{
  const int nb = p->bc_nb;
  const long b = lr / nb, off = lr % nb;
  return (b * p->bc_nprow + p->bc_myrow) * nb + off;
}

long bc_global_col(const parallel_st *p, long lc)
{
  const int nb = p->bc_nb;
  const long b = lc / nb, off = lc % nb;
  return (b * p->bc_npcol + p->bc_mycol) * nb + off;
}

/*****************************************************************************/
/* Router                                                                     */
/*****************************************************************************/
void router_init(bc_router *R, int P)
{
  R->P = P;
  R->buf = (bc_elem **)calloc((size_t)P, sizeof(bc_elem *));
  R->cnt = (long *)calloc((size_t)P, sizeof(long));
  R->cap = (long *)calloc((size_t)P, sizeof(long));
}

void router_push(bc_router *R, const parallel_st *p, long ibs, long jbs, double complex v)
{
  int owner;
  long loc;
  bc_map(p, ibs, jbs, &owner, &loc);
  if (R->cnt[owner] == R->cap[owner])
  {
    long ncap = R->cap[owner] ? R->cap[owner] * 2 : 1024;
    R->buf[owner] = (bc_elem *)realloc(R->buf[owner], (size_t)ncap * sizeof(bc_elem));
    if (!R->buf[owner])
    {
      fprintf(stderr, "OUT OF MEMORY in router_push (dest %d, %ld elems)\n", owner, ncap);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    R->cap[owner] = ncap;
  }
  bc_elem *e = &R->buf[owner][R->cnt[owner]++];
  e->loc = loc;
  e->re = creal(v);
  e->im = cimag(v);
}

void router_flush(bc_router *R, double complex *A_loc, MPI_Comm comm)
{
  const int P = R->P;
  int *sendcounts = (int *)malloc((size_t)P * sizeof(int));
  int *recvcounts = (int *)malloc((size_t)P * sizeof(int));
  int *sdispls = (int *)malloc((size_t)P * sizeof(int));
  int *rdispls = (int *)malloc((size_t)P * sizeof(int));

  long tot_send = 0;
  for (int d = 0; d < P; d++)
  {
    if (R->cnt[d] > (long)INT_MAX)
    {
      fprintf(stderr, "ERROR: router send count %ld to dest %d exceeds INT_MAX; "
                      "use more MPI ranks so each panel is smaller\n",
              R->cnt[d], d);
      MPI_Abort(comm, EXIT_FAILURE);
    }
    sendcounts[d] = (int)R->cnt[d];
    tot_send += R->cnt[d];
  }
  MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT, comm);

  long tot_recv = 0;
  for (int d = 0; d < P; d++)
    tot_recv += recvcounts[d];
  if (tot_send > (long)INT_MAX || tot_recv > (long)INT_MAX)
  {
    fprintf(stderr, "ERROR: router total count (send %ld / recv %ld) exceeds INT_MAX "
                    "displacement range; use more MPI ranks\n",
            tot_send, tot_recv);
    MPI_Abort(comm, EXIT_FAILURE);
  }

  /* one contiguous element datatype (24 bytes, no padding) */
  MPI_Datatype elem_t;
  MPI_Type_contiguous((int)sizeof(bc_elem), MPI_BYTE, &elem_t);
  MPI_Type_commit(&elem_t);

  bc_elem *sbuf = (bc_elem *)malloc((size_t)(tot_send > 0 ? tot_send : 1) * sizeof(bc_elem));
  bc_elem *rbuf = (bc_elem *)malloc((size_t)(tot_recv > 0 ? tot_recv : 1) * sizeof(bc_elem));
  if (!sbuf || !rbuf || !sendcounts || !recvcounts || !sdispls || !rdispls)
  {
    fprintf(stderr, "OUT OF MEMORY in router_flush\n");
    MPI_Abort(comm, EXIT_FAILURE);
  }

  long off = 0;
  for (int d = 0; d < P; d++)
  {
    sdispls[d] = (int)off;
    if (R->cnt[d])
      memcpy(&sbuf[off], R->buf[d], (size_t)R->cnt[d] * sizeof(bc_elem));
    off += R->cnt[d];
    free(R->buf[d]);
    R->buf[d] = NULL;
  }
  rdispls[0] = 0;
  for (int d = 1; d < P; d++)
    rdispls[d] = rdispls[d - 1] + recvcounts[d - 1];

  MPI_Alltoallv(sbuf, sendcounts, sdispls, elem_t,
                rbuf, recvcounts, rdispls, elem_t, comm);

  for (long k = 0; k < tot_recv; k++)
    A_loc[rbuf[k].loc] = rbuf[k].re + rbuf[k].im * I;

  MPI_Type_free(&elem_t);
  free(sbuf);
  free(rbuf);
  free(sendcounts);
  free(recvcounts);
  free(sdispls);
  free(rdispls);
}

void router_free(bc_router *R)
{
  if (R->buf)
    for (int d = 0; d < R->P; d++)
      free(R->buf[d]);
  free(R->buf);
  free(R->cnt);
  free(R->cap);
  R->buf = NULL;
  R->cnt = R->cap = NULL;
}

/*****************************************************************************/
/* Restart helpers (distribute-on-read).                                      */
/*****************************************************************************/
long router_load_file(bc_router *R, const parallel_st *p, const char *fileName,
                      long *a_max, long *b_max, long *i_max, long *j_max, index_st *ist)
{
  FILE *pf = fopen(fileName, "r");
  if (!pf)
  {
    fprintf(stderr, "ERROR: router_load_file could not open %s\n", fileName);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  long a, b, i, j, ibs, jbs, cnt = 0;
  double re, im;
  long ta = 0, tb = 0, ti = 0, tj = 0;
  const long max_lines = ist->n_holes * ist->n_holes * ist->n_elecs * ist->n_elecs;
  while (cnt < max_lines &&
         fscanf(pf, "%ld %ld %ld %ld %ld %ld %lg %lg", &a, &b, &i, &j, &ibs, &jbs, &re, &im) == 8)
  {
    router_push(R, p, ibs, jbs, re + im * I);
    if (a > ta) ta = a;
    if (b > tb) tb = b;
    if (i > ti) ti = i;
    if (j > tj) tj = j;
    cnt++;
  }
  fclose(pf);
  *a_max = ta;
  *b_max = tb;
  *i_max = ti;
  *j_max = tj;
  return cnt;
}

int bc_load_owned(double complex *A_loc, const char *fileName, const parallel_st *p, index_st *ist)
{
  FILE *pf = fopen(fileName, "r");
  if (!pf)
  {
    if (p->mpi_rank == 0)
      fprintf(stderr, "ERROR: bc_load_owned could not open %s\n", fileName);
    return 0;
  }
  long a, b, i, j, ibs, jbs, cnt = 0;
  double re, im;
  const long max_lines = ist->n_holes * ist->n_holes * ist->n_elecs * ist->n_elecs;
  while (cnt < max_lines &&
         fscanf(pf, "%ld %ld %ld %ld %ld %ld %lg %lg", &a, &b, &i, &j, &ibs, &jbs, &re, &im) == 8)
  {
    cnt++;
    if (ibs < jbs) /* only the stored lower triangle */
      continue;
    int owner;
    long loc;
    bc_map(p, ibs, jbs, &owner, &loc);
    if (owner == p->mpi_rank)
      A_loc[loc] = re + im * I;
  }
  fclose(pf);
  return 1;
}

/*****************************************************************************/
/* Diagonal of U^dagger A U from local block-cyclic images of Z (= U) and     */
/* AU = A*Z sharing the SAME descriptor. Local (lr,lc) contributes             */
/* conj(Z(lr,lc)) * AU(lr,lc) to global column lc's exciton. Allreduced.       */
/*****************************************************************************/
static void contract_diag(const double complex *Zloc, const double complex *AUloc,
                          const parallel_st *p, long N, MPI_Comm comm, double complex *d)
{
  const int lld = p->bc_lld, mloc = p->bc_mloc, nloc = p->bc_nloc;
  double complex *part = (double complex *)calloc((size_t)N, sizeof(double complex));
  if (!part)
  {
    fprintf(stderr, "OUT OF MEMORY in contract_diag\n");
    MPI_Abort(comm, EXIT_FAILURE);
  }
  for (int lc = 0; lc < nloc; lc++)
  {
    long gj = bc_global_col(p, lc); /* global exciton index of this local column */
    const double complex *zc = &Zloc[(size_t)lc * lld];
    const double complex *ac = &AUloc[(size_t)lc * lld];
    double complex acc = 0.0 + 0.0 * I;
    for (int lr = 0; lr < mloc; lr++)
      acc += conj(zc[lr]) * ac[lr];
    part[gj] += acc;
  }
  MPI_Allreduce(part, d, (int)(2 * N), MPI_DOUBLE, MPI_SUM, comm);
  free(part);
}

/* diag(U^H h0 U) for the diagonal h0 (h0_II = eval[lumo + I/n_ho] - eval[I%n_ho]).
 * Needs no matrix: sum_I |Z(I,gj)|^2 h0_II. Real result. */
static void contract_h0_diag(const double complex *Zloc, const parallel_st *p,
                             const double *eval, index_st *ist, long N,
                             MPI_Comm comm, double *d)
{
  const int lld = p->bc_lld, mloc = p->bc_mloc, nloc = p->bc_nloc;
  const long n_ho = ist->n_holes, lidx = ist->lumo_idx;
  double *part = (double *)calloc((size_t)N, sizeof(double));
  if (!part)
  {
    fprintf(stderr, "OUT OF MEMORY in contract_h0_diag\n");
    MPI_Abort(comm, EXIT_FAILURE);
  }
  for (int lc = 0; lc < nloc; lc++)
  {
    long gj = bc_global_col(p, lc);
    const double complex *zc = &Zloc[(size_t)lc * lld];
    double acc = 0.0;
    for (int lr = 0; lr < mloc; lr++)
    {
      long Ig = bc_global_row(p, lr);
      double h0 = eval[lidx + Ig / n_ho] - eval[Ig % n_ho];
      acc += cnorm(zc[lr]) * h0;
    }
    part[gj] += acc;
  }
  MPI_Allreduce(part, d, (int)N, MPI_DOUBLE, MPI_SUM, comm);
  free(part);
}

/*****************************************************************************/
/* Driver: diagonalize the pre-distributed block-cyclic BSE matrix in place.  */
/*****************************************************************************/
void bethe_salpeter_dist(
    double complex *H_loc,
    double complex *dir_loc,
    double complex *exc_loc,
    double complex *bs_coeff,
    double *eig_vals,
    double *xton_ene,
    index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel)
{
  (void)par;
  (void)flag;
  const int rank = parallel->mpi_rank;
  const int P = parallel->mpi_size;
  const long N = ist->n_xton;
  MPI_Comm comm = MPI_COMM_WORLD;
  const int izero = 0, i1 = 1, iN = (int)N;
  int info = 0;

  bse_setup_blockcyclic(parallel, N); /* idempotent */

  /* libsci ScaLAPACK/PBLAS are internally OpenMP-threaded (link _mp variant);
   * give every rank its full thread quota so the node is saturated. */
  omp_set_num_threads((int)ist->nthreads);

  const int ctxt2d = parallel->bc_ctxt;
  const int nb = parallel->bc_nb;
  int *descA = parallel->bc_desc;

  double t0 = omp_get_wtime();
  if (rank == 0)
  {
    printf("\n  Distributed BSE eigensolve: N = %ld, %d ranks x %ld threads (ScaLAPACK pzheevd) | %s\n",
           N, P, ist->nthreads, get_time());
    fflush(stdout);
  }

  const size_t locsz = (size_t)bc_local_size(parallel);
  double complex *Z_loc = (double complex *)malloc(locsz * sizeof(double complex));
  double complex *AU_loc = (double complex *)malloc(locsz * sizeof(double complex));
  if (!Z_loc || !AU_loc)
  {
    fprintf(stderr, "OUT OF MEMORY in bethe_salpeter_dist (rank %d)\n", rank);
    MPI_Abort(comm, EXIT_FAILURE);
  }

  /* --- diagonalize H (lower triangle only -> uplo 'L') --- */
  double complex wq;
  double rwq;
  int iwq;
  int lwork = -1, lrwork = -1, liwork = -1;
  pzheevd_("V", "L", &iN, H_loc, &i1, &i1, descA, xton_ene,
           Z_loc, &i1, &i1, descA, &wq, &lwork, &rwq, &lrwork, &iwq, &liwork, &info);
  lwork = (int)(creal(wq)) + 1;
  lrwork = (int)rwq + 1;
  liwork = iwq;
  {
    long lr_safe = 4L * N + (long)nb * nb + 16;
    if ((long)lrwork < lr_safe)
      lrwork = (int)lr_safe;
    if (liwork < 7 * (int)N + 8)
      liwork = 7 * (int)N + 8;
  }
  double complex *work = (double complex *)malloc((size_t)(lwork > 0 ? lwork : 1) * sizeof(double complex));
  double *rwork = (double *)malloc((size_t)(lrwork > 0 ? lrwork : 1) * sizeof(double));
  int *iwork = (int *)malloc((size_t)(liwork > 0 ? liwork : 1) * sizeof(int));
  if (!work || !rwork || !iwork)
  {
    fprintf(stderr, "OUT OF MEMORY: pzheevd workspace (rank %d)\n", rank);
    MPI_Abort(comm, EXIT_FAILURE);
  }
  if (rank == 0)
  {
    printf("  diagonalizing %ld x %ld BSE matrix (pzheevd, uplo=L) | %s\n", N, N, get_time());
    fflush(stdout);
  }
  pzheevd_("V", "L", &iN, H_loc, &i1, &i1, descA, xton_ene,
           Z_loc, &i1, &i1, descA, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
  free(work);
  free(rwork);
  free(iwork);
  if (info != 0)
  {
    if (rank == 0)
      fprintf(stderr, "ERROR: pzheevd returned info = %d\n", info);
    MPI_Abort(comm, EXIT_FAILURE);
  }
  MPI_Bcast(xton_ene, (int)N, MPI_DOUBLE, 0, comm);

  if (rank == 0)
  {
    printf("  done diagonalizing | %s\n", format_duration(omp_get_wtime() - t0));
    printf("\n  Ground state exciton has energy = %.5f a.u. | %.5f eV\n",
           xton_ene[0], xton_ene[0] * AUTOEV);
    fflush(stdout);
  }

  /* --- gather eigenvectors to rank 0 (bs_coeff, row-major) via a 1x1 grid --- */
  int ctxt0;
  Cblacs_get(0, 0, &ctxt0);
  Cblacs_gridinit(&ctxt0, "Row", 1, 1);
  int desc0[9];
  if (rank == 0)
    descinit_(desc0, &iN, &iN, &iN, &iN, &izero, &izero, &ctxt0, &iN, &info);
  else
  {
    for (int k = 0; k < 9; k++)
      desc0[k] = 0;
    desc0[CTXT_] = -1;
  }
  double complex *src = (rank == 0) ? (double complex *)malloc((size_t)N * N * sizeof(double complex)) : NULL;
  if (rank == 0 && !src)
  {
    fprintf(stderr, "OUT OF MEMORY: eigenvector gather buffer (rank 0)\n");
    MPI_Abort(comm, EXIT_FAILURE);
  }
  pzgemr2d_(&iN, &iN, Z_loc, &i1, &i1, descA, src, &i1, &i1, desc0, &ctxt2d);
  if (rank == 0)
  {
#pragma omp parallel for
    for (long j = 0; j < N; j++)
      for (long i = 0; i < N; i++)
        bs_coeff[i * N + j] = src[i + j * N]; /* component i of exciton j */
  }

  /* --- expectation-value diagonals for exciton.dat --- */
  double complex *d_dir = (double complex *)malloc((size_t)N * sizeof(double complex));
  double complex *d_exc = (double complex *)malloc((size_t)N * sizeof(double complex));
  double *d_h0 = (double *)malloc((size_t)N * sizeof(double));
  if (!d_dir || !d_exc || !d_h0)
  {
    fprintf(stderr, "OUT OF MEMORY: expval diagonals (rank %d)\n", rank);
    MPI_Abort(comm, EXIT_FAILURE);
  }
  const double complex one = 1.0 + 0.0 * I, zero = 0.0 + 0.0 * I;

  /* <H_dir> = diag(U^H direct U); direct is Hermitian, stored lower -> pzhemm 'L'. */
  pzhemm_("L", "L", &iN, &iN, &one, dir_loc, &i1, &i1, descA,
          Z_loc, &i1, &i1, descA, &zero, AU_loc, &i1, &i1, descA);
  contract_diag(Z_loc, AU_loc, parallel, N, comm, d_dir);

  /* <H_exc> = diag(U^H exchange U). */
  pzhemm_("L", "L", &iN, &iN, &one, exc_loc, &i1, &i1, descA,
          Z_loc, &i1, &i1, descA, &zero, AU_loc, &i1, &i1, descA);
  contract_diag(Z_loc, AU_loc, parallel, N, comm, d_exc);

  /* <H_0> = diag(U^H h0 U), h0 diagonal -> direct contraction (no matrix). */
  contract_h0_diag(Z_loc, parallel, eig_vals, ist, N, comm, d_h0);

  /* --- write exciton.dat (rank 0) --- */
  if (rank == 0)
  {
    FILE *pf = fopen("exciton.dat", "w");
    if (pf == NULL)
    {
      fprintf(stderr, "ERROR: could not open exciton.dat for writing "
                      "(is the working directory writable?)\n");
    }
    else
    {
      fprintf(pf, "#n \t E_n \t <H> \t <H_dir> \t <H_exc> \t <H_0> \t E_B (eV)\n");
      for (long i = 0; i < N; i++)
      {
        double h0d = d_h0[i];
        fprintf(pf, "%ld % .12f % .12f % .12f  % .12f  % .12f  % .12f\n", i,
                xton_ene[i], xton_ene[i], creal(d_dir[i]), creal(d_exc[i]),
                h0d, (xton_ene[i] - h0d) * AUTOEV);
      }
      fclose(pf);
      printf("  wrote exciton.dat (%ld excitons) | %s\n", N, get_time());
      fflush(stdout);
    }
  }

  free(d_dir);
  free(d_exc);
  free(d_h0);
  free(Z_loc);
  free(AU_loc);
  if (src)
    free(src);

  if (rank == 0)
    Cblacs_gridexit(ctxt0);
  bse_teardown_blockcyclic(parallel);

  MPI_Barrier(comm);
}

/*****************************************************************************/
/* Self-test (compile with -DPBSE_TEST as a standalone executable).           */
/*                                                                            */
/*   cc -DOS_Linux -DFFTW_FFT -DUSE_LIBSCI -DUSE_SCALAPACK -DPBSE_TEST -O2 \    */
/*      -I$CRAY_FFTW_PREFIX/x86_64/include pbse.c aux.o -o pbse_test.x \        */
/*      -L$CRAY_LIBSCI_PREFIX_DIR/lib -lsci_nvidia_mpi_mp -lm                   */
/*   srun -n 4 -t 2 ./pbse_test.x 400        # N = 400 (default 200)            */
/*                                                                            */
/* Builds a fixed Hermitian M on rank 0, distributes its lower triangle into   */
/* H_loc and dir_loc via the router (exc_loc = 0), then diagonalizes. Checks:   */
/*   (1) distributed eigenvalues vs serial LAPACKE_zheev reference,            */
/*   (2) eigenvector residual ||M u_j - lambda_j u_j|| (gather/transpose),      */
/*   (3) exciton.dat's <H_dir> = diag(U^H M U) equals the eigenvalues          */
/*       (router + pzhemm + diagonal-contraction path).                        */
/*****************************************************************************/
#ifdef PBSE_TEST

/* deterministic Hermitian entry M(i,j): real diagonal, M(i,j)=conj(M(j,i)). */
static double complex herm_entry(long i, long j)
{
  long a = (i < j) ? i : j, b = (i < j) ? j : i; /* order-independent hash */
  unsigned long h = (unsigned long)(a * 2654435761UL) ^ (unsigned long)(b * 40503UL + 12345UL);
  h ^= h >> 13;
  h *= 0x5bd1e995UL;
  h ^= h >> 15;
  double re = ((double)(h & 0xffff) / 65535.0) - 0.5;
  double im = ((double)((h >> 16) & 0xffff) / 65535.0) - 0.5;
  if (i == j)
    return (2.0 * re) + 0.0 * I; /* real diagonal */
  double complex v = re + im * I;
  return (i < j) ? v : conj(v); /* Hermitian */
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  parallel_st parallel = {0};
  MPI_Comm_rank(MPI_COMM_WORLD, &parallel.mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &parallel.mpi_size);
  const int rank = parallel.mpi_rank;

  long N = (argc > 1) ? strtol(argv[1], NULL, 10) : 200;

  index_st ist = {0};
  ist.n_xton = N;
  ist.n_holes = 1; /* so h0_II indexing (eval[...]) stays in range; eval = 0 */
  ist.n_elecs = N;
  ist.lumo_idx = 0;
  const char *e = getenv("OMP_NUM_THREADS");
  ist.nthreads = (e && *e) ? strtol(e, NULL, 10) : 1;
  par_st par = {0};
  flag_st flag = {0};

  bse_setup_blockcyclic(&parallel, N);
  const size_t locsz = (size_t)bc_local_size(&parallel);
  double complex *H_loc = (double complex *)calloc(locsz, sizeof(double complex));
  double complex *dir_loc = (double complex *)calloc(locsz, sizeof(double complex));
  double complex *exc_loc = (double complex *)calloc(locsz, sizeof(double complex));
  double *eig_vals = (double *)calloc((size_t)N, sizeof(double)); /* h0 = 0 */
  double *xton_ene = (double *)malloc((size_t)N * sizeof(double));

  /* rank 0 owns M and routes its lower triangle to the block-cyclic owners */
  bc_router Rh, Rd;
  router_init(&Rh, parallel.mpi_size);
  router_init(&Rd, parallel.mpi_size);
  if (rank == 0)
    for (long i = 0; i < N; i++)
      for (long j = 0; j <= i; j++) /* lower triangle */
      {
        double complex m = herm_entry(i, j);
        router_push(&Rh, &parallel, i, j, m);
        router_push(&Rd, &parallel, i, j, m);
      }
  router_flush(&Rh, H_loc, MPI_COMM_WORLD);
  router_flush(&Rd, dir_loc, MPI_COMM_WORLD);
  router_free(&Rh);
  router_free(&Rd);

  double complex *bs_coeff = (rank == 0) ? (double complex *)malloc((size_t)N * N * sizeof(double complex)) : NULL;
  double *refval = NULL, *Mref = NULL;
  if (rank == 0)
  {
    Mref = (double *)malloc((size_t)2 * N * N * sizeof(double));
    refval = (double *)malloc((size_t)N * sizeof(double));
    for (long i = 0; i < N; i++)
      for (long j = 0; j < N; j++)
        ((double complex *)Mref)[i * N + j] = herm_entry(i, j);
    LAPACKE_zheev(LAPACK_ROW_MAJOR, 'N', 'U', (lapack_int)N,
                  (lapack_complex_double *)Mref, (lapack_int)N, refval);
  }

  bethe_salpeter_dist(H_loc, dir_loc, exc_loc, bs_coeff, eig_vals, xton_ene,
                      &ist, &par, &flag, &parallel);

  if (rank == 0)
  {
    double emax = 0.0, scale = 1.0;
    for (long i = 0; i < N; i++)
      if (fabs(refval[i]) > scale)
        scale = fabs(refval[i]);
    for (long i = 0; i < N; i++)
    {
      double d = fabs(xton_ene[i] - refval[i]);
      if (d > emax)
        emax = d;
    }
    double rmax = 0.0;
    for (long j = 0; j < N; j += (N / 5 + 1))
    {
      double rr = 0.0;
      for (long i = 0; i < N; i++)
      {
        double complex acc = 0.0;
        for (long k = 0; k < N; k++)
          acc += herm_entry(i, k) * bs_coeff[k * N + j];
        acc -= xton_ene[j] * bs_coeff[i * N + j];
        rr += creal(acc) * creal(acc) + cimag(acc) * cimag(acc);
      }
      rr = sqrt(rr);
      if (rr > rmax)
        rmax = rr;
    }
    printf("\n[pbse-test] N=%ld  ranks=%d\n", N, parallel.mpi_size);
    printf("[pbse-test] max |eval_dist - eval_ref|      = %.3e (scale %.3g)\n", emax, scale);
    printf("[pbse-test] max eigenvector residual        = %.3e\n", rmax);
    printf("[pbse-test] (exciton.dat <H_dir> should match column E_n above)\n");
    printf("[pbse-test] %s\n",
           (emax < 1e-8 * scale && rmax < 1e-7 * scale) ? "PASS" : "FAIL");
    free(bs_coeff);
    free(Mref);
    free(refval);
  }
  free(H_loc);
  free(dir_loc);
  free(exc_loc);
  free(eig_vals);
  free(xton_ene);
  MPI_Finalize();
  return 0;
}

#endif

/*****************************************************************************/
