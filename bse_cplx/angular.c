/****************************************************************************/

#include "angular.h"

/****************************************************************************/

void calc_qp_spin_mtrx(
    double complex *restrict psi_qp,
    xyz_st *restrict s_mom,
    grid_st *grid,
    index_st *ist,
    par_st *par,
    parallel_st *parallel)
{

  /*******************************************************************
   * This function computes the value of the spin projection operator *
   * between the electron and hole quasiparticle states. The spin     *
   * projection operators are defined by the Pauli matrices as:       *
   * Sx = 1/2 (|up><down| + |down><up|)                               *
   * Sy = 1/2 (-i|up><down| + i|down><up|)                            *
   * Sz = 1/2 (|up><up| - |down><down|)                               *
   * and we compute the matrix elements <i|Sx|j> as                   *
   * sum_r psi_i(r, s) * Sx * psi_j(r, s) * dr                        *
   * inputs:                                                          *
   *  [s_mom] array to hold components of the spin projection elems   *
   *  [psi_qp] array holding all qp_basis states                      *
   *  [grid] grid_st instance holding values of all grid points       *
   *  [ist] ptr to counters, indices, and lengths                     *
   *  [par] ptr to par_st holding VBmin, VBmax... params              *
   *  [flag] ptr to flag_st holding job flags                         *
   *  [planfw] FFTW3 plan for executing 3D forward DFT                *
   *  [planfw] FFTW3 plan for executing 3D backwards DFT              *
   *  [fftwpsi] location to store outcome of Fourier transform        *
   * outputs: void                                                    *
   ********************************************************************/

  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  FILE *pfx;
  FILE *pfy;
  FILE *pfz;

  long i;
  long j;
  long a;
  long b;
  long i_st;
  long j_st;
  long a_st;
  long b_st;
  long jg;
  long jsg;
  long idx;

  const long n_ho = ist->n_holes;
  const long n_el = ist->n_elecs;
  const long lidx = ist->lumo_idx;
  const long ngrid = ist->ngrid;
  const long nspngr = ist->nspinngrid;

  const double dv = grid->dv;

  const int mpir = parallel->mpi_rank;
  const int mpi_size = parallel->mpi_size;
  const long mat_size = sqr(n_ho) + sqr(n_el);

  // Each rank fills a disjoint set of rows; complete with Allreduce(SUM) below.
  memset(s_mom, 0, (size_t)mat_size * sizeof(xyz_st));

  if (mpir == 0)
    write_state_dat(psi_qp, 4 * nspngr, "psi_qp_all.dat");
/************************************************************/
/*******************    CALC HOLE SPINS   *******************/
/************************************************************/

// Distribute hole rows round-robin over MPI ranks; OpenMP threads within.
#pragma omp parallel for private(i, j, i_st, j_st, jg, jsg)
  for (i = mpir; i < n_ho; i += mpi_size)
  {
    for (j = 0; j < n_ho; j++)
    {
      // nvtxRangePushA("loop ij");
      i_st = i * nspngr;
      j_st = j * nspngr;

      double complex sx;
      double complex sy;
      double complex sz;

      sx = sy = sz = 0.0 + 0.0 * I;

      jsg = ngrid;
#pragma omp simd safelen(8) aligned(psi_qp : BYTE_BOUNDARY) reduction(+ : sx, sy, sz)
      for (jg = 0; jg < ngrid; jg++)
      {
        // Handle indexing

        // S_x component
        // <j|r,dn> * <r,up|i>
        sx += conjmul(psi_qp[j_st + jsg], psi_qp[i_st + jg]);
        // printf("\n%lu\n", jg);
        // printf("psi[jdn] = %g %g psi[iup] = %g %g\n", psi_qp[j_st + jsg], psi_qp[i_st + jg] );
        // printf("conjmul  = %g %g\n", conjmul(psi_qp[j_st + jsg], psi_qp[i_st + jg]) );
        // printf("sx       = %g %g\n", sx);
        // <j|r,up> * <r,dn|i>
        sx += conjmul(psi_qp[j_st + jg], psi_qp[i_st + jsg]);
        // printf("psi[jup] = %g %g psi[idn] = %g %g\n", psi_qp[j_st + jg], psi_qp[i_st + jsg] );
        // printf("conjmul  = %g %g\n", conjmul(psi_qp[j_st + jg], psi_qp[i_st + jsg]) );
        // printf("sx       = %g %g\n", sx);

        // S_y component
        // i*<j|r,dn> * <r,up|i>
        sy += I * conjmul(psi_qp[j_st + jsg], psi_qp[i_st + jg]);
        // -i*<j|r,up> * <r,dn|i>
        sy -= I * conjmul(psi_qp[j_st + jg], psi_qp[i_st + jsg]);

        // S_z component
        // <j|r,up> * <r,up|i>
        sz += conjmul(psi_qp[j_st + jg], psi_qp[i_st + jg]);

        // - <j|r,dn> * <r,dn|i>
        sz -= conjmul(psi_qp[j_st + jsg], psi_qp[i_st + jsg]);

        jsg++;
      }

      // multiply all by (1/2)*dV and minus-conjugate bc holes are time reversed
      s_mom[i * n_ho + j].x = -0.5 * conj(sx) * dv;
      s_mom[i * n_ho + j].y = -0.5 * conj(sy) * dv;
      s_mom[i * n_ho + j].z = -0.5 * conj(sz) * dv;

      // nvtxRangePop();
    }
  }
  // nvtxRangePop();

/************************************************************/
/*******************    CALC ELEC SPINS   *******************/
/************************************************************/

// Distribute electron rows round-robin over MPI ranks; OpenMP threads within.
#pragma omp parallel for private(a, b, a_st, b_st, idx, jg, jsg)
  for (a = lidx + mpir; a < lidx + n_el; a += mpi_size)
  {
    for (b = lidx; b < lidx + n_el; b++)
    {
      // nvtxRangePop();
      a_st = a * nspngr;
      b_st = b * nspngr;

      double complex sx;
      double complex sy;
      double complex sz;

      sx = sy = sz = 0.0 + 0.0 * I;

#pragma omp simd safelen(8) aligned(psi_qp : BYTE_BOUNDARY) reduction(+ : sx, sy, sz)
      for (jg = 0; jg < ngrid; jg++)
      {
        // Handle indexing
        jsg = jg + ngrid;

        // Spin x part
        //<b|r,dn> * <r,up|a>
        sx += conjmul(psi_qp[b_st + jsg], psi_qp[a_st + jg]);
        //<b|r,up> * <r,dn|a>
        sx += conjmul(psi_qp[b_st + jg], psi_qp[a_st + jsg]);

        // Spin y part
        // i*<b|r,dn> * <r,up|a>
        sy += I * conjmul(psi_qp[b_st + jsg], psi_qp[a_st + jg]);
        //-i*<b|r,up> * <r,dn|a>
        sy -= I * conjmul(psi_qp[b_st + jg], psi_qp[a_st + jsg]);

        // Spin z part
        //<b|r,up> * <r,up|a>
        sz += conjmul(psi_qp[b_st + jg], psi_qp[a_st + jg]);
        //<b|r,dn> * <r,dn|a>
        sz -= conjmul(psi_qp[b_st + jsg], psi_qp[a_st + jsg]);
      }

      idx = sqr(n_ho) + (a - lidx) * n_el + (b - lidx);

      // Divide all by 2 and normalize
      s_mom[idx].x = 0.5 * sx * dv;
      s_mom[idx].y = 0.5 * sy * dv;
      s_mom[idx].z = 0.5 * sz * dv;
      // nvtxRangePop();
    }
  }
  // nvtxRangePop();

  // Complete the spin matrix across ranks (disjoint rows summed).
  MPI_Allreduce(MPI_IN_PLACE, s_mom, (int)(6 * mat_size), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  /************************************************************/
  /*******************      PRINT OUTPUT    *******************/
  /************************************************************/

  if (mpir == 0)
  {
    pfx = fopen("sx.dat", "w");
    pfy = fopen("sy.dat", "w");
    pfz = fopen("sz.dat", "w");

    for (i = 0; i < n_ho; i++)
      for (j = 0; j < n_ho; j++)
      {
        idx = i * n_ho + j;
        fprintf(pfx, "%ld %ld % .10g % .10g\n", i, j, creal(s_mom[idx].x), cimag(s_mom[idx].x));
        fprintf(pfy, "%ld %ld % .10g % .10g\n", i, j, creal(s_mom[idx].y), cimag(s_mom[idx].y));
        fprintf(pfz, "%ld %ld % .10g % .10g\n", i, j, creal(s_mom[idx].z), cimag(s_mom[idx].z));
      }
    for (a = lidx; a < lidx + n_el; a++)
      for (b = lidx; b < lidx + n_el; b++)
      {
        idx = sqr(n_ho) + (a - lidx) * n_el + (b - lidx);
        fprintf(pfx, "%ld %ld % .10g % .10g\n", a, b, creal(s_mom[idx].x), cimag(s_mom[idx].x));
        fprintf(pfy, "%ld %ld % .10g % .10g\n", a, b, creal(s_mom[idx].y), cimag(s_mom[idx].y));
        fprintf(pfz, "%ld %ld % .10g % .10g\n", a, b, creal(s_mom[idx].z), cimag(s_mom[idx].z));
      }
    fclose(pfx);
    fclose(pfy);
    fclose(pfz);
  }

  return;
}

// /****************************************************************************/

void calc_qp_ang_mom_mtrx(
    double complex *restrict psi_qp,
    xyz_st *restrict l_mom,
    double complex *restrict l2_mom,
    double complex *restrict ldots,
    grid_st *grid,
    index_st *ist,
    par_st *par,
    parallel_st *parallel)
{

  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  long i;
  long j;
  long a;
  long b;
  long i_st;
  long j_st;
  long a_st;
  long b_st;
  long jg;
  long jsg;
  long idx;

  const long n_ho = ist->n_holes;
  const long n_el = ist->n_elecs;
  const long lidx = ist->lumo_idx;
  const long ngrid = ist->ngrid;
  const long nspngr = ist->nspinngrid;

  const double dv = grid->dv;

  double *gx;
  double *gy;
  double *gz;
  double *g_vecs;

  double complex *Lxpsi;
  double complex *Lypsi;
  double complex *Lzpsi;
  double complex *Lx2psi;
  double complex *Ly2psi;
  double complex *Lz2psi;
  double complex *temp1;
  double complex *temp2;

  fftw_plan_loc planfw;
  fftw_plan_loc planbw;
  fftw_complex *fftwpsi;

  /************************************************************/
  /*******************    INITIALIZE FFT    *******************/
  /************************************************************/

  // Note: parallelize FFT to achieve fast performance;
  // do not parallelize over states! - Daniel C 3.13.25
  // Best strategy:
  //    - Use hybrid MPI/OpenMP to distribute initial states i/a
  //    - Compute FFT with threads
  //    - Perform integrals on GPU

  fftw_init_threads();
  fftw_plan_with_nthreads(ist->nthreads);

  // Allocate memory for multithreaded FFT
  fftwpsi = fftw_malloc(ist->ngrid * sizeof(fftw_complex));
  planfw = fftw_plan_dft_3d(grid->nz, grid->ny, grid->nx, fftwpsi, fftwpsi, FFTW_FORWARD, 0);
  planbw = fftw_plan_dft_3d(grid->nz, grid->ny, grid->nx, fftwpsi, fftwpsi, FFTW_BACKWARD, 0);

  /************************************************************/
  /*******************    INIT G vecs   *******************/
  /************************************************************/

  // G vectors for derivatives in k space
  ALLOCATE(&gx, grid->nx, "gx");
  ALLOCATE(&gy, grid->ny, "gy");
  ALLOCATE(&gz, grid->nz, "gz");
  ALLOCATE(&g_vecs, 3 * ngrid, "g_vecs");

  // Initializing the G vectors
  init_g_vecs(g_vecs, gx, gy, gz, grid, ist, par);

  /************************************************************/
  /*******************    ALLOC L OP MEM    *******************/
  /************************************************************/

  ALLOCATE(&Lxpsi, nspngr, "Lxpsi");
  ALLOCATE(&Lypsi, nspngr, "Lypsi");
  ALLOCATE(&Lzpsi, nspngr, "Lzpsi");

  ALLOCATE(&Lx2psi, nspngr, "Lx2psi");
  ALLOCATE(&Ly2psi, nspngr, "Ly2psi");
  ALLOCATE(&Lz2psi, nspngr, "Lz2psi");

  ALLOCATE(&temp1, ngrid, "temp1");
  ALLOCATE(&temp2, ngrid, "temp2");

  const int mpir = parallel->mpi_rank;
  const int mpi_size = parallel->mpi_size;
  const long mat_size = sqr(n_ho) + sqr(n_el);

  FILE *pfx = NULL, *pfy = NULL, *pfz = NULL, *pfsqr = NULL, *pfls = NULL;

  // Each rank fills a disjoint set of upper-triangle rows; the arrays are
  // completed with one Allreduce(SUM) below, then the lower triangle is filled.
  memset(l_mom, 0, (size_t)mat_size * sizeof(xyz_st));
  memset(l2_mom, 0, (size_t)mat_size * sizeof(double complex));
  memset(ldots, 0, (size_t)mat_size * sizeof(double complex));

  /************************************************************/
  /******************    CALC HOLE <L>     ********************/
  /************************************************************/

  if (mpir == 0)
  {
    printf("Hole States:\n");
    printf("i       j         Lx.re           Ly.re           Lz.re          L^2.re           L.S.re\n");
    fflush(0);
  }

  // nvtxRangePushA("Calc qp hole L elems")
  // Distribute hole states round-robin over MPI ranks (each does its L|i> FFTs).
  for (i = mpir; i < n_ho; i += mpi_size)
  {
    // nvtxRangePushA("loop over i");
    i_st = i * nspngr;

    // nvtxRangePushA("Compute L by FFTs");
    // These computations require FFT with small (max: 500^3) grids -> performed on CPU
    // 1) Compute Lx|i>, Ly|i>, Lz|i>
    // spin up part
    l_operator(&Lxpsi[0], &Lypsi[0], &Lzpsi[0], &psi_qp[i_st], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
    // spin dn part
    l_operator(&Lxpsi[ngrid], &Lypsi[ngrid], &Lzpsi[ngrid], &psi_qp[i_st + ngrid], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);

    // Compute Lx^2|i>
    // spin up part
    l_operator(&Lx2psi[0], &temp1[0], &temp2[0], &Lxpsi[0], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
    // spin dn part
    l_operator(&Lx2psi[ngrid], &temp1[0], &temp2[0], &Lxpsi[ngrid], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);

    // Compute Ly^2|i>
    // spin up part
    l_operator(&temp1[0], &Ly2psi[0], &temp2[0], &Lypsi[0], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
    // spin dn part
    l_operator(&temp1[0], &Ly2psi[ngrid], &temp2[0], &Lypsi[ngrid], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);

    // Compute Lz^2|i>
    // spin up part
    l_operator(&temp1[0], &temp2[0], &Lz2psi[0], &Lzpsi[0], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
    // spin dn part
    l_operator(&temp1[0], &temp2[0], &Lz2psi[ngrid], &Lzpsi[ngrid], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);

// nvtxRangePop();
#pragma omp parallel for private(j_st, jg, jsg)
    for (j = i; j < n_ho; j++)
    {
      // nvtxRangePushA("loop over a");
      j_st = j * nspngr;

      double complex lx;
      double complex ly;
      double complex lz;
      double complex lsqr;
      double complex ls;

      lx = ly = lz = lsqr = ls = 0.0 + 0.0 * I;

// Perform integration
#pragma omp simd
      for (jg = 0; jg < nspngr; jg++)
      {
        // Handle up spin

        // Holes are like time-reversed electrons, so -^*
        // -<j|L_x|i>^*
        lx -= conjmul(psi_qp[j_st + jg], Lxpsi[jg]);
        // -<j|L_y|i>^*
        ly -= conjmul(psi_qp[j_st + jg], Lypsi[jg]);
        // -<j|L_z|i>^*
        lz -= conjmul(psi_qp[j_st + jg], Lzpsi[jg]);

        // <j|L_x^2|i>^*
        lsqr += conjmul(psi_qp[j_st + jg], Lx2psi[jg]);
        // <j|L_y^2|i>^*
        lsqr += conjmul(psi_qp[j_st + jg], Ly2psi[jg]);
        // <j|L_z^2|i>^*
        lsqr += conjmul(psi_qp[j_st + jg], Lz2psi[jg]);
      }

      // normalize and store in allocated array
      l_mom[i * n_ho + j].x = conj(lx) * dv;
      l_mom[i * n_ho + j].y = conj(ly) * dv;
      l_mom[i * n_ho + j].z = conj(lz) * dv;
      l2_mom[i * n_ho + j] = conj(lsqr) * dv;

      jsg = ngrid;
      for (jg = 0; jg < ngrid; jg++, jsg++)
      {
        // Compute ldots, the dot product of orbital and spin angular momentum

        // sx.lx
        //<j|r,dn> * <r,up|Lxi>
        ls += conjmul(psi_qp[j_st + jsg], Lxpsi[jg]);
        //<j|r,up> * <r,dn|Lxi>
        ls += conjmul(psi_qp[j_st + jg], Lxpsi[jsg]);

        // sy.ly
        // i*<j|r,dn> * <r,up|Lyi>
        ls += I * conjmul(psi_qp[j_st + jsg], Lypsi[jg]);
        //-i*<j|r,up> * <r,dn|Lyi>
        ls -= I * conjmul(psi_qp[j_st + jg], Lypsi[jsg]);

        // sz.lz
        //<j|r,up> * <r,up|Lzi>
        ls += conjmul(psi_qp[j_st + jg], Lzpsi[jg]);
        //<j|r,dn> * <r,dn|Lzi>
        ls -= conjmul(psi_qp[j_st + jsg], Lzpsi[jsg]);
      }
      // normalize and complex conjugate ldots
      ldots[i * n_ho + j] = 0.5 * conj(ls) * dv;
      // nvtxRangePop();
    }
    // nvtxRangePop();
  }
  // nvtxRangePop();

  /************************************************************/
  /******************    CALC ELEC <L>     ********************/
  /************************************************************/

  if (mpir == 0)
  {
    printf("Electron States:\n");
    printf("a       b         Lx.re           Ly.re           Lz.re          L^2.re           L.S.re\n");
    fflush(0);
  }
  // nvtxRangePushA("Calc qp elec L elems")
  // Distribute electron states round-robin over MPI ranks.
  for (a = lidx + mpir; a < lidx + n_el; a += mpi_size)
  {
    // nvtxRangePushA("loop over a");
    a_st = a * nspngr;

    // nvtxRangePushA("Compute L by FFTs");
    // spin up part
    l_operator(Lxpsi, Lypsi, Lzpsi, &psi_qp[a_st], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
    // spin dn part
    l_operator(&Lxpsi[ngrid], &Lypsi[ngrid], &Lzpsi[ngrid], &psi_qp[a_st + ngrid], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);

    // Lxsqr parts
    // spin up part
    l_operator(Lx2psi, &temp1[0], &temp2[0], Lxpsi, g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
    // spin dn part
    l_operator(&Lx2psi[ngrid], &temp1[0], &temp2[0], &Lxpsi[ngrid], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);

    // Lysqr parts
    // spin up part
    l_operator(&temp1[0], Ly2psi, &temp2[0], Lypsi, g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
    // spin dn part
    l_operator(&temp1[0], &Ly2psi[ngrid], &temp2[0], &Lypsi[ngrid], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);

    // Lzsqr parts
    // spin up part
    l_operator(&temp1[0], &temp2[0], &Lz2psi[0], &Lzpsi[0], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
    // spin dn part
    l_operator(&temp1[0], &temp2[0], &Lz2psi[ngrid], &Lzpsi[ngrid], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
    // nvtxRangePop();

#pragma omp parallel for private(b, b_st, jg, jsg, idx)
    for (b = a; b < lidx + n_el; b++)
    {
      // nvtxRangePushA("loop over b");
      b_st = b * nspngr;

      // Declare variables on the GPU
      double complex lx;
      double complex ly;
      double complex lz;
      double complex lsqr;
      double complex ls;

      lx = ly = lz = lsqr = ls = 0.0 + 0.0 * I;

      idx = sqr(n_ho) + (a - lidx) * n_el + (b - lidx);

#pragma omp simd reduction(+ : lx, ly, lz, lsqr)
      for (jg = 0; jg < nspngr; jg++)
      {
        // Handle up spin
        // <b|L_x|a>
        lx += conjmul(psi_qp[b_st + jg], Lxpsi[jg]);
        // <b|L_y|a>
        ly += conjmul(psi_qp[b_st + jg], Lypsi[jg]);
        // <b|L_z|a>
        lz += conjmul(psi_qp[b_st + jg], Lzpsi[jg]);

        // <b|L_x^2|a>
        lsqr += conjmul(psi_qp[b_st + jg], Lx2psi[jg]);
        // <b|L_y^2|a>
        lsqr += conjmul(psi_qp[b_st + jg], Ly2psi[jg]);
        // <b|L_z^2|a>
        lsqr += conjmul(psi_qp[b_st + jg], Lz2psi[jg]);
      }

      // Normalize and store in alloc'd arrays
      l_mom[idx].x = lx * dv;
      l_mom[idx].y = ly * dv;
      l_mom[idx].z = lz * dv;
      l2_mom[idx] = lsqr * dv;

      jsg = ngrid;
      for (jg = 0; jg < ngrid; jg++)
      {
        // Compute ldots for the electron states

        // s_mom.xlx
        //<b|r,dn> * <r,up|Lxa>
        ls += conjmul(psi_qp[b_st + jsg], Lxpsi[jg]);
        //<b|r,up> * <r,dn|Lxa>
        ls += conjmul(psi_qp[b_st + jg], Lxpsi[jsg]);

        // s_mom.yly
        // i*<j|r,dn> * <r,up|Lyi>
        ls += I * conjmul(psi_qp[b_st + jsg], Lypsi[jg]);
        //-i*<j|r,up> * <r,dn|Lyi>
        ls -= I * conjmul(psi_qp[b_st + jg], Lypsi[jsg]);

        // s_mom.zlz
        //<j|r,up> * <r,up|Lzi>
        ls += conjmul(psi_qp[b_st + jg], Lzpsi[jg]);
        //<j|r,dn> * <r,dn|Lzi>
        ls -= conjmul(psi_qp[b_st + jsg], Lzpsi[jsg]);
      }

      ldots[idx] = 0.5 * ls * dv;
      // nvtxRangePop();
    }
    // nvtxRangePop();
  }
  // nvtxRangePop();

  /************************************************************/
  /*******************   REDUCE OVER RANKS  *******************/
  /************************************************************/

  // Complete the upper-triangle matrices across ranks (disjoint rows summed).
  MPI_Allreduce(MPI_IN_PLACE, l_mom, (int)(6 * mat_size), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, l2_mom, (int)(2 * mat_size), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, ldots, (int)(2 * mat_size), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  /************************************************************/
  /*******************   FILL LOWER TRIANG  *******************/
  /************************************************************/

  // Holes, then electrons: lower triangle = conj(upper triangle). Cheap; every
  // rank does it so the full matrices are available for the exciton pass.
  for (i = 0; i < n_ho; i++)
    for (j = i + 1; j < n_ho; j++)
    {
      long lt = j * n_ho + i, ut = i * n_ho + j;
      l_mom[lt].x = conj(l_mom[ut].x);
      l_mom[lt].y = conj(l_mom[ut].y);
      l_mom[lt].z = conj(l_mom[ut].z);
      l2_mom[lt] = conj(l2_mom[ut]);
      ldots[lt] = conj(ldots[ut]);
    }
  for (a = lidx; a < lidx + n_el; a++)
    for (b = a + 1; b < lidx + n_el; b++)
    {
      long lt = sqr(n_ho) + (b - lidx) * n_el + (a - lidx);
      long ut = sqr(n_ho) + (a - lidx) * n_el + (b - lidx);
      l_mom[lt].x = conj(l_mom[ut].x);
      l_mom[lt].y = conj(l_mom[ut].y);
      l_mom[lt].z = conj(l_mom[ut].z);
      l2_mom[lt] = conj(l2_mom[ut]);
      ldots[lt] = conj(ldots[ut]);
    }

  /************************************************************/
  /*******************      PRINT OUTPUT    *******************/
  /************************************************************/

  if (mpir == 0)
  {
    pfx = fopen("lx.dat", "w");
    pfy = fopen("ly.dat", "w");
    pfz = fopen("lz.dat", "w");
    pfsqr = fopen("lsqr.dat", "w");
    pfls = fopen("ls.dat", "w");

    for (i = 0; i < n_ho; i++)
      for (j = 0; j < n_ho; j++)
      {
        idx = i * n_ho + j;
        fprintf(pfx, "%ld\t%ld\t%lf\t%lf\n", i, j, creal(l_mom[idx].x), cimag(l_mom[idx].x));
        fprintf(pfy, "%ld\t%ld\t%lf\t%lf\n", i, j, creal(l_mom[idx].y), cimag(l_mom[idx].y));
        fprintf(pfz, "%ld\t%ld\t%lf\t%lf\n", i, j, creal(l_mom[idx].z), cimag(l_mom[idx].z));
        fprintf(pfsqr, "%ld\t%ld\t%lf\t%lf\n", i, j, creal(l2_mom[idx]), cimag(l2_mom[idx]));
        fprintf(pfls, "%ld\t%ld\t%lf\t%lf\n", i, j, creal(ldots[idx]), cimag(ldots[idx]));
        if (i == j)
          printf("%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\n", i, j,
                 creal(l_mom[idx].x), creal(l_mom[idx].y), creal(l_mom[idx].z),
                 creal(l2_mom[idx]), creal(ldots[idx]));
      }
    for (a = lidx; a < lidx + n_el; a++)
      for (b = lidx; b < lidx + n_el; b++)
      {
        idx = sqr(n_ho) + (a - lidx) * n_el + (b - lidx);
        fprintf(pfx, "%ld\t%ld\t%lf\t%lf\n", a, b, creal(l_mom[idx].x), cimag(l_mom[idx].x));
        fprintf(pfy, "%ld\t%ld\t%lf\t%lf\n", a, b, creal(l_mom[idx].y), cimag(l_mom[idx].y));
        fprintf(pfz, "%ld\t%ld\t%lf\t%lf\n", a, b, creal(l_mom[idx].z), cimag(l_mom[idx].z));
        fprintf(pfsqr, "%ld\t%ld\t%lf\t%lf\n", a, b, creal(l2_mom[idx]), cimag(l2_mom[idx]));
        fprintf(pfls, "%ld\t%ld\t%lf\t%lf\n", a, b, creal(ldots[idx]), cimag(ldots[idx]));
        if (a == b)
          printf("%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\n", a, b,
                 creal(l_mom[idx].x), creal(l_mom[idx].y), creal(l_mom[idx].z),
                 creal(l2_mom[idx]), creal(ldots[idx]));
      }
    fclose(pfx);
    fclose(pfy);
    fclose(pfz);
    fclose(pfsqr);
    fclose(pfls);
  }
  fflush(0);

  /************************************************************/
  /*******************   FREE DYNAMIC MEM   *******************/
  /************************************************************/

  free(gx);
  free(gy);
  free(gz);
  free(g_vecs);

  free(Lxpsi);
  free(Lypsi);
  free(Lzpsi);
  free(Lx2psi);
  free(Ly2psi);
  free(Lz2psi);
  free(temp1);
  free(temp2);

  fftw_free(fftwpsi);
  fftw_destroy_plan(planfw);
  fftw_destroy_plan(planbw);

  return;
}

/******************************************************************************************/

void calc_xton_spin_mtrx(
    double complex *bs_coeff,
    xyz_st *s_mom,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel)
{

  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  FILE *pf;

  const int mpir = parallel->mpi_rank;

  unsigned long i;
  unsigned long j;
  unsigned long a;
  unsigned long b;
  unsigned long ibs;
  unsigned long jbs;

  unsigned long n;
  unsigned long index;
  unsigned long indexba;
  unsigned long indexji;
  unsigned long *listibs;

  const long n_xton = ist->n_xton;
  const long n_ho = ist->n_holes;
  const long n_el = ist->n_elecs;
  const long lidx = ist->lumo_idx;

  const double dv = par->dv;
  (void)dv;

  const int mpi_size = parallel->mpi_size;

  xyz_st *spin_arr;
  double complex *spintot_arr;

  ALLOCATE(&listibs, n_xton, "listibs in angular spins");
  ALLOCATE(&spin_arr, n_xton, "spin_arr");
  ALLOCATE(&spintot_arr, n_xton, "spintot_arr");
  memset(spin_arr, 0, (size_t)n_xton * sizeof(xyz_st));
  memset(spintot_arr, 0, (size_t)n_xton * sizeof(double complex));

  for (ibs = 0, a = lidx; a < lidx + n_el; a++)
  {
    for (i = 0; i < n_ho; i++, ibs++)
    {
      listibs[(a - lidx) * n_ho + i] = ibs;
    }
  }

  /************************************************************/
  /*******************    CALC XTON SPIN    *******************/
  /************************************************************/

  // Distribute excitons round-robin over MPI ranks; OpenMP threads within.
#pragma omp parallel for private(a, b, i, j, ibs, jbs, index, indexba, indexji)
  for (n = mpir; n < n_xton; n += mpi_size)
  {
    xyz_st spin, stmp;
    double complex spintot, tmp;
    spin.x = spin.y = spin.z = 0.0 + 0.0 * I;
    spintot = 0.0 + 0.0 * I;

    for (a = lidx; a < lidx + n_el; a++)
    {
      for (i = 0; i < n_ho; i++)
      {
        ibs = listibs[(a - lidx) * n_ho + i];

        stmp.x = stmp.y = stmp.z = 0.0 + 0.0 * I;

        // sum over b
        for (b = lidx; b < lidx + n_el; b++)
        {
          jbs = listibs[(b - lidx) * n_ho + i];
          index = sqr(n_ho) + (a - lidx) * n_el + (b - lidx);
          // index = sqr(n_ho) + (b - lidx) * n_el + (a - lidx);

          // c_bi^* * <b| Sx |a>
          stmp.x += conjmul(bs_coeff[jbs * n_xton + n], s_mom[index].x);

          // c_bi^* * <b| Sy |a>
          stmp.y += conjmul(bs_coeff[jbs * n_xton + n], s_mom[index].y);

          // c_bi^* * <b| Sz |a>
          stmp.z += conjmul(bs_coeff[jbs * n_xton + n], s_mom[index].z);
        }

        // sum over j
        for (j = 0; j < n_ho; j++)
        {
          jbs = listibs[(a - lidx) * n_ho + j];
          index = i * n_ho + j;
          // index = j * n_ho + i;

          // c_aj^* * <j| Sx |i>
          stmp.x += conjmul(bs_coeff[jbs * n_xton + n], s_mom[index].x);

          // c_aj^* * <j| Sy |i>
          stmp.y += conjmul(bs_coeff[jbs * n_xton + n], s_mom[index].y);

          // c_aj^* * <j| Sz |i>
          stmp.z += conjmul(bs_coeff[jbs * n_xton + n], s_mom[index].z);
        }

        // multiply by the c_ai coeff
        spin.x += bs_coeff[ibs * n_xton + n] * stmp.x;
        spin.y += bs_coeff[ibs * n_xton + n] * stmp.y;
        spin.z += bs_coeff[ibs * n_xton + n] * stmp.z;

        for (b = lidx; b < lidx + n_el; b++)
        {
          for (j = 0; j < n_ho; j++)
          {
            indexba = sqr(n_ho) + (a - lidx) * n_el + (b - lidx);
            // indexba = sqr(n_ho) + (b - lidx) * n_el + (a - lidx);
            indexji = i * n_ho + j;
            // indexji = j * n_ho + i;

            jbs = listibs[(b - lidx) * n_ho + j];

            tmp = conjmul(bs_coeff[jbs * n_xton + n], bs_coeff[ibs * n_xton + n]);

            stmp.x = ((s_mom[indexba].x * s_mom[indexji].x) +
                      (s_mom[indexba].y * s_mom[indexji].y) +
                      (s_mom[indexba].z * s_mom[indexji].z));

            spintot += tmp * stmp.x;
          } // end of j
        } // end of b
      } // end of i
    } // end of a

    spin_arr[n] = spin;
    spintot_arr[n] = spintot;
  } // end of n

  // Complete the per-exciton results across ranks (disjoint excitons summed).
  MPI_Allreduce(MPI_IN_PLACE, spin_arr, (int)(6 * n_xton), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, spintot_arr, (int)(2 * n_xton), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (mpir == 0)
  {
    pf = fopen("spins.dat", "w");
    for (n = 0; n < n_xton; n++)
    {
      fprintf(pf, "%ld\t%-10.5lg\t%-10.5lg\t%-10.5lg\t%-10.5lg\t%-10.5lg\t%-10.5lg\t", n,
              creal(spin_arr[n].x), cimag(spin_arr[n].x), creal(spin_arr[n].y), cimag(spin_arr[n].y),
              creal(spin_arr[n].z), cimag(spin_arr[n].z));
      fprintf(pf, "%-10.5lg\t (%-10.5lg)\n", 1.5 + 2.0 * creal(spintot_arr[n]), 2.0 * cimag(spintot_arr[n]));
    }
    fclose(pf);
  }

  free(listibs);
  free(spin_arr);
  free(spintot_arr);

  return;
}

/******************************************************************************************/

void calc_xton_ang_mom_mtrx(
    double complex *bs_coeff,
    xyz_st *s_mom,
    xyz_st *l_mom,
    double complex *l2_mom,
    double complex *ldots,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel)
{

  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  FILE *pf;
  FILE *lspf;

  const int mpir = parallel->mpi_rank;

  unsigned long i;
  unsigned long j;
  unsigned long a;
  unsigned long b;
  unsigned long ibs;
  unsigned long jbs;

  unsigned long n;
  unsigned long index;
  unsigned long indexba;
  unsigned long indexji;

  unsigned long *listibs;

  const long n_xton = ist->n_xton;
  const long n_ho = ist->n_holes;
  const long n_el = ist->n_elecs;
  const long lidx = ist->lumo_idx;

  const double dv = par->dv;
  (void)dv;

  double complex orbittot;
  double complex lstot;
  double complex lstmp;
  double complex ctmp;

  xyz_st orbit;
  xyz_st ltmp;

  const int mpi_size = parallel->mpi_size;

  xyz_st *orbit_arr;
  double complex *orbittot_arr;
  double complex *lstot_arr;

  ALLOCATE(&listibs, n_xton, "listibs in angular spins");
  ALLOCATE(&orbit_arr, n_xton, "orbit_arr");
  ALLOCATE(&orbittot_arr, n_xton, "orbittot_arr");
  ALLOCATE(&lstot_arr, n_xton, "lstot_arr");
  memset(orbit_arr, 0, (size_t)n_xton * sizeof(xyz_st));
  memset(orbittot_arr, 0, (size_t)n_xton * sizeof(double complex));
  memset(lstot_arr, 0, (size_t)n_xton * sizeof(double complex));

  for (ibs = 0, a = lidx; a < lidx + n_el; a++)
  {
    for (i = 0; i < n_ho; i++, ibs++)
    {
      listibs[(a - lidx) * n_ho + i] = ibs;
    }
  }

  /************************************************************/
  /*******************   CALC XTON L, L2   ********************/
  /************************************************************/

  // Distribute excitons round-robin over MPI ranks; OpenMP threads within.
  // (`index` must be private -- it was shared in the original pragma, a race.)
#pragma omp parallel for private(a, b, i, j, ibs, jbs, index, indexba, indexji, ltmp, ctmp, orbit, orbittot)
  for (n = mpir; n < n_xton; n += mpi_size)
  {

    orbit.x = orbit.y = orbit.z = 0.0 + 0.0 * I;
    orbittot = 0.0 + 0.0 * I;

    for (a = lidx; a < lidx + n_el; a++)
    {
      for (i = 0; i < n_ho; i++)
      {
        ibs = listibs[(a - lidx) * n_ho + i];
        ltmp.x = ltmp.y = ltmp.z = 0.0 + 0.0 * I;

        // sum over b
        for (b = lidx; b < lidx + n_el; b++)
        {
          jbs = listibs[(b - lidx) * n_ho + i];
          index = sqr(n_ho) + (a - lidx) * n_el + (b - lidx);

          // c_bi^* * <b| Lx |a>
          ltmp.x += conjmul(bs_coeff[jbs * n_xton + n], l_mom[index].x);

          // c_bi^* * <b| Ly |a>
          ltmp.y += conjmul(bs_coeff[jbs * n_xton + n], l_mom[index].y);

          // c_bi^* * <b| Lz |a>
          ltmp.z += conjmul(bs_coeff[jbs * n_xton + n], l_mom[index].z);
        } // end of b

        // sum over j
        for (j = 0; j < n_ho; j++)
        {
          jbs = listibs[(a - lidx) * n_ho + j];
          index = i * n_ho + j;

          // c_aj^* * <j| Lx |i>
          ltmp.x += conjmul(bs_coeff[jbs * n_xton + n], l_mom[index].x);

          // c_aj^* * <j| Ly |i>
          ltmp.y += conjmul(bs_coeff[jbs * n_xton + n], l_mom[index].y);

          // c_aj^* * <j| Lz |i>
          ltmp.z += conjmul(bs_coeff[jbs * n_xton + n], l_mom[index].z);
        }

        // multiply by the c_ai coeff
        orbit.x += bs_coeff[ibs * n_xton + n] * ltmp.x;
        orbit.y += bs_coeff[ibs * n_xton + n] * ltmp.y;
        orbit.z += bs_coeff[ibs * n_xton + n] * ltmp.z;

        // Lsqr part
        for (b = lidx; b < lidx + n_el; b++)
        {
          for (j = 0; j < n_ho; j++)
          {
            indexba = sqr(n_ho) + (a - lidx) * n_el + (b - lidx);
            indexji = i * n_ho + j;

            jbs = listibs[(b - lidx) * n_ho + j];

            // c_{ai}^n * (c_{bj}^n)^*
            ctmp = conjmul(bs_coeff[jbs * n_xton + n], bs_coeff[ibs * n_xton + n]);

            ltmp.x = ((l_mom[indexba].x * l_mom[indexji].x) +
                      (l_mom[indexba].y * l_mom[indexji].y) +
                      (l_mom[indexba].z * l_mom[indexji].z));

            ltmp.x *= 2.0;

            if (i == j)
            {
              ltmp.x += l2_mom[indexba];
            }

            if (a == b)
            {
              ltmp.x += l2_mom[indexji];
            }

            orbittot += ctmp * ltmp.x;
          } // end of j
        } // end of b
      } // end of i
    } // end of a
    orbit_arr[n] = orbit;
    orbittot_arr[n] = orbittot;
  }

  /************************************************************/
  /*******************   CALC XTON LdotS   ********************/
  /************************************************************/

#pragma omp parallel for private(a, b, i, j, ibs, jbs, indexba, indexji, lstmp, ctmp, lstot)
  for (n = mpir; n < n_xton; n += mpi_size)
  {
    lstot = 0.0 + 0.0 * I;
    for (a = lidx; a < lidx + n_el; a++)
    {
      for (i = 0; i < n_ho; i++)
      {
        ibs = listibs[(a - lidx) * n_ho + i];
        for (b = lidx; b < lidx + n_el; b++)
        {
          for (j = 0; j < n_ho; j++)
          {
            jbs = listibs[(b - lidx) * n_ho + j];
            indexba = sqr(n_ho) + (a - lidx) * n_el + (b - lidx);
            indexji = i * n_ho + j;

            lstmp = 0.0 + 0.0 * I;

            // c_{ai}^n * (c_{bj}^n)^*
            ctmp = conjmul(bs_coeff[jbs * n_xton + n], bs_coeff[ibs * n_xton + n]);

            // <a|L|b>*<j|S|i>
            lstmp += ((l_mom[indexba].x * s_mom[indexji].x) +
                      (l_mom[indexba].y * s_mom[indexji].y) +
                      (l_mom[indexba].z * s_mom[indexji].z));

            // <a|S|b>*<j|L|i>
            lstmp += ((s_mom[indexba].x * l_mom[indexji].x) +
                      (s_mom[indexba].y * l_mom[indexji].y) +
                      (s_mom[indexba].z * l_mom[indexji].z));

            // delta ij part
            if (i == j)
            {
              lstmp += ldots[indexba];
            }
            // delta ab part
            if (a == b)
            {
              lstmp += ldots[indexji];
            }

            lstot += ctmp * lstmp;
          } // end of j
        } // end of b
      } // end of i
    } // end of a
    lstot_arr[n] = lstot;
  } // end of n

  // Complete the per-exciton results across ranks (disjoint excitons summed).
  MPI_Allreduce(MPI_IN_PLACE, orbit_arr, (int)(6 * n_xton), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, orbittot_arr, (int)(2 * n_xton), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, lstot_arr, (int)(2 * n_xton), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (mpir == 0)
  {
    pf = fopen("orbital.dat", "w");
    for (n = 0; n < n_xton; n++)
    {
      fprintf(pf, "%ld\t%-10.5lf\t%-10.5lf\t%-10.5lf\t", n,
              creal(orbit_arr[n].x), creal(orbit_arr[n].y), creal(orbit_arr[n].z));
      fprintf(pf, "%-10.5lf\t (%-10.5lf)\n", creal(orbittot_arr[n]), cimag(orbittot_arr[n]));
    }
    fclose(pf);

    lspf = fopen("couple.dat", "w");
    for (n = 0; n < n_xton; n++)
      fprintf(lspf, "%ld\t%-10.5lf\t (%-10.5lf)\n", n, creal(lstot_arr[n]), cimag(lstot_arr[n]));
    fclose(lspf);
  }

  free(listibs);
  free(orbit_arr);
  free(orbittot_arr);
  free(lstot_arr);

  return;
}
