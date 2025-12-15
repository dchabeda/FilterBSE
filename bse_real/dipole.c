/****************************************************************************/

#include "dipole.h"

/****************************************************************************/

void calc_elec_dipole(
    xyz_st *elec_dip,
    double *psi_qp,
    double *eig_vals,
    grid_st *grid,
    index_st *ist,
    par_st *par,
    flag_st *flag)
{
  /*******************************************************************
   * This function computes the electric transition dipole matrix     *
   * matrix elements.                                                 *
   * inputs:                                                          *
   *  [elec_dip] array to hold matrix elems in x, y, z direction  *
   *  [psi_qp] array holding all qp_basis states                      *
   *  [eig_vals] array holding the quasiparticle orbital energies     *
   *  [grid] grid_st instance holding values of all grid points       *
   *  [ist] ptr to counters, indices, and lengths                     *
   *  [par] ptr to par_st holding VBmin, VBmax... params              *
   * outputs: void                                                    *
   ********************************************************************/

  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  FILE *pf;

  //                Indices
  unsigned long i;
  unsigned long a;
  unsigned long i_st;
  unsigned long a_st;
  unsigned long idx;

  //                Scalars
  double z;
  double y;
  double x;
  double E_ia;
  double mu2;

  //                Constants
  const long nspngr = ist->nspinngrid;
  const long ngrid = ist->ngrid;
  const long lidx = ist->lumo_idx;
  const long n_el = ist->n_elecs;
  const long n_ho = ist->n_holes;

  const double dv = par->dv;

  // Output will be written to these files
  pf = fopen("OS0.dat", "w");

  fprintf(pf, "i  a   sqrt(mu2)     Ea-Ei 	  f_osc       mu_x        mu_y       mu_z\n");

  /************************************************************/
  /*******************   CALC ELEC DIPOLE   *******************/
  /************************************************************/

  double start = omp_get_wtime();

// nvtxRangePushA("Calc elec dipole");
#pragma omp parallel for private(a, i_st, a_st, x, y, z, idx)
  for (i = 0; i < n_ho; i++)
  {
    for (a = lidx; a < lidx + n_el; a++)
    {
      // // nvtxRangePushA("loop over a");
      i_st = i * nspngr;
      a_st = a * nspngr;

      long jx;
      long jy;
      long jyz;
      long jz;
      long jg;
      long jsg;

      double mu_x;
      double mu_y;
      double mu_z;
      double tmp;

      mu_x = mu_y = mu_z = 0.0;

      // If using spinors, integrate cplx psi over two spins
      // Loop over the grid
      for (jz = 0; jz < ist->nz; jz++)
      {
        z = grid->z[jz];
        for (jy = 0; jy < ist->ny; jy++)
        {
          y = grid->y[jy];
          jyz = ist->nx * (ist->ny * jz + jy);
#pragma omp simd safelen(8) aligned(psi_qp, grid : BYTE_BOUNDARY) reduction(+ : mu_x, mu_y, mu_z)
          for (jx = 0; jx < ist->nx; jx++)
          {
            x = grid->x[jx];
            jg = jyz + jx;

            // MATRIX ELEM
            tmp = psi_qp[i_st + jg] * psi_qp[a_st + jg];
            mu_x += tmp * x;
            mu_y += tmp * y;
            mu_z += tmp * z;
          }
        }
      }
      // Place these values into alloc'd memory
      idx = i * n_el + (a - lidx);

      elec_dip[idx].x = dv * mu_x;
      elec_dip[idx].y = dv * mu_y;
      elec_dip[idx].z = dv * mu_z;
    }
  }

  double end = omp_get_wtime();
  printf("Duration of computing dipole: %f s (%f s per transition)\n", end - start, (end - start) / (n_el * n_ho));
  // nvtxRangePop();

  /************************************************************/
  /*******************     PRINT VALUES     *******************/
  /************************************************************/

  for (i = 0; i < n_ho; i++)
  {
    for (a = lidx; a < lidx + n_el; a++)
    {
      idx = i * n_el + (a - lidx);

      E_ia = eig_vals[a] - eig_vals[i];
      mu2 = sqr(elec_dip[idx].x) + sqr(elec_dip[idx].y) + sqr(elec_dip[idx].z);

      fprintf(pf, "%ld % ld  %.8f %.12f % .8f % .8f % .8f % .8f\n",
              i, a, sqrt(mu2), E_ia, (2.0 / 3.0) * E_ia * mu2,
              creal(elec_dip[idx].x),
              creal(elec_dip[idx].y),
              creal(elec_dip[idx].z));
    }
  }
  fclose(pf);

  return;
}

/****************************************************************************/

void calc_mag_dipole(
    xyz_st *mag_dip,
    double *psi_qp,
    double *eig_vals,
    grid_st *grid,
    index_st *ist,
    par_st *par,
    flag_st *flag)
{
  // This function calculates the magnetic dipole matrix elements between the
  // single-particle electron (a) and hole (i) states: <psi_a|m|psi_i>
  // where m = -1/2*L = -1/2 * (r x p) where x is the cross product.

  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  //                File pointers
  FILE *pf;

  //                Indices and states
  long i;
  long a;
  long i_st;
  long a_st;
  long idx;
  long jg;

  //                Constants
  const long nspngr = ist->nspinngrid;
  const long ngrid = ist->ngrid;
  const long lidx = ist->lumo_idx;
  const long n_el = ist->n_elecs;
  const long n_ho = ist->n_holes;

  const double dv = par->dv;

  //                Double values
  double E_ia;
  double ms;

  //                k-space grid vectors
  double *gx;
  double *gy;
  double *gz;
  double *g_vecs;

  //                Complex arrays
  double complex *Lxpsi;
  double complex *Lypsi;
  double complex *Lzpsi;
  double complex *psi_tmp;

  // FFTW
  fftw_plan_loc planfw;
  fftw_plan_loc planbw;
  fftw_complex *fftwpsi;

  // Allocate memory for multithreaded FFT
  fftw_init_threads();
  fftw_plan_with_nthreads(ist->nthreads);

  fftwpsi = fftw_malloc(ngrid * sizeof(fftw_complex));
  planfw = fftw_plan_dft_3d(grid->nz, grid->ny, grid->nx, fftwpsi, fftwpsi, FFTW_FORWARD, 0);
  planbw = fftw_plan_dft_3d(grid->nz, grid->ny, grid->nx, fftwpsi, fftwpsi, FFTW_BACKWARD, 0);

  /************************************************************/
  /*********************     INIT G VECS    *******************/
  /************************************************************/

  // G vectors for derivatives in k space
  ALLOCATE(&gx, grid->nx, "gx in mag_dipole");
  ALLOCATE(&gy, grid->ny, "gy in mag_dipole");
  ALLOCATE(&gz, grid->nz, "gz in mag_dipole");
  ALLOCATE(&g_vecs, 3 * ngrid, "g_vecs in mag_dipole");

  // Initializing the G vectors
  init_g_vecs(g_vecs, gx, gy, gz, grid, ist, par);

  // Allocate memory
  ALLOCATE(&Lxpsi, nspngr, "Lxpsi in mag_dipole");
  ALLOCATE(&Lypsi, nspngr, "Lypsi in mag_dipole");
  ALLOCATE(&Lzpsi, nspngr, "Lzpsi in mag_dipole");
  ALLOCATE(&psi_tmp, nspngr, "psi_tmp in mag_dipole");

  // Output will be written to this file
  pf = fopen("M0.dat", "w");
  // fprintf(pf, "  i   a    sqrt(ms)       Ea-Ei         m_x.re      m_x.im      m_y.re      m_y.im      m_z.re      m_z.im\n");

  /************************************************************/
  /********************     CALC MAG DIP    *******************/
  /************************************************************/
  // nvtxRangePushA("Calc mag dipole");
  double start = omp_get_wtime();
  for (a = lidx; a < lidx + n_el; a++)
  {
    a_st = a * nspngr;

    real_to_complex(&psi_qp[a_st], psi_tmp, nspngr);

    l_operator(&Lxpsi[0], &Lypsi[0], &Lzpsi[0], psi_tmp, g_vecs, grid, ist, par, planfw, planbw, fftwpsi);

#pragma omp parallel for private(i, i_st, jg, idx)
    for (i = 0; i < n_ho; i++)
    {
      i_st = i * nspngr;

      double complex m_x; // <a|m|i>_x
      double complex m_y; // <a|m|i>_y
      double complex m_z; // <a|m|i>_z

      m_x = m_y = m_z = 0.0 + 0.0 * I;

#pragma omp simd safelen(8) aligned(psi_qp, Lxpsi, Lypsi, Lzpsi : BYTE_BOUNDARY) reduction(+ : m_x, m_y, m_z)
      for (jg = 0; jg < nspngr; jg++)
      {
        m_x += psi_qp[i_st + jg] * Lxpsi[jg];
        m_y += psi_qp[i_st + jg] * Lypsi[jg];
        m_z += psi_qp[i_st + jg] * Lzpsi[jg];
      }

      m_x *= 0.5 * dv;
      m_y *= 0.5 * dv;
      m_z *= 0.5 * dv;

      idx = i * n_el + (a - lidx);

      mag_dip[idx].x = m_x;
      mag_dip[idx].y = m_y;
      mag_dip[idx].z = m_z;
    }
  }

  double end = omp_get_wtime();
  printf("Duration of computing mag dipole: %f s (%f s per transition)\n", end - start, (end - start) / (n_el * n_ho));

  /************************************************************/
  /*******************     PRINT VALUES     *******************/
  /************************************************************/

  for (a = lidx; a < lidx + n_el; a++)
  {
    for (i = 0; i < n_ho; i++)
    {
      idx = i * n_el + (a - lidx);

      ms = cnorm(mag_dip[idx].x) + cnorm(mag_dip[idx].y) + cnorm(mag_dip[idx].z);

      E_ia = eig_vals[a] - eig_vals[i];

      fprintf(pf, "%3ld %3ld  %+.8f %+.12f %+.8f %+.8f %+.8f %+.8f %+.8f %+.8f %+.8f\n", i, a, sqrt(ms), E_ia, (4.0 / 3.0) * E_ia * ms,
              creal(mag_dip[idx].x), cimag(mag_dip[idx].x),
              creal(mag_dip[idx].y), cimag(mag_dip[idx].y),
              creal(mag_dip[idx].z), cimag(mag_dip[idx].z));
    }
  }
  fclose(pf);

  // Free dynamically allocated memory
  free(Lxpsi);
  free(Lypsi);
  free(Lzpsi);

  fftw_free(fftwpsi);
  fftw_destroy_plan(planfw);
  fftw_destroy_plan(planbw);

  free(gx);
  free(gy);
  free(gz);
  free(g_vecs);

  return;
}

// /****************************************************************************/

void calc_rotational_strength(
    double *rs,
    xyz_st *elec_dip,
    xyz_st *mag_dip,
    double *eval,
    index_st *ist)
{
  // RS = Im(u_ia * m_ai)
  // where m_ai = m_ia^*

  FILE *pf;
  int i, a, ia;

  pf = fopen("rs0.dat", "w");
  for (i = 0; i < ist->n_holes; i++)
  {
    for (a = ist->lumo_idx; a < ist->lumo_idx + ist->n_elecs; a++)
    {
      // R_ia = Im[(ux.re + ux.im)*(mx.re - mx.im) + ...]
      ia = i * ist->n_elecs + (a - ist->lumo_idx);
      // Compute imaginary part of u_ia \dot m_ia^*
      rs[ia] = cimag(elec_dip[ia].x * conj(mag_dip[ia].x));
      rs[ia] += cimag(elec_dip[ia].y * conj(mag_dip[ia].y));
      rs[ia] += cimag(elec_dip[ia].z * conj(mag_dip[ia].z));

      fprintf(pf, "%d %d %.12f %.12g\n", i, a, eval[a] - eval[i], rs[ia]);
    }
  }
  fclose(pf);

  return;
}
/****************************************************************************/
