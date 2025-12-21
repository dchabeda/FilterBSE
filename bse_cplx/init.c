/*****************************************************************************/

#include "init.h"

/*****************************************************************************/
/*****************************************************************************/

void init_elec_hole_kernel(
    double complex *pot_bare,
    double complex *pot_screened,
    grid_st *grid,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel)
{
  /*******************************************************************
   * This function computes the value of the Coulomb kernel at each   *
   * grid point. The actual calculation of 1/r takes place in k-space *
   * using FFT. The direct Coulomb term is screened by the average    *
   * dielectric constant, and the exchange term is unscreened (bare)  *
   * according to the derivation by Rohlfing and Louie: PRB 62 (8)    *
   * inputs:                                                          *
   *  [pot_bare] array to hold value of the direct Coulomb pot      *
   *  [pot_screened] array to hold value of the screened exchange pot *
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

  //              Integer Constants
  const int mpir = parallel->mpi_rank;
  const long ngrid = ist->ngrid;

  //              Long Integers
  long jx;
  long jy;
  long jz;
  long jyz;
  long jxyz;
  long jg;

  //              Double Precision Scalars
  double r;
  double r2;
  double x2;
  double y2;
  double z2;
  double alpha;
  double cosa;
  double sina;
  double ex_1;
  double ey_1;
  double ez_1;
  double sqrtexeyez_1;
  double sqrtaveps;
  double sqrk0;
  double xmax;
  double ymax;
  double zmax;
  double boxl;
  double boxl2;
  double minr;
  double gamma;
  double gamma2;
  double gammaeps;

  //              Double precision arrays
  double *kx2;
  double *ky2;
  double *kz2;

  double complex *potr;
  double complex *potrx;

  double complex tmp;

  /************************************************************/
  /*******************   INIT FFT PARAMS    *******************/
  /************************************************************/
  fftw_plan_loc planfw;
  fftw_plan_loc planbw;
  fftw_complex *fftwpsi;

  fftw_init_threads();
  fftw_plan_with_nthreads(ist->nthreads);

  //         Allocate memory for multithreaded FFT
  fftwpsi = fftw_malloc(ngrid * sizeof(fftw_complex));
  planfw = fftw_plan_dft_3d(grid->nz, grid->ny, grid->nx, fftwpsi, fftwpsi, FFTW_FORWARD, 0);
  planbw = fftw_plan_dft_3d(grid->nz, grid->ny, grid->nx, fftwpsi, fftwpsi, FFTW_BACKWARD, 0);

  /************************************************************/
  /******************  INITIALIZE VARIABLES  ******************/
  /************************************************************/

  // Define gamma
  xmax = fabs(grid->xmin);
  ymax = fabs(grid->ymin);
  zmax = fabs(grid->zmin);

  if ((xmax < ymax) && (xmax < zmax))
  {
    minr = xmax;
    boxl = (double)grid->nx * grid->dx;
  }

  else if ((ymax < zmax) && (ymax < xmax))
  {
    minr = ymax;
    boxl = (double)grid->ny * grid->dy;
  }

  else
  {
    minr = zmax;
    boxl = (double)grid->nz * grid->dz;
  }

  boxl2 = sqr(boxl);
  gamma = 7.0 / (2.0 * minr);
  gamma2 = sqr(gamma);
  grid->dkx = TWOPI / ((double)grid->nx * grid->dx);
  grid->dky = TWOPI / ((double)grid->ny * grid->dy);
  grid->dkz = TWOPI / ((double)grid->nz * grid->dz);

  //             Prepare scaling factors for dielectric screening
  ex_1 = 1.0 / par->epsX;
  ey_1 = 1.0 / par->epsY;
  ez_1 = 1.0 / par->epsZ;
  sqrtexeyez_1 = 1.0 / sqrt(par->epsX * par->epsY * par->epsZ);
  sqrtaveps = sqrt((par->epsX + par->epsY + par->epsZ) / 3.0);
  gammaeps = gamma * sqrtaveps;

  /*** No Yukawa screening for the exchange ***/
  sqrk0 = gamma2 * sqr(sqrtaveps);

  /************************************************************/
  /******************  ALLOC MEM FOR K-VECS  ******************/
  /************************************************************/

  ALLOCATE(&kx2, grid->nx, "kx2");
  ALLOCATE(&ky2, grid->ny, "ky2");
  ALLOCATE(&kz2, grid->nz, "kz2");

  ALLOCATE(&potr, ngrid, "potr");
  ALLOCATE(&potrx, ngrid, "potrx");

  // Initialize the arrays
  for (kx2[0] = 0.0, jx = 1; jx <= grid->nx / 2; jx++)
    kx2[jx] = (kx2[grid->nx - jx] = sqr((double)(jx)*grid->dkx));
  for (ky2[0] = 0.0, jy = 1; jy <= grid->ny / 2; jy++)
    ky2[jy] = (ky2[grid->ny - jy] = sqr((double)(jy)*grid->dky));
  for (kz2[0] = 0.0, jz = 1; jz <= grid->nz / 2; jz++)
    kz2[jz] = (kz2[grid->nz - jz] = sqr((double)(jz)*grid->dkz));

  for (jxyz = 0; jxyz < ist->ngrid; jxyz++)
    potr[jxyz] = 0.0 + 0.0 * I;
  for (jxyz = 0; jxyz < ist->ngrid; jxyz++)
    potrx[jxyz] = 0.0 + 0.0 * I;

  /************************************************************/
  /******************    CALC COULOMB POTS   ******************/
  /************************************************************/

  // Compute the screened and bare Coulomb potential at all gridpoints

  for (jz = 0; jz < grid->nz; jz++)
  {
    z2 = sqr(grid->z[jz]);
    for (jy = 0; jy < grid->ny; jy++)
    {
      y2 = sqr(grid->y[jy]);
      jyz = grid->nx * (grid->ny * jz + jy);
      for (jx = 0; jx < grid->nx; jx++)
      {
        x2 = sqr(grid->x[jx]);
        jxyz = jyz + jx;
        r2 = (x2 + y2 + z2);
        if (r2 < boxl2)
        {
          r = sqrt(x2 + y2 + z2);
          potr[jxyz] = calc_coulomb(r, gamma);
          r = sqrt(ex_1 * x2 + ey_1 * y2 + ez_1 * z2);
          potrx[jxyz] = sqrtexeyez_1 * calc_coulomb(r, gammaeps);
        }
      }
    }
  }

  // Zero out the arrays for all potentials
  for (jxyz = 0; jxyz < ngrid; jxyz++)
  {
    pot_bare[jxyz] = pot_bare[jxyz] = 0.0 + 0.0 * I;
    pot_screened[jxyz] = pot_screened[jxyz] = 0.0 + 0.0 * I;
  }

  if (mpir == 0)
    printf("\tBare Coulomb potential, v(r, r')...\n");
  memcpy(&fftwpsi[0], &potr[0], ngrid * sizeof(fftwpsi[0]));
  fftw_execute(planfw);
  memcpy(&pot_bare[0], &fftwpsi[0], ngrid * sizeof(fftwpsi[0]));

  if (mpir == 0)
    printf("\tScreened Coulomb potential, W(r, r')...\n");
  memcpy(&fftwpsi[0], &potrx[0], ngrid * sizeof(fftwpsi[0]));
  fftw_execute(planfw);
  memcpy(&pot_screened[0], &fftwpsi[0], ngrid * sizeof(fftwpsi[0]));

  for (jz = 0; jz < grid->nz; jz++)
  {
    z2 = kz2[jz];
    for (jy = 0; jy < grid->ny; jy++)
    {
      y2 = ky2[jy];
      jyz = grid->nx * (grid->ny * jz + jy);
      for (jx = 0; jx < grid->nx; jx++)
      {
        x2 = kx2[jx];

        jg = jyz + jx;

        alpha = PIE * (double)(jx + jy + jz + ngrid / 2);
        cosa = cos(alpha);
        sina = sin(alpha);

        /*** hartree term ***/
        tmp = pot_bare[jg];

        pot_bare[jg] = ((creal(tmp) * cosa - cimag(tmp) * sina) + I * (creal(tmp) * sina + cimag(tmp) * cosa)) * grid->dv;
        pot_bare[jg] += (FOURPI / (x2 + y2 + z2 + gamma2));
        pot_bare[jg] *= ist->ngrid_1;

        /*** screened Coulomb term ***/
        tmp = pot_screened[jg];

        pot_screened[jg] = ((creal(tmp) * cosa - cimag(tmp) * sina) + I * (creal(tmp) * sina + cimag(tmp) * cosa)) * grid->dv;
        pot_screened[jg] += (FOURPI * (1.0 - exp(-0.25 * (par->epsX * x2 + par->epsY * y2 + par->epsZ * z2) / sqrk0)) / (par->epsX * x2 + par->epsY * y2 + par->epsZ * z2 + EPSR));
        pot_screened[jg] *= ist->ngrid_1;
      }
    }
  }

  if (mpir == 0)
    printf("  Done generating e-h interaction potential.\n");

  free(potr);
  free(potrx);
  free(kx2);
  free(ky2);
  free(kz2);

  fftw_destroy_plan(planfw);
  fftw_destroy_plan(planbw);
  fftw_free(fftwpsi);

  return;
}

/****************************************************************************/

double calc_coulomb(double r, double gamma)
{
  /*******************************************************************
   * This function computes the 1/r potential in reciprocal space     *
   * as erf(r)/r                                                      *
   * dielectric constant, and the exchange term is unscreened (bare)  *
   * according to the derivation by Rohlfing and Louie: PRB 62 (8)    *
   * inputs:                                                          *
   *  [pot_bare] array to hold value of the direct Coulomb pot      *
   *  [pot_screened] array to hold value of the screened exchange pot *
   *  [grid] grid_st instance holding values of all grid points       *
   *  [ist] ptr to counters, indices, and lengths                     *
   *  [par] ptr to par_st holding VBmin, VBmax... params              *
   *  [flag] ptr to flag_st holding job flags                         *
   *  [planfw] FFTW3 plan for executing 3D forward DFT                *
   *  [planfw] FFTW3 plan for executing 3D backwards DFT              *
   *  [fftwpsi] location to store outcome of Fourier transform        *
   * outputs: void                                                    *
   ********************************************************************/
  if (r < EPSR)
    return (2.0 * gamma / SQRTPI);
  return (erf(gamma * r) / r);
}

/************************************************************************/
