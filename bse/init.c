/*****************************************************************************/

#include "fd.h"

/*****************************************************************************//*****************************************************************************/

void init_elec_hole_kernel(zomplex *pot_bare, zomplex *pot_screened, grid_st *grid, index_st *ist, par_st *par, fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi){
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

  long jx, jy, jz, jyz, jxyz, sx, sy, sz;
  double r, r2, x2, y2, z2, *kx2, *ky2, *kz2, alpha, cosa, sina;
  double ex_1, ey_1, ez_1, sqrtexeyez_1, sqrtaveps, sqrk0;
  double boxl, boxl2, minr;
  double gamma, gamma2, gammaeps;
  zomplex *potr, *potrx, tmp;

  // Define gamma
  if ((grid->xmax < grid->ymax) && (grid->xmax < grid->zmax)){
    minr = grid->xmax;
    boxl = (double)grid->nx * grid->dx;
  } else if ((grid->ymax < grid->zmax) && (grid->ymax < grid->xmax)){
    minr = grid->ymax;
    boxl = (double)grid->ny * grid->dy;
  } else {
    minr = grid->zmax;
    boxl = (double)grid->nz * grid->dz;
  }

  boxl2 = sqr(boxl);
  gamma = 7.0 / (2.0 * minr);
  gamma2 = sqr(gamma);
  grid->dkx = TWOPI / ((double)grid->nx * grid->dx);
  grid->dky = TWOPI / ((double)grid->ny * grid->dy);
  grid->dkz = TWOPI / ((double)grid->nz * grid->dz);
  
  // Prepare scaling factors for dielectric screening
  ex_1 = 1.0 / par->epsX;
  ey_1 = 1.0 / par->epsY;
  ez_1 = 1.0 / par->epsZ;
  sqrtexeyez_1 = 1.0 / sqrt(par->epsX * par->epsY * par->epsZ);
  sqrtaveps = sqrt((par->epsX + par->epsY + par->epsZ) / 3.0);
  gammaeps = gamma * sqrtaveps;
  
  
  /*** no yukawa screening for the exchange ***/
  sqrk0 = gamma2 * sqr(sqrtaveps);
  
  // Allocate memory to arrays for computing Coulomb kernel in k-space
  if ((kx2  = (double*)calloc(grid->nx, sizeof(double)))==NULL){
    fprintf(stderr, "ERROR: allocating memory for kx2 in init.c\n");
    exit(EXIT_FAILURE);
  }
  if ((ky2  = (double*)calloc(grid->ny, sizeof(double)))==NULL){
    fprintf(stderr, "ERROR: allocating memory for ky2 in init.c\n");
    exit(EXIT_FAILURE);
  }
  if ((kz2  = (double*)calloc(grid->nz, sizeof(double)))==NULL){
    fprintf(stderr, "ERROR: allocating memory for kz2 in init.c\n");
    exit(EXIT_FAILURE);
  }
  if ((potr = (zomplex*)calloc(ist->ngrid, sizeof(zomplex)))==NULL){
    fprintf(stderr, "ERROR: allocating memory for potr in init.c\n");
    exit(EXIT_FAILURE);
  }
  if ((potrx = (zomplex*)calloc(ist->ngrid, sizeof(zomplex)))==NULL){
    fprintf(stderr, "ERROR: allocating memory for potrx in init.c\n");
    exit(EXIT_FAILURE);
  }
  
  // Initialize the arrays
  for (kx2[0] = 0.0, jx = 1; jx <= grid->nx / 2; jx++)
    kx2[jx] = (kx2[grid->nx-jx] = sqr((double)(jx) * grid->dkx));
  for (ky2[0] = 0.0, jy = 1; jy <= grid->ny / 2; jy++)
    ky2[jy] = (ky2[grid->ny-jy] = sqr((double)(jy) * grid->dky));
  for (kz2[0] = 0.0, jz = 1; jz <= grid->nz / 2; jz++)
    kz2[jz] = (kz2[grid->nz-jz] = sqr((double)(jz) * grid->dkz));

  for (jxyz = 0; jxyz < ist->ngrid; jxyz++) potr[jxyz].re = potr[jxyz].im = 0.0;
  for (jxyz = 0; jxyz < ist->ngrid; jxyz++) potrx[jxyz].re = potrx[jxyz].im = 0.0;
  // ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

  // Compute the screened and bare Coulomb potential at all points on the grid
  for (jz = 0; jz < grid->nz; jz++) {
    z2 = sqr(grid->x[jz]);
    for (jy = 0; jy < grid->ny; jy++) {
      y2 = sqr(grid->y[jy]);
      jyz = grid->nx * (grid->ny * jz + jy);
      for (jx = 0; jx < grid->nx; jx++) {
      	x2 = sqr(grid->x[jx]);
      	jxyz = jyz + jx;
      	r2 = (x2 + y2 + z2);
      	if (r2 < boxl2) {
      	  r = sqrt(x2 + y2 + z2);
      	  potr[jxyz].re = calc_coulomb(r, gamma);
      	  r = sqrt(ex_1*x2 + ey_1*y2 + ez_1*z2);
      	  potrx[jxyz].re = sqrtexeyez_1 * calc_coulomb(r, gammaeps);
      	}
      }
    }
  }

  // Zero out the arrays for all potentials
  for (jxyz = 0; jxyz < ist->ngrid; jxyz++){
    pot_bare[jxyz].re = pot_bare[jxyz].im = pot_screened[jxyz].re = pot_screened[jxyz].im = 0.0;
  }


  printf("\tBare Coulomb potential, v(r, r')...\n");
  memcpy(&fftwpsi[0], &potr[0], ist->ngrid*sizeof(fftwpsi[0]));
  fftw_execute(planfw);
  memcpy(&pot_bare[0], &fftwpsi[0], ist->ngrid*sizeof(pot_bare[0]));

  printf("\tScreened Coulomb potential, W(r, r')...\n");
  memcpy(&fftwpsi[0], &potrx[0], ist->ngrid*sizeof(fftwpsi[0]));
  fftw_execute(planfw);
  memcpy(&pot_screened[0], &fftwpsi[0], ist->ngrid*sizeof(pot_screened[0]));
  
  
  for (jz = 0; jz < grid->nz; jz++) {
    z2 = kz2[jz];
    for (jy = 0; jy < grid->ny; jy++) {
      y2 = ky2[jy]; 
      jyz = grid->nx * (grid->ny * jz + jy);
      for (jx = 0; jx < grid->nx; jx++) {
      	x2 = kx2[jx];
      	jxyz = jyz + jx;
      	alpha = PIE * (double)(jx + jy + jz + ist->ngrid / 2);
      	cosa = cos(alpha);
      	sina = sin(alpha);
        
      	/*** hartree term ***/
      	tmp.re = pot_bare[jxyz].re;
      	tmp.im = pot_bare[jxyz].im;
      	
      	pot_bare[jxyz].re = (tmp.re * cosa - tmp.im * sina) * grid->dv;
      	pot_bare[jxyz].im = (tmp.re * sina + tmp.im * cosa) * grid->dv;
      	pot_bare[jxyz].re += (FOURPI / (x2 + y2 + z2 + gamma2));
      	pot_bare[jxyz].re *= ist->ngrid_1;
      	pot_bare[jxyz].im *= ist->ngrid_1;

      	/*** screened Coulomb term ***/
      	tmp.re = pot_screened[jxyz].re;
      	tmp.im = pot_screened[jxyz].im;

      	pot_screened[jxyz].re = (tmp.re * cosa - tmp.im * sina) * grid->dv;
      	pot_screened[jxyz].im = (tmp.re * sina + tmp.im * cosa) * grid->dv;
      	
      	pot_screened[jxyz].re += (FOURPI * (1.0 - exp(-0.25 * (par->epsX * x2 + par->epsY * y2 + par->epsZ * z2) / sqrk0)) / (par->epsX * x2 + par->epsY * y2 + par->epsZ * z2 + EPSR));
        
		    pot_screened[jxyz].re *= ist->ngrid_1;
      	pot_screened[jxyz].im *= ist->ngrid_1;
      }
    }
  }

  printf("  Done generating e-h interaction potential.\n");
  free(potr);  free(potrx); free(kx2); free(ky2); free(kz2);

  return;
}
  
/****************************************************************************/

double calc_coulomb(double r, double gamma){
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
  if (r < EPSR) return (2.0 * gamma / SQRTPI);
  return (erf(gamma * r) / r);
  
}

/************************************************************************/
