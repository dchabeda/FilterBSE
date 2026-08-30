/*****************************************************************************/

#include "init.h"

/*****************************************************************************/
/*****************************************************************************/

void init_elec_hole_kernel(
  double*           pot_bare, 
  double*           pot_screened, 
  grid_st*          grid, 
  index_st*         ist, 
  par_st*           par, 
  flag_st*          flag, 
  parallel_st*      parallel, 
  fftw_plan_loc     planfw, 
  fftw_plan_loc     planbw, 
  fftw_complex*     fftwpsi
  ){
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
  const int       mpir   = parallel->mpi_rank;
  const long      cngrid = ist->complex_idx * ist->ngrid;

  //              Long Integers
  long            jx;
  long            jy;
  long            jz;
  long            jyz;
  long            jxyz;
  long            jgr;
  long            jgi;

  //              Double Precision Scalars
  double          r;
  double          r2;
  double          x2;
  double          y2;
  double          z2;
  double          alpha;
  double          cosa;
  double          sina;
  double          ex_1;
  double          ey_1;
  double          ez_1;
  double          sqrtexeyez_1;
  double          sqrtaveps;
  double          sqrk0;
  double          boxl;
  double          boxl2;
  double          minr;
  double          gamma;
  double          gamma2;
  double          gammaeps;

  //              Double Precision Pointers
  double*         kx2;
  double*         ky2;
  double*         kz2;

  //              Complex Variables
  zomplex*        potr;
  zomplex*        potrx;
  zomplex         tmp;

  /************************************************************/
  /******************  INITIALIZE VARIABLES  ******************/
  /************************************************************/

  //              Define gamma
  double          xmax = fabs(grid->xmin);
  double          ymax = fabs(grid->ymin);
  double          zmax = fabs(grid->zmin);

  if ((xmax < ymax) && (xmax < zmax)) {
    minr = xmax;
    boxl = (double)grid->nx * grid->dx;
  }

  else if ((ymax < zmax) && (ymax < xmax)) {
    minr = ymax;
    boxl = (double)grid->ny * grid->dy;
  }

  else {
    minr = zmax;
    boxl = (double)grid->nz * grid->dz;
  }

  boxl2          =   sqr(boxl);
  gamma          =   7.0 / (2.0 * minr);
  gamma2         =   sqr(gamma);
  grid->dkx      =   TWOPI / ((double)grid->nx * grid->dx);
  grid->dky      =   TWOPI / ((double)grid->ny * grid->dy);
  grid->dkz      =   TWOPI / ((double)grid->nz * grid->dz);

  //             Prepare scaling factors for dielectric screening
  ex_1           = 1.0 / par->epsX;
  ey_1           = 1.0 / par->epsY;
  ez_1           = 1.0 / par->epsZ;
  sqrtexeyez_1   = 1.0 / sqrt(par->epsX * par->epsY * par->epsZ);
  sqrtaveps      = sqrt((par->epsX + par->epsY + par->epsZ) / 3.0);
  gammaeps       = gamma * sqrtaveps;

  /*** No Yukawa screening for the exchange ***/
  sqrk0          = gamma2 * sqr(sqrtaveps);

  /************************************************************/
  /******************  ALLOC MEM FOR K-VECS  ******************/
  /************************************************************/
  
  ALLOCATE(&kx2, grid->nx, "kx2");
  ALLOCATE(&ky2, grid->ny, "ky2");
  ALLOCATE(&kz2, grid->nz, "kz2");
  
  ALLOCATE(&potr, ist->ngrid, "potr");
  ALLOCATE(&potrx, ist->ngrid, "potrx");
  
  // Initialize the arrays
  for (kx2[0] = 0.0, jx = 1; jx <= grid->nx / 2; jx++)
    kx2[jx] = (kx2[grid->nx-jx] = sqr((double)(jx) * grid->dkx));
  for (ky2[0] = 0.0, jy = 1; jy <= grid->ny / 2; jy++)
    ky2[jy] = (ky2[grid->ny-jy] = sqr((double)(jy) * grid->dky));
  for (kz2[0] = 0.0, jz = 1; jz <= grid->nz / 2; jz++)
    kz2[jz] = (kz2[grid->nz-jz] = sqr((double)(jz) * grid->dkz));

  for (jxyz = 0; jxyz < ist->ngrid; jxyz++) potr[jxyz].re = potr[jxyz].im = 0.0;
  for (jxyz = 0; jxyz < ist->ngrid; jxyz++) potrx[jxyz].re = potrx[jxyz].im = 0.0;
  
  /************************************************************/
  /******************    CALC COULOMB POTS   ******************/
  /************************************************************/

  // Compute the screened and bare Coulomb potential at all gridpoints

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
  for (jxyz = 0; jxyz < cngrid; jxyz++){
    pot_bare[jxyz]     = pot_bare[jxyz]     = 0.0;
    pot_screened[jxyz] = pot_screened[jxyz] = 0.0;
  }

  if (mpir == 0) printf("\tBare Coulomb potential, v(r, r')...\n");
  memcpy(&fftwpsi[0], &potr[0], ist->ngrid*sizeof(fftwpsi[0]));
  fftw_execute(planfw);
  memcpy(&pot_bare[0], &fftwpsi[0], ist->ngrid*sizeof(fftwpsi[0]));

  if (mpir == 0) printf("\tScreened Coulomb potential, W(r, r')...\n");
  memcpy(&fftwpsi[0], &potrx[0], ist->ngrid*sizeof(fftwpsi[0]));
  fftw_execute(planfw);
  memcpy(&pot_screened[0], &fftwpsi[0], ist->ngrid*sizeof(fftwpsi[0]));
  
  
  for (jz = 0; jz < grid->nz; jz++) {
    z2 = kz2[jz];
    for (jy = 0; jy < grid->ny; jy++) {
      y2 = ky2[jy]; 
      jyz = grid->nx * (grid->ny * jz + jy);
      for (jx = 0; jx < grid->nx; jx++) {
      	x2 = kx2[jx];

      	jxyz = jyz + jx;
        jgr = ist->complex_idx * jxyz;
        jgi = jgr + 1;

      	alpha = PIE * (double)(jx + jy + jz + ist->ngrid / 2);
      	cosa = cos(alpha);
      	sina = sin(alpha);
        
      	/*** hartree term ***/
      	tmp.re = pot_bare[jgr];
      	tmp.im = pot_bare[jgi];
      	
      	pot_bare[jgr] = (tmp.re * cosa - tmp.im * sina) * grid->dv;
      	pot_bare[jgi] = (tmp.re * sina + tmp.im * cosa) * grid->dv;
      	pot_bare[jgr] += (FOURPI / (x2 + y2 + z2 + gamma2));
      	pot_bare[jgr] *= ist->ngrid_1;
      	pot_bare[jgi] *= ist->ngrid_1;

      	/*** screened Coulomb term ***/
      	tmp.re = pot_screened[jgr];
      	tmp.im = pot_screened[jgi];

      	pot_screened[jgr] = (tmp.re * cosa - tmp.im * sina) * grid->dv;
      	pot_screened[jgi] = (tmp.re * sina + tmp.im * cosa) * grid->dv;
      	
      	pot_screened[jgr] += (
          FOURPI * (1.0 - exp(-0.25 * (par->epsX * x2 + par->epsY * y2 + par->epsZ * z2) / sqrk0)) 
          / (par->epsX * x2 + par->epsY * y2 + par->epsZ * z2 + EPSR)
        );
        
		    pot_screened[jgr] *= ist->ngrid_1;
      	pot_screened[jgi] *= ist->ngrid_1;
      }
    }
  }

  if (mpir == 0) printf("  Done generating e-h interaction potential.\n");


  free(potr);
  free(potrx);
  free(kx2);
  free(ky2);
  free(kz2);

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
