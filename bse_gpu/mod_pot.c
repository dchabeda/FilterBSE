#include "mod_pot.h"

void mod_pot(
  double**       pot_bare,
  double**       pot_screened,
  grid_st*       grid,
  index_st*      ist,
  par_st*        par,
  flag_st*       flag,
  parallel_st*   parallel
  ){

  /************************************************************/
	/*******************  DECLARE VARIABLES   *******************/
	/************************************************************/

  const int              mpir = parallel->mpi_rank;

  const unsigned long    cngrid = ist->complex_idx * ist->ngrid;
  
  fftw_plan_loc          planfw;
  fftw_plan_loc          planbw;
  fftw_complex*          fftwpsi;
  
  /************************************************************/
	/*******************   INIT FFT PARAMS    *******************/
	/************************************************************/

  fftw_init_threads();
	fftw_plan_with_nthreads(ist->nthreads);

	//         Allocate memory for multithreaded FFT
	fftwpsi    = fftw_malloc(ist->ngrid * sizeof(fftw_complex));
  planfw     = fftw_plan_dft_3d(grid->nz, grid->ny, grid->nx, fftwpsi, fftwpsi, FFTW_FORWARD, 0);
	planbw     = fftw_plan_dft_3d(grid->nz, grid->ny, grid->nx, fftwpsi, fftwpsi, FFTW_BACKWARD, 0);
 

  /************************************************************/
	/*******************  ALLOC MEM FOR POTS  *******************/
	/************************************************************/
  if (mpir == 0){ 
    write_separation(stdout, "T");
    printf("\n2.\tCOMPUTING ELECTRON-HOLE INTERACTION POTENTIALS | %s\n", get_time());
    write_separation(stdout, "B"); fflush(stdout);
  }

  // Allocate memory for the hartree, screened coulomb and bare exchange pots
  if (mpir == 0) printf("Allocating memory for hartree, direct, and exchange potentials..."); fflush(stdout);
  
  ALLOCATE(pot_bare, cngrid, "pot_bare");
  ALLOCATE(pot_screened, cngrid, "pot_screened");
  
  if (mpir == 0) printf(" done.\n"); fflush(stdout);

  /************************************************************/
	/*******************  CALC E-H INT. POT   *******************/
	/************************************************************/
  if (mpir == 0) printf("Computing interaction potential on grid...\n"); fflush(stdout);
  
  init_elec_hole_kernel(*pot_bare, *pot_screened, grid, ist, par, flag, parallel, planfw, planbw, &fftwpsi[0]);


  fftw_destroy_plan(planfw);
  fftw_destroy_plan(planbw);
  fftw_free(fftwpsi);
  
  return;
}
