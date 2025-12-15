#include "mod_pot.h"

void mod_pot(
  double complex** pot_bare,
  double complex** pot_screened,
  grid_st*         grid,
  index_st*        ist,
  par_st*          par,
  flag_st*         flag,
  parallel_st*     parallel
  ){

  /************************************************************/
	/*******************  DECLARE VARIABLES   *******************/
	/************************************************************/

  const int              mpir = parallel->mpi_rank;

  const unsigned long    ngrid = ist->ngrid;
  
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
  
  ALLOCATE(pot_bare,     ngrid, "pot_bare");
  ALLOCATE(pot_screened, ngrid, "pot_screened");
  
  if (mpir == 0) printf(" done.\n"); fflush(stdout);

  /************************************************************/
	/*******************  CALC E-H INT. POT   *******************/
	/************************************************************/
  if (mpir == 0) printf("Computing interaction potential on grid...\n"); fflush(stdout);
  
  init_elec_hole_kernel(*pot_bare, *pot_screened, grid, ist, par, flag, parallel);
  
  return;
}
