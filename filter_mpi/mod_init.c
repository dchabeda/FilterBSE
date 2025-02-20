#include "mod_init.h"

/************************************************************/
/************************************************************/
/************************************************************/

// This module initializes the Filter job

void mod_init(
  grid_st*      grid,
  double**      gridx,
  double**      gridy,
  double**      gridz,
  xyz_st**      R,
  atom_info**   atom,
  double**      ene_targets,
  double**      ksqr,
  index_st*     ist,
  par_st*       par,
  flag_st*      flag,
  parallel_st*  parallel){


  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  const int mpir = parallel->mpi_rank;
  char *top;
  char *bottom;
  top    = malloc(2 * sizeof(top[0]));
  bottom = malloc(2 * sizeof(bottom[0]));
  strcpy(top, "T\0");
  strcpy(bottom, "B\0");


  /************************************************************/
  /*******************     READ INPUT       *******************/
  /************************************************************/


  if (mpir == 0) write_separation(stdout, top);
  if (mpir == 0) printf("\n1.\tINITIALIZING JOB | %s\n", get_time());
  if (mpir == 0) write_separation(stdout, bottom); fflush(stdout);


  if (mpir == 0) printf("\nReading job specs from input.par:\n");
  read_input(flag, grid, ist, par, parallel);


  /************************************************************/
  /*******************  ALLOCATING MEMORY   *******************/
  /************************************************************/


  ALLOCATE(R, ist->natoms, "R array");
  ALLOCATE(atom, ist->natoms, "atom struct");
  ALLOCATE(ene_targets, ist->m_states_per_filter, "ene_targets");


  /************************************************************/
  /*******************  INIT ATOMS & GRID   *******************/
  /************************************************************/
  
  
  if (mpir == 0) printf("\nReading atomic coords from conf.par:\n");
  read_conf(*R, *atom, ist, par, flag, parallel);


  if (mpir == 0) printf("\nInitializing the grid parameters:\n");
  init_grid_params(grid, *R, ist, par, flag, parallel);

  /************************************************************/
  /*******************      BUILD GRID      *******************/
  /************************************************************/

  ALLOCATE(gridx, grid->nx, "grid.x");
  ALLOCATE(gridy, grid->ny, "grid.y");
  ALLOCATE(gridz, grid->nz, "grid.z");
  ALLOCATE(ksqr, ist->ngrid, "ksqr");


  if (mpir == 0) printf("\nBuilding the r- and k-space grids:\n");
  build_grid_ksqr(*ksqr, *R, grid, ist, par, flag, parallel);
  
  /************************************************************/
  /*******************   SET ENE TARGETS    *******************/
  /************************************************************/
  
  if (mpir == 0) printf("\nSetting the filter energy targets:\n");
  set_ene_targets(*ene_targets, ist, par, flag, parallel);
  
  free(top); free(bottom);
  
  return;
}
