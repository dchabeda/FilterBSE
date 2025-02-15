#include "fd.h"

/************************************************************/
/************************************************************/
/************************************************************/

// This module is the driver for initializing the Filter job

void mod_init(
  grid_st      *grid,
  xyz_st       *R,
  atom_info    *atom,
  double       *ene_targets, 
  double       *ksqr,
  index_st     *ist,
  par_st       *par,
  flag_st      *flag, 
  parallel_st  *parallel){


  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  const int mpir = parallel->mpi_rank;
  char *top; 
  char *bottom; 
  top = malloc(2*sizeof(top[0])); 
  bottom = malloc(2*sizeof(bottom[0]));
  strcpy(top, "T\0"); 
  strcpy(bottom, "B\0");

  /************************************************************/
  /*******************     READ INPUT       *******************/
  /************************************************************/

  if (mpir == 0) write_separation(stdout, top);
  if (mpir == 0) printf("\n1.\tINITIALIZING JOB | %s\n", get_time());
  if (mpir == 0) write_separation(stdout, bottom); fflush(stdout);

  /*** read initial setup from input.par ***/
  //
  //
  if (mpir == 0) printf("\nReading job specs from input.par:\n");
  read_input(&flag, grid, &ist, &par, &parallel);
  //
  //

  /************************************************************/
  /*******************  ALLOCATING MEMORY   *******************/
  /************************************************************/

  // the positions of the atoms in the x, y, and z directions 
  ALLOCATE(R, ist->natoms, "R array");
  
  // the atom specific information 
  ALLOCATE(atom, ist->natoms, "atom struct");
  
  // the energies of each energy target in filter
  ALLOCATE(ene_targets, ist->m_states_per_filter, "ene_targets");
  
  // the number of atom types
  ALLOCATE(&(ist->atom_types), N_MAX_ATOM_TYPES, "ist->atom_types");

  /************************************************************/
  /*******************  INIT ATOMS & GRID   *******************/
  /************************************************************/

  /*** read the nanocrystal configuration ***/
  //
  //
  if (mpir == 0) printf("\nReading atomic coords from conf.par:\n");
  read_conf(*R, *atom, ist, par, flag, parallel);
  //
  //

  /*** initialize parameters for the grid ***/
  //
  //
  if (mpir == 0) printf("\nInitializing the grid parameters:\n");
  init_grid_params(grid, R, &ist, &par, &flag, &parallel);
  //
  //

  // Allocate memory for the grid in the x, y, and z directions ***/
  ALLOCATE(&(grid->x), grid->nx, "grid.x");
  ALLOCATE(&(grid->y), grid->ny, "grid.y");
  ALLOCATE(&(grid->z), grid->nz, "grid.z");
  
  // the kinetic energy stored on the grid
  ALLOCATE(ksqr, ist->ngrid, "ksqr");

  /*** build the real- and k-space grids ***/
  if (mpir == 0) printf("\nBuilding the real-space and k-space grids:\n");
  build_grid_ksqr(ksqr, R, grid, &ist, &par, &flag, &parallel);
  
  /*** set the energy targets ***/
  if (mpir == 0) printf("\nSetting the filter energy targets:\n");
  set_ene_targets(ene_targets, &ist, &par, &flag, &parallel);

  return;
}