#include "mod_pot.h"

void mod_pot(
  double*       pot_local,
  pot_st*       pot,
  xyz_st*       R,
  atom_info*    atom,
  grid_st*      grid,
  zomplex*      LS,
  nlc_st*       nlc,
  long*         nl,
  double*       SO_projs,
  index_st*     ist,
  par_st*       par,
  flag_st*      flag,
  parallel_st*  parallel){

  
  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  const int atyp_tot = ist->ngeoms * ist->n_atom_types;
  const int potfl_tot = atyp_tot * ist->max_pot_file_len;
  const int mpir = parallel->mpi_rank;
  
  if (mpir == 0) printf("\nInitializing potentials...\n"); 
  
  /************************************************************/
  /*******************    LOCAL POTENTIAL   *******************/
  /************************************************************/

  if (mpir == 0) printf("\nLocal pseudopotential:\n");
  
  // Alloc mem for reading the atomic potentials ***/
  
  ALLOCATE(&(pot->dr), atyp_tot, "pot->dr");
  ALLOCATE(&(pot->r), potfl_tot, "pot->r");
  ALLOCATE(&(pot->pseudo), potfl_tot, "pot->pseudo");
  if (1.0 != par->scale_surface_Cs){
    // allocate mem for separate LR potentials 
    // if the surface atoms will be charge balanced
    ALLOCATE(&(pot->r_LR), potfl_tot, "pot->r_LR");
    ALLOCATE(&(pot->pseudo), potfl_tot, "pot->pseudo_LR");
  }
  ALLOCATE(&(pot->file_lens), atyp_tot, "pot->file_lens");


  build_local_pot(pot_local, pot, R, atom, grid, ist, par, flag, parallel);
  write_cube_file(pot_local, grid, "local-pot.cube");
  
  /************************************************************/
  /*******************    SPIN-ORBIT POT    *******************/
  /************************************************************/

  if(flag->SO==1) {
    if (mpir == 0) printf("\nSpin-orbit pseudopotential:\n");
    
    init_SO_projectors(SO_projs, grid, R, atom, ist, par, flag, parallel);
    
    def_LS(LS, ist, par);

    if (mpir == 0) printf("\tSO projectors generated.\n");
  }
  
  /************************************************************/
  /*******************     NON-LOCAL POT    *******************/
  /************************************************************/

  if (flag->NL == 1){
    if (mpir == 0) printf("\nNon-local pseudopotential:\n"); fflush(0);
    
    init_NL_projectors(nlc, nl, SO_projs, grid, R, atom, ist, par, flag, parallel);
    
    if (mpir == 0) printf("\tNL projectors generated.\n");
  }
  
  /************************************************************/
  /*******************    FREE POT MEMORY   *******************/
  /************************************************************/

  free(pot->r); pot->r = NULL; 
  free(pot->pseudo); pot->pseudo = NULL; 
  free(pot->dr); pot->dr = NULL; 
  free(pot->file_lens); pot->file_lens = NULL;
  if (1.0 != par->scale_surface_Cs){
    free(pot->r_LR); pot->r_LR = NULL;
    free(pot->pseudo_LR); pot->pseudo_LR = NULL; 
  }   
  // SO/NL pots
  // free memory allocated to SO_projectors
  if ( (flag->SO == 1) || (flag->NL == 1) ){
    free(SO_projs); SO_projs = NULL;
  }

  return;
}