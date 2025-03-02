#include "mod_ortho.h"

void mod_ortho(
  double*       psitot,
  double*       pot_local,
  grid_st*      grid,
  nlc_st*       nlc,
  long*         nl,
  zomplex*      an,
  double*       zn,
  double*       ene_targets,
  double*       ksqr,
  index_st*     ist,
  par_st*       par,
  flag_st*      flag,
  parallel_st*  parallel
  ){

  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/
  const int mpir = parallel->mpi_rank;
  
  long tot_sz = ist->complex_idx*ist->nspinngrid*ist->mn_states_tot;
  
  double init_clock;
  double init_wall;
  
  
  /************************************************************/
  /*******************  RESTART FROM CHKPT  *******************/
  /************************************************************/
  
  if ( (1 == flag->restartFromCheckpoint) && (1 != flag->restartFromOrtho) ){
    par->checkpoint_id = flag->restartFromCheckpoint;
    // 
    restart_from_save(
      "checkpoint_1.dat", par->checkpoint_id, psitot, pot_local, ksqr, an, zn,
      ene_targets, nl, nlc, grid, ist, par, flag, parallel);
    // 
  }
  par->checkpoint_id++;
  
  // Time reverse spinors

  if ( (1 == flag->useSpinors) && (1 != flag->noTimeRev) ){
    printf("\nTime-reversing filtered states (2x no. of orthogonal states)"); fflush(0);
    time_reverse_all(&psitot[0], &psitot[tot_sz], ist, parallel);
  }

  /************************************************************/
  /*******************    ORTHOGONALIZE     *******************/
  /************************************************************/

  if (mpir == 0){
    write_separation(stdout, "T");
    printf("\n5. ORTHOGONALIZATING FILTERED STATES | %s\n", get_time());
    write_separation(stdout, "B");

    printf("\nmn_states_tot before ortho = %ld\n", ist->mn_states_tot);
    fflush(stdout);
  }

  init_clock = (double)clock(); 
  init_wall = (double)time(NULL);

  //
  //
  if (1 == flag->isComplex){
    ist->mn_states_tot = ortho_cplx((MKL_Complex16*)psitot, grid->dv, ist, par, flag, parallel);      
  } 
  else {
    ist->mn_states_tot = ortho_real(psitot, grid->dv, ist, par, flag, parallel);
  }
  //
  //

  if (mpir == 0) printf("mn_states_tot after ortho = %ld\n", ist->mn_states_tot); fflush(0);
  // psitot = (double *) realloc(psitot, ist->mn_states_tot * ist->nspinngrid * ist->complex_idx * sizeof(psitot[0]) );
  

  // Normalize the orthogonalized states

  normalize_all(psitot, ist->mn_states_tot, ist, par, flag, parallel);
        
  if (mpir == 0){
    printf(
      "\ndone calculating ortho, CPU time (sec) %g, wall run time (sec) %g\n",
      ((double)clock()-init_clock)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-init_wall); 
    fflush(stdout);
  } 

  return;
}
