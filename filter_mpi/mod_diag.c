#include "mod_diag.h"

void mod_diag(
  double*       psitot,
  double*       pot_local,
  double*       eig_vals,
  double*       sigma_E,
  grid_st*      grid,
  zomplex*      LS,
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

  double init_clock;
  double init_wall;

  /************************************************************/
  /*******************  RESTART FROM DIAG   *******************/
  /************************************************************/


  if (flag->restartFromCheckpoint == 2){
    par->checkpoint_id = flag->restartFromCheckpoint;
    restart_from_save(
      "checkpoint_2.dat", par->checkpoint_id, psitot, pot_local, ksqr, an, zn,
      ene_targets, nl, nlc, grid, ist, par, flag, parallel);
  }
  par->checkpoint_id++;


  /************************************************************/
  /*******************   DIAGONALIZE HAM    *******************/
  /************************************************************/
  /*** diagonalize hamiltonian in the subspace spanned by   ***/
  /*** orthogonal filtered states, generating eigenstates   ***/
  /*** of the hamiltonian within the desired energy range   ***/
  /************************************************************/
  
  if (mpir == 0){
    write_separation(stdout, "T");
    printf("\n6. DIAGONALIZING HAMILTONIAN | %s\n", get_time()); 
    write_separation(stdout, "B"); fflush(stdout);
  }

  init_clock = (double)clock(); 
  init_wall = (double)time(NULL);

  diag_H(psitot, pot_local, LS, nlc, nl, ksqr, eig_vals, ist, par, flag, parallel);
  normalize_all(&psitot[0],ist->mn_states_tot, ist, par, flag, parallel);
  
  
  if (mpir == 0) printf("\ndone calculating Hmat, CPU time (sec) %g, wall run time (sec) %g\n",
              ((double)clock()-init_clock)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-init_wall);
  fflush(stdout);

  /************************************************************/
  /*******************     CALC SIGMA E     *******************/
  /************************************************************/
  /************************************************************/
  /*** calculate the standard deviation of these states     ***/
  /*** this is used to check if there are ghost states      ***/
  /************************************************************/

  if (mpir == 0){
    write_separation(stdout, "T");
    printf("\n7. CALCULATING VARIANCE OF EIGENVALUES | %s\n", get_time());
    write_separation(stdout, "B"); fflush(stdout);
  }

  calc_sigma_E(psitot, pot_local, LS, nlc, nl, ksqr, sigma_E, ist, par, flag);
      
  return;
}