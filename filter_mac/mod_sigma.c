#include "mod_sigma.h"

/*****************************************************************************/


void mod_sigma(
  double*       psitot,
  double*       pot_local,
  double*       eig_vals,
  double*       sigma_E,
  grid_st*      grid,
  zomplex*      LS,
  nlc_st*       nlc,
  long*         nl,
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
  
  unsigned long long psitot_sz;
  psitot_sz = ist->mn_states_tot * ist->nspinngrid * ist->complex_idx * sizeof(double);
  
  unsigned long long st_sz = (long long) ist->complex_idx * ist->nspinngrid * sizeof(double);
  unsigned long long mem_thresh;
  mem_thresh = 380ULL * 1024 * 1024 * 1024; // 380GiB RAM threshold
  
  
  /************************************************************/
  /*******************    CALC SIGMA E     ********************/
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

  if (psitot_sz + (ist->nthreads * st_sz) < mem_thresh){
    // This function parallelizes over states and allocates an additional
    // nthreads*state_mem to heap. Can be ~40GB overhead!
    printf("Total mem < Large mem threshold. Parallelizing over states\n"); fflush(0);
    calc_sigma_E(psitot, pot_local, LS, nlc, nl, ksqr, sigma_E, ist, par, flag);
  } else{
    // The large memory compatible function parallelizes the Hamiltonian to
    // reduce the memory footprint
    printf("Total mem exceeds large mem threshold. Parallelizing Hamiltonian\n"); fflush(0);
    calc_sigma_E_lg_mem(psitot, pot_local, LS, nlc, nl, ksqr, sigma_E, ist, par, flag);
  }

  return;
}

/*****************************************************************************/

void restart_from_sigma(
  double**      psitot,
  double*       pot_local,
  double*       eig_vals,
  grid_st*      grid,
  zomplex*      LS,
  nlc_st*       nlc,
  long*         nl,
  double*       ksqr,
  index_st*     ist,
  par_st*       par,
  flag_st*      flag,
  parallel_st*  parallel
  ){

  FILE *pf;

  long long tot_sz = ist->complex_idx * ist->nspinngrid * ist->mn_states_tot;
  
  /************************************************************/
  /*******************     READ IN PSITOT   *******************/
  /************************************************************/

  write_separation(stdout, "T");
  printf("****  RESTART FROM SIGMA  **  RESTART FROM SIGMA  **  RESTART FROM SIGMA  ****");
  write_separation(stdout, "B"); fflush(stdout);
  
  
  printf("\nNo. states for Sigma, <dH^2> = %ld\n", ist->mn_states_tot);
  printf("Size of psitot array = %.2g GiB\n", 
    sizeof(double) * (double)tot_sz/1024/1024/1024
  ); 
  fflush(stdout);
  
  // Allocate memory for psitot
  ALLOCATE(psitot, tot_sz, "psitot in restart_from_sigma"); 
  
  pf = fopen("psi-diag.dat", "r");
  if (pf == NULL){
    printf("ERROR: could not open psi-diag.dat\n");
    fprintf(stderr, "ERROR: could not open psi-diag.dat\n");
    exit(EXIT_FAILURE);
  }

  fread(*psitot, sizeof(double) , tot_sz, pf);

  fclose(pf);

  /************************************************************/
  /*******************   CALC ALL ENERGIES   ******************/
  /************************************************************/

  // printf("Computing energies of all states | %s", get_time()); fflush(0);
  // energy_all(*psitot, ist->mn_states_tot, pot_local, LS, nlc, nl, ksqr, eig_vals, ist, par, flag, parallel);

  return;
}

