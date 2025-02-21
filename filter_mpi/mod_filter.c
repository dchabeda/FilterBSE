#include "mod_filter.h"


/**************************************************************/

/**************************************************************/
/********************  RUN FILTER MODULE   ********************/
/**************************************************************/

void mod_filter(
  double*       psi_rank,
  zomplex*      psi,
  zomplex*      phi,
  double*       pot_local,
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
  parallel_st*  parallel){

  
  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  const int cidx = ist->complex_idx;
  const int nspngrd = ist->nspinngrid;
  const int ntVB = par->n_targets_VB;
  const int ntCB = par->n_targets_CB;
  const int mpir = parallel->mpi_rank;

  long rand_seed;

  char    *top;
  char    *bottom; 
  top    = malloc(2*sizeof(top[0]));
  bottom = malloc(2*sizeof(bottom[0]));
  strcpy(top, "T\0");
  strcpy(bottom, "B\0");
  

  /************************************************************/
  /*******************   CALC ENERGY RANGE  *******************/
  /************************************************************/


  if (mpir == 0){
    write_separation(stdout, top);
    printf("\n2. CALCULATING HAMILTONIAN ENERGY RANGE | %s\n", get_time());
    write_separation(stdout, bottom); fflush(stdout);
  }

  time_t init_clock = (double) clock(); 
  time_t init_wall = (double) time(NULL);
  
  //
  //
  get_energy_range(
    psi, phi, pot_local, grid, LS, nlc, nl, ksqr, ist, par, flag, parallel);
  //
  //

  if (mpir == 0){
    printf("\nDone w calc energy range, CPU time (sec) %g, wall run time (sec) %g\n",
    ((double)clock()-init_clock)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-init_wall); 
    fflush(stdout);
  }

  
  /************************************************************/
  /*******************   GEN CHEBY COEFFS   *******************/
  /************************************************************/


  if (mpir == 0){
    write_separation(stdout, top);
    printf("\n3. GENERATING COEFFICIENTS | %s\n", get_time());
    write_separation(stdout, bottom); fflush(stdout);
  }

  par->dt = sqr((double)(ist->ncheby) / (2.5*par->dE));

  // 
  // 
  gen_newton_coeff(an, zn, ene_targets, ist, par, parallel);
  // 
  // 

  if (mpir == 0){
    printf("\n  ncheby = %ld dt = %g dE = %g\n", ist->ncheby, par->dt, par->dE);
    printf("  Energy width, sigma, of filter function = %.6g a.u.\n", sqrt(1 / (2*par->dt)));
    printf("  Suggested max span of spectrum for filtering = %.6g a.u.\n", ist->m_states_per_filter * sqrt(1 / (2*par->dt)));
    printf("  Requested span of spectrum to filter = %.6g a.u.\n", 
      (ene_targets[ntVB-1] - ene_targets[0]) + (ene_targets[ntVB+ntCB-1] - ene_targets[ntVB])); 
    fflush(stdout);
  }
  

  /************************************************************/
  /*******************  INIT FILTER STATES  *******************/
  /************************************************************/
  

  if (mpir == 0) printf("\nInitializing random filter states\n");

  // Generate the random seed for initial states
  // or read it from input.par if preset (for debugging)
  if (0 == flag->setSeed){
    Randomize();  rand_seed = -random() + mpir;
  } else {
    rand_seed = - par->rand_seed;
  }
  
  // Gen initial random wavefunctions for each filter cycle
  // Out: array of len n_filter_cycles * m_states_per_filter
  // every block of len [m_states] has the same random psi
  
  // 
  // 
  init_filter_states(psi_rank, psi, grid, &rand_seed, ist, par, flag, parallel);
  // 
  // 


  /************************************************************/
  /*******************   TIME HAMILTONIAN   *******************/
  /************************************************************/
  

  if (1 == flag->timeHamiltonian){
    if (mpir == 0) {
      printf("\nTiming Hamiltonian operator... | %s\n", get_time()); 
      fflush(0);
    }
    
    // Initialize state on which to act Hamiltonian
    memcpy(&phi[0], &psi_rank[0], cidx * nspngrd * sizeof(psi_rank[0]));
  
    //
    //
    time_hamiltonian(phi, psi, pot_local, LS, nlc, nl, ksqr, ist, par, flag, parallel);
    //
    //

    if (mpir == 0) printf("Done timing Hamiltonian | %s\n", get_time()); fflush(0);
  }


  /************************************************************/
  /*******************   RUN FILTER CYCLE   *******************/
  /************************************************************/

  /************************************************************/
  /*******  For [n_filter_cycles] random initial states *******/
  /*******   perform filtering at [m_states] E targets  *******/
  

  if (mpir == 0){
    write_separation(stdout, top);
    printf("\n4. RUN FILTER CYCLE | %s\n", get_time());
    write_separation(stdout, bottom);
    printf("\n  4.1 Running filter cycle\n");
  }
  
  init_clock = (double)clock(); 
  init_wall = (double)time(NULL);

  // 
  //
  run_filter_cycle(
    psi_rank, pot_local, LS, nlc, nl, ksqr, an, zn, 
    ene_targets, grid, ist, par, flag, parallel);
  // 
  // 

  // Ensure all ranks synchronize here
  // MPI_Barrier(MPI_COMM_WORLD);

  if (mpir == 0){
    printf("\ndone calculating filter, CPU time (sec) %g, wall run time (sec) %g\n",
    ((double)clock()-init_clock)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-init_wall); 
    fflush(stdout);
  }

  free(top); free(bottom);

  
  return;
}
