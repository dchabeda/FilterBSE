#include "mod_mem.h"

void mod_mem_alloc(
  double**      psi_rank,
  zomplex**     psi,
  zomplex**     phi,
  double**      pot_local,
  zomplex**     LS,
  nlc_st**      nlc,
  long**        nl,
  double**      SO_projectors,
  zomplex**     an,
  double**      zn,
  double**      eig_vals,
  double**      sigma_E,
  index_st*     ist,
  par_st*       par,
  flag_st*      flag,
  parallel_st*  parallel){

  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/
  
  const int mpir = parallel->mpi_rank;
  
  const int nj = ist->n_j_ang_mom;
  
  ist->psi_rank_size = ist->n_states_per_rank*ist->nspinngrid*ist->complex_idx;
  
  
  /************************************************************/
  /*******************    ALLOCATE MEMORY   *******************/
  /************************************************************/

  if (mpir == 0) printf("\nAllocating memory for pot, psi, eig_vals...\n");
  
  // memory allocation for the spin-orbit/NL potential 
  if (flag->SO == 1){
    ALLOCATE(SO_projectors, PROJ_LEN * ist->nproj, "SO_projectors");
    ALLOCATE(LS, nj * nj, "LS");
    
    if (flag->NL == 1){
      ALLOCATE(nlc, ist->n_NL_atoms*ist->n_NL_gridpts, "nlc");
      ALLOCATE(nl, ist->natoms, "nl");
    }
  }
  
  // Wavefunction-type objects
  // The value of isComplex is 1 if wavefunctions are real valued, 2 if functions are complex valued
  // The factor of par->t_rev_factor (2 w spinors, 1 w/o spinors) in the psitot memory allocation
  // is because we time reverse the spinors to get double the orthogonal states
  // Spinor calculations are 8 times more memory intensive than scalar calculations
  
  ALLOCATE(psi_rank, ist->psi_rank_size, "psi_rank");
  ALLOCATE(psi, ist->nspinngrid, "psi in mod_mem");
  ALLOCATE(phi, ist->nspinngrid, "phi in mod_mem");
  ALLOCATE(pot_local, ist->ngrid, "pot_local in mod_mem");
  
  // For Newton interpolation coefficients
  ALLOCATE(an, ist->ncheby*ist->m_states_per_filter, "an cheby");
  ALLOCATE(zn, ist->ncheby, "zn cheby");
  
  // the quasiparticle energies and standard deviations
  ALLOCATE(eig_vals, par->t_rev_factor*ist->mn_states_tot, "eig_vals");
  ALLOCATE(sigma_E, par->t_rev_factor*ist->mn_states_tot, "sigma_E");
  
  
  if (mpir == 0) printf("\tdone allocating memory.\n"); fflush(stdout);

}

