#include "mod_mem.h"

void mod_mem_alloc(
    double **psi_rank,
    zomplex **psi,
    zomplex **phi,
    double **pot_local,
    zomplex **LS,
    nlc_st **nlc,
    long **nl,
    double **SO_projectors,
    zomplex **an,
    double **zn,
    double **eig_vals,
    double **sigma_E,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel)
{

  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  const int mpir = parallel->mpi_rank;

  const int nj = ist->n_j_ang_mom;

  long eig_sz = par->t_rev_factor * ist->mn_states_tot;

  if (1 == flag->periodic)
  {
    // psi_rank_size was already computed for the grouped layout in
    // setup_k_communicators (it includes the per-group n_k_local factor).
    // eig_vals/sigma_E on a group master hold the per-k upper bound across its
    // local k-points: n_k_local * t_rev_factor * mn_states_per_k.
    eig_sz = (long)parallel->n_k_local * par->t_rev_factor * ist->mn_states_per_k;
  }
  else
  {
    ist->psi_rank_size = ist->n_states_per_rank * ist->nspinngrid * ist->complex_idx;
  }

  /************************************************************/
  /*******************    ALLOCATE MEMORY   *******************/
  /************************************************************/

  if (mpir == 0)
    printf("\nAllocating memory for pot, psi, eig_vals...\n");

  // memory allocation for the spin-orbit/NL potential
  if ((flag->SO == 1) || (flag->NL == 1))
  {
    ALLOCATE(SO_projectors, PROJ_LEN * ist->nproj, "SO_projectors");
    ALLOCATE(LS, nj * nj, "LS");

    ALLOCATE(nlc, ist->n_NL_atoms * ist->n_NL_gridpts, "nlc");
    ALLOCATE(nl, ist->natoms, "nl");
  }

  // Wavefunction-type objects
  // The value of isComplex is 1 if wavefunctions are real valued, 2 if functions are complex valued
  // The factor of par->t_rev_factor (2 w spinors, 1 w/o spinors) in the psitot memory allocation
  // is because we time reverse the spinors to get double the orthogonal states
  // Spinor calculations are 8 times more memory intensive than scalar calculations

  ALLOCATE(psi_rank, ist->psi_rank_size, "psi_rank");
  ALLOCATE(psi, ist->nspinngrid, "psi in mod_mem");
  ALLOCATE(phi, ist->nspinngrid, "phi in mod_mem");
  ALLOCATE(pot_local, 2 * ist->ngrid, "pot_local in mod_mem");

  // For Newton interpolation coefficients.
  // In the periodic case the spectrum range is k-dependent, so a separate set of
  // coefficients is stored per k-point. The arrays hold n_k_pts contiguous blocks:
  // an has block ncheby*m_states_per_filter and zn has block ncheby, ordered by
  // global k index (kpt 0, kpt 1, ...). gen_newton_coeff_k populates each block.
  long n_coeff_blocks = (1 == flag->periodic) ? (long)ist->n_k_pts : 1;
  ALLOCATE(an, n_coeff_blocks * ist->ncheby * ist->m_states_per_filter, "an cheby");
  ALLOCATE(zn, n_coeff_blocks * ist->ncheby, "zn cheby");

  // Per-k Hamiltonian spectrum bounds, one entry per global k-point.
  if (1 == flag->periodic)
  {
    ALLOCATE(&par->Emin, ist->n_k_pts, "par Emin");
    ALLOCATE(&par->Emax, ist->n_k_pts, "par Emax");
  }

  // the quasiparticle energies and standard deviations
  ALLOCATE(eig_vals, eig_sz, "eig_vals");
  ALLOCATE(sigma_E, eig_sz, "sigma_E");

  if (mpir == 0)
    printf("\tdone allocating memory.\n");
  fflush(stdout);
}
