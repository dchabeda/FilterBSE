#include "mod_diag.h"

/*****************************************************************************/

void run_periodic_postfilter(
    double *psi_rank,
    double *pot_local,
    xyz_st *R,
    grid_st *grid,
    vector *G_vecs,
    vector *k_vecs,
    zomplex *LS,
    nlc_st *nlc,
    long *nl,
    double *ksqr,
    double *eig_vals,
    double *sigma_E,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel)
{
  /*******************************************************************
   * Periodic post-filter pipeline. Each k-group processes the        *
   * k-points it owns one at a time: gather the group's filtered      *
   * states to the master (k_rank == 0), then on the master           *
   * orthogonalize (with optional time reversal), diagonalize and     *
   * compute sigma_E at that k-point, and write eval-[k].dat /         *
   * psi-[k].dat. The k-point matching par->diag_k_idx also writes     *
   * output.dat for downstream BSE.                                    *
   *                                                                  *
   * The master's per-k buffer psitot_k holds t_rev_factor *           *
   * mn_states_per_k states (room for the time-reversed copies); the   *
   * SVD then reduces this to an adaptive `cutoff` per k-point.        *
   ********************************************************************/

  const int master = (parallel->k_rank == 0);
  const long stlen = ist->complex_idx * ist->nspinngrid;
  const long prs = ist->n_states_per_rank * stlen;     // per-rank elems for one k
  const long n0 = ist->mn_states_per_k;                // gathered states per k (across group)
  const long kbuf = (long)par->t_rev_factor * n0 * stlen; // master buffer per k

  const int do_time_rev = ((1 == flag->useSpinors) && (1 != flag->noTimeRev));

  double init_clock = (double)clock();
  double init_wall = (double)time(NULL);

  if (parallel->mpi_rank == 0)
  {
    write_separation(stdout, "T");
    printf("\n5-7. PERIODIC ORTHO / DIAG / SIGMA (per k-point) | %s\n", get_time());
    write_separation(stdout, "B");
    fflush(stdout);
  }

  double *psitot_k = NULL;
  if (master)
    ALLOCATE(&psitot_k, kbuf, "psitot_k");

  for (int ikl = 0; ikl < parallel->n_k_local; ikl++)
  {
    const int ik_global = parallel->k_local_to_global[ikl];

    // 1. Gather this k-point's filtered states from the group into psitot_k
    gather_mpi_filt_k(&psi_rank[ikl * prs], psitot_k, prs, stlen,
                      ist->n_states_per_rank, parallel->k_comm, 0);

    if (!master)
      continue;

    // 2. Orthogonalize (optionally time-reverse to double the basis first)
    long cols = n0;
    if (do_time_rev)
    {
      time_reverse_all(psitot_k, &psitot_k[n0 * stlen], n0, ist, parallel);
      cols = (long)par->t_rev_factor * n0;
    }

    long cutoff;
    if (1 == flag->isComplex)
      cutoff = ortho_cplx((MKL_Complex16 *)psitot_k, grid->dv, cols, ist, par, flag, parallel);
    else
      cutoff = ortho_real(psitot_k, grid->dv, cols, ist, par, flag, parallel);

    normalize_all(psitot_k, cutoff, ist, par, flag, parallel);

    printf("  ik = %d: %ld filtered -> %ld orthogonal states\n", ik_global, cols, cutoff);
    fflush(stdout);

    // 3. Diagonalize the Hamiltonian at this k-point (writes eig_vals[0, cutoff))
    diag_H(psitot_k, pot_local, G_vecs, k_vecs, LS, nlc, nl, ksqr, eig_vals,
           cutoff, ik_global, ist, par, flag, parallel);

    normalize_all(psitot_k, cutoff, ist, par, flag, parallel);

    // 4. Eigenvalue variance (ghost-state diagnostic), sigma_E[0, cutoff)
    calc_sigma_E_k(psitot_k, pot_local, G_vecs, k_vecs[ik_global], grid, LS, nlc, nl,
                   sigma_E, cutoff, ist, par, flag);

    // 5. Per-k output: eval-[k].dat (always), psi-[k].dat, and output.dat for the
    //    BSE reference k-point (par->diag_k_idx).
    char name[64];
    sprintf(name, "eval-%d.dat", ik_global);
    write_eval_dat(eig_vals, sigma_E, cutoff, name);

    if ((1 == flag->getAllStates) || (1 == flag->printPsiDiag))
    {
      sprintf(name, "psi-%d.dat", ik_global);
      write_psi_dat(psitot_k, cutoff * stlen, name);
    }

    if ((ik_global == par->diag_k_idx) && (0 != flag->saveOutput))
    {
      // save_output writes ist->mn_states_tot states; set it to this k's count
      ist->mn_states_tot = cutoff;
      save_output("output.dat", psitot_k, eig_vals, sigma_E, R, grid, ist, par, flag, parallel);
    }
  } // end local-k loop

  if (master)
    free(psitot_k);

  if (parallel->mpi_rank == 0)
  {
    printf("\ndone with periodic post-filter, CPU time (sec) %g, wall run time (sec) %g\n",
           ((double)clock() - init_clock) / (double)(CLOCKS_PER_SEC), (double)time(NULL) - init_wall);
    fflush(stdout);
  }

  return;
}

/*****************************************************************************/

void mod_diag(
    double *psitot,
    double *pot_local,
    double *eig_vals,
    double *sigma_E,
    grid_st *grid,
    vector *G_vecs,
    vector *k_vecs,
    zomplex *LS,
    nlc_st *nlc,
    long *nl,
    zomplex *an,
    double *zn,
    double *ene_targets,
    double *ksqr,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel)
{

  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  FILE *pf;

  const int mpir = parallel->mpi_rank;

  double init_clock;
  double init_wall;

  /************************************************************/
  /*******************  RESTART FROM DIAG   *******************/
  /************************************************************/

  if (flag->restartFromCheckpoint == 2)
  {
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

  if (mpir == 0)
  {
    write_separation(stdout, "T");
    printf("\n6. DIAGONALIZING HAMILTONIAN | %s\n", get_time());
    write_separation(stdout, "B");
    fflush(stdout);
  }

  init_clock = (double)clock();
  init_wall = (double)time(NULL);

  if (1 == flag->periodic)
  {
    /*** Periodic path: serial (rank-0) construction & diagonalization of the ***/
    /*** k-dependent Hamiltonian, followed by the eigenvalue variance using   ***/
    /*** the k-aware operator. (calc_sigma_E_k replaces the mod_sigma module, ***/
    /*** which main.c skips when flag->periodic is set.) psitot only lives on ***/
    /*** rank 0 after the gather, so this whole block is guarded by mpir==0.  ***/
    if (mpir == 0)
    {
      diag_H(psitot, pot_local, G_vecs, k_vecs, LS, nlc, nl, ksqr, eig_vals, ist->mn_states_tot, par->diag_k_idx, ist, par, flag, parallel);
      normalize_all(&psitot[0], ist->mn_states_tot, ist, par, flag, parallel);

      printf("\ndone calculating Hmat, CPU time (sec) %g, wall run time (sec) %g\n",
             ((double)clock() - init_clock) / (double)(CLOCKS_PER_SEC), (double)time(NULL) - init_wall);
      fflush(stdout);

      /************************************************************/
      /*******************     CALC SIGMA E     *******************/
      /*** standard deviation of the eigenvalues, used to check ***/
      /*** for ghost states                                     ***/
      /************************************************************/

      write_separation(stdout, "T");
      printf("\n7. CALCULATING VARIANCE OF EIGENVALUES | %s\n", get_time());
      write_separation(stdout, "B");
      fflush(stdout);

      calc_sigma_E_k(psitot, pot_local, G_vecs, k_vecs[par->diag_k_idx], grid, LS, nlc, nl, sigma_E, ist->mn_states_tot, ist, par, flag);
    }
  }
  else
  {
    // Construct Hamiltonian on single rank
    if ((0 == flag->MPIDiag) && (0 == mpir))
    {
      diag_H(psitot, pot_local, G_vecs, k_vecs, LS, nlc, nl, ksqr, eig_vals, ist->mn_states_tot, par->diag_k_idx, ist, par, flag, parallel);
    }
    // Construct and diagonalize Hamiltonian with distributed implementation
    else if (1 == flag->MPIDiag)
    {
      diag_H_mpi(psitot, pot_local, LS, nlc, nl, ksqr, eig_vals, ist, par, flag, parallel);
    }

    if (0 == mpir)
    {
      normalize_all(&psitot[0], ist->mn_states_tot, ist, par, flag, parallel);

      printf("\ndone calculating Hmat, CPU time (sec) %g, wall run time (sec) %g\n",
             ((double)clock() - init_clock) / (double)(CLOCKS_PER_SEC), (double)time(NULL) - init_wall);
      fflush(stdout);

      if (1 == flag->printPsiDiag)
      {
        pf = fopen("psi-diag.dat", "w");
        fwrite(psitot, ist->mn_states_tot * ist->complex_idx, ist->nspinngrid * sizeof(double), pf);
        fclose(pf);
      }
    }
  }

  return;
}
