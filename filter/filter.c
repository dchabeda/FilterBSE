#include "fd.h"
#include "aux.h"
/*****************************************************************************/

void run_filter_cycles(
    double *psitot,
    double *pot_local,
    atom_info *atom,
    nlc_st *nlc,
    long *nl,
    double *ksqr,
    zomplex *an,
    double *zn,
    double *ene_targets,
    grid_st *grid,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel)
{
  /*******************************************************************
   * This function runs a filter cycle on m ene_targets with one of   *
   * the initial random states as the starting point.                 *
   * inputs:                                                          *
   *  [psims] arr that will hold all ms filtered wavefunctions        *
   *  [pot_local] ngrid-long arr holding the value of the local pot   *
   *  [nlc] nlc struct holding values for computing SO and NL pots    *
   *  [nl] natom-long arr holding the number of NL gridpts per atom   *
   *  [ksqr] ngrid-long arr holding the values of k^2 for KE calc     *
   *  [an] the Newton interpolation coefficients for the filter func  *
   *  [zn] Chebyshev polynomial support points                        *
   *  [ene_targets] target energies where filter funcs are centered   *
   *  [thread_id] filtering of each random state occurs on one thread *
   *  [jns] index of filter cycle                                     *
   *  [grid] struct holding the grid and grid parameters (dx, xmin...)*
   *  [ist] ptr to counters, indices, and lengths                     *
   *  [par] ptr to par_st holding VBmin, VBmax... params              *
   *  [flag] ptr to flag_st holding job flags                         *
   *  [parallel] holds options for parallelization                    *
   * outputs: void                                                    *
   ********************************************************************/

  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  const long ms = ist->m_states_per_filter;
  const long stlen = ist->nspinngrid * ist->complex_idx;

  //              Counters and array indexing
  long jns;
  long jms;
  long jstate;
  long ns_block;

  //              File I/O
  FILE *pf;
  char fileName[100];

  //              Filtered state energies
  double *ene_filters;

  //              Arrays for Hamiltonian evaluation
  zomplex *psi;
  zomplex *phi;

  //              FFT plans and arrays
  // fftw_init_threads();
  // fftw_plan_with_nthreads(par->nthreads);
  fftw_plan_loc *planfw;
  fftw_plan_loc *planbw;
  fftw_complex *fftwpsi;
  long fft_flags = FFTW_MEASURE;

  /************************************************************/
  /*******************  INIT ARRAYS & MEM   *******************/
  /************************************************************/

  // Create FFT structs and plans for Fourier transform
  fftwpsi = fftw_malloc(parallel->nthreads * sizeof(fftw_complex) * ist->ngrid);
  planfw = (fftw_plan_loc *)malloc(parallel->nthreads * sizeof(fftw_plan_loc));
  planbw = (fftw_plan_loc *)malloc(parallel->nthreads * sizeof(fftw_plan_loc));
  for (jns = 0; jns < parallel->nthreads; jns++)
  {
    jstate = jns * ist->ngrid;
    planfw[jns] = fftw_plan_dft_3d(ist->nz, ist->ny, ist->nx, &fftwpsi[jstate], &fftwpsi[jstate], FFTW_FORWARD, fft_flags);
    planbw[jns] = fftw_plan_dft_3d(ist->nz, ist->ny, ist->nx, &fftwpsi[jstate], &fftwpsi[jstate], FFTW_BACKWARD, fft_flags);
  }

  ALLOCATE(&psi, parallel->nthreads * ist->nspinngrid, "psi in filter_mpi");
  ALLOCATE(&phi, parallel->nthreads * ist->nspinngrid, "phi in filter_mpi");
  ALLOCATE(&ene_filters, ist->mn_states_tot, "ene_filters");

  /************************************************************/
  /*******************   BEGIN FILTERING    *******************/
  /************************************************************/
  /* In MPI implementation, every MPI rank has n_filt/n_ranks */
  /* initial random states to filter. OpenMP parallelization  */
  /* used to accelerate Hamiltonian application as well as    */
  /* propagating at each energy target                        */
  /************************************************************/

  char str[40];
  int n_loop_trips = (int)(ist->n_filter_cycles / parallel->nthreads);
  int loop_count;
  int thread_id;
  loop_count = 0;
  thread_id = 0;
  // Loop over all of the random states handled by this mpi-rank
  // *
  // *

  omp_set_num_threads(parallel->nthreads);
#pragma omp parallel for private(jns, thread_id)
  for (jns = 0; jns < ist->n_filter_cycles; jns++)
  {
    FILE *pft;
    thread_id = omp_get_thread_num();
    // Keep track of how many filter iterations have taken place
    if (thread_id == 0)
    {
      printf("\n  Filter cycle %d / %d | %s\n", loop_count + 1, n_loop_trips, get_time());
      loop_count++;
      fflush(0);
    }

    // Run the filter cycle (apply Hamiltonian ncheby times)
    filter_cycle(
        psitot, jns, &psi[thread_id * ist->nspinngrid], &phi[thread_id * ist->nspinngrid], pot_local, atom, nlc, nl, ksqr, an, zn,
        ene_targets, grid, ist, par, flag, parallel, planfw[thread_id], planbw[thread_id], &fftwpsi[thread_id * ist->ngrid], thread_id);

    // Print out all ms states from initial random state jns
    if (1 == flag->printPsiFilt)
    {
      // double *rho;
      // char cubstr[20];
      long jgrid, jgrid_real, jgrid_imag;
      // rho = calloc(ist->ngrid, sizeof(double));
      ns_block = jns * (ms * stlen);

      sprintf(str, "psi-filt-%ld-%d.dat", jns, thread_id);
      pft = fopen(str, "w");

      for (jms = 0; jms < ms; jms++)
      {
        jstate = ns_block + jms * stlen;
        fwrite(&psitot[jstate], sizeof(double), stlen, pft);
      }

      fclose(pft);

      // free(rho);
    } // end of printPsiFilt
  } // end of jns

  free(psi);
  free(phi);
  fftw_free(fftwpsi);
  for (jns = 0; jns < parallel->nthreads; jns++)
  {
    fftw_destroy_plan(planfw[jns]);
    fftw_destroy_plan(planbw[jns]);
  }

  /***********************************************************************/

  // normalize the states and get their energies
  printf("\n  4.3 Normalizing filtered states\n");
  fflush(stdout);

  normalize_all(psitot, ist, par, flag, parallel);

  // Get the energy of all the filtered states
  printf("\n  4.4 Computing the energies of all filtered states\n");
  fflush(stdout);

  energy_all(
      psitot, pot_local, atom,
      nlc, nl, ksqr, ene_filters, ist, par, flag, parallel);

  // Print out energies of all filters

  for (jns = 0; jns < ist->n_filter_cycles; jns++)
  {

    sprintf(fileName, "ene-filt-jns-%ld.dat", jns);
    pf = fopen(fileName, "w");

    for (jms = 0; jms < ms; jms++)
    {
      jstate = jns * ms + jms;
      fprintf(pf, "%ld %.16g %.16g\n", jms, ene_filters[jstate], ene_targets[jms]);
    }

    fclose(pf);
  }

  free(ene_filters);

  return;
}

/*****************************************************************************/

void filter_cycle(
    double *psitot,
    long jns,
    zomplex *psi,
    zomplex *phi,
    double *pot_local,
    atom_info *atom,
    nlc_st *nlc,
    long *nl,
    double *ksqr,
    zomplex *an,
    double *zn,
    double *ene_targets,
    grid_st *grid,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel,
    fftw_plan_loc planfw,
    fftw_plan_loc planbw,
    fftw_complex *fftwpsi,
    int thread_id)
{
  /*******************************************************************
   * This function runs a filter cycle on m ene_targets with one of   *
   * the initial random states as the starting point.                 *
   * inputs:                                                          *
   *  [psims] arr that will hold all ms filtered wavefunctions        *
   *  [pot_local] ngrid-long arr holding the value of the local pot   *
   *  [nlc] nlc struct holding values for computing SO and NL pots    *
   *  [nl] natom-long arr holding the number of NL gridpts per atom   *
   *  [ksqr] ngrid-long arr holding the values of k^2 for KE calc     *
   *  [an] the Newton interpolation coefficients for the filter func  *
   *  [zn] Chebyshev polynomial support points                        *
   *  [ene_targets] target energies where filter funcs are centered   *
   *  [thread_id] filtering of each random state occurs on one thread *
   *  [jns] index of filter cycle                                     *
   *  [grid] struct holding the grid and grid parameters (dx, xmin...)*
   *  [ist] ptr to counters, indices, and lengths                     *
   *  [par] ptr to par_st holding VBmin, VBmax... params              *
   *  [flag] ptr to flag_st holding job flags                         *
   *  [parallel] holds options for parallelization                    *
   * outputs: void                                                    *
   ********************************************************************/

  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  //             Array indexing
  long jms;
  long ns_block;
  long jstate;
  long jgrid;
  long jgrid_real;
  long jgrid_imag;
  long jc;
  long ncjms;

  double zn_jc;

  const int is_cmplx = flag->isComplex;
  const long ms = ist->m_states_per_filter;
  const long stlen = ist->nspinngrid * ist->complex_idx;

  //

  ns_block = jns * (ms * stlen);

  /************************************************************/
  /*******************   READ RANDOM PSI    *******************/
  /************************************************************/

  // Populate psi with the random initial wavefunc, jns
  if (1 == is_cmplx)
  {
    for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++)
    {
      jgrid_real = jgrid * ist->complex_idx;
      jgrid_imag = jgrid_real + 1;

      psi[jgrid].re = psitot[ns_block + jgrid_real];
      psi[jgrid].im = psitot[ns_block + jgrid_imag];
    }
  }
  else
  {
    for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++)
    {
      psi[jgrid].re = psitot[ns_block + jgrid];
      psi[jgrid].im = 0.0;
    }
  }

  /************************************************************/
  /*******************   TERM 0 OF CHEBY    *******************/
  /************************************************************/

  for (jms = 0; jms < ms; jms++)
  {
    jstate = ns_block + jms * stlen;
    ncjms = ist->ncheby * jms;

    if (1 == is_cmplx)
    {
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++)
      {
        jgrid_real = ist->complex_idx * jgrid;
        jgrid_imag = jgrid_real + 1;

        // Real part
        psitot[jstate + jgrid_real] =
            an[ncjms + 0].re * psi[jgrid].re - an[ncjms + 0].im * psi[jgrid].im;
        // Imag part
        psitot[jstate + jgrid_imag] =
            an[ncjms + 0].re * psi[jgrid].im + an[ncjms + 0].im * psi[jgrid].re;
      }
    }
    else
    {
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++)
      {
        psitot[jstate + jgrid] =
            an[ncjms + 0].re * psi[jgrid].re - an[ncjms + 0].im * psi[jgrid].im;
      }
    }
  }

  /************************************************************/
  /*******************   CHEBY EXPANSION    *******************/
  /************************************************************/

  for (jc = 1; jc < ist->ncheby; jc++)
  {
    zn_jc = zn[jc - 1];

    // Apply Hamiltonian and renormalize eigs for Cheby stability
    hnorm(
        psi, phi, pot_local, atom, nlc, nl, ksqr, zn_jc,
        ist, par, flag, parallel, planfw, planbw, fftwpsi);

    // Propagate the state at all energy targets based on coeffs
    // #pragma omp parallel for private(jms, ncjms, jstate, jgrid, jgrid_real, jgrid_imag)
    for (jms = 0; jms < ms; jms++)
    {
      ncjms = ist->ncheby * jms;
      jstate = ns_block + jms * stlen;

      if (1 == is_cmplx)
      {
        for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++)
        {
          jgrid_real = jgrid * ist->complex_idx;
          jgrid_imag = jgrid_real + 1;

          // Real part
          // #pragma omp atomic
          psitot[jstate + jgrid_real] +=
              (an[ncjms + jc].re * psi[jgrid].re - an[ncjms + jc].im * psi[jgrid].im);
          // Imag part
          // #pragma omp atomic
          psitot[jstate + jgrid_imag] +=
              (an[ncjms + jc].re * psi[jgrid].im + an[ncjms + jc].im * psi[jgrid].re);
        }
      }
      else
      {
        for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++)
        {
          // #pragma omp atomic
          psitot[jstate + jgrid] += an[ncjms + jc].re * psi[jgrid].re - an[ncjms + jc].im * psi[jgrid].im;
        }
      }
    }

    if (thread_id == 0)
    {
      if ((1 == jc) || (0 == (jc % ((long)(ist->ncheby / 4 + 1)))) || ((ist->ncheby - 1) == jc))
      {
        print_progress_bar(jc, ist->ncheby);
      }
    }
  }

  return;
}

/*******************************************************************/

void hnorm(
    zomplex *psi_out,
    zomplex *psi_tmp,
    double *pot_local,
    atom_info *atom,
    nlc_st *nlc,
    long *nl,
    double *ksqr,
    double zn,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel,
    fftw_plan_loc planfw,
    fftw_plan_loc planbw,
    fftw_complex *fftwpsi)
{

  long jgrid;

  memcpy(&psi_tmp[0], &psi_out[0], ist->nspinngrid * sizeof(psi_tmp[0]));

  hamiltonian(
      psi_out, psi_tmp, pot_local, atom, nlc, nl, ksqr, ist,
      par, flag, planfw, planbw, fftwpsi);

  // Scale eigenvalues to remain within disk of convergence of expansion
  for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++)
  {
    /*** par->dE_1 = 4.0 / par->dE and therefore I don't multiply by 4 ***/
    psi_out[jgrid].re = par->dE_1 * psi_out[jgrid].re - (2.0 + zn + par->Vmin * par->dE_1) * psi_tmp[jgrid].re;
    psi_out[jgrid].im = par->dE_1 * psi_out[jgrid].im - (2.0 + zn + par->Vmin * par->dE_1) * psi_tmp[jgrid].im;
  }

  return;
}

/*****************************************************************************/

// void time_hamiltonian(
//     zomplex *psi_out,
//     zomplex *psi_tmp,
//     double *pot_local,
//     zomplex *LS,
//     nlc_st *nlc,
//     long *nl,
//     double *ksqr,
//     index_st *ist,
//     par_st *par,
//     flag_st *flag,
//     parallel_st *parallel)
// {
//   /*******************************************************************
//    * This function applies the Hamiltonian onto a state               *
//    * inputs:                                                          *
//    *  [psi_tmp] nspinngrid-long arr of orig. wavefnc                  *
//    *  [psi_out] nspinngrid-long arr to hold |psi_out> = H|psi_tmp>    *
//    *  [pot_local] ngrid-long arr holding the value of the local pot   *
//    *  [nlc] nlc struct holding values for computing SO and NL pots    *
//    *  [nl] natom-long arr holding the number of NL gridpts per atom   *
//    *  [ksqr] ngrid-long arr holding the values of k^2 for KE calc     *
//    *  [ist] ptr to counters, indices, and lengths                     *
//    *  [par] ptr to par_st holding VBmin, VBmax... params              *
//    *  [flag] ptr to flag_st holding job flags                         *
//    *  [planfw] FFTW3 plan for executing 3D forward DFT                *
//    *  [planfw] FFTW3 plan for executing 3D backwards DFT              *
//    *  [fftwpsi] location to store outcome of Fourier transform        *
//    * outputs: void                                                    *
//    ********************************************************************/

//   const int mpir = parallel->mpi_rank;

//   struct timespec start, end;
//   int jspin, j, jtmp;
//   int n_iter = 150;

//   omp_set_num_threads(par->ham_threads);

//   fftw_init_threads();
//   fftw_plan_with_nthreads(par->ham_threads);
//   fftw_plan_loc planfw, planbw;
//   fftw_complex *fftwpsi;
//   long fft_flags = 0;

//   fftwpsi = fftw_malloc(sizeof(fftw_complex) * ist->ngrid);
//   /*** initialization for the fast Fourier transform ***/
//   planfw = fftw_plan_dft_3d(ist->nz, ist->ny, ist->nx, fftwpsi, fftwpsi, FFTW_FORWARD, fft_flags);
//   planbw = fftw_plan_dft_3d(ist->nz, ist->ny, ist->nx, fftwpsi, fftwpsi, FFTW_BACKWARD, fft_flags);

//   // Copy psi_out into psi_tmp
//   memcpy(&psi_tmp[0], &psi_out[0], ist->nspinngrid * sizeof(psi_tmp[0]));

//   // Calculate the action of the kinetic energy part of the Hamiltonian on psi_tmp: |psi_out> = T|psi_tmp>
//   // Warmup runs to avoid including caching time, optimizations, innitial overhead etc.
//   for (j = 0; j < 10; j++)
//   {
//     for (jspin = 0; jspin < ist->nspin; jspin++)
//     {
//       kinetic(&psi_out[jspin * ist->ngrid], ksqr, planfw, planbw, fftwpsi, ist); // spin up/down
//     }
//   }
//   clock_gettime(CLOCK_MONOTONIC, &start);
//   for (j = 0; j < n_iter; j++)
//   {
//     for (jspin = 0; jspin < ist->nspin; jspin++)
//     {
//       kinetic(&psi_out[jspin * ist->ngrid], ksqr, planfw, planbw, fftwpsi, ist); // spin up/down
//     }
//   }
//   clock_gettime(CLOCK_MONOTONIC, &end);
//   double elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
//   if (mpir == 0)
//     printf("\tKinetic energy: %.4g (msec)\n", (elapsed_seconds * 1000.0) / (double)n_iter);
//   fflush(0);

//   // Calculate the action of the potential operator on the wavefunction: |psi_out> = V|psi_tmp>

//   if (flag->SO == 1)
//   {
//     // Calculate |psi_out> = V_SO|psi_tmp>

//     // Warmup runs
//     for (j = 0; j < 3; j++)
//     {
//       p_spin_orbit_proj_pot(psi_out, psi_tmp, LS, nlc, nl, ist, par, par->ham_threads);
//     }
//     clock_gettime(CLOCK_MONOTONIC, &start);
//     for (j = 0; j < n_iter; j++)
//     {
//       p_spin_orbit_proj_pot(psi_out, psi_tmp, LS, nlc, nl, ist, par, par->ham_threads);
//     }
//     clock_gettime(CLOCK_MONOTONIC, &end);
//     double elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
//     if (mpir == 0)
//       printf("\tSpin-Orbit potential: %.4g (msec)\n", (elapsed_seconds * 1000.0) / (double)n_iter);
//   }

//   if (flag->NL == 1)
//   {
//     // Calculate |psi_out> += V_NL|psi_tmp>
//     for (j = 0; j < 3; j++)
//     {
//       p_nonlocal_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par, par->ham_threads);
//     }
//     clock_gettime(CLOCK_MONOTONIC, &start);
//     for (j = 0; j < n_iter; j++)
//     {
//       p_nonlocal_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par, par->ham_threads);
//     }
//     clock_gettime(CLOCK_MONOTONIC, &end);
//     elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
//     if (mpir == 0)
//       printf("\tNon-local potential: %.4g (msec)\n", (elapsed_seconds * 1000.0) / (double)n_iter);
//   }

//   // Calculate the action of the local potential energy part of the Hamiltonian on psi_tmp
//   clock_gettime(CLOCK_MONOTONIC, &start);
//   for (int i = 0; i < n_iter; i++)
//   {
//     if (1 == flag->useSpinors)
//     {
//       for (jspin = 0; jspin < 2; jspin++)
//       {
//         for (j = 0; j < ist->ngrid; j++)
//         {
//           jtmp = jspin * ist->ngrid + j;
//           psi_out[jtmp].re += (pot_local[j] * psi_tmp[jtmp].re);
//           psi_out[jtmp].im += (pot_local[j] * psi_tmp[jtmp].im);
//           // // handle spin down component

//           // psi_out[jtmp].re += (pot_local[j] * psi_tmp[jtmp].re);
//           // psi_out[jtmp].im += (pot_local[j] * psi_tmp[jtmp].im);
//         }
//       }
//     }
//     else if (0 == flag->useSpinors)
//     {
//       for (j = 0; j < ist->ngrid; j++)
//       {
//         psi_out[j].re += (pot_local[j] * psi_tmp[j].re);
//       }
//     }
//     else
//     {
//       printf("ERROR: unrecognized value for flag->useSpinors: %d\n", flag->useSpinors);
//       fprintf(stderr, "ERROR: unrecognized value for flag->useSpinors: %d\n", flag->useSpinors);
//       exit(EXIT_FAILURE);
//     }
//   }
//   clock_gettime(CLOCK_MONOTONIC, &end);
//   elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
//   if (mpir == 0)
//     printf("\tLocal potential: %.4g (msec)\n", (elapsed_seconds * 1000.0) / (double)n_iter);

//   return;
// }

/*****************************************************************************/
