#include "filter.h"

/*****************************************************************************/

void run_filter_cycles_k(
  double*       psi_rank, 
  double*       pot_local,
  vector*       G_vecs,
  vector*       k_vecs,
  grid_st*      grid,
  zomplex*      LS,
  nlc_st*       nlc, 
  long*         nl,
  zomplex*      an, 
  double*       zn, 
  double*       ene_targets, 
  index_st*     ist, 
  par_st*       par, 
  flag_st*      flag, 
  parallel_st*  parallel
  ){
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

  const int       mpir  = parallel->mpi_rank;
  const long      ms    = ist->m_states_per_filter;
  const long      stlen = ist->nspinngrid * ist->complex_idx;

  //              Counters and array indexing
  long            ik;
  long            jns;
  long            jms;
  long            jstate;
  long            ns_block;
  long            ik_block;
  long            mn_sz = ist->n_filters_per_rank * ms * stlen;

  //              File I/O
  FILE*           pf; 
  char            fileName[100];

  //              Filtered state energies
  double*         ene_filters;

  //              Arrays for Hamiltonian evaluation
  zomplex*        psi;
  zomplex*        phi;
  
  //              FFT plans and arrays
  fftw_plan_loc   planfw;
  fftw_plan_loc   planbw;
  fftw_complex*   fftwpsi;
  long            fft_flags = FFTW_MEASURE;

  vector          k;

  /************************************************************/
  /*******************  INIT ARRAYS & MEM   *******************/
  /************************************************************/

  fftw_init_threads();
  fftw_plan_with_nthreads(par->ham_threads);
  
  // Create FFT structs and plans for Fourier transform
  fftwpsi = fftw_malloc(sizeof (fftw_complex) * ist->ngrid);
  planfw = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_FORWARD,fft_flags);
  planbw = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_BACKWARD,fft_flags);
  
  ALLOCATE(&psi, ist->nspinngrid, "psi in filter_mpi"); 
  ALLOCATE(&phi, ist->nspinngrid, "phi in filter_mpi");
  ALLOCATE(&ene_filters, ist->n_states_per_rank, "ene_filters");
  
  /************************************************************/
  /*******************   BEGIN FILTERING    *******************/
  /************************************************************/
  /* In MPI implementation, every MPI rank has n_filt/n_ranks */ 
  /* initial random states to filter. OpenMP parallelization  */ 
  /* used to accelerate Hamiltonian application as well as    */ 
  /* propagating at each energy target                        */
  /************************************************************/
  
  char str[40];

  // Set the number of OMP threads for Hamiltonian eval
  omp_set_num_threads(par->ham_threads);

  // Loop over all of the random states handled by this mpi-rank
  // *
  // *
  for (ik = 0; ik < ist->n_k_pts; ik++){
    
    k = k_vecs[ik];
    printf("ik = %ld, k = %.4lg %.4lg %.4lg\n", ik, k.x, k.y, k.z);

    ik_block = ik * mn_sz;

    for (jns = 0; jns < ist->n_filters_per_rank; jns++) {
      
      // Keep track of how many filter iterations have taken place
      if (mpir == 0) {
        printf("  Random psi %ld / %ld | %s\n", jns + 1, ist->n_filters_per_rank, get_time());
        fflush(0);
      }

      // Run the filter cycle (apply Hamiltonian ncheby times)
      filter_cycle_k(
        &psi_rank[ik_block], jns, psi, phi, pot_local, G_vecs, k, grid, LS, nlc, nl, an, zn,
        ene_targets, ist, par, flag, parallel, planfw, planbw, fftwpsi
      );

      // Print out all ms states from initial random state jns
      if (1 == flag->printPsiFilt){
        ns_block = jns * (ms * stlen);

        sprintf(str, "psi-filt-%ld-%ld-%d.dat", ik, jns, mpir);
        pf = fopen(str, "w");

        for (jms = 0; jms < ms; jms++){
          jstate = ik_block + ns_block + jms * stlen;
          fwrite(&psi_rank[jstate], sizeof(double), stlen, pf);
        }

        fclose(pf);
      } // end of printPsiFilt
    } // end of jns
  
    /***********************************************************************/

    // normalize the states and get their energies
    if (mpir == 0) printf("\n  4.2 Normalizing filtered states\n"); fflush(stdout);
    
    normalize_all(&psi_rank[ik_block], ist->n_states_per_rank, ist, par, flag, parallel);

    // Get the energy of all the filtered states
    if (mpir == 0) printf("\n  4.3 Computing the energies of all filtered states\n"); fflush(stdout);
    
    energy_all_k(
      &psi_rank[ik_block], ist->n_states_per_rank, pot_local, G_vecs, k, grid, LS,
      nlc, nl, ene_filters, ist, par, flag, parallel
    );

    // Print out energies of all filters

    for (jns = 0; jns < ist->n_filters_per_rank; jns++){

      sprintf (fileName,"ene-filt-jns-%ld-%ld-%d.dat", ik, jns, mpir);
      pf = fopen(fileName , "w");

      for (jms = 0; jms < ms; jms++){
        jstate = jns * ms + jms;
        fprintf (pf,"%ld %.16g %.16g\n", jms, ene_filters[jstate], ene_targets[jms]);
      }

      fclose(pf);
    }


  } // end of ik

  free(psi); free(phi);
  fftw_destroy_plan(planfw);
  fftw_destroy_plan(planbw);
  fftw_free(fftwpsi);
  free(ene_filters);

  return;
}



/*****************************************************************************/

void filter_cycle_k(
  double*       psi_rank,
  long          jns,
  zomplex*      psi,
  zomplex*      phi,
  double*       pot_local,
  vector*       G_vecs,
  vector        k,
  grid_st*      grid,
  zomplex*      LS,
  nlc_st*       nlc, 
  long*         nl,
  zomplex*      an, 
  double*       zn, 
  double*       ene_targets,
  index_st*     ist, 
  par_st*       par, 
  flag_st*      flag, 
  parallel_st*  parallel,
  fftw_plan_loc planfw,
  fftw_plan_loc planbw,
  fftw_complex* fftwpsi
  ){
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
  long           jms; 
  long           ns_block;
  long           jstate;
  long           jgrid;
  long           jgrid_real;
  long           jgrid_imag;
  long           jc;
  long           ncjms;

  double zn_jc;

  const int      is_cmplx = flag->isComplex;
  const long     ms = ist->m_states_per_filter;
  const long     stlen = ist->nspinngrid * ist->complex_idx;

  const long     mpir = mpir;
  // 

  ns_block = jns * (ms * stlen);
  
  /************************************************************/
  /*******************   READ RANDOM PSI    *******************/
  /************************************************************/

  // Populate psi with the random initial wavefunc, jns
  if (1 == is_cmplx){
    for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
      jgrid_real = jgrid * ist->complex_idx;
      jgrid_imag = jgrid_real + 1;

      psi[jgrid].re = psi_rank[ns_block + jgrid_real];
      psi[jgrid].im = psi_rank[ns_block + jgrid_imag];
    }
  } else{
    for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
      psi[jgrid].re = psi_rank[ns_block + jgrid];
    }
  }

  /************************************************************/
  /*******************   TERM 0 OF CHEBY    *******************/
  /************************************************************/

  for (jms = 0; jms < ms; jms++){
    jstate = ns_block + jms * stlen;
    ncjms = ist->ncheby * jms;

    if (1 == is_cmplx){
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
        jgrid_real = ist->complex_idx * jgrid;
        jgrid_imag = jgrid_real + 1;

        // Real part
        psi_rank[jstate + jgrid_real] = \
          an[ncjms+0].re * psi[jgrid].re - an[ncjms+0].im * psi[jgrid].im;
        // Imag part
        psi_rank[jstate + jgrid_imag] = \
          an[ncjms+0].re * psi[jgrid].im + an[ncjms+0].im * psi[jgrid].re;
      }
    } else{
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
        psi_rank[jstate + jgrid] = \
          an[ncjms+0].re * psi[jgrid].re - an[ncjms+0].im * psi[jgrid].im;
      }
    }
  }

  /************************************************************/
  /*******************   CHEBY EXPANSION    *******************/
  /************************************************************/

  for (jc = 1; jc < ist->ncheby; jc++){
    zn_jc = zn[jc - 1];

    // Apply Hamiltonian and renormalize eigs for Cheby stability
    p_hnorm_k(
      psi, phi, pot_local, G_vecs, k, grid, LS, nlc, nl, zn_jc, 
      ist, par, flag, parallel, planfw, planbw, fftwpsi
    );

    // Propagate the state at all energy targets based on coeffs
    #pragma omp parallel for private(jms, ncjms, jstate, jgrid, jgrid_real, jgrid_imag)
    for (jms = 0; jms < ms; jms++){
      ncjms = ist->ncheby * jms;
      jstate = ns_block + jms * stlen;

      if (1 == is_cmplx){
        for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
          jgrid_real =  jgrid * ist->complex_idx;
          jgrid_imag = jgrid_real + 1;

          // Real part
          psi_rank[jstate + jgrid_real] += \
            (an[ncjms+jc].re * psi[jgrid].re - an[ncjms+jc].im * psi[jgrid].im);
          // Imag part
          psi_rank[jstate + jgrid_imag] += \
            (an[ncjms+jc].re * psi[jgrid].im + an[ncjms+jc].im * psi[jgrid].re);
        }
      }
      else{
        for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
          psi_rank[jstate + jgrid] += an[ncjms+jc].re * psi[jgrid].re;
        }
      }
    }
  }

  return;
}


/*******************************************************************/

void p_hnorm_k(
  zomplex*      psi_out, 
  zomplex*      psi_tmp, 
  double*       pot_local,
  vector*       G_vecs,
  vector        k,
  grid_st*      grid,
  zomplex*      LS,
  nlc_st*       nlc, 
  long*         nl,
  double        zn,
  index_st*     ist, 
  par_st*       par, 
  flag_st*      flag, 
  parallel_st*  parallel,
  fftw_plan_loc planfw,
  fftw_plan_loc planbw,
  fftw_complex* fftwpsi
  ){
  
  long jgrid; 

  memcpy(&psi_tmp[0], &psi_out[0], ist->nspinngrid * sizeof(psi_tmp[0]));

  p_hamiltonian_k(
    psi_out, psi_tmp, pot_local, G_vecs, k, grid, LS, nlc, nl, ist, 
    par, flag, planfw, planbw, fftwpsi, par->ham_threads
  );
  
  for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
    /*** par->dE_1 = 4.0 / par->dE and therefore I don't multiply by 4 ***/
    psi_out[jgrid].re = par->dE_1 * psi_out[jgrid].re - (2.0 + zn + par->Vmin * par->dE_1) * psi_tmp[jgrid].re;
    psi_out[jgrid].im = par->dE_1 * psi_out[jgrid].im - (2.0 + zn + par->Vmin * par->dE_1) * psi_tmp[jgrid].im;
  }

  return;
}

/*****************************************************************************/

void time_hamiltonian_k(
  zomplex*      psi_out, 
  zomplex*      psi_tmp, 
  double*       pot_local, 
  vector*       G_vecs,
  vector        k,
  grid_st*      grid,
  zomplex*      LS,
  nlc_st*       nlc, 
  long*         nl, 
  index_st*     ist, 
  par_st*       par, 
  flag_st*      flag, 
  parallel_st*  parallel){
  /*******************************************************************
  * This function applies the Hamiltonian onto a state               *
  * inputs:                                                          *
  *  [psi_tmp] nspinngrid-long arr of orig. wavefnc                  *
  *  [psi_out] nspinngrid-long arr to hold |psi_out> = H|psi_tmp>    *
  *  [pot_local] ngrid-long arr holding the value of the local pot   *
  *  [nlc] nlc struct holding values for computing SO and NL pots    *
  *  [nl] natom-long arr holding the number of NL gridpts per atom   *
  *  [ksqr] ngrid-long arr holding the values of k^2 for KE calc     *
  *  [ist] ptr to counters, indices, and lengths                     *
  *  [par] ptr to par_st holding VBmin, VBmax... params              *
  *  [flag] ptr to flag_st holding job flags                         *
  *  [planfw] FFTW3 plan for executing 3D forward DFT                *
  *  [planfw] FFTW3 plan for executing 3D backwards DFT              *
  *  [fftwpsi] location to store outcome of Fourier transform        *
  * outputs: void                                                    *
  ********************************************************************/
  
  const int       mpir  = parallel->mpi_rank;

  struct timespec start, end;
  int jspin, j, jtmp;
  int n_iter = 20;

  omp_set_num_threads(par->ham_threads);

  fftw_init_threads();
  fftw_plan_with_nthreads(par->ham_threads);
  fftw_plan_loc planfw, planbw; fftw_complex *fftwpsi; 
  long fft_flags=0;

  fftwpsi = fftw_malloc(sizeof(fftw_complex)*ist->ngrid);
  /*** initialization for the fast Fourier transform ***/
  planfw = fftw_plan_dft_3d(ist->nz, ist->ny, ist->nx, fftwpsi, fftwpsi, FFTW_FORWARD, fft_flags);
  planbw = fftw_plan_dft_3d(ist->nz, ist->ny, ist->nx, fftwpsi, fftwpsi, FFTW_BACKWARD, fft_flags);
  
  // Copy psi_out into psi_tmp
  memcpy(&psi_tmp[0], &psi_out[0], ist->nspinngrid*sizeof(psi_tmp[0]));
  
  
  // Calculate the action of the kinetic energy part of the Hamiltonian on psi_tmp: |psi_out> = T|psi_tmp>
  // Warmup runs to avoid including caching time, optimizations, innitial overhead etc.
  for (j = 0; j < 10; j++){
    for (jspin = 0; jspin < ist->nspin; jspin++){
      kinetic_k(&psi_out[jspin*ist->ngrid], G_vecs, k, planfw, planbw, fftwpsi, ist); //spin up/down
    } 
  }
  clock_gettime(CLOCK_MONOTONIC, &start); 
  for (j = 0; j < n_iter; j++){
    for (jspin = 0; jspin < ist->nspin; jspin++){
      kinetic_k(&psi_out[jspin*ist->ngrid], G_vecs, k, planfw, planbw, fftwpsi, ist); //spin up/down
    } 
  }
  clock_gettime(CLOCK_MONOTONIC, &end);
  double elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
  if (mpir == 0) printf("\tKinetic energy: %.4g (msec)\n", (elapsed_seconds*1000.0)/(double)n_iter ); fflush(0);
  
  
  // Calculate the action of the potential operator on the wavefunction: |psi_out> = V|psi_tmp>
  //
  if (1 == flag->periodic){
    e_ikr(psi_tmp, k, grid, ist, par, flag);
  }
  //
  if(flag->SO==1){
    // Calculate |psi_out> = V_SO|psi_tmp>

    // Warmup runs
    for (j = 0; j < 3; j++){
      p_spin_orbit_proj_pot(psi_out, psi_tmp, LS, nlc, nl, ist, par, par->ham_threads);
    }
    clock_gettime(CLOCK_MONOTONIC, &start); 
    for (j = 0; j < n_iter; j++){
      p_spin_orbit_proj_pot(psi_out, psi_tmp, LS, nlc, nl, ist, par, par->ham_threads);
    }
    clock_gettime(CLOCK_MONOTONIC, &end); 
    double elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    if (mpir == 0) printf("\tSpin-Orbit potential: %.4g (msec)\n", (elapsed_seconds*1000.0)/(double)n_iter);
  }

  
  if (flag->NL == 1){
    // Calculate |psi_out> += V_NL|psi_tmp>
    for (j = 0; j < 3; j++){
      p_nonlocal_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par, par->ham_threads);
    }
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (j = 0; j < n_iter; j++){
      p_nonlocal_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par, par->ham_threads);
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    if (mpir == 0) printf("\tNon-local potential: %.4g (msec)\n", (elapsed_seconds*1000.0)/(double)n_iter );
    
  }
  
  
  // Calculate the action of the local potential energy part of the Hamiltonian on psi_tmp
  clock_gettime(CLOCK_MONOTONIC, &start);
  for (int i = 0; i < n_iter; i++){
    if (1 == flag->useSpinors){
      for (j = 0; j < ist->ngrid; j++) {
        psi_out[j].re += (pot_local[j] * psi_tmp[j].re);
        psi_out[j].im += (pot_local[j] * psi_tmp[j].im);
        // handle spin down component
        jtmp = ist->ngrid + j;
        psi_out[jtmp].re += (pot_local[j] * psi_tmp[jtmp].re);
        psi_out[jtmp].im += (pot_local[j] * psi_tmp[jtmp].im);
      }
    } 
    else if (0 == flag->useSpinors){
      for (j = 0; j < ist->ngrid; j++) {
        psi_out[j].re += (pot_local[j] * psi_tmp[j].re);
      }
    } else {
      printf("ERROR: unrecognized value for flag->useSpinors: %d\n", flag->useSpinors);
      fprintf(stderr, "ERROR: unrecognized value for flag->useSpinors: %d\n", flag->useSpinors);
      exit(EXIT_FAILURE);
    }
  }
  clock_gettime(CLOCK_MONOTONIC, &end);
  elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
  if (mpir == 0) printf("\tLocal potential: %.4g (msec)\n", (elapsed_seconds*1000.0)/(double)n_iter );

  return;
}
