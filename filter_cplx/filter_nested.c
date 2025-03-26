#include "fd.h"

/*****************************************************************************/

void run_filter_cycle(double *psi_rank, double *pot_local, nlc_st *nlc, long *nl, 
  double *ksqr, zomplex *an, double *zn, double *ene_targets, grid_st *grid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){
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
  // Array indexing
  long jns, jms, jmn, cntr = 0;
  int i, n_outer_threads = parallel->nthreads, n_inner_threads = 1;
  // File I/O
  FILE *pf; char fileName[100];
  // energy of the filtered states
  double *ene_filters;
  FILE *ptime; char str[100];
  ptime = fopen("time.dat", "w");
  
  
  //************************************************************************
  //************************************************************************
  // BEGIN FILTERING HERE
  // In MPI implementation, every MPI task will get ~n_filter_cycles/nranks 
  // initial random states to filter. The filtering of each energy target 
  // using this random state happens serially and OpenMP parallelization is 
  // used to accelerate evaluation of the Hamiltonian operation
  //************************************************************************
  //************************************************************************
  if (parallel->mpi_rank == 0) printf("Starting filtering loop\n\n"); fflush(0);
  
  if (1 == parallel->nestedOMP){
    omp_set_max_active_levels(2);
    if (parallel->mpi_rank == 0) printf("Nested OMP parallel region enabled\n");
    
    if (0.0 != (parallel->nthreads % ist->m_states_per_filter) ){
      if (parallel->mpi_rank == 0) fprintf(stderr, "Error: number of total threads not a factor of 4!\n");
      exit(EXIT_FAILURE);
    }
    if (parallel->mpi_rank == 0) printf("  n_outer_threads = %d\n", parallel->n_outer_threads);
    if (parallel->mpi_rank == 0) printf("  n_inner_threads = %d\n", parallel->n_inner_threads);

    n_outer_threads = parallel->n_outer_threads;
    n_inner_threads = parallel->n_inner_threads;
  }

  // Allocate thread-local memory for the outer loop
  zomplex *psi, *phi;
  if ((psi = (zomplex*)calloc(n_outer_threads * ist->nspinngrid,sizeof(zomplex)))==NULL){ 
    fprintf(stderr, "\nOUT OF MEMORY: psi in run_filter_cycle\n\n"); exit(EXIT_FAILURE);
  }
  if ((phi = (zomplex*)calloc(n_outer_threads * ist->nspinngrid,sizeof(zomplex)))==NULL){ 
    fprintf(stderr, "\nOUT OF MEMORY: phi in run_filter_cycle\n\n"); exit(EXIT_FAILURE);
  }
  
  // // FFT
  fftw_init_threads();
  // fftw_init_threads();
  // fftw_plan_with_nthreads(parallel->n_inner_threads);
  // fftw_complex *fftwpsi;
  // long fft_flags = FFTW_MEASURE;
  // fftw_plan_loc *planfw, *planbw;
  // // Create FFT structs and plans for Fourier transform
  // fftwpsi = fftw_malloc(n_outer_threads * ist->ngrid * sizeof (fftw_complex));
  // planfw = (fftw_plan_loc*) calloc(n_outer_threads, sizeof(planfw[0]));
  // planbw = (fftw_plan_loc*) calloc(n_outer_threads, sizeof(planfw[0]));

  // for (i = 0; i < n_outer_threads; i++){
  //   planfw[i] = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_FORWARD,fft_flags);
  //   planbw[i] = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_BACKWARD,fft_flags);
  // }
  long t_nsng;
  // *
  // *
  // Loop over all of the random states handled by this mpi-rank
  // *
  // *
  for (jns = 0; jns < ist->n_filters_per_rank; jns++) {
     // Keep track of how many filter iterations have taken place
    if (parallel->mpi_rank == 0) {
      cntr++;
      time_t current_time = time(NULL);
      char* c_time_string = ctime(&current_time);
      if (parallel->mpi_rank == 0) printf("  Random psi %ld of %ld on node %d | %s\n", cntr, ist->n_filters_per_rank, parallel->mpi_rank, c_time_string); fflush(0);
    }

     
    #pragma omp parallel for num_threads(n_outer_threads) private(jns, jms, t_nsng)
    for (jms = 0; jms < ist->m_states_per_filter; jms++){
      if (parallel->mpi_rank == 0) {
        if (parallel->mpi_rank == 0) printf("\tEnergy target %ld centered at %.3lg a.u.\n", jms, ene_targets[jms]); fflush(0);
      }
      // Arrays for hamiltonian evaluation
      
      struct timespec start, end;
      clock_gettime(CLOCK_MONOTONIC, &start);

      int t_id = omp_get_thread_num();
      // Array indexing
      long jstate, ns_block; // across states
      long jgrid, jgrid_real, jgrid_imag; // within a single state
      long jc, ncjms; // chebyshev coefficients
      t_nsng = t_id * ist->nspinngrid;

      // Set index values for this iteration of the loop
      ns_block = (jns * ist->m_states_per_filter * ist->complex_idx * ist->nspinngrid); 
      jstate = ns_block + jms*ist->nspinngrid*ist->complex_idx;
      ncjms = ist->ncheby * jms;
      
      // FFT
      
      fftw_plan_with_nthreads(parallel->n_inner_threads);
      fftw_complex *fftwpsi;
      long fft_flags = FFTW_MEASURE;
      fftw_plan_loc planfw, planbw;
      // Create FFT structs and plans for Fourier transform
      fftwpsi = fftw_malloc(ist->ngrid * sizeof (fftw_complex));
      planfw = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_FORWARD,fft_flags);
      planbw = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_BACKWARD,fft_flags);
      
      
      // Populate psi with the random initial wavefunction corresponding to filter cycle jns
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
        jgrid_real = ist->complex_idx * jgrid;
        jgrid_imag = ist->complex_idx * jgrid + 1;

        psi[t_nsng + jgrid].re = psi_rank[ns_block + jgrid_real];
        if (1 == flag->isComplex){
          psi[t_nsng + jgrid].im = psi_rank[ns_block + jgrid_imag];
        }
      }
      
      // For each energy target, perform a filter procedure to reshape psi
      // into a state with energy close to the target energy
      
      // Calculate term 0 of the expansion
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
        jgrid_real = ist->complex_idx * jgrid;
        jgrid_imag = ist->complex_idx * jgrid + 1;

        psi_rank[jstate + jgrid_real] = an[ncjms+0].re * psi[t_nsng + jgrid].re - an[ncjms+0].im * psi[t_nsng + jgrid].im;
        if (1 == flag->isComplex){
          psi_rank[jstate + jgrid_imag] = an[ncjms+0].re * psi[t_nsng + jgrid].im + an[ncjms+0].im * psi[t_nsng + jgrid].re;
        }
      }
      
      // Calculate the subsequent terms of the expansion
      for (jc = 1; jc < ist->ncheby; jc++){

        memcpy(&phi[t_nsng], &psi[t_nsng], ist->nspinngrid * sizeof(phi[0]));
        // write_state_dat(phi, ist->nspinngrid, "phi_filter.dat");
        p_hamiltonian(&psi[t_nsng],&phi[t_nsng],pot_local,nlc,nl,ksqr,ist,par,flag,planfw,planbw,fftwpsi,n_inner_threads);
        //exit(0);
        for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
          /*** par->dE_1 = 4.0 / par->dE and therefore I don't multiply by 4 ***/
          psi[t_nsng + jgrid].re = par->dE_1 * psi[t_nsng + jgrid].re - (2.0 + zn[jc-1] + par->Vmin * par->dE_1) * phi[t_nsng + jgrid].re;
          psi[t_nsng + jgrid].im = par->dE_1 * psi[t_nsng + jgrid].im - (2.0 + zn[jc-1] + par->Vmin * par->dE_1) * phi[t_nsng + jgrid].im;
        }

        for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
          jgrid_real = ist->complex_idx * jgrid;
          jgrid_imag = ist->complex_idx * jgrid + 1;

          psi_rank[jstate + jgrid_real] += (an[ncjms+jc].re * psi[t_nsng + jgrid].re - an[ncjms+jc].im * psi[t_nsng + jgrid].im);
          if (1 == flag->isComplex){
            psi_rank[jstate + jgrid_imag] += (an[ncjms+jc].re * psi[t_nsng + jgrid].im + an[ncjms+jc].im * psi[t_nsng + jgrid].re);
          }
        }
      }
    
      // Print out the filtered state for debugging
      // if (1 == flag->printPsiFilt){
        // for (jgrid = 0; jgrid < ist->ngrid; jgrid++){
        //   jgrid_real = ist->complex_idx * jgrid;
        //   jgrid_imag = ist->complex_idx * jgrid + 1;

        //   rho[jgrid] = sqr(psi_rank[jstate + jgrid_real]);
        //   sgn_val = psi_rank[jstate + jgrid_real];
        //   if (1 == flag->isComplex){
        //     rho[jgrid] += sqr(psi_rank[jstate + jgrid_imag]);
        //     if (sgn_val < psi_rank[jstate + jgrid_imag]) sgn_val = psi_rank[jstate + jgrid_imag];
        //   }
        //   rho[jgrid] *= sign(sgn_val);
        // }
        // sprintf(str, "psi-filt-%ld-%ld-%d.cube", jns, jms, parallel->mpi_rank);
        // write_cube_file(rho, grid, str);
      //   sprintf(str, "psi-filt-%ld-%ld-%d.dat", jns, jms, parallel->mpi_rank);
      //   pf = fopen(str, "w");
      //   fwrite(&psi_rank[jstate], sizeof(psi_rank[0]), ist->nspinngrid*ist->complex_idx, pf);
      //   fclose(pf);
      // }

      // for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
      //   jgrid_real = ist->complex_idx * jgrid;
      //   jgrid_imag = ist->complex_idx * jgrid + 1;
      //   psi_rank[jstate + jgrid_real] = psi_loc[t_id * ist->nspinngrid + jgrid].re;
      //   if (1 == flag->isComplex){
      //     psi_rank[jstate + jgrid_imag] = psi_loc[t_id * ist->nspinngrid + jgrid].im;
      //   }
      // }
      fftw_destroy_plan(planfw);
      fftw_destroy_plan(planbw);
      fftw_free(fftwpsi);

      clock_gettime(CLOCK_MONOTONIC, &end);
      double elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
      if (parallel->mpi_rank == 0) fprintf(ptime, "\tOuter loop %ld time ham_threads (): %.4g (msec)\n", jms, (elapsed_seconds*1000.0)/10 );
      
    }
    fclose(ptime);
  }

  free(psi); free(phi);
  // for (i = 0; i < n_outer_threads; i++){
  //   fftw_destroy_plan(planfw[i]);
  //   fftw_destroy_plan(planbw[i]);
  // }
  // fftw_free(fftwpsi);
  
  
  /***********************************************************************/
  /*** normalize the states and get their energies***/
  if (parallel->mpi_rank == 0) printf("Normalizing filtered states\n"); fflush(stdout);
  normalize_all(psi_rank,ist->n_states_per_rank,ist,par,flag,parallel);
  
  // double *rho, sgn_val;
  // long jgrid, jgrid_imag, jgrid_real, jstate;
  // char str[100];
  // if ((rho = (double *) calloc(ist->ngrid, sizeof(rho[0]))) == NULL){
  //   fprintf(stderr, "\nOUT OF MEMORY: filter rho\n\n"); exit(EXIT_FAILURE);
  // }

  // for (jns = 0; jns < ist->n_filters_per_rank; jns++){
  //   for (jms = 0; jms < ist->m_states_per_filter; jms++){
  //     for (jgrid = 0; jgrid < ist->ngrid; jgrid++){
  //       jstate = ist->complex_idx * (jns*ist->nspinngrid*ist->m_states_per_filter + jms*ist->nspinngrid);
  //       jgrid_real = ist->complex_idx * jgrid;
  //       jgrid_imag = ist->complex_idx * jgrid + 1;

  //       rho[jgrid] = sqr(psi_rank[jstate + jgrid_real]);
  //       sgn_val = psi_rank[jstate + jgrid_real];
  //       if (1 == flag->isComplex){
  //         rho[jgrid] += sqr(psi_rank[jstate + jgrid_imag]);
  //         if (sgn_val < psi_rank[jstate + jgrid_imag]) sgn_val = psi_rank[jstate + jgrid_imag];
  //       }
  //       rho[jgrid] *= sign(sgn_val);
  //     }
  //     sprintf(str, "psi-filt-%ld-%ld-%d.cube", jns, jms, parallel->mpi_rank);
  //     write_cube_file(rho, grid, str);
  //   }
  // }

  // Get the energy of all the filtered states
  if ((ene_filters = (double*)calloc(ist->n_filters_per_rank*ist->m_states_per_filter,sizeof(double)))==NULL)nerror("ene_filters");
  
  if (parallel->mpi_rank == 0) printf("Computing the energies of all filtered states\n"); fflush(stdout);
  /*** calculate and print the energy of the filtered states ***/
  energy_all(psi_rank,ist->n_states_per_rank,pot_local,nlc,nl,ksqr,ene_filters,ist,par,flag,parallel);

  for (jns = 0; jns < ist->n_filters_per_rank; jns++){
    sprintf (fileName,"ene-filt-jns-%ld-%d.dat", jns, parallel->mpi_rank);
    pf = fopen(fileName , "w");
    for (jms = 0; jms < ist->m_states_per_filter; jms++){
      jmn = jns*ist->m_states_per_filter + jms;
      fprintf (pf,"%ld %.16g %.16g\n", jms, ene_filters[jmn], ene_targets[jms]);
    }
    fclose(pf);
  }

  free(ene_filters); 
  
  return;
}


/*****************************************************************************/
int sign(float x) {
    return (int)copysign(1.0, x);  // copysign gives the sign of x
}


/*****************************************************************************/

void time_hamiltonian(zomplex *psi_out, zomplex *psi_tmp, double *pot_local, nlc_st *nlc, long *nl, double *ksqr,
  index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){
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
  struct timespec start, end;
  int jspin, j, jtmp; 

  fftw_init_threads();
  fftw_plan_with_nthreads(parallel->n_inner_threads);
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
      kinetic(&psi_out[jspin*ist->ngrid], ksqr, planfw, planbw, fftwpsi, ist); //spin up/down
    } 
  }
  clock_gettime(CLOCK_MONOTONIC, &start); 
  for (j = 0; j < 10; j++){
    for (jspin = 0; jspin < ist->nspin; jspin++){
      kinetic(&psi_out[jspin*ist->ngrid], ksqr, planfw, planbw, fftwpsi, ist); //spin up/down
    } 
  }
  clock_gettime(CLOCK_MONOTONIC, &end);
  double elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
  if (parallel->mpi_rank == 0) printf("\tKinetic energy: %.4g (msec)\n", (elapsed_seconds*1000.0)/10 );
  
  
  // Calculate the action of the potential operator on the wavefunction: |psi_out> = V|psi_tmp>
  
  if(flag->SO==1){
    // Calculate |psi_out> = V_SO|psi_tmp>
    // Warmup runs
    for (j = 0; j < 3; j++){
      p_spin_orbit_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par, parallel->n_inner_threads);
    }
    clock_gettime(CLOCK_MONOTONIC, &start); 
    for (j = 0; j < 10; j++){
      p_spin_orbit_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par, parallel->n_inner_threads);
    }
    clock_gettime(CLOCK_MONOTONIC, &end); 
    double elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    if (parallel->mpi_rank == 0) printf("\tSpin-Orbit potential: %.4g (msec)\n", (elapsed_seconds*1000.0)/10.0 );
  }

  if (flag->NL == 1){
    // Calculate |psi_out> += V_NL|psi_tmp>
    for (j = 0; j < 3; j++){
      p_nonlocal_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par, parallel->n_inner_threads);
    }
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (j = 0; j < 10; j++){
      p_nonlocal_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par, parallel->n_inner_threads);
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    if (parallel->mpi_rank == 0) printf("\tNon-local potential: %.4g (msec)\n", (elapsed_seconds*1000.0)/10.0 );
    
  }
  
  // Calculate the action of the local potential energy part of the Hamiltonian on psi_tmp
  clock_gettime(CLOCK_MONOTONIC, &start);
  for (int i = 0; i < 10; i++){
    for (jspin = 0; jspin < ist->nspin; jspin++){
      for (j = 0; j < ist->ngrid; j++) {
        jtmp = ist->ngrid * jspin + j ; // generalized indexing to handle spinors or spinless wavefuncs

        psi_out[jtmp].re += (pot_local[j] * psi_tmp[jtmp].re);
        psi_out[jtmp].im += (pot_local[j] * psi_tmp[jtmp].im);
      }
    }
  }
  clock_gettime(CLOCK_MONOTONIC, &end);
  elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
  if (parallel->mpi_rank == 0) printf("\tLocal potential: %.4g (msec)\n", (elapsed_seconds*1000.0)/10.0 );

  
  return;
}