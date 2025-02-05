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
  long jmn, jns, jms, jstate, cntr = 0; // across states
  long jgrid, jgrid_real, jgrid_imag; // within a single state
  long jc, ncjms; // chebyshev coefficients
  int i;
  long MAX_SIZE = ist->complex_idx*ist->nspinngrid*(parallel->mpi_rank*ist->n_states_per_rank  + ist->n_states_per_rank);
  
  // File I/O
  FILE *pf; 
  char fileName[100];
  // energy of the filtered states
  double *ene_filters;
  
  // Arrays for hamiltonian evaluation
  zomplex *psi, *phi;  
  // Allocate memory for arrays
  if ((psi = (zomplex*)calloc(ist->nspinngrid,sizeof(zomplex)))==NULL){ 
    fprintf(stderr, "\nOUT OF MEMORY: psi in run_filter_cycle\n\n"); exit(EXIT_FAILURE);
  }
  if ((phi = (zomplex*)calloc(ist->nspinngrid,sizeof(zomplex)))==NULL){ 
    fprintf(stderr, "\nOUT OF MEMORY: phi in run_filter_cycle\n\n"); exit(EXIT_FAILURE);
  }
  // FFT
  fftw_init_threads();
  fftw_plan_with_nthreads(par->ham_threads);
  long fft_flags = FFTW_MEASURE;
  fftw_plan_loc planfw, planbw;
  fftw_complex *fftwpsi;
  // Create FFT structs and plans for Fourier transform
  fftwpsi = fftw_malloc(sizeof (fftw_complex) * ist->ngrid);
  planfw = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_FORWARD,fft_flags);
  planbw = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_BACKWARD,fft_flags);
  
  for (jns = 0; jns < ist->n_filter_cycles; jns++){
    for (jms = 0; jms < ist->m_states_per_filter; jms++){
      jmn = jns * ist->m_states_per_filter + jms;
      parallel->jns[jmn] = jns;
      parallel->jms[jmn] = jms;
    }
  }
  
  // Check that mn_states is equal to ist->mn_states_tot
  if ( (jmn + 1) != ist->mn_states_tot){
    fprintf(stderr, "ERROR: parallelizing filter jmn %ld NOT equal to ist->mn_states_tot %ld\n", jmn, ist->mn_states_tot);
    exit(EXIT_FAILURE);
  }


  //************************************************************************
  //************************************************************************
  // BEGIN FILTERING HERE
  // In MPI implementation, every MPI task will get ~n_filter_cycles/nranks 
  // initial random states to filter. The filtering of each energy target 
  // using this random state happens serially and OpenMP parallelization is 
  // used to accelerate evaluation of the Hamiltonian operation
  //************************************************************************
  //************************************************************************
  
  int start = parallel->mpi_rank * ist->n_states_per_rank;
  int end = (parallel->mpi_rank == parallel->mpi_size - 1) ? parallel->mpi_size*ist->n_states_per_rank : start + ist->n_states_per_rank;
  

  // Loop over all of the states handled by this mpi-rank 
  for (jmn = start; jmn < end; jmn++){
    // Keep track of how many filter iterations have taken place
    
    if (parallel->mpi_rank == 0){
      printf("\tCurrently working on iteration %ld/%ld of filtering cycle\n", cntr+1, ist->n_states_per_rank); fflush(0);
      printf("\t  (~%ld states)\n\n", (cntr+1) * parallel->mpi_size); fflush(0);
    }
    
    // Set index values for this iteration of the loop    
    jms = parallel->jms[jmn]; 
    jstate = cntr * ist->complex_idx * ist->nspinngrid; 
    ncjms = ist->ncheby * jms; 
    
    if (jstate >= MAX_SIZE || jstate < 0) {
    fprintf(stderr, "jstate index out of bounds %ld!\n", jstate);
    exit(EXIT_FAILURE);
    }
    
    // Populate psi with the random initial wavefunction of filter cycle jns
    if (1 == flag->isComplex){
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
        jgrid_real = ist->complex_idx * jgrid;
        jgrid_imag = ist->complex_idx * jgrid + 1;

        psi[jgrid].re = psi_rank[jstate + jgrid_real];
        psi[jgrid].im = psi_rank[jstate + jgrid_imag];
      }
    } else{
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
        psi[jgrid].re = psi_rank[jstate + jgrid];
      }
    }
    
    // For each energy target, perform a filter procedure to reshape psi
    // into a state with energy close to the target energy
    
    // Calculate term 0 of the expansion
    if (1 == flag->isComplex){
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
        jgrid_real = ist->complex_idx * jgrid;
        jgrid_imag = ist->complex_idx * jgrid + 1;

        psi_rank[jstate + jgrid_real] = an[ncjms+0].re * psi[jgrid].re - an[ncjms+0].im * psi[jgrid].im;
        psi_rank[jstate + jgrid_imag] = an[ncjms+0].re * psi[jgrid].im + an[ncjms+0].im * psi[jgrid].re;
      }
    } else{
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
        psi_rank[jstate + jgrid_real] = an[ncjms+0].re * psi[jgrid].re - an[ncjms+0].im * psi[jgrid].im;
      }
    }
    
    // Calculate the subsequent terms of the expansion
    for (jc = 1; jc < ist->ncheby; jc++){
      
      memcpy(&phi[0], &psi[0], ist->nspinngrid * sizeof(phi[0]));
      
      // omp_set_num_threads(par->ham_threads);
      p_hamiltonian(psi,phi,pot_local,nlc,nl,ksqr,ist,par,flag,planfw,planbw,fftwpsi,par->ham_threads);
      
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
        /*** par->dE_1 = 4.0 / par->dE and therefore I don't multiply by 4 ***/
        psi[jgrid].re = par->dE_1 * psi[jgrid].re - (2.0 + zn[jc-1] + par->Vmin * par->dE_1) * phi[jgrid].re;
        psi[jgrid].im = par->dE_1 * psi[jgrid].im - (2.0 + zn[jc-1] + par->Vmin * par->dE_1) * phi[jgrid].im;
      }
      
      if (1 == flag->isComplex){
        for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
          jgrid_real = ist->complex_idx * jgrid;
          jgrid_imag = ist->complex_idx * jgrid + 1;

          psi_rank[jstate + jgrid_real] += (an[ncjms+jc].re * psi[jgrid].re - an[ncjms+jc].im * psi[jgrid].im);
          psi_rank[jstate + jgrid_imag] += (an[ncjms+jc].re * psi[jgrid].im + an[ncjms+jc].im * psi[jgrid].re);
        }
      } else{
        for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
          psi_rank[jstate + jgrid] += an[ncjms+jc].re * psi[jgrid].re ; 
        }
      }

    }

    cntr++;
  }
  free(psi); free(phi);
  fftw_destroy_plan(planfw);
  fftw_destroy_plan(planbw);
  fftw_free(fftwpsi);

  /*****************************************************************************/
  // Copy the filtered wavefunctions back into psi_rank
  /*****************************************************************************/

  if (1 == flag->printPsiFilt){
    // Print the normalized filtered states to disk
    sprintf(fileName, "psi-filt-%d.dat", parallel->mpi_rank);
    pf = fopen(fileName, "wb");
    for (jmn = 0; jmn < ist->n_states_per_rank; jmn++){
      fwrite(&psi_rank[jmn*ist->complex_idx*ist->nspinngrid], sizeof(double), ist->complex_idx*ist->nspinngrid, pf);
    }
    fclose(pf);
  }

  /***********************************************************************/
  /*** normalize the states and get their energies***/
  if (parallel->mpi_rank == 0) printf("\n  4.3 Normalizing filtered states\n"); fflush(stdout);
  // printf("Normalizing states\n"); fflush(0);
  normalize_all(psi_rank, ist->n_states_per_rank, ist, par, flag, parallel);
  

  // Get the energy of all the filtered states
  if ((ene_filters = (double*)calloc(ist->n_states_per_rank,sizeof(double)))==NULL)nerror("ene_filters");
  
  if (parallel->mpi_rank == 0) printf("\n  4.4 Computing the energies of all filtered states\n"); fflush(stdout);
  /*** calculate and print the energy of the filtered states ***/
  energy_all(psi_rank,ist->n_states_per_rank,pot_local,nlc,nl,ksqr,ene_filters,ist,par,flag,parallel);
  // printf("exited energy_all\n"); fflush(0);
  sprintf(fileName,"ene-filt-rank-%d.dat", parallel->mpi_rank);
  pf = fopen(fileName , "w");
  cntr = 0;
  for (jmn = start; jmn < end; jmn++){
    jns = parallel->jns[jmn];
    jms = parallel->jms[jmn];
    
    fprintf (pf,"%ld %ld %.16g %.16g\n", jns, jms, ene_filters[cntr], ene_targets[jms]);
    cntr++;
  }
  fclose(pf);
  
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
      p_spin_orbit_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par, par->ham_threads);
    }
    clock_gettime(CLOCK_MONOTONIC, &start); 
    for (j = 0; j < 10; j++){
      p_spin_orbit_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par, par->ham_threads);
    }
    clock_gettime(CLOCK_MONOTONIC, &end); 
    double elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    if (parallel->mpi_rank == 0) printf("\tSpin-Orbit potential: %.4g (msec)\n", (elapsed_seconds*1000.0)/10.0 );
  }

  if (flag->NL == 1){
    // Calculate |psi_out> += V_NL|psi_tmp>
    for (j = 0; j < 3; j++){
      p_nonlocal_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par, par->ham_threads);
    }
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (j = 0; j < 10; j++){
      p_nonlocal_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par, par->ham_threads);
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    if (parallel->mpi_rank == 0) printf("\tNon-local potential: %.4g (msec)\n", (elapsed_seconds*1000.0)/10.0 );
    
  }
  
  // Calculate the action of the local potential energy part of the Hamiltonian on psi_tmp
  clock_gettime(CLOCK_MONOTONIC, &start);
  for (int i = 0; i < 10; i++){
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
  if (parallel->mpi_rank == 0) printf("\tLocal potential: %.4g (msec)\n", (elapsed_seconds*1000.0)/10.0 );

  
  return;
}

