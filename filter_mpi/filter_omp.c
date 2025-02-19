#include "filter.h"


/*****************************************************************************/

void run_filter_cycle(
  double*       psi_rank, 
  double*       pot_local,
  zomplex*      LS,
  nlc_st*       nlc, 
  long*         nl, 
  double*       ksqr, 
  zomplex*      an, 
  double*       zn, 
  double*       ene_targets, 
  grid_st*      grid, 
  index_st*     ist, 
  par_st*       par, 
  flag_st*      flag, 
  parallel_st*  parallel){
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
  long MAX_SIZE = ist->complex_idx*ist->nspinngrid*((ist->n_filters_per_rank - 1)*ist->m_states_per_filter + (ist->m_states_per_filter - 1)) + ist->nspinngrid*ist->complex_idx - 1;
  
  // File I/O
  FILE *pf; char fileName[100];
  // energy of the filtered states
  double *ene_filters;
  
  ALLOCATE(&(parallel->jns), ist->mn_states_tot, "parallel->jns");
  ALLOCATE(&(parallel->jms), ist->mn_states_tot, "parallel->jms");
  
  for (jns = 0; jns < ist->n_filters_per_rank; jns++){
    for (jms = 0; jms < ist->m_states_per_filter; jms++){
      jmn = jns * ist->m_states_per_filter + jms;
      parallel->jns[jmn] = jns;
      parallel->jms[jmn] = jms;
    }
  }
  
  // Check that mn_states is equal to ist->mn_states_tot
  if ( (jmn + 1) != ist->n_states_per_rank){
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
  
  // Loop over all of the random states handled by this mpi-rank
  
  omp_set_max_active_levels(1);
  omp_set_num_threads(parallel->nthreads);
  #pragma omp parallel for private(jmn) 
  for (jmn = 0; jmn < ist->n_states_per_rank; jmn++) {
     // Keep track of how many filter iterations have taken place
    if (0 == omp_get_thread_num()){
      cntr++;
      if (parallel->mpi_rank == 0) printf("\tCurrently working on iteration %ld of filtering cycle\n", cntr); fflush(0);
      if (parallel->mpi_rank == 0) printf("\t  (~%ld states)\n\n", cntr * omp_get_num_threads()); fflush(0);
      
    }
    // All variables declared within this parallel region are private for each thread
    // thread tracking

    // File I/O
    FILE *pf; char str[100];
    // Arrays for hamiltonian evaluation
    zomplex *psi, *phi;  
    // Array indexing
    long jns, jms, jns_ms, ns_block; // across states
    long jgrid, jgrid_real, jgrid_imag; // within a single state
    long jc, ncjms; // chebyshev coefficients
    // FFT
    long fft_flags = FFTW_MEASURE;
    fftw_plan_loc planfw, planbw;
    fftw_complex *fftwpsi;
    // visualization
    // double *rho, sgn_val;
    // if ((rho = (double *) calloc(ist->ngrid, sizeof(rho[0]))) == NULL){
    //   fprintf(stderr, "\nOUT OF MEMORY: filter rho\n\n"); exit(EXIT_FAILURE);
    // }

    // Allocate memory for arrays
    if ((psi = (zomplex*)calloc(ist->nspinngrid,sizeof(zomplex)))==NULL){ 
      fprintf(stderr, "\nOUT OF MEMORY: psi in run_filter_cycle\n\n"); exit(EXIT_FAILURE);
    }
    if ((phi = (zomplex*)calloc(ist->nspinngrid,sizeof(zomplex)))==NULL){ 
      fprintf(stderr, "\nOUT OF MEMORY: phi in run_filter_cycle\n\n"); exit(EXIT_FAILURE);
    }
  
    // Create FFT structs and plans for Fourier transform
    fftwpsi = fftw_malloc(sizeof (fftw_complex) * ist->ngrid);
    planfw = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_FORWARD,fft_flags);
    planbw = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_BACKWARD,fft_flags);
  
    // Set index values for this iteration of the loop
    jns = parallel->jns[jmn]; 
    jms = parallel->jms[jmn]; 
    
    ns_block = (jns * ist->m_states_per_filter * ist->complex_idx * ist->nspinngrid); 
    jns_ms = ns_block + (jms * ist->complex_idx * ist->nspinngrid); 
    ncjms = ist->ncheby * jms; 
    
    if (jns_ms >= MAX_SIZE || jns_ms < 0) {
    fprintf(stderr, "jns_ms index out of bounds %ld!\n", jns_ms);
    exit(EXIT_FAILURE);
    }
    
    // Populate psi with the random initial wavefunction of filter cycle jns
    if (1 == flag->isComplex){
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
        jgrid_real = ist->complex_idx * jgrid;
        jgrid_imag = ist->complex_idx * jgrid + 1;

        psi[jgrid].re = psi_rank[ns_block + jgrid_real];
        psi[jgrid].im = psi_rank[ns_block + jgrid_imag];
      }
    } else{
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
        psi[jgrid].re = psi_rank[ns_block + jgrid];
      }
    }
    
    // For each energy target, perform a filter procedure to reshape psi
    // into a state with energy close to the target energy
    
    // Calculate term 0 of the expansion
    if (1 == flag->isComplex){
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
	      jgrid_real = ist->complex_idx * jgrid;
        jgrid_imag = ist->complex_idx * jgrid + 1;

        psi_rank[jns_ms + jgrid_real] = an[ncjms+0].re * psi[jgrid].re - an[ncjms+0].im * psi[jgrid].im;
        psi_rank[jns_ms + jgrid_imag] = an[ncjms+0].re * psi[jgrid].im + an[ncjms+0].im * psi[jgrid].re;
      }
    } 
    else {
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
        psi_rank[jns_ms + jgrid] = an[ncjms+0].re * psi[jgrid].re - an[ncjms+0].im * psi[jgrid].im;
      }
    }
    
    // Calculate the subsequent terms of the expansion

    for (jc = 1; jc < ist->ncheby; jc++){
      
      memcpy(&phi[0], &psi[0], ist->nspinngrid * sizeof(phi[0]));
      hamiltonian(psi,phi,pot_local,LS,nlc,nl,ksqr,ist,par,flag,planfw,planbw,fftwpsi);
      // printf("Finished Hamiltonian\n"); fflush(0);
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
        /*** par->dE_1 = 4.0 / par->dE and therefore I don't multiply by 4 ***/
        psi[jgrid].re = par->dE_1 * psi[jgrid].re - (2.0 + zn[jc-1] + par->Vmin * par->dE_1) * phi[jgrid].re;
        psi[jgrid].im = par->dE_1 * psi[jgrid].im - (2.0 + zn[jc-1] + par->Vmin * par->dE_1) * phi[jgrid].im;
      }

      if (1 == flag->isComplex){
        for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
          jgrid_real = ist->complex_idx * jgrid;
          jgrid_imag = ist->complex_idx * jgrid + 1;

          // #pragma omp atomic
          psi_rank[jns_ms + jgrid_real] += (an[ncjms+jc].re * psi[jgrid].re - an[ncjms+jc].im * psi[jgrid].im);
          // #pragma omp atomic
          psi_rank[jns_ms + jgrid_imag] += (an[ncjms+jc].re * psi[jgrid].im + an[ncjms+jc].im * psi[jgrid].re);
        } 
      }
      else{
        for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
          // #pragma omp atomic
          psi_rank[jns_ms + jgrid] += (an[ncjms+jc].re * psi[jgrid].re - an[ncjms+jc].im * psi[jgrid].im);
        }
      }

      if ((0 == omp_get_thread_num()) && (parallel->mpi_rank == 0)){
        if ( (0 == (jc % ((long) (ist->ncheby / 4)) )) || (1 == jc) || ( (ist->ncheby - 1) == jc) ) {
          // print the filtering progress to the output file
          print_progress_bar(jc, ist->ncheby);
        }
      }
    }
    
    // fclose(pf);
    
    /*****************************************************************************/
    // Copy the filtered wavefunctions back into psi_rank
    /*****************************************************************************/

    // printf("Copying data\n"); fflush(0);
    // for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
    //   jgrid_real = ist->complex_idx * jgrid;
    //   jgrid_imag = ist->complex_idx * jgrid + 1;
    //   psi_rank[jns_ms + jgrid_real] = psi_out[jgrid].re;
    //   if (1 == flag->isComplex){
    //     psi_rank[jns_ms + jgrid_imag] = psi_out[jgrid].im;
    //   }
    // }
    
    if (1 == flag->printPsiFilt){
      // Print the normalized filtered states to disk
      sprintf(str, "psi-filt-%ld-%d.dat", jmn, parallel->mpi_rank);
      pf = fopen(str, "w");
      fwrite(&psi_rank[jmn*ist->complex_idx*ist->nspinngrid], sizeof(double), ist->complex_idx*ist->nspinngrid, pf);
    }
    
    free(psi); free(phi);
    fftw_destroy_plan(planfw);
    fftw_destroy_plan(planbw);
    fftw_free(fftwpsi);
  }
  
  /***********************************************************************/
  /*** normalize the states and get their energies***/
  if (parallel->mpi_rank == 0) printf("\n  4.3 Normalizing filtered states\n"); fflush(stdout);
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
  
  if (parallel->mpi_rank == 0) printf("\n  4.4 Computing the energies of all filtered states\n"); fflush(stdout);
  /*** calculate and print the energy of the filtered states ***/
  energy_all(psi_rank,ist->n_states_per_rank,pot_local,LS,nlc,nl,ksqr,ene_filters,ist,par,flag,parallel);

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

void time_hamiltonian(
  zomplex*      psi_out, 
  zomplex*      psi_tmp, 
  double*       pot_local, 
  zomplex*      LS,
  nlc_st*       nlc, 
  long*         nl, 
  double*       ksqr,
  index_st*     ist, 
  par_st*       par, 
  flag_st*      flag, 
  parallel_st*  parallel
  ){
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
      spin_orbit_proj_pot(psi_out, psi_tmp, LS, nlc, nl, ist, par);
    }
    clock_gettime(CLOCK_MONOTONIC, &start); 
    for (j = 0; j < 10; j++){
      spin_orbit_proj_pot(psi_out, psi_tmp, LS, nlc, nl, ist, par);
    }
    clock_gettime(CLOCK_MONOTONIC, &end); 
    double elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    if (parallel->mpi_rank == 0) printf("\tSpin-Orbit potential: %.4g (msec)\n", (elapsed_seconds*1000.0)/10.0 );
  }

  if (flag->NL == 1){
    // Calculate |psi_out> += V_NL|psi_tmp>
    for (j = 0; j < 3; j++){
      nonlocal_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par);
    }
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (j = 0; j < 10; j++){
      nonlocal_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par);
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

