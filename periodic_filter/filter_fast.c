#include "fd.h"

/*****************************************************************************/

void run_filter_cycle_k(double *psi_rank, double *pot_local, vector *G_vecs, vector *k_vecs, nlc_st *nlc, long *nl, 
  zomplex *an, double *zn, double *ene_targets, grid_st *grid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){
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
  int i, ik;
  vector k;
  // File I/O
  FILE *pf; char fileName[100];
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
  // Array indexing
  long jstate, ns_block; // across states
  long jgrid, jgrid_real, jgrid_imag; // within a single state
  long jc, ncjms; // chebyshev coefficients
  
  // FFT
  fftw_init_threads();
  fftw_plan_with_nthreads(par->ham_threads);
  fftw_complex *fftwpsi;
  long fft_flags = FFTW_MEASURE;
  fftw_plan_loc planfw, planbw;
  // Create FFT structs and plans for Fourier transform
  fftwpsi = fftw_malloc(sizeof (fftw_complex) * ist->ngrid);
  planfw = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_FORWARD,fft_flags);
  planbw = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_BACKWARD,fft_flags);
  
  //************************************************************************
  //************************************************************************
  // BEGIN FILTERING HERE
  // In MPI implementation, every MPI task will get ~n_filter_cycles/nranks 
  // initial random states to filter. The filtering of each energy target 
  // using this random state happens serially and OpenMP parallelization is 
  // used to accelerate evaluation of the Hamiltonian operation
  //************************************************************************
  //************************************************************************
  
  for (ik = 0; ik < ist->n_k_pts; ik++){
    k = k_vecs[ik];
  
    // *
    // *
    // Loop over all of the random states handled by this mpi-rank
    // *
    // *

    double *rho = (double*)calloc(ist->ngrid, sizeof(double)) ;
    double sgn_val;
    char str[40];
    
    for (jns = 0; jns < ist->n_filters_per_rank; jns++) {
      // Keep track of how many filter iterations have taken place
      if (parallel->mpi_rank == 0) {
        cntr++;
        time_t current_time = time(NULL);
        char* c_time_string = ctime(&current_time);
        if (parallel->mpi_rank == 0) printf("  Random psi %ld of %ld on rank %d | %s\n", cntr, ist->n_filters_per_rank, parallel->mpi_rank, c_time_string); fflush(0);
      }

      for (jms = 0; jms < ist->m_states_per_filter; jms++){
        if (parallel->mpi_rank == 0) {
          if (parallel->mpi_rank == 0) printf("\tEnergy target %ld centered at %.3lg a.u.\n", jms, ene_targets[jms]); fflush(0);
        }
        
        // Set index values for this iteration of the loop
        ns_block = (jns * ist->m_states_per_filter * ist->complex_idx * ist->nspinngrid); 
        jstate = ns_block + jms*ist->nspinngrid*ist->complex_idx;
        ncjms = ist->ncheby * jms;
        
        // Populate psi with the random initial wavefunction corresponding to filter cycle jns
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

            psi_rank[jstate + jgrid_real] = an[ncjms+0].re * psi[jgrid].re - an[ncjms+0].im * psi[jgrid].im;
            psi_rank[jstate + jgrid_imag] = an[ncjms+0].re * psi[jgrid].im + an[ncjms+0].im * psi[jgrid].re;
          }
        } else{
          for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
            psi_rank[jstate + jgrid] = an[ncjms+0].re * psi[jgrid].re - an[ncjms+0].im * psi[jgrid].im;
          }
        }
        
        // Calculate the subsequent terms of the expansion
        for (jc = 1; jc < ist->ncheby; jc++){

          memcpy(&phi[0], &psi[0], ist->nspinngrid * sizeof(phi[0]));
          p_hamiltonian_k(psi,phi,pot_local,G_vecs,k,grid,nlc,nl,ist,par,flag,planfw,planbw,fftwpsi,par->ham_threads);
          
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
          }
          else{
            for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
              psi_rank[jstate + jgrid] += an[ncjms+jc].re * psi[jgrid].re;
            }
          }
        }
      
        // Print out the filtered state
        if (1 == flag->printPsiFilt){
          sprintf(str, "psi-filt-%ld-%ld-%d.dat", jns, jms, parallel->mpi_rank);
          pf = fopen(str, "w");
          fwrite(&psi_rank[jstate], sizeof(double), ist->complex_idx*ist->nspinngrid, pf);
          fclose(pf);
        }
        
      }
    }
    free(psi); free(phi);
    fftw_destroy_plan(planfw);
    fftw_destroy_plan(planbw);
    fftw_free(fftwpsi);
    
    /***********************************************************************/
    /*** normalize the states and get their energies***/
    if (parallel->mpi_rank == 0) printf("\n  4.3 Normalizing filtered states\n"); fflush(stdout);
    normalize_all(psi_rank,ist->n_states_per_rank,ist,par,flag,parallel);

    // Get the energy of all the filtered states
    if ((ene_filters = (double*)calloc(ist->n_states_per_rank, sizeof(double)))==NULL)nerror("ene_filters");
    
    if (parallel->mpi_rank == 0) printf("\n  4.4 Computing the energies of all filtered states\n"); fflush(stdout);
    /*** calculate and print the energy of the filtered states ***/
    energy_all_k(psi_rank,ist->n_states_per_rank,pot_local,G_vecs,k,grid,nlc,nl,ene_filters,ist,par,flag,parallel);


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
  }

  return;
}


/*****************************************************************************/
int sign(float x) {
    return (int)copysign(1.0, x);  // copysign gives the sign of x
}


/*****************************************************************************/

void time_hamiltonian_k(zomplex *psi_out, zomplex *psi_tmp, double *pot_local, vector *G_vecs, vector k, grid_st *grid, nlc_st *nlc, long *nl,
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
  

  int **t_jatom;
  int n_atm_p_thr = (int)(ist->natoms / parallel->nthreads) + 1;
  if ((t_jatom = (int **)calloc(parallel->nthreads, sizeof(t_jatom[0]))) == NULL){
    if (parallel->mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: t_jatom\n\n"); exit(EXIT_FAILURE);
  }
  for (j = 0; j < parallel->nthreads; j++){
    t_jatom[j] = (int*) calloc(n_atm_p_thr, sizeof(int));
  }
  
  // Calculate the action of the kinetic energy part of the Hamiltonian on psi_tmp: |psi_out> = T|psi_tmp>
  // Warmup runs to avoid including caching time, optimizations, innitial overhead etc.
  for (j = 0; j < 10; j++){
    for (jspin = 0; jspin < ist->nspin; jspin++){
      kinetic_k(&psi_out[jspin*ist->ngrid], G_vecs, k, planfw, planbw, fftwpsi, ist); //spin up/down
    } 
  }
  clock_gettime(CLOCK_MONOTONIC, &start); 
  for (j = 0; j < 10; j++){
    for (jspin = 0; jspin < ist->nspin; jspin++){
      kinetic_k(&psi_out[jspin*ist->ngrid], G_vecs, k, planfw, planbw, fftwpsi, ist); //spin up/down
    } 
  }
  clock_gettime(CLOCK_MONOTONIC, &end);
  double elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
  if (parallel->mpi_rank == 0) printf("\tKinetic energy: %.4g (msec)\n", (elapsed_seconds*1000.0)/10 ); fflush(0);
  
  
  // Calculate the action of the potential operator on the wavefunction: |psi_out> = V|psi_tmp>
  if (1 == flag->periodic){
    e_ikr(psi_tmp, k, grid, ist, par, flag);
  }

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
    if (parallel->mpi_rank == 0) printf("\tSpin-Orbit potential: %.4g (msec)\n", (elapsed_seconds*1000.0)/10.0 ); fflush(0);
  }

  
  if (flag->NL == 1){
    // Calculate |psi_out> += V_NL|psi_tmp>
    for (j = 0; j < 3; j++){
      p_nonlocal_proj_pot(psi_out, psi_tmp,nlc, nl, ist, par, par->ham_threads);
    }
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (j = 0; j < 10; j++){
      p_nonlocal_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par, par->ham_threads);
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    if (parallel->mpi_rank == 0) printf("\tNon-local potential: %.4g (msec)\n", (elapsed_seconds*1000.0)/10.0 ); fflush(0);
    
  }
  
  // Calculate the action of the local potential energy part of the Hamiltonian on psi_tmp
  clock_gettime(CLOCK_MONOTONIC, &start);
  for (int i = 0; i < 10; i++){
    if (1 == flag->periodic){
      e_ikr(psi_tmp, k, grid, ist, par, flag);
    }
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
  if (parallel->mpi_rank == 0) printf("\tLocal potential: %.4g (msec)\n", (elapsed_seconds*1000.0)/10.0 ); fflush(0);
  
  return;
}
