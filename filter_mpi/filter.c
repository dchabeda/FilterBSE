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
  int i;
  // File I/O
  FILE *pf; char fileName[100];
  // energy of the filtered states
  double *ene_filters;
  
  // FFT
  fftw_plan_with_nthreads(omp_get_max_threads());
  long fft_flags = FFTW_MEASURE;
  fftw_plan_loc planfw, planbw;
  fftw_complex *fftwpsi;
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
  if (parallel->mpi_rank == 0) printf("Starting filtering loop\n\n"); fflush(0);
  //long nthreads = (long) (parallel->nthreads / 2);
  // int chunk_size = ist->n_filter_cycles / parallel->mpi_size;
  // int start = parallel->mpi_rank * ist->n_filters_per_rank;
  // int end = (parallel->mpi_rank == parallel->mpi_size - 1) ? ist->n_filter_cycles : start + ist->n_filters_per_rank;
  
  // Loop over all of the random states handled by this mpi-rank
  for (jns = 0; jns < ist->n_filters_per_rank; jns++) {
     // Keep track of how many filter iterations have taken place
    if (parallel->mpi_rank == 0) {
      cntr++;
      time_t current_time = time(NULL);
      char* c_time_string = ctime(&current_time);
      if (parallel->mpi_rank == 0) printf("  Random psi %ld of %ld on node %d | %s\n", cntr, ist->n_filters_per_rank, parallel->mpi_rank, c_time_string); fflush(0);
    }
    for (jms = 0; jms < ist->m_states_per_filter; jms++){
      if (parallel->mpi_rank == 0) {
        if (parallel->mpi_rank == 0) printf("\tEnergy target %ld centered at %.3lg a.u.\n", jms, ene_targets[jms]); fflush(0);
      }
      // All variables declared within this parallel region are private for each thread
      // File I/O
      FILE *pf; char str[100];
      // Arrays for hamiltonian evaluation
      zomplex *psi, *phi;  
      // Array indexing
      long jstate, ns_block; // across states
      long jgrid, jgrid_real, jgrid_imag; // within a single state
      long jc, ncjms; // chebyshev coefficients
      
      // visualization
      double *rho, sgn_val;
      if ((rho = (double *) calloc(ist->ngrid, sizeof(rho[0]))) == NULL){
        fprintf(stderr, "\nOUT OF MEMORY: filter rho\n\n"); exit(EXIT_FAILURE);
      }

      // Allocate memory for arrays
      if ((psi = (zomplex*)calloc(ist->nspinngrid,sizeof(zomplex)))==NULL){ 
        fprintf(stderr, "\nOUT OF MEMORY: psi in run_filter_cycle\n\n"); exit(EXIT_FAILURE);
      }
      if ((phi = (zomplex*)calloc(ist->nspinngrid,sizeof(zomplex)))==NULL){ 
        fprintf(stderr, "\nOUT OF MEMORY: phi in run_filter_cycle\n\n"); exit(EXIT_FAILURE);
      }
      
      // Set index values for this iteration of the loop
      ns_block = (jns * ist->m_states_per_filter * ist->complex_idx * ist->nspinngrid); 
      jstate = ns_block + jms*ist->nspinngrid*ist->complex_idx;
      ncjms = ist->ncheby * jms; 
      
      
      // Populate psi with the random initial wavefunction corresponding to filter cycle jns
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
        jgrid_real = ist->complex_idx * jgrid;
        jgrid_imag = ist->complex_idx * jgrid + 1;

        psi[jgrid].re = psi_rank[ns_block + jgrid_real];
        if (1 == flag->isComplex){
          psi[jgrid].im = psi_rank[ns_block + jgrid_imag];
        }
      }
      
      // For each energy target, perform a filter procedure to reshape psi
      // into a state with energy close to the target energy
      
      // Calculate term 0 of the expansion
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
        jgrid_real = ist->complex_idx * jgrid;
        jgrid_imag = ist->complex_idx * jgrid + 1;

        psi_rank[jstate + jgrid_real] = an[ncjms+0].re * psi[jgrid].re - an[ncjms+0].im * psi[jgrid].im;
        if (1 == flag->isComplex){
          psi_rank[jstate + jgrid_imag] = an[ncjms+0].re * psi[jgrid].im + an[ncjms+0].im * psi[jgrid].re;
        }
      }
      
      // Calculate the subsequent terms of the expansion
      // sprintf (str, "prop-%ld-%ld-%ld.dat", thread_id, jns, jms); // for debugging parallelized loops
      // pf = fopen(str , "w");

      for (jc = 1; jc < ist->ncheby; jc++){

        memcpy(&phi[0], &psi[0], ist->nspinngrid * sizeof(phi[0]));
        p_hamiltonian(psi,phi,pot_local,nlc,nl,ksqr,ist,par,flag,planfw,planbw,fftwpsi);

        for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
          /*** par->dE_1 = 4.0 / par->dE and therefore I don't multiply by 4 ***/
          psi[jgrid].re = par->dE_1 * psi[jgrid].re - (2.0 + zn[jc-1] + par->Vmin * par->dE_1) * phi[jgrid].re;
          psi[jgrid].im = par->dE_1 * psi[jgrid].im - (2.0 + zn[jc-1] + par->Vmin * par->dE_1) * phi[jgrid].im;
        }

        for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
          jgrid_real = ist->complex_idx * jgrid;
          jgrid_imag = ist->complex_idx * jgrid + 1;

          psi_rank[jstate + jgrid_real] += (an[ncjms+jc].re * psi[jgrid].re - an[ncjms+jc].im * psi[jgrid].im);
          if (1 == flag->isComplex){
            psi_rank[jstate + jgrid_imag] += (an[ncjms+jc].re * psi[jgrid].im + an[ncjms+jc].im * psi[jgrid].re);
          }
        }
      }
    
      // Print out the filtered state for debugging
      if (1 == flag->printPsiFilt){
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
        sprintf(str, "psi-filt-%ld-%ld-%d.dat", jns, jms, parallel->mpi_rank);
        pf = fopen(str, "w");
        fwrite(&psi_rank[jstate], sizeof(psi_rank[0]), ist->nspinngrid*ist->complex_idx, pf);
        fclose(pf);
      }
      
      
      free(psi); free(phi);
    }
  }
  fftw_destroy_plan(planfw);
  fftw_destroy_plan(planbw);
  fftw_free(fftwpsi);
  
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