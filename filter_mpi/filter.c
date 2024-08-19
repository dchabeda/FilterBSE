#include "fd.h"

/*****************************************************************************/

void run_filter_cycle(double *psitot, double *pot_local, nlc_st *nlc, long *nl, 
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
  long jmn, jns, jms, cntr = 0;
  // File I/O
  FILE *pf; char fileName[100];
  // energy of the filtered states
  double *ene_filters;
  long MAX_SIZE = ist->complex_idx*ist->nspinngrid*((ist->n_filter_cycles - 1)*ist->m_states_per_filter + (ist->m_states_per_filter - 1)) + ist->nspinngrid*ist->complex_idx - 1;
  
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
  
  // Send necessary arrays to each node
  if (mpi_rank == 0) {
  // Fill data array with values
  for (int i = 1; i < mpi_size; i++) {
      MPI_Send(data, 100, MPI_INT, i, 0, MPI_COMM_WORLD);
  }
  } else {
      // Other processes receive the data
      MPI_Recv(data, 100, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  //************************************************************************
  //************************************************************************
  // BEGIN FILTERING HERE
  //************************************************************************
  //************************************************************************
  printf("Starting filtering loop\n\n"); fflush(0);
  //long nthreads = (long) (parallel->nthreads / 2);
  int chunk_size = ist->mn_states_tot / parallel->mpi_size;
  int start = mpi_rank * chunk_size;
  int end = (mpi_rank == parallel->mpi_size - 1) ? ist->mn_states_tot : start + chunk_size;
  
  for (jmn = start; jmn < end; jmn++) {
     // Keep track of how many filter iterations have taken place
    if (mpi_rank == 0) {
      cntr++;
      printf("\tCurrently working on iteration %ld of filtering cycle\n", cntr); fflush(0);
      printf("\t  (~%ld states)\n", cntr * nthreads); fflush(0);
    }
    // All variables declared within this parallel region are private for each thread
    // thread tracking
    long thread_id;
    // File I/O
    FILE *pf; char str[100];
    // Arrays for hamiltonian evaluation
    zomplex *psi, *phi, *psi_out;  
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

    thread_id = omp_get_thread_num();	
    // Allocate memory for arrays
    if ((psi = (zomplex*)calloc(ist->nspinngrid,sizeof(zomplex)))==NULL){ 
      fprintf(stderr, "\nOUT OF MEMORY: psi in run_filter_cycle\n\n"); exit(EXIT_FAILURE);
    }
    if ((phi = (zomplex*)calloc(ist->nspinngrid,sizeof(zomplex)))==NULL){ 
      fprintf(stderr, "\nOUT OF MEMORY: phi in run_filter_cycle\n\n"); exit(EXIT_FAILURE);
    }
    if ((psi_out = (zomplex*)calloc(ist->nspinngrid,sizeof(zomplex)))==NULL){ 
      fprintf(stderr, "\nOUT OF MEMORY: psi_out in run_filter_cycle\n\n"); exit(EXIT_FAILURE);
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
    
    // Populate psi with the random initial wavefunction corresponding to filter cycle jns
    for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
      jgrid_real = ist->complex_idx * jgrid;
      jgrid_imag = ist->complex_idx * jgrid + 1;

      psi[jgrid].re = psitot[ns_block + jgrid_real];
      if (1 == flag->isComplex){
        psi[jgrid].im = psitot[ns_block + jgrid_imag];
      }
    }
    
    // For each energy target, perform a filter procedure to reshape psi
    // into a state with energy close to the target energy
    
    // Calculate term 0 of the expansion
    for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
      psi_out[jgrid].re = an[ncjms+0].re * psi[jgrid].re - an[ncjms+0].im * psi[jgrid].im;
      if (1 == flag->isComplex){
        psi_out[jgrid].im = an[ncjms+0].re * psi[jgrid].im + an[ncjms+0].im * psi[jgrid].re;
      }
    }
    
    // Calculate the subsequent terms of the expansion
    sprintf (str, "prop-%ld-%ld-%ld.dat", thread_id, jns, jms); // for debugging parallelized loops
    pf = fopen(str , "w");

    for (jc = 1; jc < ist->ncheby; jc++){

      memcpy(&phi[0], &psi[0], ist->nspinngrid * sizeof(phi[0]));
      hamiltonian(psi,phi,pot_local,nlc,nl,ksqr,ist,par,flag,planfw,planbw,fftwpsi);

      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
        /*** par->dE_1 = 4.0 / par->dE and therefore I don't multiply by 4 ***/
        psi[jgrid].re = par->dE_1 * psi[jgrid].re - (2.0 + zn[jc-1] + par->Vmin * par->dE_1) * phi[jgrid].re;
        psi[jgrid].im = par->dE_1 * psi[jgrid].im - (2.0 + zn[jc-1] + par->Vmin * par->dE_1) * phi[jgrid].im;
      }

      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
        psi_out[jgrid].re += (an[ncjms+jc].re * psi[jgrid].re - an[ncjms+jc].im * psi[jgrid].im);
        if (1 == flag->isComplex){
          psi_out[jgrid].im += (an[ncjms+jc].re * psi[jgrid].im + an[ncjms+jc].im * psi[jgrid].re);
        }
      }

      if (0 == (jc % 100)) {
        fprintf (pf,"%ld %ld %ld\n", jns, jms, jc); fflush(pf);
      }

      // if (jc < 20){
      //   // Print out the evolution of the filter states through the cheby iterations
      //   for (jgrid = 0; jgrid < ist->ngrid; jgrid++){
      //     jgrid_real = ist->complex_idx * jgrid;
      //     jgrid_imag = ist->complex_idx * jgrid + 1;

      //     rho[jgrid] = sqr(psitot[jns_ms + jgrid_real]);
      //     sgn_val = psitot[jns_ms + jgrid_real];
      //     if (1 == flag->isComplex){
      //       rho[jgrid] += sqr(psitot[jns_ms + jgrid_imag]);
      //       if (sgn_val < psitot[jns_ms + jgrid_imag]) sgn_val = psitot[jns_ms + jgrid_imag];
      //     }
      //     rho[jgrid] *= sign(sgn_val);
      //   }
      //   sprintf(str, "psi-filt-%ld-%ld-%ld.cube", jns, jms, jc);
      //   write_cube_file(rho, grid, str);
      // }
    }
    
    fclose(pf);
    
    /*****************************************************************************/
    // Copy the filtered wavefunctions back into psitot
    /*****************************************************************************/


    for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
      jgrid_real = ist->complex_idx * jgrid;
      jgrid_imag = ist->complex_idx * jgrid + 1;
      psitot[jns_ms + jgrid_real] = psi_out[jgrid].re;
      if (1 == flag->isComplex){
        psitot[jns_ms + jgrid_imag] = psi_out[jgrid].im;
      }
    }
    
    free(psi); free(phi); free(psi_out);
    fftw_destroy_plan(planfw);
    fftw_destroy_plan(planbw);
    fftw_free(fftwpsi);
  }
  
  /***********************************************************************/
  /*** normalize the states and get their energies***/
  printf("Normalizing filtered states\n"); fflush(stdout);
  normalize_all(psitot,ist,par,flag,parallel);
  // exit(0);
  // Get the energy of all the filtered states
  if ((ene_filters = (double*)calloc(ist->mn_states_tot,sizeof(double)))==NULL)nerror("ene_filters");
  
  printf("Computing the energies of all filtered states\n"); fflush(stdout);
  /*** calculate and print the energy of the filtered states ***/
  energy_all(psitot,pot_local,nlc,nl,ksqr,ene_filters,ist,par,flag,parallel);

  for (jns = 0; jns < ist->n_filter_cycles; jns++){
    sprintf (fileName,"ene-filt-jns-%ld.dat", jns);
    pf = fopen(fileName , "w");
    for (jms = 0; jms < ist->m_states_per_filter; jms++){
      jmn = jns * ist->m_states_per_filter + jms;
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