/*****************************************************************************/
// File contains functions that are involved with normalizing wavefunctions

#include "fd.h"

/*****************************************************************************/
// Returns the norm of a single state of length ngrid and grid volume of dv
// works with complex a wavefunction (complex spinors as well)

double calc_norm(zomplex *psi, double dv, long ngrid, long nthreads){
  /*******************************************************************
  * This function integrates the square mod of the wavefunc          *
  *     N = sum^{ngrid}_{0} |Psi(gridpt)|^2 * dv
  * inputs:                                                          *
  *  [psi] ngrid-long array of double/zomplex for storing wavefnc    *
  *  [dv] grid volume element dx * dy * dz                           *
  *  [ngrid] number of grid points                                   *
  *  [nthreads] number of OMP threads for parallelization            *
  * outputs: [double] the norm                                       *
  ********************************************************************/

  long i;
  double norm = 0.0;
  
  omp_set_dynamic(0);
  omp_set_num_threads(nthreads);
#pragma omp parallel for reduction(+:norm)
  for (i = 0; i < ngrid; i++){
    norm += (sqr(psi[i].re) + sqr(psi[i].im));
  }
  
  return (sqrt(norm * dv));
}

/*****************************************************************************/
// Normalizes a single state of length ngrid and grid volume of dv
// works with complex a wavefunction (complex spinors as well)

double normalize(zomplex *psi, long ngrid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){
  /*******************************************************************
  * This function integrates the square mod of the wavefunc          *
  *     N = sum^{ngrid}_{0} |Psi(gridpt)|^2 * dv                     *
  * inputs:                                                          *
  *  [psi] ngrid-long array of double/zomplex for storing wavefnc    *
  *  [dv] grid volume element dx * dy * dz                           *
  *  [ngrid] number of grid points                                   *
  *  [nthreads] number of OMP threads for parallelization            *
  * outputs: [double] the norm                                       *
  ********************************************************************/

  long k;
  double N = calc_norm(psi,par->dv,ngrid,parallel->nthreads);

  // if (parallel->mpi_rank == 0) printf("norm in normalize = %g\n", N);
  // if (flag->printNorm == 1) {
  //     FILE *pf;
  //     char str[30];
  //     sprintf(str, "norm-%ld.dat", thread_id);
  //     pf = fopen(str, "w");
  //     fprintf(pf, "norm = %g\n", N); fflush(pf);
  //     fclose(pf);
  //   }
  omp_set_dynamic(0);
  omp_set_num_threads(parallel->nthreads);
#pragma omp parallel for private(k)
  for (k = 0; k < ngrid; k++){
    psi[k].re /= N;
    psi[k].im /= N;
  }
  
  return (N);
}

/*****************************************************************************/
// Normalizes to 1 all numStates states in psi each of which is of length ngrid 
// with a grid volume of dv works with a complex wavefunctions (complex spinors as well)

void normalize_all(double *psitot, long n_states, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){
  /*******************************************************************
  * This function normalizes to 1.0 all states in psitot             *
  * inputs:                                                          *
  *  [psitot] ngrid-long array of double/zomplex holding wavefncs    *
  *  [dv] grid volume element dx * dy * dz                           *
  *  [ms] number of states held in psitot                            *
  *  [ngrid] number of grid points                                   *
  *  [nthreads] number of OMP threads for parallelization            *
  *  [printNorm] flag to print wavefunc norms (debugging)            *
  * outputs: [double] the norm                                       *
  ********************************************************************/

  // NOTE: If your calculation fails with error "norm is 0/NaN!", then  
  // check what the calculated energy range of the Hamiltonian is in your calculation.
  // If the range is very small, then it likely means that the energy min was not
  // found by the algorithm. You can track the iterative convergence of E_min in the file
  // Emin-init.dat. To better ensure convergence of the E_min, decrease the value of the tau
  // parameter in get_energy_range energy.c
  // if (parallel->mpi_rank == 0) printf("Entered normalize_all on node %d\n", parallel->mpi_rank);
  // if (parallel->mpi_rank == 0) printf("n_states = %ld node %d\n", n_states, parallel->mpi_rank); fflush(0);
  long jns;
  //if (parallel->mpi_rank == 0) printf("This is dv coming into normalize: %g\n", dv);
  //if (parallel->mpi_rank == 0) printf("This is ngrid coming into normalize: %ld\n", ngrid);
  // FILE *pf;
  //     char str[30];
  //     sprintf(str, "norm-all-%ld.dat", thread_id);
  //     pf = fopen(str, "w");
  // omp_set_dynamic(0);
  // omp_set_num_threads(parallel->nthreads);
// #pragma omp parallel for private(jns)
  for (jns = 0; jns < n_states; jns++){
    // if (parallel->mpi_rank == 0) printf("jns = %ld on node %d\n", jns, parallel->mpi_rank); fflush(0);
    long j_state;
    long jgrid, jgrid_real, jgrid_imag;
    double norm;
    // Loop over all states
    // complex_idx is 2 if psi is complex (1 if real), this accounts for storing real and imag components
    j_state = jns * ist->complex_idx * ist->nspinngrid; 
    // Calculate the norm integral for each state
    norm = 0.0;
    for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
      jgrid_real = ist->complex_idx * jgrid; // indexes the real valued wavefunc OR real component of complex val func
      jgrid_imag = ist->complex_idx * jgrid + 1; // indexes complex component

      norm += sqr(psitot[j_state+jgrid_real]) ;
      // only add imag part to norm if psi is complex
      if (1 == flag->isComplex){ 
        norm += sqr(psitot[j_state+jgrid_imag]); 
      }
    }
    norm *= par->dv;
    
    if (flag->printNorm == 1) {
      if (parallel->mpi_rank == 0) printf("state %ld norm = %g\n", jns, norm); fflush(0);
    }
    
    // Handle exception behavior. If norm is zero or nan, something failed
    if (norm != 0.0) {
      norm = 1.0 / sqrt(norm); // 1/srt(N)
    } else if (isnan(norm)) {
      fprintf(stderr, "Norm is NAN! exiting!\n"); exit(EXIT_FAILURE); 
    } else {
      fprintf(stderr, "Norm is %.2g! exiting!\n", norm); exit(EXIT_FAILURE); 
    }
    // Rescale the wavefunctions by the normalization constant to normalize them
    for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
      jgrid_real = ist->complex_idx * jgrid;
      jgrid_imag = ist->complex_idx * jgrid + 1;

      psitot[j_state+jgrid_real] *= norm;
      if (1 == flag->isComplex){
        psitot[j_state+jgrid_imag] *= norm;
      }
    }
  }
  // fclose(pf);
  return;
}

/*****************************************************************************/
