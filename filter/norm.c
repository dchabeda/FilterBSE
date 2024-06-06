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
  for (i = 0; i < ngrid; i++) norm += (sqr(psi[i].re) + sqr(psi[i].im));
  
  return (sqrt(norm * dv));
}

/*****************************************************************************/
// Normalizes a single state of length ngrid and grid volume of dv
// works with complex a wavefunction (complex spinors as well)

double normalize(zomplex *psi,double dv,long ngrid,long nthreads){
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
  double N = calc_norm(psi,dv,ngrid,nthreads);
  //printf("norm in normalize = %g\n", N);
  omp_set_dynamic(0);
  omp_set_num_threads(nthreads);
#pragma omp parallel for 
  for (k = 0; k < ngrid; k++){
    psi[k].re /= N;
    psi[k].im /= N;
  }

  return (N);
}

/*****************************************************************************/
// Normalizes to 1 all numStates states in psi each of which is of length ngrid 
// with a grid volume of dv works with a complex wavefunctions (complex spinors as well)

void normalize_all(double *psitot, double dv, long ms, long ngrid, long nthreads, int complex_idx, int printNorm){
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

  long jgrid, jgrid_real, jgrid_imag, ie, ieg;
  double norm;
  //printf("This is dv coming into normalize: %g\n", dv);
  //printf("This is ngrid coming into normalize: %ld\n", ngrid);
  
  omp_set_dynamic(0);
  omp_set_num_threads(nthreads);
#pragma omp parallel for private(ie, norm, ieg, jgrid)
  for (ie = 0; ie < ms; ie++){
    // Loop over all states
    // complex_idx is 2 if psi is complex (1 if real), this accounts for storing real and imag components
    ieg = complex_idx * ie * ngrid; 
    
    // Calculate the norm integral for each state
    norm = 0.0;
    for (jgrid = 0; jgrid < ngrid; jgrid++){
      jgrid_real = complex_idx * jgrid; // indexes the real valued wavefunc OR real component of complex val func
      jgrid_imag = complex_idx * jgrid + 1; // indexes complex component

      norm += sqr(psitot[ieg+jgrid_real]) ;
      // only add imag part to norm if psi is complex
      if (2 == complex_idx){ 
        norm += sqr(psitot[ieg+jgrid_imag]); 
      }
    }
    norm *= dv;
    
    if (printNorm == 1) {
      printf("state %ld norm = %g\n", ie, norm); fflush(0);
    }
    
    // Handle exception behavior. If norm is zero or nan, something failed
    if (norm != 0.0) {
      norm = 1.0 / sqrt(norm);
    } else if (isnan(norm)) {
      fprintf(stderr, "Norm is NAN! exiting!\n"); exit(EXIT_FAILURE); 
    } else {
      fprintf(stderr, "Norm is %.2g! exiting!\n", norm); exit(EXIT_FAILURE); 
    }

    for (jgrid = 0; jgrid < ngrid; jgrid++) {
      jgrid_real = complex_idx * jgrid;
      jgrid_imag = complex_idx * jgrid + 1;

      psitot[ieg+jgrid_real] *= norm;
      if (2 == complex_idx){
        psitot[ieg+jgrid_imag] *= norm;
      }
    }
  }

  return;
}

/*****************************************************************************/
