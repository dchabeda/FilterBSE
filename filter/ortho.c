#include "fd.h"

/*****************************************************************************/

long ortho(double *psitot, double dv, index_st *ist, par_st *par, flag_st *flag){
  /*******************************************************************
  * This function computes the singular value decomposition of a     *
  * real or complex M x N matrix, A. For us, A is the matrix of all  *
  * filtered states, with M = nspinngrid and N = t_rev*mn_states_tot *
  * The matrix is tall and skinny. We use the MKL LAPACK code        *
  * REAL VERSION:                                                    *
  *   https://www.netlib.org/lapack/explore-3.1.1-html/dgesvd.f.html *
  * COMPLEX VERSION:                                                 *
  *   https://www.netlib.org/lapack/explore-3.1.1-html/zgesvd.f.html *
  * inputs:                                                          *
  *  [psitot] arr holding all filtered wavefunctions                 *
  *  [dv] grid volume element dx*dy*dz                               *
  *  [ist] ptr to counters, indices, and lengths                     *
  *  [flag] ptr to flag_st holding job flags                         *
  * outputs: [long] number of orthogonal states after SVD            *
  ********************************************************************/
  

  long long lwork;  long long info, one=1, i, cutoff;
  long long ngrid = (long long)(ist->nspinngrid), mn_states_tot = (long long)(par->t_rev_factor*ist->mn_states_tot);
  double *S; 
  double *work; MKL_Complex16 *work_z;
  double *rwork;
  
  lwork = 5*(long long)(mn_states_tot*mn_states_tot+ngrid);
  S = (double*) malloc(mn_states_tot * sizeof(S[0]));
  
  if (0 == flag->isComplex) {
    work = (double*) malloc(lwork * sizeof(work[0]));

    // Do real-valued SVD
    dgesvd_("O","N",&(ngrid),&(mn_states_tot),&(psitot[0]),&(ngrid),&(S[0]),
	    NULL,&one,NULL,&one,&(work[0]),&(lwork),&info);
    
    if (info != 0) {
      fprintf(stderr, "error in dgesvd(2) %lld, exiting\n",info); exit(EXIT_FAILURE);
  }
  }
  if (1 == flag->isComplex) {
    work_z = (MKL_Complex16*) malloc(lwork * sizeof(MKL_Complex16));
    rwork = (double*) malloc(5*mn_states_tot * sizeof(rwork[0]));

    // Do complex-valued SVD
    zgesvd_("O","N",&(ngrid),&(mn_states_tot),&(psitot[0]),&(ngrid),&(S[0]),
	      NULL,&one,NULL,&one,&(work_z[0]),&(lwork),&(rwork[0]),&info);
    
    if (info != 0) {
      fprintf(stderr, "error in zgesvd(1) %lld, exiting\n",info); exit(EXIT_FAILURE);
    }

  }
  /*
  if (ngrid*mn_states_tot<2147483647){
    mkl_set_dynamic(0);
    mkl_set_num_threads(ist->nthreads);
    omp_set_nested(1);
  }
  else {
    mkl_set_dynamic(0);
    mkl_set_num_threads(1);
    omp_set_nested(1);
  }*/
  
  for (cutoff = ist->mn_states_tot, i = 0; i<ist->mn_states_tot; i++) {
    if ((S[i] / S[0]) < SVDEPS) {
      cutoff = i;
      break;
    }
  }
  printf("\nSVD cutoff (no. orthogonal vectors) is %lld\n",cutoff); fflush(0);
  
  // Free memory
  free(S);
  if (0 == flag->isComplex){
    free(work); 
  }
  if (1 == flag->isComplex) {
    free(rwork);
    free(work_z); 
  }

  return (cutoff);
}

/*****************************************************************************/

