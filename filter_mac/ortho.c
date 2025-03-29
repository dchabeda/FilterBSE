#include "ortho.h"

/*****************************************************************************/

long ortho_real(double *psitot, double dv, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){
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
  

  int info;
  int i, cutoff;
  const int one=1;
  const int ngrid = (ist->nspinngrid);
  const int mn_states_tot = (par->t_rev_factor * ist->mn_states_tot);
  double *S; 
  double *work;
  const int mpir = parallel->mpi_rank;

  const int lwork = 5 * (mn_states_tot * mn_states_tot+ngrid);

  S = (double*) malloc(mn_states_tot * sizeof(S[0]));
  
  if (mpir == 0) printf("Scalar wavefunctions used. Allocating \"work\" for real SVD\n"); fflush(0);
  ALLOCATE(&work, lwork, "ortho work");
  // if ((work = (double*) malloc(lwork * sizeof(work[0]))) == NULL){
  //   fprintf(stderr, "\nOUT OF MEMORY: ortho work\n\n"); exit(EXIT_FAILURE);
  // }
  
  if (mpir == 0) printf("Computing real-valued SVD...\n"); fflush(0);
  
  // Do real-valued SVD
  char jobu = 'O'; 
  char jobvt = 'N';
  dgesvd_(&jobu, &jobvt, &(ngrid), &(mn_states_tot), &(psitot[0]), &(ngrid), &(S[0]),
    NULL,&one,NULL,&one,&(work[0]),&(lwork),&info , 1, 1);
  
  if (info != 0) {
    fprintf(stderr, "error in dgesvd(2) %d, exiting\n",info); fflush(stderr); 
    exit(EXIT_FAILURE);
  }
  
  if (mpir == 0) printf("Done with SVD!\n"); fflush(0);


  for (cutoff = mn_states_tot, i = 0; i < mn_states_tot; i++) {
    if ((S[i] / S[0]) < SVDEPS) {
      cutoff = i;
      break;
    }
  }

  if (mpir == 0) printf("\nSVD cutoff (no. orthogonal vectors) is %d\n", cutoff); fflush(0);
  
  // Free memory
  free(S);
  
  if (mpir == 0) printf("Scalar wavefunctions used. Freeing work\n"); fflush(0);
  free(work); 
  
  return (cutoff);
}

/*****************************************************************************/

long ortho_cplx(zomplex *psitot, double dv, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){
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
  
  const int mpir = parallel->mpi_rank;
  int i, cutoff;
  long long lwork;
  long long info, one=1;
  long long ngrid = (long long)(ist->nspinngrid);
  long long mn_states_tot = (long long)(par->t_rev_factor * ist->mn_states_tot);

  double *S = malloc(MIN(ngrid, mn_states_tot) * sizeof(double)); 
  lapack_complex_double *work_z;
  double *rwork;

  // Query optimal work size
  lapack_complex_double work_query[1];
  lwork = -1;

  char jobu = 'O';
  char jobvt = 'N';
  zgesvd_(&jobu, &jobvt, &ngrid, &mn_states_tot, 
          (lapack_complex_double*)psitot, &ngrid, 
          S, 
          NULL, &one, 
          NULL, &one, 
          work_query, &lwork, 
          rwork, 
          &info 
          FORTRAN_STRLEN_PASS);

  // Allocate work arrays based on queried size
  lwork = (long long)work_query[0];
  work_z = malloc(lwork * sizeof(lapack_complex_double));
  rwork = malloc(lwork * sizeof(double));

  // Perform SVD
  zgesvd_(&jobu, &jobvt, &ngrid, &mn_states_tot, 
          (lapack_complex_double*)psitot, &ngrid, 
          S, 
          NULL, &one, 
          NULL, &one, 
          work_z, &lwork, 
          rwork, 
          &info 
          FORTRAN_STRLEN_PASS);

  // Check for errors
  if (info != 0) {
      fprintf(stderr, "ZGESVD error: %lld\n", info);
      // Handle error
  }
  
  for (cutoff = mn_states_tot, i = 0; i < mn_states_tot; i++) {
    if ((S[i] / S[0]) < SVDEPS) {
      cutoff = i;
      break;
    }
  }

  if (mpir == 0) printf("Done with SVD\n");
  if (mpir == 0) printf("\nSVD cutoff (no. orthogonal vectors) is %lld\n",cutoff); fflush(0);
  
  // Free memory
  free(S);
  
  if (mpir == 0) printf("Spinor wavefunctions used. Freeing rwork and work_z\n"); fflush(0);
  free(rwork);
  free(work_z); 
  

  return (cutoff);
}

/*****************************************************************************/


void restart_from_ortho(
  double*       psitot,
  index_st*     ist,
  par_st*       par,
  flag_st*      flag,
  parallel_st*  parallel){
  
  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  FILE *ppsi;
  long tot_sz = ist->complex_idx*ist->nspinngrid*ist->mn_states_tot;
  
  /************************************************************/
  /*******************     READ IN PSITOT   *******************/
  /************************************************************/

  write_separation(stdout, "T");
  printf("****  RESTART FROM ORTHO  **  RESTART FROM ORTHO  **  RESTART FROM ORTHO  ****");
  write_separation(stdout, "B"); fflush(stdout);
  
  
  printf("\nNo. states for orthogonalization = %ld\n", ist->mn_states_tot);
  printf("Size of psitot array = %.2g GB\n", 
    sizeof(double) * (double)tot_sz/1024/1024/1024
  ); 
  fflush(stdout);
  
  // Read in the states from psi-filt.dat file
  
  ppsi = fopen("psi-filt.dat", "r");
  
  if (ppsi != NULL){
    printf("Reading psi-filt.dat\n"); fflush(stdout);
    fread(psitot, sizeof(double), tot_sz, ppsi);
    fclose(ppsi);
  } else{
    fprintf(stderr, "ERROR: psi-filt.dat could not be opened\n");
    exit(EXIT_FAILURE);
  }

  printf("psitot[max] = %lg\n", psitot[(tot_sz) - 1]);
  
  
  printf("\nNormalizing filtered states (for safety)\n"); fflush(stdout);

  normalize_all(psitot, ist->mn_states_tot, ist, par, flag, parallel);
  
  printf("Done normalizing\n"); fflush(stdout);
  
  return;
}
