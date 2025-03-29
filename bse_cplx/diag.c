/*****************************************************************************/

#include "fd.h"

/*****************************************************************************/
// Uses LAPACK zheev_ subroutine to diagonalize the cplx, Hermitian matrix mat
// the eigenvalues are stored in eval and eigenvectors go into mat if jobz='V'

void diag(const int n, int nthreads, double complex *mat, double *eval) {
  // n is order of the matrix mat (also leading dimension of mat)
  int info; // signifies successful exit or not -> used to handle errors
  
  // Set what will be computed along with if matrix is upper or lower traingular
  const char jobz = 'V'; // Compute eigenvalues and eigenvectors ('N' for eigenvalues only) 
  const char uplo = 'U'; // mat matrix is upper triangular

  const int N = n;
  const int lwork = 3 * N;
  
  double complex *work_z; // dimension (MAX(1,lwork))
  double *rwork;          //dimension (MAX(1,3*n-2))

  // Allocate memory for work array
  ALLOCATE(&rwork, 3*N - 2, "rwork in diag");
  ALLOCATE(&work_z, lwork, "work_z in diag");
  
  // Diagonalize the complex, hermitian matrix

  zheev_(&jobz, &uplo, &N, &mat[0], &N, &eval[0], &work_z[0], &lwork, &rwork[0], &info);
  
  // Error handling
  if (info){
    fprintf(stderr, "error in zheev_ H\n"); exit(EXIT_FAILURE);
  }

  // Free dynamically allocated memory that was only internally used by LAPACK 
  free(work_z);
  free(rwork);

  return;
}

/*****************************************************************************/
