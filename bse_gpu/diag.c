/*****************************************************************************/

#include "fd.h"

/*****************************************************************************/
// Uses LAPACK dysev_ subroutine to diagonalize the real, symmetric matrix mat
// the eigenvalues are stored in eval and eigenvectors go into mat if jobz='V'

void diag(const int n, int nthreads, double _Complex *mat, double *eval) {
  // n is order of the matrix mat (also leading dimension of mat)
  //MKL_INT lwork = 2*n-1; // length of the array work. lwork >= max(1,3*n-1)
  int info; // signifies successful exit or not -> used to handle errors
  
  // Set what will be computed along with if matrix is upper or lower traingular
  const char jobz = 'V'; // Compute eigenvalues and eigenvectors ('N' for eigenvalues only) 
  const char uplo = 'U'; // mat matrix is upper triangular

  const int N = n;
  const int lwork = 3 * N;
  // double _Complex *mat_conv = (double _Complex *) malloc(N * N * sizeof(double _Complex));
  double _Complex *work_z; // dimension (MAX(1,lwork))
  double *rwork; //dimension (MAX(1,3*n-2))

  // Use multiple threads
  //mkl_set_dynamic(0);
  //omp_set_nested(1);


  // Allocate memory for work array
  rwork = (double *) calloc(3*N - 2, sizeof(double));
  work_z = (double _Complex *) calloc(lwork, sizeof(double _Complex));
  
  // Diagonalize the complex, hermitian matrix

  // Convert from column major to row major
  // for (int i = 0; i < N * N; i++) {
  //   mat_conv[i] = mat[i].re + mat[i].im * I;
  // }

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
