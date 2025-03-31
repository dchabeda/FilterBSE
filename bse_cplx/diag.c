/*****************************************************************************/

#include "fd.h"

/*****************************************************************************/
// Uses LAPACK zheev_ subroutine to diagonalize the cplx, Hermitian matrix mat
// the eigenvalues are stored in eval and eigenvectors go into mat if jobz='V'

void diag(const int n, int nthreads, double complex *mat, double *eval) {
  // n is order of the matrix mat (also leading dimension of mat)
  //MKL_INT lwork = 2*n-1; // length of the array work. lwork >= max(1,3*n-1)
  lapack_int info; // signifies successful exit or not -> used to handle errors
  char jobz; 
  char uplo;
  lapack_int N = n;
  //MKL_Complex16 *work; // dimension (MAX(1,lwork))
  double *rwork; //dimension (MAX(1,3*n-2))

  // Use multiple threads
  //mkl_set_dynamic(0);
  mkl_set_num_threads(1);
  //omp_set_nested(1);

  // Allocate memmory for work array
  //rwork = (double *) calloc(3*n-2,sizeof(double));
  //work = (MKL_Complex16 *) calloc(lwork, sizeof(MKL_Complex16));
  
  // Set what will be computed along with if matrix is upper or lower traingular
  jobz = 'V'; // Compute eigenvalues and eigenvectors ('N' for eigenvalues only) 
  uplo = 'U'; // mat matrix is upper triangular

  
  // Diagonalize the complex, hermitian matrix

  //printf("zheev( jobz= %c, uplo= %c, n= %d, &mat= %x, n= %d, &eval= %x, &work= %x, lwork= %d, &rwork= %x, info= %d)\n",
  //  jobz, uplo, n, (MKL_Complex16*)&mat[0], n, &eval[0], &work[0], lwork, &rwork[0], info);

  info = LAPACKE_zheev(LAPACK_ROW_MAJOR, jobz, uplo, N, (lapack_complex_double*)&mat[0], N, &eval[0]);

  // Error handling
  if (info < 0) {
    printf("Argument %lld of zysev_ in diag.c had an illegal argument\n", (-info));
    printf("The program is exiting due to this error!!!\n");
    exit(EXIT_FAILURE);
  } else if (info > 0) {
    printf("zysev_ failed to converge -> the program is exiting!\n");
    exit(EXIT_FAILURE);
  }

  // Free dynamically allocated memory that was only internally used by MKL 
  //free(work);
  //free(rwork);

  return;
}

// void diag(const int n, int nthreads, double complex *mat, double *eval) {
//   // n is order of the matrix mat (also leading dimension of mat)
//   int info; // signifies successful exit or not -> used to handle errors
  
//   // Set what will be computed along with if matrix is upper or lower traingular
//   const char jobz = 'V'; // Compute eigenvalues and eigenvectors ('N' for eigenvalues only) 
//   const char uplo = 'U'; // mat matrix is upper triangular

//   const int N = n;
//   const int lwork = 4 * N;
  
//   // double complex *work_z; // dimension (MAX(1,lwork))
//   MKL_Complex16 *work_z;
//   double *rwork;          //dimension (MAX(1,3*n-2))

//   // Allocate memory for work array
//   ALLOCATE(&rwork, lwork, "rwork in diag");
//   ALLOCATE(&work_z, lwork, "work_z in diag");
  
//   // Diagonalize the complex, hermitian matrix

//   zheev_(&jobz, &uplo, &N, (MKL_Complex16*)&mat[0], &N, &eval[0],&work_z[0], &lwork, &rwork[0], &info);
  
//   // Error handling
//   if (info){
//     fprintf(stderr, "error in zheev_ H\n"); exit(EXIT_FAILURE);
//   }

//   // Free dynamically allocated memory that was only internally used by LAPACK 
//   free(work_z);
//   free(rwork);

//   return;
// }

/*****************************************************************************/
