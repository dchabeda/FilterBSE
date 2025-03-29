#include "aux.h"

/*****************************************************************************/
char* format_duration(double elapsed_seconds) {
    // Calculate hours, minutes, and seconds
    int hours = (int)(elapsed_seconds / 3600);
    int minutes = (int)((elapsed_seconds - (hours * 3600)) / 60);
    int seconds = (int)(elapsed_seconds - (hours * 3600) - (minutes * 60));
    
    // Allocate memory for the string
    char* result = (char*)malloc(16 * sizeof(char)); // "hh:mm:ss" + null terminator

    // Format the string
    if (result != NULL) {
        snprintf(result, 16, "%02d:%02d:%02d", hours, minutes, seconds);
    }

    return result;
}

/*****************************************************************************/
char* get_time(void){
    time_t current_time;
    char* c_time_string;

    current_time = time(NULL);
    c_time_string = ctime(&current_time);
    return c_time_string;
}

/*****************************************************************************/
void allocate_memory(void **ptr, size_t length, size_t type_size, char* message) {
    *ptr = calloc(length, type_size);
    if (*ptr == NULL) {
        fprintf(stderr, "Memory allocation failed: %s\n", message);
        exit(EXIT_FAILURE);
    }
    return;
}

/*****************************************************************************/
void matmul(int M, int N, int K, double *A, double *B, double *X){
    // Function to multiply two matrices together:
    // [X] = [A] * [B]
    // where A (M,K) and B (K,N)

    int i, j, k;
    double sum;

    for (i = 0; i < M; i++){
        for (j = 0; j < N; j++){
            sum = 0.0;
            for (k = 0; k < K; k++){
                sum += A[i * M + k] * B[k * N + j];
            }
            X[i * M + j] = sum;
        }
    }

    return;
}

/*****************************************************************************/

void trans_mat(int N, double *U, double *A, double *Ap){
    // Function that transforms matrix A according to matrix X:
    // Ap = X^T * A * X
    // A must be a square matrix, so M = K = N.
    int i, j, k, q;
    double sum;

    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            sum = 0.0;
            for (k = 0; k < N; k++){
                for (q = 0; q < N; q++){
                    sum += U[k * N + i] * A[k * N + q] * U[q * N + j];
                }
            }
            Ap[i * N + j] = sum;
        }
    }
}


/*****************************************************************************/

void diag_mat(double *mat, double *eigv, int n_dim){
    /*******************************************************************
    * This function calculates eigenvalues and vectors of the real or  *
    * complex valued matrix, H, where H_ij = <psi_i|H|psi_j>           *
    * It utilizes the MKL LAPACK routines to do the diagonalization    *
    * REAL VALUED                                                      *
    *  https://www.netlib.org/lapack/explore-html-3.6.1/d2/d8a/group__double_s_yeigen_ga442c43fca5493590f8f26cf42fed4044.html#ga442c43fca5493590f8f26cf42fed4044
    * COMPLEX VALUED                                                   *
    *   https://www.netlib.org/lapack/explore-html-3.6.1/d6/dee/zheev_8f_a70c041fd19635ff621cfd5d804bd7a30.html
    */
    long long info; 
    long long lwk = (long long) 3 * n_dim;
    long long N = (long long) n_dim;
    double *work;
  
    // Allocate memory for scratch work
    work = (double *) calloc(lwk, sizeof(double));
    
    //
    //
    char jobu = 'V';
    char jobvt = 'U';
    dsyev_(&jobu, &jobvt, &N, &mat[0], &N, &eigv[0], &work[0], &lwk, &info , 1, 1);
    //
    //
  
    free(work);
}

/*****************************************************************************/

// void mpi_print(const char *message) {
//     int rank;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     if (rank == 0) {
//         printf("%s", message);
//     }
// }

/*****************************************************************************/

void print_progress_bar(int cur, int tot){
    // print the filtering progress to the output file
    int barWidth = 16; // Width of the progress bar
    float percent = (float)cur / tot * 100;
    int pos = barWidth * cur / tot;

    printf("\t  [");
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) printf("#");
        else printf(" ");
    }
    printf("] %3.0f%% | %s\n", percent, get_time());
    fflush(stdout);

    return;
}

/*****************************************************************************/

int sign(float x) {
    return (int)copysign(1.0, x);  // copysign gives the sign of x
}

/*****************************************************************************/

