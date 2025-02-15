#include "fd.h"

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
char* get_time(){
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

    // Obtain the current time
    time_t current_time = time(NULL);
    // Convert to local time format and print
    char* c_time_string = ctime(&current_time);
    printf("\t  [");
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) printf("#");
        else printf(" ");
    }
    printf("] %3.0f%% | %s\n", percent, c_time_string);
    fflush(stdout);

    return;
}

/*****************************************************************************/
