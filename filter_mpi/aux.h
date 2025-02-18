#include "fd.h"

void allocate_memory(void **ptr, size_t length, size_t type_size, char* message);
void matmul(int M, int N, int K, double *A, double *B, double *X);
void trans_mat(int N, double *U, double *A, double *Ap);
void diag_mat(double *mat, double *eigv, int n_dim);
void print_progress_bar(int cur, int tot);
char* format_duration(double elapsed_seconds);
char* get_time();