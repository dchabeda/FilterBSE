#include "fd.h"

void allocate_memory(void **ptr, size_t length, size_t type_size, char *message);
void matmul(int M, int N, int K, double *A, double *B, double *X);
void trans_mat(int N, double *U, double *A, double *Ap);
void diag_mat(double *mat, double *eigv, int n_dim);
void print_progress_bar(int cur, int tot);
void read_psi_dat(FILE *pf, double *psi, long start, long end, long len, int n_every);
char *format_duration(double elapsed_seconds);
char *get_time(void);
int sign(float x);
void load_eval(double *eig_vals, double *sigma_E, index_st *ist, par_st *par);
void load_psitot(double *psitot, index_st *ist, par_st *par);
