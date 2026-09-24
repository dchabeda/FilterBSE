#include "fd.h"

void write_cube_file(double *rho, grid_st *grid, char *fileName);
void write_current_time(FILE *pf);
void write_separation(FILE *pf, char *top_bttm);
void write_state_dat(double complex *psi, long n_elems, char* fileName);
void write_BSE_coeff(double complex *bs_coeff, long n_xton, char *fileName);
double density_phase_sign(double re, double im);
void write_exciton_densities(double complex *bs_coeff, double complex *psi_qp,
                             grid_st *grid, index_st *ist, par_st *par,
                             flag_st *flag, parallel_st *parallel);
