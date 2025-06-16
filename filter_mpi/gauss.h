#include "fd.h"
#include "aux.h"
#include "integral.h"

void mod_gauss(
    gauss_st *gauss, double *pot_local, double * eig_vals, xyz_st *R, atom_info *atom, grid_st *grid,
    index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel
);

void init_gauss_params( 
    gauss_st     *gauss,   xyz_st *R, 
    atom_info    *atom,    index_st *ist, 
    par_st       *par,     flag_st *flag, 
    parallel_st  *parallel
);

void build_gauss_hamiltonian(double *H_mat, double *T_mat, double *V_mat, index_st *ist, par_st *par);

void proj_gauss_on_grid(gauss_st *gauss, grid_st *grid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);

void calc_X_canonical(double *S_mat, double *X, double *U, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);

void transform_H(double *H_mat, double *X, double *eig_vals, MO_st *MO, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);

void proj_gauss_on_grid_MO(MO_st *MO, gauss_st *gauss, grid_st *grid, int start, int end, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
