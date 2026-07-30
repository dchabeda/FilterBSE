#include "fd.h"
#include "angular.h"

void mod_optional_output(
    double *psitot,
    grid_st *grid,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel);

void print_eigstate_densities(
    double *psitot,
    grid_st *grid,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel);

// Periodic path: write cube files of the frontier orbitals at one k-point,
// visualizing the full Bloch function e^{ik.r} u_nk(r) on an auto-sized supercell.
int bloch_supercell_factor(double frac);

void print_bloch_cubes_k(
    double *psitot_k,
    double *eig_vals,
    double *sigma_E,
    long n_states,
    int ik,
    vector k,
    grid_st *grid,
    index_st *ist,
    par_st *par,
    flag_st *flag);