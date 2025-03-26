#include "fd.h"
#include "hamiltonian.h"

double energy(
    zomplex*        psi, 
    zomplex*        phi, 
    double*         pot_local, 
    zomplex*        LS, 
    nlc_st*         nlc, 
    long*           nl, 
    double*         ksqr, 
    index_st*       ist,
    par_st*         par, 
    flag_st*        flag, 
    fftw_plan_loc   planfw, 
    fftw_plan_loc   planbw, 
    fftw_complex*   fftwpsi
);

void energy_all(
    double*         psitot, 
    long            n_states, 
    double*         pot_local, 
    zomplex*        LS, 
    nlc_st*         nlc, 
    long*           nl, 
    double*         ksqr,
    double*         ene_filters, 
    index_st*       ist, 
    par_st*         par, 
    flag_st*        flag, 
    parallel_st*    parallel
);

void get_energy_range(
    zomplex*        psi, 
    zomplex*        phi, 
    double*         pot_local, 
    grid_st*        grid, 
    zomplex*        LS, 
    nlc_st*         nlc, 
    long*           nl, 
    double*         ksqr,
    index_st*       ist, 
    par_st*         par, 
    flag_st*        flag, 
    parallel_st*    parallel
);

void calc_sigma_E(
    double*         psitot, 
    double*         pot_local, 
    zomplex*        LS, 
    nlc_st*         nlc, 
    long*           nl, 
    double*         ksqr,
    double*         sigma_E, 
    index_st*       ist, 
    par_st*         par, 
    flag_st*        flag
);

