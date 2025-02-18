#include "fd.h"

void hamiltonian(
    zomplex*       psi_out, 
    zomplex*       psi_tmp, 
    double*        pot_local, 
    zomplex*       LS, 
    nlc_st*        nlc, 
    long*          nl, 
    double*        ksqr,
    index_st*      ist, 
    par_st*        par, 
    flag_st*       flag, 
    fftw_plan_loc  planfw, 
    fftw_plan_loc  planbw, 
    fftw_complex*  fftwpsi);

void kinetic(
    zomplex*         psi_out, 
    double*          ksqr, 
    fftw_plan_loc    planfw, 
    fftw_plan_loc    planbw, 
    fftw_complex*    fftwpsi, 
    index_st*        ist
);

void potential(
    zomplex*        psi_out, 
    zomplex*        psi_tmp, 
    double*         pot_local, 
    zomplex*        LS, 
    nlc_st*         nlc, 
    long*           nl, 
    index_st*       ist,
    par_st*         par, 
    flag_st*        flag
);

void spin_orbit_proj_pot(
    zomplex*        psi_out, 
    zomplex*        psi_tmp,
    zomplex*        LS,
    nlc_st*         nlc, 
    long*           nl, 
    index_st*       ist, 
    par_st*         par
);

void nonlocal_proj_pot(
    zomplex*        psi_out, 
    zomplex*        psi_tmp, 
    nlc_st*         nlc, 
    long*           nl, 
    index_st*       ist, 
    par_st*         par
);

void time_reverse_all(
    double*         psitot, 
    double*         dest, 
    index_st*       ist, 
    parallel_st*    parallel
);
    
// OMP parallelized Hamiltonian

void p_hamiltonian(
    zomplex*       psi_out, 
    zomplex*       psi_tmp, 
    double*        pot_local, 
    zomplex*       LS, 
    nlc_st*        nlc, 
    long*          nl, 
    double*        ksqr,
    index_st*      ist, 
    par_st*        par, 
    flag_st*       flag, 
    fftw_plan_loc  planfw, 
    fftw_plan_loc  planbw, 
    fftw_complex*  fftwpsi,
    int            ham_threads
);

void p_potential(
    zomplex*        psi_out, 
    zomplex*        psi_tmp, 
    double*         pot_local, 
    zomplex*        LS, 
    nlc_st*         nlc, 
    long*           nl, 
    index_st*       ist,
    par_st*         par, 
    flag_st*        flag,
    int             ham_threads
);

void p_spin_orbit_proj_pot(
    zomplex*        psi_out, 
    zomplex*        psi_tmp,
    zomplex*        LS,
    nlc_st*         nlc, 
    long*           nl, 
    index_st*       ist, 
    par_st*         par,
    int             ham_threads
);

void p_nonlocal_proj_pot(
    zomplex*        psi_out, 
    zomplex*        psi_tmp, 
    nlc_st*         nlc, 
    long*           nl, 
    index_st*       ist, 
    par_st*         par,
    int             ham_threads
);

void def_LS(
    zomplex*      LS, 
    index_st*     ist, 
    par_st*       par);

