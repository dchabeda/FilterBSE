#include "fd.h"
#include "hamiltonian.h"
#include "energy.h"
#include "aux.h"

void run_filter_cycle(
    double*       psi_rank, 
    double*       pot_local,
    zomplex*      LS,
    nlc_st*       nlc, 
    long*         nl, 
    double*       ksqr, 
    zomplex*      an, 
    double*       zn, 
    double*       ene_targets, 
    grid_st*      grid, 
    index_st*     ist, 
    par_st*       par, 
    flag_st*      flag, 
    parallel_st*  parallel);

void filter_cycle(
    double*       psi_rank,
    long          jns,
    zomplex*      psi,
    zomplex*      phi,
    double*       pot_local,
    zomplex*      projs,
    zomplex*      LS,
    nlc_st*       nlc, 
    long*         nl, 
    double*       ksqr, 
    zomplex*      an, 
    double*       zn, 
    double*       ene_targets, 
    grid_st*      grid, 
    index_st*     ist, 
    par_st*       par, 
    flag_st*      flag, 
    parallel_st*  parallel,
    fftw_plan_loc planfw,
    fftw_plan_loc planbw,
    fftw_complex* fftwpsi
);

void p_hnorm(
    zomplex*      psi_out, 
    zomplex*      psi_tmp, 
    double*       pot_local,
    zomplex*      projs, 
    zomplex*      LS,
    nlc_st*       nlc, 
    long*         nl, 
    double*       ksqr,
    double        zn,
    index_st*     ist, 
    par_st*       par, 
    flag_st*      flag, 
    parallel_st*  parallel,
    fftw_plan_loc planfw,
    fftw_plan_loc planbw,
    fftw_complex* fftwpsi
);

void time_hamiltonian(
    zomplex*      phi, 
    zomplex*      psi, 
    double*       pot_local, 
    zomplex*      LS,
    nlc_st*       nlc, 
    long*         nl, 
    double*       ksqr,
    index_st*     ist, 
    par_st*       par, 
    flag_st*      flag, 
    parallel_st*  parallel
);
