#include "fd.h"
#include "hamiltonian.h"
#include "energy.h"
#include "aux.h"

void run_filter_cycles(
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

void gather_mpi_filt(
    double*       psi_rank,
    double**      psitot,
    index_st*     ist,
    par_st*       par,
    flag_st*      flag,
    parallel_st*  parallel
);


/*************************************************************/

void run_filter_cycles_k(
    double*       psi_rank, 
    double*       pot_local,
    vector*       G_vecs,
    vector*       k_vecs,
    grid_st*      grid,
    zomplex*      LS,
    nlc_st*       nlc, 
    long*         nl,
    zomplex*      an, 
    double*       zn, 
    double*       ene_targets, 
    index_st*     ist, 
    par_st*       par, 
    flag_st*      flag, 
    parallel_st*  parallel
);


void filter_cycle_k(
    double*       psi_rank,
    long          jns,
    zomplex*      psi,
    zomplex*      phi,
    double*       pot_local,
    vector*       G_vecs,
    vector        k,
    grid_st*      grid, 
    zomplex*      LS,
    nlc_st*       nlc, 
    long*         nl, 
    zomplex*      an, 
    double*       zn, 
    double*       ene_targets, 
    index_st*     ist, 
    par_st*       par, 
    flag_st*      flag, 
    parallel_st*  parallel,
    fftw_plan_loc planfw,
    fftw_plan_loc planbw,
    fftw_complex* fftwpsi
);

void p_hnorm_k(
    zomplex*      psi_out, 
    zomplex*      psi_tmp, 
    double*       pot_local,
    vector*       G_vecs,
    vector        k,
    grid_st*      grid,
    zomplex*      LS,
    nlc_st*       nlc, 
    long*         nl,
    double        zn,
    index_st*     ist, 
    par_st*       par, 
    flag_st*      flag, 
    parallel_st*  parallel,
    fftw_plan_loc planfw,
    fftw_plan_loc planbw,
    fftw_complex* fftwpsi
);

void time_hamiltonian_k(
    zomplex*      psi_out, 
    zomplex*      psi_tmp, 
    double*       pot_local, 
    vector*       G_vecs,
    vector        k,
    grid_st*      grid,
    zomplex*      LS,
    nlc_st*       nlc, 
    long*         nl, 
    index_st*     ist, 
    par_st*       par, 
    flag_st*      flag, 
    parallel_st*  parallel
);
