#include "fd.h"
#include "hamiltonian.h"

void hamiltonian_k(
    zomplex       *psi_out,  zomplex      *psi_tmp,    double *pot_local,    
    vector        *G_vecs,   vector        k,          grid_st *grid, zomplex *LS, 
    nlc_st        *nlc,      long         *nl,         index_st *ist, 
    par_st        *par,      flag_st      *flag,       fftw_plan_loc planfw, 
    fftw_plan_loc planbw,    fftw_complex *fftwpsi
);

void p_hamiltonian_k(zomplex *psi_out, zomplex *psi_tmp, double *pot_local, vector *G_vecs, vector k, grid_st *grid, zomplex* LS, nlc_st *nlc, long *nl,
    index_st *ist, par_st *par, flag_st *flag, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi, int ham_threads
);

void kinetic_k(zomplex *psi_out, vector *G_vecs, vector k, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi, index_st *ist
);

void e_ikr(zomplex *psi, vector k, grid_st *grid, index_st *ist, par_st *par, flag_st *flag
);