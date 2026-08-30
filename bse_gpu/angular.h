#include "fd.h"
#include "aux.h"

void calc_spin_mtrx(
	xyz_st*      s_mom, 
	double*      psi_qp, 
	grid_st*     grid, 
	index_st*    ist, 
	par_st*      par
);

void calc_ang_mom_mtrx(
	xyz_st*        l_mom, 
	zomplex*       l2_mom,
	zomplex*       LdotS,
	double*        psi_qp,
	grid_st*       grid,
	index_st*      ist,
	par_st*        par
);

void l_operator(
    zomplex* Lxpsi, 
    zomplex* Lypsi, 
    zomplex* Lzpsi, 
    zomplex* psi_qp, 
    double* g_vecs,
	grid_st *grid, 
    index_st *ist, 
    par_st *par, 
    fftw_plan_loc planfw, 
    fftw_plan_loc planbw, 
    fftw_complex* fftwpsi
);

void init_g_vecs(
    double *kindex, 
    double *kx, 
    double *ky, 
    double *kz, 
    grid_st *grid, 
    index_st *ist, 
    par_st *par
);

void p_operator(
    char* direc, 
    double *kindex, 
    zomplex *psi, 
    zomplex *Lpsi, 
    grid_st *grid, 
    index_st *ist, 
    par_st *par, 
    fftw_plan_loc planfw, 
    fftw_plan_loc planbw, 
    fftw_complex *fftwpsi
);
