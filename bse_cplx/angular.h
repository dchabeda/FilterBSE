#include "fd.h"

void calc_qp_spin_mtrx(
  double complex* restrict psi_qp,
	xyz_st *restrict         s_mom,  
	grid_st*                 grid, 
	index_st*                ist, 
	par_st*                  par
);

void calc_qp_ang_mom_mtrx(
  double complex* restrict psi_qp, 
  xyz_st* restrict         l_mom,
	double complex* restrict l2_mom,
	double complex* restrict ldots,
	grid_st*                 grid,
	index_st*                ist,
	par_st*                  par
);

void l_operator(
  double complex* Lxpsi, 
  double complex* Lypsi, 
  double complex* Lzpsi, 
  double complex* psi_qp, 
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
  double complex *psi, 
  double complex *Lpsi, 
  grid_st *grid, 
  index_st *ist, 
  par_st *par, 
  fftw_plan_loc planfw, 
  fftw_plan_loc planbw, 
  fftw_complex *fftwpsi
);
