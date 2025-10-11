#include "fd.h"
#include "init.h"

void calc_angular_exp(
    double *psitot,
    grid_st *grid,
    int start,
    int stop,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel);

void apply_J_op(
    zomplex *Jxpsi,
    zomplex *Jypsi,
    zomplex *Jzpsi,
    zomplex *psi,
    grid_st *grid,
    fftw_plan_loc planfw,
    fftw_plan_loc planbw,
    fftw_complex *fftwpsi,
    index_st *ist,
    par_st *par,
    parallel_st *parallel);

void apply_L_op(
    zomplex *Lxpsi,
    zomplex *Lypsi,
    zomplex *Lzpsi,
    zomplex *psi,
    grid_st *grid,
    fftw_plan_loc planfw,
    fftw_plan_loc planbw,
    fftw_complex *fftwpsi,
    index_st *ist,
    par_st *par,
    parallel_st *parallel);

void low_pass_filter(
    zomplex *psi,
    grid_st *grid,
    fftw_plan_loc planfw,
    fftw_plan_loc planbw,
    fftw_complex *fftwpsi,
    index_st *ist,
    par_st *par,
    parallel_st *parallel);
