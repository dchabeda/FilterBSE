#include "fd.h"

void gen_newton_coeff(
    zomplex *an,
    double *samp,
    double *ene_targets,
    index_st *ist,
    par_st *par,
    parallel_st *parallel);

void chebyshev_reordered(double *, double, double, long);

double samp_points_ashkenazy(zomplex *point, double min, double max, long ncheby);

void check_function(zomplex *an, zomplex *samp, index_st *ist, par_st *par, double ene_targets);

void read_newton_coeff(zomplex *an, double *samp, index_st *ist, par_st *par);
