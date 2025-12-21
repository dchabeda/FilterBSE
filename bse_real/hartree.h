#include "fd.h"

void hartree(
    double complex *rho,
    double complex *potq,
    double *poth,
    index_st *ist,
    fftw_plan_loc planfw,
    fftw_plan_loc planbw,
    fftw_complex *fftwpsi);
