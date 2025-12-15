#include "fd.h"

void bethe_salpeter(
  double complex*  direct, 
  double complex*  exchange,
  double complex*  bsmat, 
  double complex*  bs_coeff, 
  double*          h0mat, 
  double*          xton_ene,
  grid_st*         grid, 
  index_st*        ist, 
  par_st*          par,
  flag_st*         flag,
  parallel_st*     parallel
);

void diag(
  const int n, 
  int nthreads, 
  double complex *mat, 
  double *eval
);
