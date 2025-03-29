#include "fd.h"
#include "aux.h"

void qp_spin_frac(
  double complex*  psi_qp,
  double*          eig_vals,
  grid_st*         grid,
  index_st*        ist,
  par_st*          par,
  flag_st*         flag,
  parallel_st*     parallel
);