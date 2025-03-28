#include "fd.h"
#include "coulomb.h"

void mod_kernel(
  double complex*   psi_qp,
  double complex**  direct,
  double complex**  exchange,
  double complex*   pot_bare,
  double complex*   pot_screened,
  index_st*         ist,
  par_st*           par,
  flag_st*          flag,
  parallel_st*      parallel
);
