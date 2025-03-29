#include "fd.h"
#include "init.h"
#include "aux.h"

void mod_pot(
  double complex** pot_bare,
  double complex** pot_screened,
  grid_st*         grid,
  index_st*        ist,
  par_st*          par,
  flag_st*         flag,
  parallel_st*     parallel
);