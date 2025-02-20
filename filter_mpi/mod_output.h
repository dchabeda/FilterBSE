#include "fd.h"

void mod_output(
  double*       psitot,
  xyz_st*       R,
  double*       eig_vals,
  double*       sigma_E,
  grid_st*      grid,
  index_st*     ist,
  par_st*       par,
  flag_st*      flag,
  parallel_st*  parallel
);

int select_conv_states(
  double*       sigma_E, 
  int*          arr,
  double        secut,
  int           n_elems
);

