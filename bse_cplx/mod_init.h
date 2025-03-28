#include "fd.h"
#include "basis.h"
#include "aux.h"
#include "read.h"

void mod_init(
  double complex**  psitot,
  double complex**  psi_qp,
  double**          eig_vals,
  double**          sigma_E,
  xyz_st**          R,
  grid_st*          grid,
  double**          gridx,
  double**          gridy,
  double**          gridz,
  index_st*         ist,
  par_st*           par,
  flag_st*          flag,
  parallel_st*      parallel
);