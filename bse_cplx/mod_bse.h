#include "fd.h"
#include "bethe-salpeter.h"

void mod_bse(
  double complex*  psi_qp,
  double complex*  direct, 
  double complex*  exchange,
  double complex** bsmat, 
  double complex** bs_coeff, 
  double**         h0mat, 
  double**         xton_ene,
  double*          eig_vals,
	xyz_st*          s_mom, 
  xyz_st*          l_mom, 
  double complex*  l2_mom, 
  double complex*  ldots, 
  grid_st*         grid, 
  index_st*        ist, 
  par_st*          par,
  flag_st*         flag,
  parallel_st*     parallel
);

void build_h0_mat(
  double *h0mat, 
  double *eval, 
  index_st* ist
);

void build_BSE_mat(
  double complex* bsmat, 
  double complex* direct, 
  double complex* exchange, 
  index_st*       ist
);
