#include "fd.h"
#include "aux.h"
#include "hartree.h"
#include "write.h"

void calc_eh_kernel_cplx(
	double complex* psi_qp, 
	double complex* pot_bare,
	double complex* pot_screened,
	double complex* direct,
	double complex* exchange,
	index_st*       ist,
	par_st*         par,
	flag_st*        flag,
	parallel_st*    parallel
);

int load_coulomb_mat(
  double complex* mat, 
  char*           fileName, 
  long*           a,
  long*           b,
  long*           i,
  long*           j,
  index_st*       ist
);
