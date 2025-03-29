#include "fd.h"
#include <float.h>
#include "aux.h"

void calc_eh_kernel_cplx(
	double        *psi_qp, 
	double        *pot_bare,
	double        *pot_screened,
	zomplex       *direct,
	zomplex       *exchange,
	index_st*     ist,
	par_st        *par,
	flag_st       *flag,
	parallel_st   *parallel
);