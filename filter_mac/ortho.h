#include "fd.h"
#include "hamiltonian.h"


long ortho_real(
  double *psitot, 
  double dv, 
  index_st *ist, 
  par_st *par, 
  flag_st *flag, 
  parallel_st *parallel
);

long ortho_cplx(
  zomplex *psitot, 
  double dv, 
  index_st *ist, 
  par_st *par, 
  flag_st *flag, 
  parallel_st *parallel
);

void restart_from_ortho(
    double*       psitot,
    index_st*     ist,
    par_st*       par,
    flag_st*      flag,
    parallel_st*  parallel
);
