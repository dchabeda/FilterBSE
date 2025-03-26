#include "fd.h"
#include "hamiltonian.h"

void diag_H(
    double*        psitot, 
    double*        pot_local, 
    zomplex*       LS, 
    nlc_st*        nlc, 
    long*          nl, 
    double*        ksqr, 
    double*        eval, 
    index_st*      ist, 
    par_st*        par, 
    flag_st*       flag, 
    parallel_st*   parallel
);
