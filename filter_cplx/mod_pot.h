#include "fd.h"
#include "init.h"
#include "hamiltonian.h"

void mod_pot(
    double*       pot_local,
    pot_st*       pot,
    xyz_st*       R,
    atom_info*    atom,
    grid_st*      grid,
    zomplex*      LS,
    nlc_st*       nlc,
    long*         nl,
    double*       SO_projectors,
    index_st*     ist,
    par_st*       par,
    flag_st*      flag,
    parallel_st*  parallel);



