#include "fd.h"
#include "aux.h"

void init_elec_hole_kernel(
    double complex*   pot_bare, 
    double complex*   pot_screened, 
    grid_st*          grid, 
    index_st*         ist, 
    par_st*           par, 
    flag_st*          flag, 
    parallel_st*      parallel
);

double calc_coulomb(double r, double gamma);