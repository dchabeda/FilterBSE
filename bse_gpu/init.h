#include "fd.h"
#include "aux.h"

void init_elec_hole_kernel(
    double*           pot_bare, 
    double*           pot_screened, 
    grid_st*          grid, 
    index_st*         ist, 
    par_st*           par, 
    flag_st*          flag, 
    parallel_st*      parallel, 
    fftw_plan_loc     planfw, 
    fftw_plan_loc     planbw, 
    fftw_complex*     fftwpsi
);

double calc_coulomb(double r, double gamma);