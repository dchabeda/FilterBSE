#include "fd.h"
#include "ortho.h"
#include "Hmat.h"
#include "aux.h"
#include "mod_ortho.h"
#include "mod_sigma.h"
#include "mod_output.h"

void mod_portho(
    double**       psi_rank,
    double**       psitot,
    double*        pot_local,
    double*        eig_vals,
    double*        sigma_E,
    xyz_st*        R,
    zomplex*       LS,
    nlc_st*        nlc, 
    long*          nl, 
    double*        ksqr,
    grid_st*       grid,
    index_st*      ist,
    par_st*        par,
    flag_st*       flag,
    parallel_st*   parallel
);


void read_psi_from_disk(
    double*             psibuf, 
    unsigned long       start, 
    unsigned long       end, 
    unsigned long       stlen, 
    char*               filename
);

