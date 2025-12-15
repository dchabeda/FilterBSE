#include "fd.h"
#include "aux.h"
#include "dipole.h"

void mod_dipole(
    double *psi_qp,
    double *eig_vals,
    grid_st *grid,
    xyz_st **elec_dip,
    xyz_st **mag_dip,
    double **rot_strength,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel);