#include "fd.h"
#include "aux.h"
#include "angular.h"

/**************************************************************/

void calc_elec_dipole(
    xyz_st *elec_dip,
    double complex *psi_qp,
    double *eig_vals,
    grid_st *grid,
    index_st *ist,
    par_st *par,
    flag_st *flag);

/**************************************************************/

void calc_mag_dipole(
    xyz_st *mag_dip,
    double complex *psi_qp,
    double *eig_vals,
    grid_st *grid,
    index_st *ist,
    par_st *par,
    flag_st *flag);

/**************************************************************/

void calc_rotational_strength(
    double *rs,
    xyz_st *elec_dip,
    xyz_st *mag_dip,
    double *eval,
    index_st *ist);
