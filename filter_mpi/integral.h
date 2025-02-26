#include "fd.h"
#include <math.h> // For pow() and exp()

void overlap_gauss(double *S_mat, gauss_st *gauss, atom_info *atom, index_st *ist, par_st *par, flag_st *flag);

double S(
    xyz_st R_a, xyz_st R_b, 
    double a, double b, 
    int ix, int iy, int iz,
    int jx, int jy, int jz
);

void kinetic_gauss(double *T_mat, gauss_st *gauss, atom_info *atom, index_st *ist, par_st *par, flag_st *flag);

double T(
    xyz_st R_a, xyz_st R_b, 
    double a, double b, 
    int ix, int iy, int iz,
    int jx, int jy, int jz
);

void potential_gauss(double *V_mat, double *pot_local, gauss_st *gauss, grid_st *grid, atom_info *atom, index_st *ist, par_st *par, flag_st *flag);

double E_hermite_gauss(int i, int j, int t, double Qx, double a, double b);

void get_gauss_polarization(int polz, int* ix, int *iy, int *iz);