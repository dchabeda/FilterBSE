#include "mod_dipole.h"

void mod_dipole(
    double*      psi_qp,
    double*      eig_vals,
    grid_st*     grid,
    xyz_st**     elec_dip,
    xyz_st**     mag_dip,
    double**     rot_strength,
    index_st*    ist,
    par_st*      par,
    flag_st*     flag,
    parallel_st*    parallel
    ){

    /************************************************************/
	/*******************  DECLARE VARIABLES   *******************/
	/************************************************************/

    const int          mpir = parallel->mpi_rank;

    unsigned long      n_el = ist->n_elecs;
    unsigned long      n_ho = ist->n_holes;
    
    /************************************************************/
	/*******************  DECLARE VARIABLES   *******************/
	/************************************************************/

    write_separation(stdout, "T");
    printf("\n4.\tCOMPUTING SINGLE-PARTICLE PROPERTIES | %s\n", get_time());
    write_separation(stdout, "B"); fflush(stdout);

    /************************************************************/
	/*****************   ALLOC MEM MTRX ELEMS   *****************/
	/************************************************************/

    if (mpir == 0) printf("Allocating memory for single particle matrix elements... "); fflush(stdout);
    
    ALLOCATE(elec_dip,     n_el * n_ho, "elec_dip");
    ALLOCATE(mag_dip,      n_el * n_ho, "mag_dip");
    ALLOCATE(rot_strength, n_el* n_ho,  "rot_strength");
        
    if (mpir == 0) printf("done\n"); fflush(stdout);

    /************************************************************/
	/****************  CALC TRANSITION DIPOLES   ****************/
	/************************************************************/

    if (mpir == 0) printf("\nElectric transition dipole moment...\n");
    calc_elec_dipole(*elec_dip, psi_qp, eig_vals, grid, ist, par, flag);

    if (mpir == 0) printf("\nMagnetic transition dipole moment...\n");
    calc_mag_dipole(*mag_dip, psi_qp, eig_vals, grid, ist, par, flag);

    // calc_rot_strength(rs, mux, muy, muz, mx, my, mz, eig_vals, &ist);

    return;
}