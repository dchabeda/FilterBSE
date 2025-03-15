#include "mod_init.h"

void mod_init(
    double**        psitot,
    double**        psi_qp,
    double**        eig_vals,
    double**        sigma_E,
    xyz_st**        R,
    grid_st*        grid,
    double**        gridx,
    double**        gridy,
    double**        gridz,
    index_st*       ist,
    par_st*         par,
    flag_st*        flag,
    parallel_st*    parallel
    ){

    /************************************************************/
	/*******************  DECLARE VARIABLES   *******************/
	/************************************************************/

    long             stlen;
    const int        mpir = parallel->mpi_rank;

    /************************************************************/
	/********************     READ INPUT     ********************/
	/************************************************************/

    if (mpir == 0){
        write_separation(stdout, "T");
        printf("\n1.\tINITIALIZING JOB | %s\n", get_time());
        write_separation(stdout, "B"); fflush(stdout);
    } 
    
        
    if (0 == flag->initUnsafe){
        if (mpir == 0) printf("\nReading filter output from output.dat:\n"); fflush(stdout);
        read_filter_output("output.dat", psitot, eig_vals, sigma_E, R, grid, gridx, gridy, gridz, ist, par, flag);

    } else if (1 == flag->initUnsafe){
        if (mpir == 0) printf("\nReading UNSAFE input from unsafe_input.par:\n"); fflush(stdout);
        read_unsafe_input(
            psitot, eig_vals, sigma_E, R, grid, gridx, gridy, gridz, 
            ist, par, flag, parallel
        );
    }

    
    /*** Read initial setup from input.par ***/
    if (mpir == 0) printf("\nReading BSE job specifications from input.par:\n"); fflush(stdout);
    
    read_input(flag, grid, ist, par, parallel);
    
    stlen = ist->complex_idx * ist->nspinngrid;

    /************************************************************/
	/*******************   GET QP BASIS IDXS   ******************/
	/************************************************************/

    if (mpir == 0) printf("\nSetting quasiparticle basis indices:\n"); fflush(stdout);
    
    // Allocate memory for the lists of the indices of eigenstates
    // The maximum possible number of hole states is mn_states_tot from filter.
    // We allocate an entire block of that size for both elecs and holes 
    // because we will reallocate after the get_qp_basis_indices function.
    
    ALLOCATE(&(ist->eval_hole_idxs), ist->mn_states_tot, "eval_hole_idxs");
    ALLOCATE(&(ist->eval_elec_idxs), ist->mn_states_tot, "eval_elec_idxs");
    
    get_qp_basis_indices(*eig_vals, *sigma_E, &ist->eval_hole_idxs, &ist->eval_elec_idxs, ist, par, flag, parallel);
    
    
    /************************************************************/
	/********************    BUILD QP BASIS    ******************/
	/************************************************************/

    if (mpir == 0) printf("\nReading quasiparticle basis wavefunctions:\n"); fflush(stdout);
    printf("n_qp = %ld stlen = %ld\n", ist->n_qp, stlen);
    ALLOCATE(psi_qp, ist->n_qp * stlen, "psi_qp");
    
    get_qp_basis(*psi_qp, *psitot, *eig_vals, *sigma_E, ist, par, flag);
    
    ist->n_xton = ist->n_elecs * ist->n_holes;
    
    free(*psitot); free(*sigma_E);
    
    return;
}