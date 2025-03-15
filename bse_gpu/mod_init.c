#include "mod_init.h"

void mod_init(
    double**        psitot,
    ){

    /************************************************************/
	/*******************  DECLARE VARIABLES   *******************/
	/************************************************************/


    /************************************************************/
	/********************     READ INPUT     ********************/
	/************************************************************/

    if (mpir == 0){
        write_separation(stdout, top);
        printf("\n1.\tINITIALIZING JOB | %s\n", get_time());
        write_separation(stdout, bottom); fflush(stdout);
    } 
    
        
    if (0 == init_unsafe){
        if (mpir == 0) printf("\nReading filter output from output.dat:\n"); fflush(stdout);
        read_filter_output("output.dat", &psitot, &eig_vals, &sigma_E, &R, &grid, &(grid.x), &(grid.y), &(grid.z), &ist, &par, &flag);

    } else if (1 == init_unsafe){
        if (mpir == 0) printf("\nReading UNSAFE input from unsafe_input.par:\n"); fflush(stdout);
        read_unsafe_input(&psitot, &eig_vals, &sigma_E, &R, &grid, &(grid.x), &(grid.y), &(grid.z), &ist, &par, &flag, &parallel);
    }

    /*************************************************************************/
    /*** Read initial setup from input.par ***/
    if (mpir == 0) printf("\nReading BSE job specifications from input.par:\n"); fflush(stdout);
    // ******
    // ******
    read_input(&flag, &grid, &ist, &par, &parallel);
    // ******
    // ******
    fflush(stdout);

    /*************************************************************************/
    /*** Determine the configuration of the quasiparticle basis ***/
    if (mpir == 0) printf("\nSetting quasiparticle basis indices:\n"); fflush(stdout);
    
    // Allocate memory for the lists of the indices of eigenstates
    // The maximum possible number of hole states is mn_states_tot from filter.
    // We allocate an entire block of that size for both elecs and holes 
    // because we will reallocate after the get_qp_basis_indices function.
    if((ist.eval_hole_idxs = (long*) malloc(ist.mn_states_tot * sizeof(ist.eval_hole_idxs[0]))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for ist.eval_hole_idxs in main.c\n");
        exit(EXIT_FAILURE);
    }
    if((ist.eval_elec_idxs = (long*) malloc(ist.mn_states_tot * sizeof(ist.eval_elec_idxs[0]))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for ist.eval_elec_idxs in main.c\n");
        exit(EXIT_FAILURE);
    }
    
    // ******
    // ******
    get_qp_basis_indices(eig_vals, sigma_E, &ist.eval_hole_idxs, &ist.eval_elec_idxs, &ist, &par, &flag, &parallel);
    // ******
    // ******
    
    // Reallocate the eig_vals and sigma_E arrays to only contain the n_qp states
    // eig_vals = realloc(eig_vals, ist.n_qp * sizeof(eig_vals[0]));
    // sigma_E = realloc(sigma_E, ist.n_qp * sizeof(sigma_E[0]));
    
    // Resize the index arrays to tightly contain the indices of eigenstates in VB/CB
    // The conditional checks whether n_holes saturated the memory block (if so, no need to resize)
    // or if they are not equal, then n_holes is smaller and we should shrink the array
    // if (ist.n_holes != ist.max_hole_states){
    //     ist.eval_hole_idxs = realloc(ist.eval_hole_idxs, ist.n_holes * sizeof(long));
    // }
    // if (ist.n_elecs != ist.max_elec_states){
    //     ist.eval_elec_idxs = realloc(ist.eval_elec_idxs, ist.n_elecs * sizeof(long));
    // }
    // fprintf(pmem, "alloc ist.eval_hole_idxs %ld B\n", ist.mn_states_tot * sizeof(ist.eval_hole_idxs[0])); mem += ist.mn_states_tot * sizeof(ist.eval_hole_idxs[0]);
    // fprintf(pmem, "alloc ist.eval_elec_idxs %ld B\n", ist.mn_states_tot * sizeof(ist.eval_elec_idxs[0])); mem += ist.mn_states_tot * sizeof(ist.eval_elec_idxs[0]);
    

    /*************************************************************************/
    // Allocate memory for the quasiparticle wavefunctions
    
    if ((psi_qp = (double *) malloc( ist.complex_idx * ist.nspinngrid * (ist.n_qp) * sizeof(double))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for psi_qp in main.c\n");
        exit(EXIT_FAILURE);
    }
    
    fprintf(pmem, "alloc psi_qp %ld B", ist.complex_idx * ist.nspinngrid * (ist.n_holes + ist.n_elecs) * sizeof(double)); mem += ist.complex_idx * ist.nspinngrid * (ist.n_holes + ist.n_elecs) * sizeof(double);
    if (mpir == 0) write_separation(pmem, top);
    fprintf(pmem, "\ntotal mem usage %ld MB\n", mem / 1000000 );
    if (mpir == 0) write_separation(pmem, bottom); fflush(pmem);

    // if (mpir == 0){
    //     printf("\nTHE HOLES:\n");
    //     for (i = 0; i < ist.n_holes; i++){
    //         printf("%ld\n", ist.eval_hole_idxs[i]);
    //     }
        
    //     printf("\nTHE ELECS:\n");
    //     for (i = 0; i < ist.n_elecs; i++){
    //         printf("%ld\n", ist.eval_elec_idxs[i]);
    //     }
    // }
    
    // ******
    // ******
    if (mpir == 0) printf("\nReading quasiparticle basis wavefunctions:\n"); fflush(stdout);
    get_qp_basis(psi_qp, psitot, eig_vals, sigma_E, &ist, &par, &flag);
    ist.n_xton = ist.n_elecs * ist.n_holes;
    // ******
    // ******
    free(psitot); free(sigma_E); //free(psi_hole); free(psi_elec);
    fprintf(pmem, "free psitot %ld B\n", ist.complex_idx * ist.nspinngrid * ist.mn_states_tot * sizeof(double)); mem -= ist.complex_idx * ist.nspinngrid * ist.mn_states_tot * sizeof(double);
    // fprintf(pmem, "free psi_hole %ld B\n", ist.complex_idx * ist.nspinngrid * ist.n_holes * sizeof(psi_hole[0])); mem -= ist.complex_idx*ist.nspinngrid*ist.n_holes*sizeof(psi_hole[0]);
    // fprintf(pmem, "free psi_elec %ld B\n", ist.complex_idx * ist.nspinngrid * ist.n_elecs * sizeof(psi_elec[0])); mem -= ist.complex_idx*ist.nspinngrid*ist.n_elecs*sizeof(psi_elec[0]);
    fprintf(pmem, "free sigma_E %ld B", ist.mn_states_tot * sizeof(double)); mem -= ist.mn_states_tot * sizeof(double);
    if (mpir == 0) write_separation(pmem, top);
    fprintf(pmem, "\ntotal mem usage %ld MB\n", mem / 1000000 );
    if (mpir == 0) write_separation(pmem, bottom); fflush(pmem);

    return;
}