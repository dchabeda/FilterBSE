/*****************************************************************************/

#include "pot_coupling.h"

/*****************************************************************************/

int main(int argc, char *argv[]){
    /*****************************************************************
    * This is the main function for pot_coupling.x. It is the driver *
    * for computing integrals of <psi_i|deltaU(r;R)|psi_j> used in   *
    * the derivation of dynamics expressions (see Dipti Jasrasaria)  *
    * thesis                                                         *
    ******************************************************************/ 

    // DECLARE VARIABLES AND STRUCTS
    // file pointers
    FILE *pf, *pmem;
    pmem = fopen("mem.dat", "w");
    // zomplex types
    zomplex *psi, *phi;
    // custom structs 
    flag_st flag; index_st ist; par_st par; 
    atom_info *atom; atom_info *atom_equil;
    pot_st pot;
    grid_st grid;
    xyz_st *R_equil; xyz_st *R; 
    nlc_st *nlc_equil = NULL; nlc_st *nlc = NULL;  
    parallel_st parallel; 
    // double arrays
    double *pot_local_equil, *pot_local;
    double *SO_projectors_equil; double *SO_projectors; 
    double *gridx = NULL, *gridy = NULL, *gridz = NULL;
    double *psitot = NULL, *psi_hole = NULL, *psi_elec = NULL, *psi_qp;
    double *eig_vals = NULL, *sigma_E = NULL;
    double *gridx = NULL, *gridy = NULL, *gridz = NULL;
    double *rho;
    // long int arrays and counters
    long *nl_equil = NULL; long *nl = NULL;
    long i, j, state;
    long jgrid, jgrid_real, jgrid_imag;
    long mem = 0;
    ist.atom_types = malloc(N_MAX_ATOM_TYPES*sizeof(ist.atom_types[0]));
    fprintf(pmem, "alloc ist.atom_types %ld B\n", N_MAX_ATOM_TYPES*sizeof(ist.atom_types[0])); mem += N_MAX_ATOM_TYPES*sizeof(ist.atom_types[0]);
    // Clock/Wall time output and stdout formatting
    time_t start_time = time(NULL); // Get the actual time for total wall runtime
    time_t start_clock = clock(); // Get the starting CPU clock time for total CPU runtime
    char *top; top = malloc(2*sizeof(top[0])); 
    fprintf(pmem, "alloc top %ld B\n", 2*sizeof(top[0])); mem += 2*sizeof(top[0]);
    char *bottom; bottom = malloc(2*sizeof(bottom[0]));
    fprintf(pmem, "alloc bottom %ld B\n", 2*sizeof(bottom[0])); mem += 2*sizeof(bottom[0]);
    strcpy(top, "T\0"); strcpy(bottom, "B\0");


    /*************************************************************************/
    fprintf(stdout, "******************************************************************************\n");
    printf("\nRUNNING PROGRAM: POT COUPLING INTEGRALS\n");
    printf("This calculation began at: %s", ctime(&start_time)); 
    write_separation(stdout, bottom);
    fflush(stdout);

    /*************************************************************************/
    /*************************************************************************/
    // READ IN EQUILIBRIUM WAVEFUNCTIONS
    /*************************************************************************/
    /*************************************************************************/
    
    /*************************************************************************/
    // 1. Initialize job from input file
    
    write_separation(stdout, top);
    printf("\n1.\tREADING EQUILIBRIUM EIGENFUNCTIONS\n");
    write_separation(stdout, bottom); fflush(stdout);

    /*************************************************************************/
    /*** Read output from filter output.par ***/
    printf("\nReading filter output from output.dat:\n"); fflush(stdout);
    // ******
    // ******
    read_filter_output("output.dat", &psitot, &eig_vals, &sigma_E, &R, &grid_par, &gridx, &gridy, &gridz, &ist, &par, &flag);
    // ******
    // ******

    // Move the grid values into the grid struct arrays.
    // The memory allocation got too thorny, so we do this simple transfer
    
    if ((grid.x = malloc( grid.nx * sizeof(grid.x[0]))) == NULL ){
        fprintf(stderr, "ERROR: allocating memory for grid.x in main\n");
        exit(EXIT_FAILURE);
    }
    if ((grid.y = malloc(grid.ny * sizeof(grid.y[0]))) == NULL ){
            fprintf(stderr, "ERROR: allocating memory for grid.y in main\n");
            exit(EXIT_FAILURE);
        }
    if ((grid.z = malloc( grid.nz * sizeof(grid.z[0]))) == NULL ){
        fprintf(stderr, "ERROR: allocating memory for grid.z in main\n");
        exit(EXIT_FAILURE);
    }
    fprintf(pmem, "alloc grid.x %ld B\n", grid.nx * sizeof(grid.x[0])); mem += grid.nx * sizeof(grid.x[0]);
    fprintf(pmem, "alloc grid.y %ld B\n", grid.ny * sizeof(grid.y[0])); mem += grid.ny * sizeof(grid.y[0]);
    fprintf(pmem, "alloc grid.z %ld B\n", grid.nz * sizeof(grid.z[0])); mem += grid.nz * sizeof(grid.z[0]);

    for (i = 0; i < grid.nx; i++){
        grid.x[i] = gridx[i];
    }
    for (i = 0; i < grid.ny; i++){
        grid.y[i] = gridy[i];
    }
    for (i = 0; i < grid.nz; i++){
        grid.z[i] = gridz[i];
    }

    /*************************************************************************/
    /*** Read initial setup from input.par ***/
    printf("\nReading BSE job specifications from input.par:\n"); fflush(stdout);
    // ******
    // ******
    read_input(&flag, &grid, &ist, &par, &parallel);
    // ******
    // ******
    fflush(stdout);

    /*************************************************************************/
    /*** Determine the configuration of the quasiparticle basis ***/
    printf("\nSetting quasiparticle basis indices:\n"); fflush(stdout);
    
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
    get_qp_basis_indices(eig_vals, sigma_E, &ist.eval_hole_idxs, &ist.eval_elec_idxs, &ist, &par, &flag);
    // ******
    // ******
    
    // Reallocate the eig_vals and sigma_E arrays to only contain the n_qp states
    eig_vals = realloc(eig_vals, ist.n_qp * sizeof(eig_vals[0]));
    sigma_E = realloc(sigma_E, ist.n_qp * sizeof(sigma_E[0]));
    
    // Resize the index arrays to tightly contain the indices of eigenstates in VB/CB
    // The conditional checks whether n_holes saturated the memory block (if so, no need to resize)
    // or if they are not equal, then n_holes is smaller and we should shrink the array
    if (ist.n_holes != ist.max_hole_states){
        ist.eval_hole_idxs = realloc(ist.eval_hole_idxs, ist.n_holes * sizeof(long));
    }
    if (ist.n_elecs != ist.max_elec_states){
        ist.eval_elec_idxs = realloc(ist.eval_elec_idxs, ist.n_elecs * sizeof(long));
    }
    fprintf(pmem, "alloc ist.eval_hole_idxs %ld B\n", ist.mn_states_tot * sizeof(ist.eval_hole_idxs[0])); mem += ist.mn_states_tot * sizeof(ist.eval_hole_idxs[0]);
    fprintf(pmem, "alloc ist.eval_elec_idxs %ld B\n", ist.mn_states_tot * sizeof(ist.eval_elec_idxs[0])); mem += ist.mn_states_tot * sizeof(ist.eval_elec_idxs[0]);
    

    /*************************************************************************/
    // Allocate memory for the electron and hole wavefunctions
    if ((psi_hole = (double *) malloc( ist.complex_idx * ist.nspinngrid * ist.n_holes * sizeof(psi_hole[0]))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for psi_hole in main.c\n");
        exit(EXIT_FAILURE);
    }
    if ((psi_elec = (double *) malloc( ist.complex_idx * ist.nspinngrid * ist.n_elecs * sizeof(psi_elec[0]))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for psi_elec in main.c\n");
        exit(EXIT_FAILURE);
    }
    if ((psi_qp = (double *) malloc( ist.complex_idx * ist.nspinngrid * (ist.n_holes + ist.n_elecs) * sizeof(double))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for psi_elec in main.c\n");
        exit(EXIT_FAILURE);
    }
    fprintf(pmem, "alloc psi_hole %ld B\n", ist.complex_idx * ist.nspinngrid * ist.n_holes * sizeof(psi_hole[0])); mem += ist.complex_idx*ist.nspinngrid*ist.n_holes*sizeof(psi_hole[0]);
    fprintf(pmem, "alloc psi_elec %ld B\n", ist.complex_idx * ist.nspinngrid * ist.n_elecs * sizeof(psi_elec[0])); mem += ist.complex_idx*ist.nspinngrid*ist.n_elecs*sizeof(psi_elec[0]);
    fprintf(pmem, "alloc psi_qp %ld B", ist.complex_idx * ist.nspinngrid * (ist.n_holes + ist.n_elecs) * sizeof(double)); mem += ist.complex_idx * ist.nspinngrid * (ist.n_holes + ist.n_elecs) * sizeof(double);
    write_separation(pmem, top);
    fprintf(pmem, "\ntotal mem usage %ld MB\n", mem / 1000000 );
    write_separation(pmem, bottom); fflush(pmem);

    printf("\nTHE HOLES:\n");
    for (i = 0; i < ist.n_holes; i++){
        printf("%ld\n", ist.eval_hole_idxs[i]);
    }
    printf("\nTHE ELECS:\n");
    for (i = 0; i < ist.n_elecs; i++){
        printf("%ld\n", ist.eval_elec_idxs[i]);
    }
    
    // ******
    // ******
    get_qp_basis(psi_qp, psitot, psi_hole, psi_elec, eig_vals, sigma_E, &ist, &par, &flag);
    // ******
    // ******
    free(psitot); free(psi_hole); free(psi_elec); free(sigma_E);
    fprintf(pmem, "free psitot %ld B\n", ist.complex_idx * ist.nspinngrid * ist.mn_states_tot * sizeof(double)); mem -= ist.complex_idx * ist.nspinngrid * ist.mn_states_tot * sizeof(double);
    fprintf(pmem, "free psi_hole %ld B\n", ist.complex_idx * ist.nspinngrid * ist.n_holes * sizeof(psi_hole[0])); mem -= ist.complex_idx*ist.nspinngrid*ist.n_holes*sizeof(psi_hole[0]);
    fprintf(pmem, "free psi_elec %ld B\n", ist.complex_idx * ist.nspinngrid * ist.n_elecs * sizeof(psi_elec[0])); mem -= ist.complex_idx*ist.nspinngrid*ist.n_elecs*sizeof(psi_elec[0]);
    fprintf(pmem, "free sigma_E %ld B", ist.mn_states_tot * sizeof(double)); mem -= ist.mn_states_tot * sizeof(double);
    write_separation(pmem, top);
    fprintf(pmem, "\ntotal mem usage %ld MB\n", mem / 1000000 );
    write_separation(pmem, bottom); fflush(pmem);

    /*************************************************************************/
    /*****************************************************************************/
    // READ IN EQUILIBRIUM AND DISTORTED GEOMETRIES
    /*************************************************************************/
    /*****************************************************************************/

    /*************************************************************************/
    // Initialize job from input file
    
    write_separation(stdout, top);
    printf("\n2.\tREADING EQUILIBRIUM & DISTORTED GEOMETRIES\n");
    write_separation(stdout, bottom); fflush(stdout);

    /*** read initial setup from input.par ***/
    /*** allocating memory ***/
    // the positions of the atoms in the x, y, and z directions 
    if ((R_equil = (xyz_st *) calloc(ist.natoms, sizeof(xyz_st))) == NULL) {
        fprintf(stderr, "\nOUT OF MEMORY: R array\n\n"); exit(EXIT_FAILURE);
    }
    // the atom specific information 
    if ((atom_equil = (atom_info *) calloc(ist.natoms, sizeof(atom_info))) == NULL){
        fprintf(stderr, "\nOUT OF MEMORY: atom struct\n\n"); exit(EXIT_FAILURE);
    }
    
    // the distorted positions of the atoms in the x, y, and z directions 
    if ((R = (xyz_st *) calloc(ist.natoms, sizeof(xyz_st))) == NULL) {
        fprintf(stderr, "\nOUT OF MEMORY: R array\n\n"); exit(EXIT_FAILURE);
    }
    // the atom specific information 
    if ((atom = (atom_info *) calloc(ist.natoms, sizeof(atom_info))) == NULL){
        fprintf(stderr, "\nOUT OF MEMORY: atom struct\n\n"); exit(EXIT_FAILURE);
    }
    
    /*** create strings to differentiate files ***/
    char str[100];
    char *file_name_equil; file_name_equil = malloc(16*sizeof(file_name_equil[0]));
    strcpy(file_name_equil, "conf_equil.par");
    char *file_name; file_name = malloc(9*sizeof(file_name[0]));
    strcpy(file_name, "conf.par");

    /*** read the equilibrium configuration ***/
    sprintf(str, "\nReading atomic configuration from %s:\n", file_name_equil);
    printf("%s", str);
    read_conf(file_name_equil, R_equil, atom_equil, &ist, &par, &flag);

    /*** read the distorted configuration ***/
    sprintf(str, "\nReading atomic configuration from %s:\n", file_name_equil);
    printf("%s", str);
    read_conf(file_name, R, atom, &ist, &par, &flag);

    /*************************************************************************/
    /*** allocating memory for the rest of the program ***/
    
    // For reading the atomic potentials ***/
    pot.dr = (double *) calloc(ist.ngeoms * ist.n_atom_types, sizeof(double));
    pot.r = (double *) calloc(ist.ngeoms * ist.max_pot_file_len * ist.n_atom_types, sizeof(double));
    if (1.0 != par.scale_surface_Cs){
        // allocate memory for separately read LR potentials if the surface atoms will be charge balanced
        pot.r_LR = (double *) calloc(ist.ngeoms * ist.max_pot_file_len * ist.n_atom_types, sizeof(double));
    }
    pot.pseudo = (double *) calloc(ist.ngeoms * ist.max_pot_file_len * ist.n_atom_types, sizeof(double));
    if (1.0 != par.scale_surface_Cs){
        pot.pseudo_LR = (double *) calloc(ist.ngeoms * ist.max_pot_file_len * ist.n_atom_types, sizeof(double));
    }
    pot.file_lens = (long *) calloc(ist.ngeoms * ist.n_atom_types, sizeof(long));
    
    // Wavefunction-type objects
    if ((psi = (zomplex *)calloc(ist.nspinngrid, sizeof(zomplex))) == NULL){
        fprintf(stderr, "\nOUT OF MEMORY: psi\n\n"); exit(EXIT_FAILURE);
    }
    if ((phi = (zomplex *)calloc(ist.nspinngrid, sizeof(zomplex))) == NULL){
        fprintf(stderr, "\nOUT OF MEMORY: phi\n\n"); exit(EXIT_FAILURE);
    }
    if ((pot_local_equil = (double *) calloc(ist.ngrid, sizeof(double))) == NULL){
        fprintf(stderr, "\nOUT OF MEMORY: pot_local_equil\n\n"); exit(EXIT_FAILURE);
    }
    if ((pot_local = (double *) calloc(ist.ngrid, sizeof(double))) == NULL){
        fprintf(stderr, "\nOUT OF MEMORY: pot_local\n\n"); exit(EXIT_FAILURE);
    }
    
    // memory allocation for the spin-orbit potential 
    if ( (flag.SO == 1) || (flag.NL == 1) ){
        if ((SO_projectors_equil = (double*) calloc(PROJ_LEN * ist.nproj, sizeof(double)))==NULL){nerror("mem_SO_projector");}  
        if ((SO_projectors = (double*) calloc(PROJ_LEN * ist.nproj, sizeof(double)))==NULL){nerror("mem_SO_projector");}  
    }
    // memory allocation for the non-local potential
    if (1 == flag.NL){
        if ((nlc_equil = (nlc_st *) calloc(ist.n_NL_atoms*ist.n_NL_gridpts, sizeof(nlc_st))) == NULL){ 
        fprintf(stderr, "\nOUT OF MEMORY: nlc\n\n"); exit(EXIT_FAILURE);
        }
        if ((nl_equil = (long *) calloc(ist.natoms, sizeof(nl[0]))) == NULL) {
        fprintf(stderr, "\nOUT OF MEMORY: nl\n\n"); exit(EXIT_FAILURE);
        }

        if ((nlc = (nlc_st *) calloc(ist.n_NL_atoms*ist.n_NL_gridpts, sizeof(nlc_st))) == NULL){ 
        fprintf(stderr, "\nOUT OF MEMORY: nlc\n\n"); exit(EXIT_FAILURE);
        }
        if ((nl = (long *) calloc(ist.natoms, sizeof(nl[0]))) == NULL) {
        fprintf(stderr, "\nOUT OF MEMORY: nl\n\n"); exit(EXIT_FAILURE);
        }
    }
    
    printf("\tdone allocating memory.\n"); fflush(stdout);

    /*************************************************************************/
    /*****************************************************************************/
    // CALCULATE POTENTIALS IN EQUILIBRIUM AND DISTORTED GEOMETRIES
    /*************************************************************************/
    /*****************************************************************************/

    printf("\nInitializing potentials...\n");
    
    /**************************************************************************/
    
    printf("\nEquilibrium local pseudopotential:\n");
    build_local_pot(pot_local_equil, &pot, R_equil, ksqr, atom_equil, grid, &grid_par, &ist, &par, &flag, &parallel);
    
    free(pot.r); pot.r = NULL; 
    free(pot.pseudo); pot.pseudo = NULL; 
    free(pot.dr); pot.dr = NULL; 
    free(pot.file_lens); pot.file_lens = NULL;
    if (1.0 != par.scale_surface_Cs){
    free(pot.r_LR); pot.r_LR = NULL;
    free(pot.pseudo_LR); pot.pseudo_LR = NULL; 
    }
    
    write_cube_file(pot_local_equil, &grid_par, "local_pot_equil.cube");
    
    if(flag.SO==1) {
    printf("\nEquilibrium spin-orbit pseudopotential:\n");
    init_SO_projectors(SO_projectors_equil, R_equil, atom_equil, &grid_par, &ist, &par);
    }
    /*** initialization for the non-local potential ***/
    if (flag.NL == 1){
    printf("\nEquilibrium non-local pseudopotential:\n"); fflush(0);
    init_NL_projectors(nlc_equil, nl_equil, SO_projectors_equil, R_equil, atom_equil, grid, &grid_par, &ist, &par, &flag);
    }
    // free memory allocated to SO_projectors
    if ( (flag.SO == 1) || (flag.NL == 1) ){
    free(SO_projectors_equil); SO_projectors_equil = NULL;
    }

    /**************************************************************************/
    
    // Allocate memory for the distorted potentials
    pot.dr = (double *) calloc(ist.ngeoms * ist.n_atom_types, sizeof(double));
    pot.r = (double *) calloc(ist.ngeoms * ist.max_pot_file_len * ist.n_atom_types, sizeof(double));
    if (1.0 != par.scale_surface_Cs){
        // allocate memory for separately read LR potentials if the surface atoms will be charge balanced
        pot.r_LR = (double *) calloc(ist.ngeoms * ist.max_pot_file_len * ist.n_atom_types, sizeof(double));
    }
    pot.pseudo = (double *) calloc(ist.ngeoms * ist.max_pot_file_len * ist.n_atom_types, sizeof(double));
    if (1.0 != par.scale_surface_Cs){
        pot.pseudo_LR = (double *) calloc(ist.ngeoms * ist.max_pot_file_len * ist.n_atom_types, sizeof(double));
    }
    pot.file_lens = (long *) calloc(ist.ngeoms * ist.n_atom_types, sizeof(long));
    
    // Calculate the distorted local potential
    build_local_pot(pot_local, &pot, R, atom, grid, &grid_par, &ist, &par, &flag, &parallel);
    
    write_cube_file(pot_local, &grid_par, "local_pot_distorted.cube");
    
    free(pot.r); pot.r = NULL; 
    free(pot.pseudo); pot.pseudo = NULL; 
    free(pot.dr); pot.dr = NULL; 
    free(pot.file_lens); pot.file_lens = NULL;
    if (1.0 != par.scale_surface_Cs){
    free(pot.r_LR); pot.r_LR = NULL;
    free(pot.pseudo_LR); pot.pseudo_LR = NULL; 
    }

    // Calculate the distorted nonlocal potentials

    if(flag.SO==1) {
    printf("\nSpin-orbit pseudopotential:\n");
    init_SO_projectors(SO_projectors, R, atom, &grid_par, &ist, &par);
    }
    /*** initialization for the non-local potential ***/
    if (flag.NL == 1){
    printf("\nNon-local pseudopotential:\n"); fflush(0);
    init_NL_projectors(nlc, nl, SO_projectors, R, atom, grid, &grid_par, &ist, &par, &flag);
    }
    // free memory allocated to SO_projectors
    if ( (flag.SO == 1) || (flag.NL == 1) ){
    free(SO_projectors); SO_projectors = NULL;
    }
    // ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 
    // ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
    // Potential mat elems
    // This code computes dU = U(r; R) - U(r; R_equil)
    // the difference between the potential energy surfaces of a Hamiltonian
    // that parametrically depends on the nuclear coordinates vs one that
    // only depends on the equilibrium nuclear configuration.
    // We calculate the variation in the energies as a perturbation expansion
    // E_0(R) = E_0(R_equil) + <phi_0|dU|phi_0> + sum_{i neq 0} |<phi_0|dU|phi_i>|^2/dE + ...
    // We verify that the first order term is sufficient to describe the changes
    // by computing its magnitude in comparison to the energy differences between
    // states in the electronic manifold at the equilibrium geometry.
    // ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 
    // ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    

    // 2. Calculate the distorted local potential
    if ((pot_local = (double *) calloc(ist.ngrid, sizeof(double))) == NULL){
        fprintf(stderr, "\nOUT OF MEMORY: pot_local_equil\n\n"); exit(EXIT_FAILURE);
    }
    pot.dr = (double *) calloc(ist.ngeoms * ist.n_atom_types, sizeof(double));
    pot.r = (double *) calloc(ist.ngeoms * ist.max_pot_file_len * ist.n_atom_types, sizeof(double));
    if (1.0 != par.scale_surface_Cs){
        // allocate memory for separately read LR potentials if the surface atoms will be charge balanced
        pot.r_LR = (double *) calloc(ist.ngeoms * ist.max_pot_file_len * ist.n_atom_types, sizeof(double));
    }
    pot.pseudo = (double *) calloc(ist.ngeoms * ist.max_pot_file_len * ist.n_atom_types, sizeof(double));
    if (1.0 != par.scale_surface_Cs){
        pot.pseudo_LR = (double *) calloc(ist.ngeoms * ist.max_pot_file_len * ist.n_atom_types, sizeof(double));
    }
    pot.file_lens = (long *) calloc(ist.ngeoms * ist.n_atom_types, sizeof(long));
    
    build_local_pot(pot_local, &pot, R, atom, grid, &grid_par, &ist, &par, &flag, &parallel);
    
    write_cube_file(pot_local, &grid_par, "local_pot.cube");
    
    free(pot.r); pot.r = NULL; 
    free(pot.pseudo); pot.pseudo = NULL; 
    free(pot.dr); pot.dr = NULL; 
    free(pot.file_lens); pot.file_lens = NULL;
    if (1.0 != par.scale_surface_Cs){
    free(pot.r_LR); pot.r_LR = NULL;
    free(pot.pseudo_LR); pot.pseudo_LR = NULL; 
    }


    // 3. Calculate the distorted nonlocal potentials

    if(flag.SO==1) {
    printf("\nSpin-orbit pseudopotential:\n");
    init_SO_projectors(SO_projectors, R, atom, &grid_par, &ist, &par);
    }
    /*** initialization for the non-local potential ***/
    if (flag.NL == 1){
    printf("\nNon-local pseudopotential:\n"); fflush(0);
    init_NL_projectors(nlc, nl, SO_projectors, R, atom, grid, &grid_par, &ist, &par, &flag);
    }
    // free memory allocated to SO_projectors
    if ( (flag.SO == 1) || (flag.NL == 1) ){
    free(SO_projectors); SO_projectors = NULL;
    }

    // 4. Read in the equilibrium wavefunction.

    //count number of states found
    jms = countlines("eval.dat");
    printf("%ld total states in psi.dat\n", jms); fflush(0);
  
    // allocate memory for all equilibrium wavefunctions, psi
    if ((psitot = (double *) calloc(ist.complex_idx * jms * ist.nspinngrid, sizeof(psitot[0]))) == NULL) nerror("psitot");

    //read psi from file
	ppsi = fopen("psi.dat" , "r");
    fread(&psi[0], sizeof(psitot[0]), ist.complex_idx*jms*ist.nspinngrid, ppsi);
	printf("Successfully read wavefunction\n"); fflush(0);

    // Determine which states are hole states and which states are electron states
    long a, ieof, nval;
    double evalloc, deloc;

    ist.homo_idx = ist.lumo_idx = 0;
    pf = fopen("eval.dat" , "r");
    for (i = ieof = 0; ieof != EOF; i++){
        ieof = fscanf(pf, "%ld %lg %lg", &a, &evalloc, &deloc);
        if (deloc < par.sigma_E_cut && evalloc < par.fermi_E) ist.homo_idx = i;
    }
    fclose(pf);

    nval = i - 1;
    pf = fopen("eval.dat" , "r");
    for (i = 0; i <= ist.homo_idx; i++){
        fscanf(pf, "%ld %lg %lg", &a, &evalloc, &deloc);
    }
    for (i = ist.homo_idx+1; i < nval; i++) {
        fscanf(pf, "%ld %lg %lg", &a, &evalloc, &deloc);
        if (deloc < par.sigma_E_cut) {
            ist.lumo_idx = i;
            break;
        }
    }
    fclose(pf);

    printf("index of homo, homo_idx = %ld; index of lumo, nlumo = %ld\n", ist.homo_idx, ist.lumo_idx); fflush(0);
    // Set the total number of electron and hole states in order to calculate the potential overlap integrals
    ist.total_homo = ist.homo_idx + 1; ist.total_lumo = ist.mn_states_tot - ist.total_homo;
    printf("total_homo = %ld total_lumo = %ld\n", ist.total_homo, ist.total_lumo); fflush(0);

    // Select which states should be used for the calculation. Must be converged.


    // // 4. Calculate matrix elements of the equilibrium wavefunction with U(r;R_equil)
    
    // long jgrid, jgrid_real, jgrid_imag, ims;
    // double eval, eval2;
    // // Loop over all M*N states
    // for (ims = 0; ims < ist.mn_states_tot; ims++) {
    //     // select the current state to compute sigma_E for
    //     for (jgrid = 0; jgrid < ist.nspinngrid; jgrid++) {
    //     jgrid_real = ist.complex_idx * jgrid;
    //     jgrid_imag = ist.complex_idx * jgrid + 1;

    //     psi[jgrid].re = psitot[ist.complex_idx*ims*ist.nspinngrid+jgrid_real];
    //     if (1 == flag.isComplex){
    //         psi[jgrid].im = psitot[ist.complex_idx*ims*ist.nspinngrid+jgrid_imag];
    //     } else if (0 == flag.isComplex){
    //         psi[jgrid].im = 0.0;
    //     }
        
    //     }
    //     memcpy(&phi[0],&psi[0],ist.nspinngrid*sizeof(phi[0]));
    //     // Apply the potential to |psi>: |phi> = V|psi>
    //     potential(phi, psi, pot_local, nlc, nl, &ist, &par, &flag);
    //     // Calculate the expectation value of H for wavefunc psi: <psi|H|psi> = <psi|phi> = sum_{jgrid} psi[jgrid] * phi[jgrid] * dv
    //     for (eval = 0.0, jgrid = 0; jgrid < ist.nspinngrid; jgrid++) {
    //     jgrid_real = ist.complex_idx*jgrid;
    //     jgrid_imag = ist.complex_idx*jgrid + 1;

    //     eval += psitot[ist.complex_idx*ims*ist.nspinngrid+jgrid_real] * phi[jgrid].re;
    //     if (1 == flag.isComplex){
    //         eval += psitot[ist.complex_idx*ims*ist.nspinngrid+jgrid_imag] * phi[jgrid].im;
    //     }
    //     }
    //     eval *= par.dv;

    // }
    // // 5. Calculate matrix elements of the equilibrium wavefunction with U(r;R)
    

    return 0;
}


/*****************************************************************************/
int countlines(char *filename){
  FILE* fp = fopen(filename,"r");
  int lines = 0;
  int ch;
  while(1){
    ch = fgetc(fp);
    if (feof(fp)){ break; }
    if(ch == '\n')
    {
      lines++;
    }
  }
  fclose(fp);
  return lines;
}

