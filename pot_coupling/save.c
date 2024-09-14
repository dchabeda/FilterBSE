#include "pot_coupling.h"

/*****************************************************************************/

void print_input_state(FILE *pf, flag_st *flag, grid_st *grid, par_st *par, index_st *ist, parallel_st *parallel){
    // ****** ****** ****** ****** ****** ****** 
    // Print grid parameters
    // ****** ****** ****** ****** ****** ****** 
    fprintf(pf, "\n\tGrid parameters from filter job:\n");
    fprintf(pf, "\t----------------------------------\n");
    fprintf(pf, "\tnx = %ld, ", grid->nx);
    fprintf(pf, "ny = %ld, ", grid->ny);
    fprintf(pf, "nz = %ld\n", grid->nz);
    fprintf(pf, "\tdx = %lg, ", grid->dx);
    fprintf(pf, "dy = %lg, ", grid->dy);
    fprintf(pf, "dz = %lg\n", grid->dz);
    fprintf(pf, "\tngrid = %ld nspin = %d\n", ist->ngrid, ist->nspin);

    // ****** ****** ****** ****** ****** ****** 
    // Set parameters & counters for BSE algorithm
    // ****** ****** ****** ****** ****** ****** 
    fprintf(pf, "\n\tParameters & counters for BSE algorithm:\n");
    fprintf(pf, "\t-------------------------------------------\n");
    if (-1 == ist->max_hole_states){
        fprintf(pf, "\tmaxHoleStates (# h+ in basis) not set, using all h+ states\n");
    } else if (ist->max_hole_states > 0){
        fprintf(pf, "\tmaxHoleStates (# h+ in basis) = %ld\n", ist->max_hole_states);
    } else {
        fprintf(pf, "ERROR: invalid maxHoleStates!\n");
        exit(EXIT_FAILURE);
    }

    if (-1 == ist->max_elec_states){
        fprintf(pf, "\tmaxElecStates (# e- in basis) not set, using all e- states\n");
    } else if (ist->max_elec_states > 0){
        fprintf(pf, "\tmaxElecStates (# e- in basis) = %ld\n", ist->max_elec_states);
    } else {
        fprintf(pf, "ERROR: invalid maxElecStates!\n");
        exit(EXIT_FAILURE);
    }
    
    fprintf(pf, "\tsigmaECut = %lg a.u.\n", par->sigma_E_cut);
    fprintf(pf, "\tDeltaE_hole = %lg a.u.", par->delta_E_hole);
    fprintf(pf, "  DeltaE_elec = %lg a.u.\n", par->delta_E_elec);
    fprintf(pf, "\tKEmax = %lg a.u.\n", par->KE_max);
    fprintf(pf, "\tfermiEnergy = %lg a.u.\n", par->fermi_E); 

    // ****** ****** ****** ****** ****** ****** 
    // Set options for dielectric properties
    // ****** ****** ****** ****** ****** ****** 
    fprintf(pf, "\n\tDielectric screening options:\n");
    fprintf(pf, "\t-------------------------------\n");
    fprintf(pf, "\tepsX = %lg\n", par->epsX);
    fprintf(pf, "\tepsY = %lg\n", par->epsY);
    fprintf(pf, "\tepsZ = %lg\n", par->epsZ);
    if (1 == flag->LR) fprintf(pf, "\tPseudopotentials were long ranged.\n");
    else fprintf(pf, "\tPseudopotentials were short ranged.\n");
     
    // ****** ****** ****** ****** ****** ****** 
    // Set options for spinor representation
    // ****** ****** ****** ****** ****** ****** 
    fprintf(pf, "\n\tParameters for two-component representation:\n");
    fprintf(pf, "\t----------------------------------------------\n");

    if (flag->useSpinors == 1) fprintf(pf, "\tSpinor wavefunctions turned ON!\n");
    else fprintf(pf, "\tSpinor wavefunctions turned OFF!\n");

    if (flag->SO == 1) fprintf(pf, "\tSpin-orbit coupling is ON!\n");
    else fprintf(pf, "\tSpin-orbit coupling is OFF!\n");

    if (flag->NL == 1) fprintf(pf, "\tNon-local potential is ON!\n");
    else fprintf(pf, "\tNon-local potential is OFF!\n");

    if (0 == flag->isComplex) fprintf(pf, "\tWavefunctions are REAL valued!\n\t  complex_idx = %d\n", ist->complex_idx);
    else if (1 == flag->isComplex) fprintf(pf, "\tWavefunctions are COMPLEX valued! complex_idx = %d\n", ist->complex_idx); 

    // ****** ****** ****** ****** ****** ****** 
    // Print optional output flags
    // ****** ****** ****** ****** ****** ******
    fprintf(pf, "\n\tFlags for optional output:\n");
    fprintf(pf, "\t----------------------------\n");

    if (flag->calcDarkStates == 1) fprintf(pf, "\tExchange operator not computed -> dark states\n");
    else fprintf(pf, "\tComputing spin-allowed sector of BSE matrix -> bright states\n");

    if (flag->timingSpecs == 1) fprintf(pf, "\tTiming specs will be printed\n");
    else fprintf(pf, "\tTiming specs will NOT be printed\n");

    if (flag->printFPDensity == 1) fprintf(pf, "\tFixed point densities of %d excitons will be printed\n", ist->n_FP_density);
    else fprintf(pf, "\tFixed point densities will NOT be printed\n");

    if (flag->calcSpinAngStat == 1) fprintf(pf, "\tSpin and Ang. Mom. statistics will be computed\n");
    else fprintf(pf, "\tSpin and Ang. Mom. statistics will NOT be computed\n");
    
    // ****** ****** ****** ****** ****** ****** 
    // Set options for parallelization
    // ****** ****** ****** ****** ****** ****** 
    fprintf(pf, "\n\tParameters for parallelization:\n");
    fprintf(pf, "\t-------------------------------\n");

    fprintf(pf, "\tnThreads (# OMP threads) = %ld\n", parallel->nthreads);
    
    // ****** ****** ****** ****** ****** ****** 
    // Restart flags
    // ****** ****** ****** ****** ****** ******
    fprintf(pf, "\n\tFlags for restarting computation:\n");
    fprintf(pf, "\t-----------------------------------\n");
    if (flag->saveCheckpoints == 1) fprintf(pf, "\tFilter will save checkpoints along the job run\n");
    else fprintf(pf, "\tNo checkpoint saves requested\n");
    
    if (flag->restartFromCheckpoint > -1) fprintf(pf, "\tFilter will restart from checkpoint %d\n", flag->restartFromCheckpoint);
    else fprintf(pf, "\tNo checkpoint specified for restart. Job will run in normal sequence.\n");
    
    return;
}

/*****************************************************************************/

void read_filter_output(char *file_name, double **psitot, double **eig_vals, double **sigma_E, xyz_st **R, grid_st *grid, double **gridx, double **gridy, double **gridz, index_st *ist, par_st *par, flag_st *flag){

    FILE *pf;
    long j;
    long output_tag;
    char *end_buffer, *eof; 
    eof = malloc(4*sizeof(eof[0])); end_buffer = malloc(4*sizeof(end_buffer[0]));

    strcpy(eof, "EOF");

    if( access(file_name, F_OK) == -1 ){
        printf("ERROR: no checkpoint file %s exists in directory\n", file_name);
        fprintf(stderr, "ERROR: no checkpoint file %s exists in directory\n", file_name);
        exit(EXIT_FAILURE);
    }

    pf = fopen(file_name, "r");
    // Read the output tag to confirm that the output.dat file is from the intended filter run
    // This is not enforced, but you can check that this was the intended run by looking at the 
    // number printed at the bottom of the filter "run.dat"
    fscanf(pf, "%ld\n", &output_tag);
    printf("\n\tOutput tag = %ld\n", output_tag);
    printf("\tCheck that this output tag matches your filter run.dat!\n");

    // Read ist
    printf("\n\tindex_st from filter...\n"); fflush(stdout);
    fscanf(pf, "%ld %ld", &ist->ngrid, &ist->nspinngrid); 
    fscanf(pf, "%ld", &ist->mn_states_tot); 
    fscanf(pf, "%ld %ld", &ist->natoms, &ist->n_atom_types); 
    for (j = 0; j < ist->n_atom_types; j++){ fscanf(pf, "%ld ", &ist->atom_types[j]); }
    fscanf(pf, "%d", &ist->nspin); 
    fscanf(pf, "%d", &ist->complex_idx); 
    /*
    printf("%ld %ld\n", ist->ngrid, ist->nspinngrid);
    printf("%ld", ist->mn_states_tot);
    printf("%ld %ld\n", ist->natoms, ist->n_atom_types);
    printf("\n%d\n", ist->nspin);
    printf("%d\n", ist->complex_idx); */
    // Read par
    printf("\tpar_st from filter...\n"); fflush(stdout);
    fscanf(pf, "%lg %lg", &par->KE_max, &par->fermi_E);
    printf("%lg %lg\n", par->KE_max, par->fermi_E); fflush(0);
    // Read flags
    printf("\tflag_st from filter...\n"); fflush(stdout);
    fscanf(pf, "%d %d %d %d %d", &flag->SO, &flag->NL, &flag->LR, &flag->useSpinors, &flag->isComplex);fflush(0);
    //printf("%d %d %d %d %d\n", flag->SO, flag->NL, flag->LR, flag->useSpinors, flag->isComplex);fflush(0);
    // Read conf
    printf("\tconf from filter...\n"); fflush(stdout);
    if(( *R = malloc(ist->natoms * sizeof(xyz_st))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for R in read_filter_output\n");
        exit(EXIT_FAILURE);
    }
    
    for (j = 0; j < ist->natoms; j++){
        fscanf(pf, "%lg %lg %lg", &((*R)[j].x), &((*R)[j].y), &((*R)[j].z));
        //printf("%lg %lg %lg\n", (*R)[j].x, (*R)[j].y, (*R)[j].z); 
    }

    // Read grid
    printf("\tgrid from filter...\n"); fflush(stdout);
    fscanf(pf, "%lg %lg %lg %lg %lg %lg %lg %lg", &grid->dx, &grid->dy, &grid->dz, &grid->dr, &grid->dv, &grid->dkx, &grid->dky, &grid->dkz);
    fscanf(pf, "%lg %lg %lg %lg %lg %lg", &grid->xmin, &grid->xmax, &grid->ymin, &grid->ymax, &grid->zmin, &grid->zmax);
    fscanf(pf, "%ld %ld %ld", &grid->nx, &grid->ny, &grid->nz); //printf("%ld %ld %ld\n", grid->nx, grid->ny, grid->nz);fflush(0);
    fscanf(pf, "%lg %lg %lg", &grid->nx_1, &grid->ny_1, &grid->nz_1);
    fscanf(pf, "%7ld", &grid->ngrid); //printf("\ngrid->ngrid = %ld\n", grid->ngrid);

    if((*gridx = malloc(grid->nx * sizeof(double))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for grid.x in read_filter_output\n");
        exit(EXIT_FAILURE);
    }
    if((*gridy = malloc(grid->ny * sizeof(double))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for grid.y in read_filter_output\n");
        exit(EXIT_FAILURE);
    }
    if((*gridz = malloc(grid->nz * sizeof(double))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for grid.x in read_filter_output\n");
        exit(EXIT_FAILURE);
    }

    //fseek(pf, 1 , SEEK_CUR);
    fread(*gridx, sizeof(double), grid->nx, pf);
    fread(*gridy, sizeof(double), grid->ny, pf);
    fread(*gridz, sizeof(double), grid->nz, pf);
    printf("gridz\n");
    //for (j = 0; j< grid->nx; j++){
        //printf("%lg\n", *(gridz)[j]);
    //}
    // Read eig_vals and sigma_E
    if ((*eig_vals = malloc(ist->mn_states_tot * sizeof(*eig_vals[0]))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for eig_vals in read_filter_output\n");
        exit(EXIT_FAILURE);
    }
    if ((*sigma_E = malloc(ist->mn_states_tot * sizeof(*sigma_E[0]))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for eig_vals in read_filter_output\n");
        exit(EXIT_FAILURE);
    }
    printf("\teig_vals from filter...\n"); fflush(stdout);
    fread(*eig_vals, sizeof(*eig_vals[0]), ist->mn_states_tot, pf);
    
    printf("\tsigma_E from filter...\n"); fflush(stdout);
    fread(*sigma_E, sizeof(*sigma_E[0]), ist->mn_states_tot, pf);
    for (j = 0; j< ist->mn_states_tot; j++){
        printf("%lg %lg\n", (*eig_vals)[j], (*sigma_E)[j]);
    }
    // Read psitot
    if ((*psitot = malloc(ist->complex_idx * ist->mn_states_tot * ist->nspinngrid * sizeof(psitot[0]))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for psitot in read_filter_output\n");
        exit(EXIT_FAILURE);
    }

    printf("\tpsitot from filter...\n"); fflush(stdout);
    fread(*psitot, sizeof(psitot[0]), ist->mn_states_tot * ist->nspinngrid * ist->complex_idx, pf);
    // The psitot will not be read in yet because there is no allocated memory for it.
    
    
    fseek(pf, 1 , SEEK_CUR);
    fscanf(pf, "%3s", end_buffer); 
    fclose(pf);
    printf("FOUND END BUFFER: %s\n", end_buffer);
    // printf(" The %s end buffer: %s\n", file_name, end_buffer); 
    if (strcmp((const char *) end_buffer, (const char *) eof) != 0){
        fprintf(stderr, "ERROR: restarting from %s failed. Bad END.\n", file_name);
        exit(EXIT_FAILURE);
    }

    return;
}
