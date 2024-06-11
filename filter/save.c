#include "fd.h"

/*****************************************************************************/

void print_input_state(FILE *pf, flag_st *flag, grid_st *grid, par_st *par, index_st *ist, parallel_st *parallel){
    // ****** ****** ****** ****** ****** ****** 
    // Print grid parameters
    // ****** ****** ****** ****** ****** ****** 
    fprintf(pf, "\n\tGrid parameters (preliminary):\n");
    fprintf(pf, "\t------------------------------\n");
    fprintf(pf, "\tnx = %ld, ", grid->nx);
    fprintf(pf, "ny = %ld, ", grid->ny);
    fprintf(pf, "nz = %ld\n", grid->nz);
    fprintf(pf, "\tdx = %lg, ", grid->dx);
    fprintf(pf, "dy = %lg, ", grid->dy);
    fprintf(pf, "dz = %lg\n", grid->dz);
    fprintf(pf, "\tngrid = %ld nspin = %d\n", ist->ngrid, ist->nspin);

    if (flag->centerConf == 1){
        fprintf(pf, "\tNanocrystal COM will be centered at the origin\n");
    } else if (flag->centerConf == 0){
        fprintf(pf, "\tNanocrystal configuration will not be centered.\n");
    } else {fprintf(stderr, "\tInvalid centerConf parameter!\n"); exit(EXIT_FAILURE);}
    // ****** ****** ****** ****** ****** ****** 
    // Set options for pseudopotential
    // ****** ****** ****** ****** ****** ****** 
    fprintf(pf, "\n\tPseudopotential options:\n");
    fprintf(pf, "\t--------------------------\n");
    if (1 == flag->interpolatePot) fprintf(pf, "\tPseudopotential parameters will be interpolated based on NC geometry\n");
    else fprintf(pf, "\tPseudopotential parameters will not be interpolated\n");
    if (1 == flag->useStrain) fprintf(pf, "\tStrain dependent pseudopotential will be calculated\n");
    else fprintf(pf, "\tStrain dependent pseudopotential will NOT be calculated\n");
    fprintf(pf, "\tCrystal structure: %s\n", par->crystal_structure);
    fprintf(pf, "\tOutmost material: %s\n", par->outmost_material);
    if (par->scale_surface_Cs != 1.0) fprintf(pf, "\tLong range component of surface Cs atom potentials scaled by %lg\n", par->scale_surface_Cs);
    else fprintf(pf, "\tSurface Cs atoms NOT rescaled for charge balancing\n");

    // ****** ****** ****** ****** ****** ****** 
    // Set parameters & counters for filter algorithm
    // ****** ****** ****** ****** ****** ****** 
    fprintf(pf, "\n\tParameters & counters for filter algorithm:\n");
    fprintf(pf, "\t-------------------------------------------\n");
    fprintf(pf, "\tstatesPerFilter (# energy targets) = %ld\n", ist->m_states_per_filter);
    fprintf(pf, "\tnFilterCycles (# random initial states) = %ld\n", ist->n_filter_cycles);
    fprintf(pf, "\tnCheby = %ld\n", ist->ncheby);
    fprintf(pf, "\tVBmin = %lg, ", par->VBmin);
    fprintf(pf, "VBmax = %lg\n", par->VBmax);
    fprintf(pf, "\tCBmin = %lg, ", par->CBmin);
    fprintf(pf, "CBmax = %lg\n", par->CBmax);
    fprintf(pf, "\tKEmax = %lg a.u.\n", par->KE_max);
    fprintf(pf, "\tfermiEnergy = %lg a.u.\n", par->fermi_E);

    if (flag->setTargets == 1){
        fprintf(pf, "\tsetNumEnergyTargets flag on: n_targets_VB = %ld n_targets_CB = %ld\n", par->n_targets_VB, par->n_targets_CB);
    } else {
        fprintf(pf, "\tDefault setting of energy targets will be used.\n");
    } 

    if (flag->setSeed == 1){
            fscanf(pf, "%ld", &par->rand_seed);
            fprintf(pf,"\tSetting initial filter random seed to -%ld\n", par->rand_seed);
    } else {fprintf(pf,"\tRandom seed will be generated based on clock at runtime\n");}
    
    if (flag->printPsiFilt == 1){
        fprintf(pf, "\tprintPsiFilt is on. psi-filt.dat and psi-filt.cube files will be printed\n");
    } else {
        fprintf(pf, "\tIntermediate filter wavefunctions will not be output to disk\n");
    }
    
    if (flag->printOrtho == 1){
        fprintf(pf, "\tprintOrtho is on. ortho.dat files will be printed\n");
    } else {
        fprintf(pf, "\tOrthogonalized wavefunctions will not be output to disk\n");
    }
    // ****** ****** ****** ****** ****** ****** 
    // Set options for parallelization
    // ****** ****** ****** ****** ****** ****** 
    fprintf(pf, "\n\tParameters for parallelization:\n");
    fprintf(pf, "\t-------------------------------\n");

    fprintf(pf, "\tnThreads (# OMP threads) = %ld\n", parallel->nthreads);
    
    // ****** ****** ****** ****** ****** ****** 
    // Set options for spin-orbit calculation
    // ****** ****** ****** ****** ****** ****** 
    fprintf(pf, "\n\tParameters for spin-orbit and non-local:\n");
    fprintf(pf, "\t----------------------------------------\n");

    if (flag->useSpinors == 1) fprintf(pf, "\tSpinor wavefunctions turned ON!\n");
    else fprintf(pf, "\tSpinor wavefunctions turned OFF!\n");

    if (flag->SO == 1) fprintf(pf, "\tSpin-orbit coupling is ON!\n");
    else fprintf(pf, "\tSpin-orbit coupling is OFF!\n");

    if (flag->NL == 1) fprintf(pf, "\tNon-local potential is ON!\n");
    else fprintf(pf, "\tNon-local potential is OFF!\n");

    if (flag->SO == 1) fprintf(pf, "\tRnlcut = %g Bohr\n", sqrt(par->R_NLcut2));

    if (0 == flag->isComplex) fprintf(pf, "\tWavefunctions are REAL valued!\n\t  complex_idx = %d\n", ist->complex_idx);
    else if (1 == flag->isComplex) fprintf(pf, "\tWavefunctions are COMPLEX valued! complex_idx = %d\n", ist->complex_idx); 

    // ****** ****** ****** ****** ****** ****** 
    // Print optional output flags
    // ****** ****** ****** ****** ****** ******
    fprintf(pf, "\n\tFlags for optional output:\n");
    fprintf(pf, "\t--------------------------\n");

    if (flag->calcPotOverlap == 1) fprintf(pf, "\tMatrix elements of pseudopotential will be calculated\n");
    else fprintf(pf, "\tMatrix elements of the potential will not be calculated\n");

    if (flag->getAllStates == 1) fprintf(pf, "\tAll eigenstates will be saved to disk regardless of variance\n");
    else fprintf(pf, "\tEigenstates with variance above %lg a.u. will not be saved to disk\n",par->sigma_E_cut);

    if (flag->timeHamiltonian == 1) fprintf(pf, "\tHamiltonian timing specs will be printed\n");
    else fprintf(pf, "\tHamiltonian timing specs will not be printed\n");

    if (flag->printCubes == 1) fprintf(pf, "\tCube files of top %d states (sigma_E_cut < %lg) in VB & CB will be printed\n", ist->ncubes, par->sigma_E_cut);
    else fprintf(pf, "\tCube files will not be printed\n");

    if (flag->calcSpinAngStat == 1) fprintf(pf, "\tSpin and Ang. Mom. statistics will be computed\n");
    else fprintf(pf, "\tSpin and Ang. Mom. statistics will not be computed\n");
    
    // ****** ****** ****** ****** ****** ****** 
    // Restart flags
    // ****** ****** ****** ****** ****** ******
    fprintf(pf, "\n\tFlags for restarting computation:\n");
    fprintf(pf, "\t---------------------------------\n");

    if (flag->retryFilter == 1) fprintf(pf, "\tFilter will restart if no eigenstates acquired after diag\n");
    else fprintf(pf, "\tFilter will NOT restart if no eigenstates obtained after diag.\n");
    if (flag->saveCheckpoints == 1) fprintf(pf, "\tFilter will save checkpoints along the job run\n");
    else fprintf(pf, "\tNo checkpoint saves requested\n");
    if (flag->restartFromCheckpoint > -1) fprintf(pf, "\tFilter will restart from checkpoint %d\n", flag->restartFromCheckpoint);
    else fprintf(pf, "\tNo checkpoint specified for restart. Job will run in normal sequence.\n");
    
    return;
}

/*****************************************************************************/

void save_job_state(char *file_name, int checkpoint_id, double *psitot, double *pot_local, double *ksqr, zomplex *an, double *zn, double *ene_targets, long *nl,\
    nlc_st *nlc, grid_st *grid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){

    // This function writes an output file that saves all the data about the job state
    // So that filter can be restarted

    FILE *pf;
    long j;
    long checkpoint_tag = random();

    printf("Save state tag: %ld\n", checkpoint_tag);

    pf = fopen(file_name, "w");
    fprintf(pf, "%ld\n", checkpoint_tag);

    fprintf(pf, "%ld %ld %ld\n", ist->m_states_per_filter, ist->n_filter_cycles, ist->mn_states_tot);
    fprintf(pf, "%ld %ld %ld %ld\n", ist->homo_idx, ist->lumo_idx, ist->total_homo, ist->total_lumo);
    fprintf(pf, "%ld %ld %ld\n", ist->ngrid, ist->nspinngrid, ist->ncheby);
    fprintf(pf, "%ld %ld %ld\n", ist->natoms, ist->n_atom_types, ist->n_max_atom_types);
    for (j = 0; j < ist->n_atom_types; j++){ fprintf(pf, "%ld ", ist->atom_types[j]);}
    fprintf(pf, "\n%ld %ld %ld %ld\n", ist->max_pot_file_len, ist->n_NL_gridpts, ist->n_NL_atoms, ist->nproj);
    fprintf(pf, "%d %d %d\n", ist->nspin, ist->ncubes, ist->ngeoms);
    fprintf(pf, "%d\n", ist->complex_idx);
    fprintf(pf, "%d %d\n", ist->crystal_structure_int, ist->outmost_material_int);
    fprintf(pf, "%ld %ld %ld\n", ist->nx, ist->ny, ist->nz);
    fprintf(pf, "%ld\n", ist->nthreads);
    
    fprintf(pf, "%lg %lg %lg %lg %lg %lg\n", par->Vmin, par->Vmax, par->VBmin, par->VBmax, par->CBmin, par->CBmax);
    fprintf(pf, "%lg\n", par->scale_surface_Cs);
    fprintf(pf, "%ld %ld\n", par->n_targets_VB, par->n_targets_CB);
    fprintf(pf, "%ld\n", par->rand_seed);
    fprintf(pf, "%lg %lg %lg %lg %lg\n", par->KE_max, par->fermi_E, par->dt, par->dE, par->dE_1);
    fprintf(pf, "%lg %lg\n", par->R_NLcut2, par->sigma_E_cut);
    fprintf(pf, "%d\n", par->t_rev_factor);
    fprintf(pf, "%d\n", par->checkpoint_id);
    fprintf(pf, "%s %s\n", par->crystal_structure, par->outmost_material);
    fprintf(pf, "%lg\n", par->dv);
    
    fprintf(pf, "%d %d %d %d %d\n", flag->centerConf, flag->setTargets, flag->setSeed, flag->interpolatePot, flag->useStrain);
    fprintf(pf, "%d %d %d %d\n", flag->SO, flag->NL, flag->useSpinors, flag->isComplex);
    fprintf(pf, "%d %d %d %d\n", flag->calcPotOverlap, flag->getAllStates, flag->timeHamiltonian, flag->calcSpinAngStat);
    fprintf(pf, "%d %d %d\n", flag->retryFilter, flag->alreadyTried, flag->saveCheckpoints);
    
    fprintf(pf, "%ld\n", parallel->nthreads);

    for (j = 0; j < ist->n_NL_atoms * ist-> n_NL_gridpts; j++){
        fprintf(pf, "%ld\n", nlc[j].jxyz);
        fprintf(pf, "%lg %lg %lg %lg %lg %lg\n", nlc[j].y1[0].re, nlc[j].y1[0].im, nlc[j].y1[1].re, nlc[j].y1[1].im, nlc[j].y1[2].re, nlc[j].y1[2].im);
        fprintf(pf, "%lg %lg %lg %lg %lg\n", nlc[j].proj[0], nlc[j].proj[1], nlc[j].proj[2], nlc[j].proj[3], nlc[j].proj[4]);
        fprintf(pf, "%lg %lg %lg %lg %lg\n", nlc[j].NL_proj[0], nlc[j].NL_proj[1], nlc[j].NL_proj[2], nlc[j].NL_proj[3], nlc[j].NL_proj[4]);
        fprintf(pf, "%d %d %d %d %d\n", nlc[j].NL_proj_sign[0], nlc[j].NL_proj_sign[1], nlc[j].NL_proj_sign[2], nlc[j].NL_proj_sign[3], nlc[j].NL_proj_sign[4]);
        fprintf(pf, "%lg %lg %lg %lg", nlc[j].r, nlc[j].r2_1, nlc[j].r2, nlc[j].Vr);
    }
    
    fprintf(pf, "%lg %lg %lg %lg %lg %lg %lg %lg\n", grid->dx, grid->dy, grid->dz, grid->dr, grid->dv, grid->dkx, grid->dky, grid->dkz);
    fprintf(pf, "%lg %lg %lg %lg %lg %lg\n", grid->xmin, grid->xmax, grid->ymin, grid->ymax, grid->zmin, grid->zmax);
    fprintf(pf, "%ld %ld %ld\n", grid->nx, grid->ny, grid->nz);
    fprintf(pf, "%lg %lg %lg\n", grid->nx_1, grid->ny_1, grid->nz_1);
    fprintf(pf, "%ld", grid->ngrid);
    fwrite(grid->x, sizeof(grid->x[0]), grid->nx, pf);
    fwrite(grid->y, sizeof(grid->x[0]), grid->ny, pf);
    fwrite(grid->z, sizeof(grid->x[0]), grid->nz, pf);

    fwrite(ksqr, sizeof(double), ist->ngrid, pf); 
    fwrite(pot_local, sizeof(double), ist->ngrid, pf);
    fwrite(an, sizeof(an[0]), ist->ncheby * ist->m_states_per_filter, pf);
    fwrite(zn, sizeof(zn[0]), ist->ncheby, pf);
    fwrite(ene_targets, sizeof(ene_targets[0]), ist->m_states_per_filter, pf);
    fwrite(nl, sizeof(nl[0]), ist->natoms, pf);

    if ((checkpoint_id > 0) && checkpoint_id <= 1){
        fwrite(psitot, sizeof(psitot[0]), ist->complex_idx * par->t_rev_factor * ist->nspinngrid * ist->mn_states_tot, pf);
    } else if (checkpoint_id > 1){
        fwrite(psitot, sizeof(psitot[0]), ist->complex_idx * ist->nspinngrid * ist->mn_states_tot, pf);
    }
    fprintf(pf, "\nEOF\n");

    fclose(pf);
    
    return;
}


/*****************************************************************************/

void restart_from_save(char *file_name, int checkpoint_id, double *psitot, double *pot_local, double *ksqr, zomplex *an, double *zn, double *ene_targets, long *nl,\
    nlc_st *nlc, grid_st *grid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){

    FILE *pf;
    long j;
    long checkpoint_tag;
    char *end_buffer, *eof; 
    eof = malloc(4*sizeof(eof[0])); end_buffer = malloc(4*sizeof(eof[0]));

    strcpy(eof, "EOF");

    // Create the checkpoint file name header from the checkpoint_id
    //sprintf(str, "checkpoint_%d", checkpoint_id);
    //sprintf(chkfile_name, "%s_ist.dat", str);

    if( access(file_name, F_OK) == -1 ){
        printf("ERROR: no checkpoint file %s exists in directory\n", file_name);
        fprintf(stderr, "ERROR: no checkpoint file %s exists in directory\n", file_name);
        exit(EXIT_FAILURE);
    }

    pf = fopen(file_name, "r");
    // Read the checkpoint tag to confirm that all checkpoint files are from the same job run
    fscanf(pf, "%ld", &checkpoint_tag);
    printf("Save state tag: %ld\n", checkpoint_tag);

    fscanf(pf, "%ld %ld %ld", &ist->m_states_per_filter, &ist->n_filter_cycles, &ist->mn_states_tot);
    fscanf(pf, "%ld %ld %ld %ld", &ist->homo_idx, &ist->lumo_idx, &ist->total_homo, &ist->total_lumo);
    fscanf(pf, "%ld %ld %ld", &ist->ngrid, &ist->nspinngrid, &ist->ncheby);
    fscanf(pf, "%ld %ld %ld", &ist->natoms, &ist->n_atom_types, &ist->n_max_atom_types);
    for (j = 0; j < ist->n_atom_types; j++){fscanf(pf, "%ld ", &ist->atom_types[j]);}
    fscanf(pf, "%ld %ld %ld %ld", &ist->max_pot_file_len, &ist->n_NL_gridpts, &ist->n_NL_atoms, &ist->nproj);
    fscanf(pf, "%d %d %d", &ist->nspin, &ist->ncubes, &ist->ngeoms);
    fscanf(pf, "%d", &ist->complex_idx);
    fscanf(pf, "%d %d", &ist->crystal_structure_int, &ist->outmost_material_int);
    fscanf(pf, "%ld %ld %ld", &ist->nx, &ist->ny, &ist->nz);
    fscanf(pf, "%ld", &ist->nthreads);
    
    fscanf(pf, "%lg %lg %lg %lg %lg %lg", &par->Vmin, &par->Vmax, &par->VBmin, &par->VBmax, &par->CBmin, &par->CBmax);
    fscanf(pf, "%lg", &par->scale_surface_Cs);
    fscanf(pf, "%ld %ld", &par->n_targets_VB, &par->n_targets_CB);
    fscanf(pf, "%ld", &par->rand_seed);
    fscanf(pf, "%lg %lg %lg %lg %lg", &par->KE_max, &par->fermi_E, &par->dt, &par->dE, &par->dE_1);
    fscanf(pf, "%lg %lg", &par->R_NLcut2, &par->sigma_E_cut);
    fscanf(pf, "%d", &par->t_rev_factor);
    fscanf(pf, "%d", &par->checkpoint_id);
    fscanf(pf, "%s %s", par->crystal_structure, par->outmost_material);
    fscanf(pf, "%lg", &par->dv);
    
    fscanf(pf, "%d %d %d %d %d", &flag->centerConf, &flag->setTargets, &flag->setSeed, &flag->interpolatePot, &flag->useStrain);
    fscanf(pf, "%d %d %d %d", &flag->SO, &flag->NL, &flag->useSpinors, &flag->isComplex);
    fscanf(pf, "%d %d %d %d", &flag->calcPotOverlap, &flag->getAllStates, &flag->timeHamiltonian, &flag->calcSpinAngStat);
    fscanf(pf, "%d %d %d", &flag->retryFilter, &flag->alreadyTried, &flag->saveCheckpoints);
    
    fscanf(pf, "%ld", &parallel->nthreads);

    for (j = 0; j < ist->n_NL_atoms * ist-> n_NL_gridpts; j++){
        fscanf(pf, "%ld", &nlc[j].jxyz);
        fscanf(pf, "%lg %lg %lg %lg %lg %lg", &nlc[j].y1[0].re, &nlc[j].y1[0].im, &nlc[j].y1[1].re, &nlc[j].y1[1].im, &nlc[j].y1[2].re, &nlc[j].y1[2].im);
        fscanf(pf, "%lg %lg %lg %lg %lg", &nlc[j].proj[0], &nlc[j].proj[1], &nlc[j].proj[2], &nlc[j].proj[3], &nlc[j].proj[4]);
        fscanf(pf, "%lg %lg %lg %lg %lg", &nlc[j].NL_proj[0], &nlc[j].NL_proj[1], &nlc[j].NL_proj[2], &nlc[j].NL_proj[3], &nlc[j].NL_proj[4]);
        fscanf(pf, "%d %d %d %d %d", &nlc[j].NL_proj_sign[0], &nlc[j].NL_proj_sign[1], &nlc[j].NL_proj_sign[2], &nlc[j].NL_proj_sign[3], &nlc[j].NL_proj_sign[4]);
        fscanf(pf, "%lg %lg %lg %lg", &nlc[j].r, &nlc[j].r2_1, &nlc[j].r2, &nlc[j].Vr);
    }
    
    fscanf(pf, "%lg %lg %lg %lg %lg %lg %lg %lg", &grid->dx, &grid->dy, &grid->dz, &grid->dr, &grid->dv, &grid->dkx, &grid->dky, &grid->dkz);
    fscanf(pf, "%lg %lg %lg %lg %lg %lg", &grid->xmin, &grid->xmax, &grid->ymin, &grid->ymax, &grid->zmin, &grid->zmax);
    fscanf(pf, "%ld %ld %ld", &grid->nx, &grid->ny, &grid->nz);
    fscanf(pf, "%lg %lg %lg", &grid->nx_1, &grid->ny_1, &grid->nz_1);
    fscanf(pf, "%ld", &grid->ngrid);
    fread(grid->x, sizeof(grid->x[0]), ist->nx, pf);
    fread(grid->y, sizeof(grid->x[0]), ist->ny, pf);
    fread(grid->z, sizeof(grid->x[0]), ist->nz, pf);

    fread(ksqr, sizeof(double), ist->ngrid, pf);
    fread(pot_local, sizeof(double), ist->ngrid, pf);
    fread(an, sizeof(an[0]), ist->ncheby * ist->m_states_per_filter, pf);
    fread(zn, sizeof(zn[0]), ist->ncheby, pf);
    fread(ene_targets, sizeof(ene_targets[0]), ist->m_states_per_filter, pf);
    fread(nl, sizeof(nl[0]), ist->natoms, pf);
    
    if ((checkpoint_id > 0) && (checkpoint_id <= 1)) {
        fread(psitot, sizeof(psitot[0]), ist->complex_idx * par->t_rev_factor * ist->nspinngrid * ist->mn_states_tot, pf);
    } else if (checkpoint_id > 1) {
        fread(psitot, sizeof(psitot[0]), ist->complex_idx * ist->nspinngrid * ist->mn_states_tot, pf);
    }

    fseek(pf, 1, SEEK_CUR);
    fscanf(pf, "%3s", end_buffer); 
    fclose(pf);

    // printf(" The %s end buffer: %s\n", file_name, end_buffer); 
    if (strcmp((const char *) end_buffer, (const char *) eof) != 0){
        fprintf(stderr, "ERROR: restarting from %s failed. Bad END.\n", file_name);
        exit(EXIT_FAILURE);
    }

    return;

}


