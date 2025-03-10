/****************************************************************************/

#include "fd.h"

/*****************************************************************************/

void read_input(flag_st *flag, grid_st *grid, index_st *ist, par_st *par, parallel_st *parallel){
  /*******************************************************************
  * This function reads the input.par file and initializes the job   *
  * inputs:                                                          *
  *  [flag_st *par] ptr to flag_st holding job flags                 *
  *  [grid_st *grid] pointer to grid struct, holds the grid          *
  *  [index_st *ist] pointer to counters, indices, and lengths       *
  *  [par_st *par] ptr to par_st holding VBmin, VBmax... params      *
  *  [parallel] struct holding parallelization options (for norm)    *
  * outputs: void                                                    *
  ********************************************************************/

  FILE *pf;
  int i = 0;
  char field[1000], tmp[1000], *endptr;
  // long tmp = 0;

  // ****** ****** ****** ****** ****** ****** 
  // Setting the default behavior 
  // ****** ****** ****** ****** ****** ****** 

  // BSE algorithm parameters
  ist->max_hole_states = -1;
  ist->max_elec_states = -1;
  par->KE_max = 20.0; //NOTE: lowered this to 10 for tests
  par->checkpoint_id = 0;
  // Pseudopotential parameters
  flag->LR = 0; // Long range flag. By default, pseudopotentials are short ranged.
  
  // Spin-orbit and non-local terms
  flag->useSpinors = 0; // default is to use non-spinor wavefunctions
  flag->isComplex = 0;
  ist->nspin = 1; // turning on the useSpinors flag will set nspin to 2.
  flag->SO = 0; // computes the spin-orbit terms in the Hamiltonian
  flag->NL = 0; // computes the non-local terms in the Hamiltonian; automatically on if SO flag on
  flag->restartCoulomb = 0;
  flag->coulombDone = 0;
  flag->calcCoulombOnly = 0;
  // Optional output flags
  flag->calcDarkStates = 0;
  flag->calcSpinAngStat = 1; // Are angular momentum statistics computed. Only available with SO coupling
  flag->timingSpecs = 0; // print timing info for computing the Hamiltonian
  par->fermi_E = -0.18; // This default value is not good for LR potentials.
  par->sigma_E_cut = 0.0001;
  flag->saveCheckpoints = 0; // by default do not save job checkpoints
  flag->restartFromCheckpoint = -1; // by default do not restart from checkpoint
  flag->printFPDensity=0;
  ist->n_FP_density=0;
  

  // Parse the input file
  if( access( "input.par", F_OK) != -1 ) {
    pf = fopen("input.par", "r");

    while (fscanf(pf, "%s = %s", field, tmp) != EOF && i < 100) {
      // ****** ****** ****** ****** ****** ****** 
      // Set parameters&counters for BSE algorithm
      // ****** ****** ****** ****** ****** ****** 
      if (!strcmp(field, "maxHoleStates")) {
          ist->max_hole_states = strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "maxElecStates")) {
          ist->max_elec_states = strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "sigmaECut")) {
          par->sigma_E_cut = strtod(tmp, &endptr);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to double.\n"); exit(EXIT_FAILURE);}
          if (0 > par->sigma_E_cut) {fprintf(stderr, "Error: sigma_E_cut is negative.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "deltaEhole")) {
          par->delta_E_hole = strtod(tmp, &endptr);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to double.\n"); exit(EXIT_FAILURE);}
          if (0 > par->delta_E_hole) {fprintf(stderr, "Error: delta_E_hole is negative.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "deltaEelec")) {
          par->delta_E_elec = strtod(tmp, &endptr);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to double.\n"); exit(EXIT_FAILURE);}
          if (0 > par->delta_E_elec) {fprintf(stderr, "Error: delta_E_elec is negative.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "KEmax")) {
          par->KE_max = strtod(tmp, &endptr);
          if (*endptr != '\0') {printf("Error converting string to double.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "fermiEnergy")) {
          par->fermi_E = strtod(tmp, &endptr);
          if (*endptr != '\0') {printf("Error converting string to double.\n"); exit(EXIT_FAILURE);}
      }
      // ****** ****** ****** ****** ****** ****** 
      // Set options for dielectric properties
      // ****** ****** ****** ****** ****** ******
      else if (!strcmp(field, "epsX")) {
          par->epsX = strtod(tmp, &endptr);
          if (*endptr != '\0') {printf("Error converting string to double.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "epsY")) {
          par->epsY = strtod(tmp, &endptr);
          if (*endptr != '\0') {printf("Error converting string to double.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "epsZ")) {
          par->epsZ = strtod(tmp, &endptr);
          if (*endptr != '\0') {printf("Error converting string to double.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "longRange")) {
          flag->LR = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {printf("Error converting string to double.\n"); exit(EXIT_FAILURE);}
      }
      // ****** ****** ****** ****** ****** ****** 
      // Set options for parallelization
      // ****** ****** ****** ****** ****** ****** 
        else if (!strcmp(field, "nThreads")) {
          parallel->nthreads = strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } 
      // ****** ****** ****** ****** ****** ****** 
      // Set options for spin-orbit calculation
      // ****** ****** ****** ****** ****** ****** 
      else if (!strcmp(field, "useSpinors")) {
          flag->useSpinors = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "spinOrbit")) {
          flag->SO = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "NonLocal")) {
          flag->NL = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      }
      // ****** ****** ****** ****** ****** ****** 
      // Set options for additional output
      // ****** ****** ****** ****** ****** ******
      else if (!strcmp(field, "timingSpecs")) {
          flag->timingSpecs = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "printFPDensity")) {
          flag->printFPDensity = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
          if (flag->printFPDensity == 1){fscanf(pf, "%d", &ist->n_FP_density);}
      } else if (!strcmp(field, "calcDarkStates")) {
          flag->calcSpinAngStat = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "calcSpinAngStat")) {
          flag->calcSpinAngStat = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "saveCheckpoints")) {
          flag->saveCheckpoints = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "restartFromCheckpoint")) {
          flag->restartFromCheckpoint = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "restartCoulomb")) {
        flag->restartCoulomb = (int) strtol(tmp, &endptr, 10);
        if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "coulombDone")) {
        flag->coulombDone = (int) strtol(tmp, &endptr, 10);
        if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "calcCoulombOnly")) {
        flag->calcCoulombOnly = (int) strtol(tmp, &endptr, 10);
        if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      }
      // ****** ****** ****** ****** ****** ****** 
      // Handle exceptions
      // ****** ****** ****** ****** ****** ******
      else {
          printf("\nFIELD NOT RECOGNIZED: %s\n", field);
          printf("\nInvalid input field and/ or format -> equal sign required after each field\n");
          printf("Only allowed fields are (case-sensitive):\n\n");
          printf("maxHoleStates = int (max no. of hole states to include in exciton basis)\n");
          printf("maxElecStates = int (max no. of electron states to include in exciton basis)\n");
          printf("sigmaECut = double, the cutoff std. dev allowed for filtered quasiparticle states\n");
          printf("deltaEhole = double, the desired energy span of h+ basis states\n");
          printf("deltaEelec = double, the desired energy span of e- basis states\n");
          printf("epsX = double (dielectric constant along X direction)\n");
          printf("epsY = double (dielectric constant along Y direction)\n");
          printf("epsZ = double (dielectric constant along Z direction)\n");
          printf("longRange = int (if 1, the pseudopots have long range terms; no truncate)\n");
          printf("nThreads = 1-16 (number of compute threads to parallelize over)\n");
          printf("fermi_E = double (fermi_E of the system)\n");
          printf("KEmax = double (maximum kinetic energy value considered)\n");
          printf("spinOrbit = int (0 for no spinOrbit, 1 for spinOrbit)\n");
          printf("NonLocal = int (0 for no non-local, 1 for non-local potential)\n");
          printf("timingSpecs = int (if 1 print timing specs)\n");
          printf("printFPDensity = int (1 to print fixed point exciton densities, 0 to not)\n");
          printf("If printCubes = 1, the next entry MUST specify the number of excitons to print\n");
          printf("calcDarkStates = int (if 1, calculate dark (spin-forbidden) manifold of excitons)\n");
          printf("calcSpinAngStat = int (if 1, calculate spin and ang. mom. statistics for eigenstates)\n");
          printf("saveCheckpoints = int, if 1 then save states will be generated along the job run.\n");
          printf("restartFromCheckpoint = int, value is the ID of the checkpoint that the job should restart from.\n");
          printf("restartCoulomb = int, (if 1, Coulomb matrix elements will be read to restart job).\n");
          printf("coulombDone = int, (if 1, then Kernel has already been calc'd, will compute BSE).\n");
          
          fflush(stdout);
          exit(EXIT_FAILURE);
      }
      i++;
    } 
  } else{
      printf("PROGRAM EXITING: input.par does not exist in directory\n");
      fprintf(stderr, "PROGRAM EXITING: input.par does not exist in directory\n");
      exit(EXIT_FAILURE);
  }
  fclose(pf);

  // Sanity check the input.par parameters
  if ((par->epsX < 0.0) || (par->epsY < 0.0) || (par->epsZ < 0.0)){
    fprintf(stderr, "ERROR: one of the dielectric constants is negative!");
    exit(EXIT_FAILURE);
  }

  // ****** ****** ****** ****** ****** ****** 
  // Set dependencies based on input parameters
  // ****** ****** ****** ****** ****** ****** 
  if (flag->SO == 1) {
    flag->useSpinors = 1;
    flag->NL = 1; // SO automatically switches on NL
    // par->R_NLcut2 = 1.5 +  6.0 * log(10.0) + 3.0 * grid->dx;
  }
  if (flag->useSpinors == 1) {
    flag->isComplex = 1;
    ist->nspin = 2; // generate spinor wavefunctions
  }

  // Set parameters that depend on these inputs
  ist->ngrid_1 = 1.0 / (double)(ist->ngrid); // for rescaling FFT
  ist->nthreads = parallel->nthreads;
  ist->complex_idx = flag->isComplex + 1;
  par->dv = grid->dv;
  ist->nx = grid->nx; ist->ny = grid->ny; ist->nz = grid->nz;

  // Using input parameters, print the current job state
  print_input_state(stdout, flag, grid, par, ist, parallel);

  return;
}

/****************************************************************************/

void read_unsafe_input(
    double** psitot,
    double** eig_vals,
    double** sigma_E,
    xyz_st** R, 
    grid_st *grid,
    double** gridx,
    double** gridy,
    double** gridz,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel
    ){

    // Read in psitot, set up the grid, etc.
    // but in an unsafe way without using output.dat
    // Wavefunctions might not be aligned with the grid,
    // leading to uncontrolled errors
    // USE AT YOUR OWN RISK
    FILE *pf;
    int i = 0;
    char field[1000], tmp[1000], *endptr;

    unsigned long   stlen;
    unsigned long   j;
    unsigned long   cntr;
    unsigned long   itmp;
    unsigned long   homo_idx;
    unsigned long   vb_mindex;
    unsigned long   cb_maxdex;
    unsigned long   nst_tot;
    unsigned long   nholes = 0;
    unsigned long   nelecs = 0;

    double          dbltmp;
    
    if( access( "unsafe_input.par", F_OK) != -1 ) {
        pf = fopen("unsafe_input.par", "r");
    
        while (fscanf(pf, "%s = %s", field, tmp) != EOF && i < 100) {
            // ****** ****** ****** ****** ****** ****** 
            // Set parameters&counters for BSE algorithm
            // ****** ****** ****** ****** ****** ****** 
            if (!strcmp(field, "mnStatesTot")) {
                nst_tot = strtol(tmp, &endptr, 10);
                if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
            } else if (!strcmp(field, "nHoles")) {
                nholes = strtol(tmp, &endptr, 10);
                if (*endptr != '\0') {fprintf(stderr, "Error converting string to double.\n"); exit(EXIT_FAILURE);}
            } else if (!strcmp(field, "nElecs")) {
                nelecs = strtol(tmp, &endptr, 10);
                if (*endptr != '\0') {fprintf(stderr, "Error converting string to double.\n"); exit(EXIT_FAILURE);}
            } else if (!strcmp(field, "nAtoms")) {
                ist->natoms = strtol(tmp, &endptr, 10);
                if (*endptr != '\0') {fprintf(stderr, "Error converting string to double.\n"); exit(EXIT_FAILURE);}
            } else if (!strcmp(field, "nx")) {
                grid->nx = strtod(tmp, &endptr);
                if (*endptr != '\0') {fprintf(stderr, "Error converting string to double.\n"); exit(EXIT_FAILURE);}
            } else if (!strcmp(field, "ny")) {
                grid->ny = strtod(tmp, &endptr);
                if (*endptr != '\0') {fprintf(stderr, "Error converting string to double.\n"); exit(EXIT_FAILURE);}
            } else if (!strcmp(field, "nz")) {
                grid->nz = strtod(tmp, &endptr);
                if (*endptr != '\0') {fprintf(stderr, "Error converting string to double.\n"); exit(EXIT_FAILURE);}
            } else if (!strcmp(field, "xmin")) {
                grid->xmin = strtod(tmp, &endptr);
                if (*endptr != '\0') {fprintf(stderr, "Error converting string to double.\n"); exit(EXIT_FAILURE);}
            }  else if (!strcmp(field, "ymin")) {
                grid->ymin = strtod(tmp, &endptr);
                if (*endptr != '\0') {fprintf(stderr, "Error converting string to double.\n"); exit(EXIT_FAILURE);}
            } else if (!strcmp(field, "zmin")) {
                grid->zmin = strtod(tmp, &endptr);
                if (*endptr != '\0') {fprintf(stderr, "Error converting string to double.\n"); exit(EXIT_FAILURE);}
            } else if (!strcmp(field, "dx")) {
                grid->dx = strtod(tmp, &endptr);
                if (*endptr != '\0') {fprintf(stderr, "Error converting string to double.\n"); exit(EXIT_FAILURE);}
            } else if (!strcmp(field, "dy")) {
                grid->dy = strtod(tmp, &endptr);
                if (*endptr != '\0') {fprintf(stderr, "Error converting string to double.\n"); exit(EXIT_FAILURE);}
            } else if (!strcmp(field, "dz")) {
                grid->dz = strtod(tmp, &endptr);
                if (*endptr != '\0') {fprintf(stderr, "Error converting string to double.\n"); exit(EXIT_FAILURE);}
            }
            // ****** ****** ****** ****** ****** ****** 
            // Handle exceptions
            // ****** ****** ****** ****** ****** ******
            else {
                printf("\nFIELD NOT RECOGNIZED in unsafe_input.par: %s\n", field);
                printf("\nInvalid input field and/ or format -> equal sign required after each field\n");
                
                fflush(stdout);
                exit(EXIT_FAILURE);
            }
            i++;
        } 
        fclose(pf);
    } 
    else{
        printf("PROGRAM EXITING: unsafe_input.par does not exist in directory\n");
        fprintf(stderr, "PROGRAM EXITING: unsafe_input.par does not exist in directory\n");
        exit(EXIT_FAILURE);
    }
    
    // Assign and construct system size variables
    if ((nholes != 0) && (nelecs != 0) ){
        ist->mn_states_tot = nholes + nelecs;
    } else{
        printf("ERROR: unsafe_input has 0 nholes or 0 nelecs!\n");
        exit(EXIT_FAILURE);
    }
    

    ist->nx = grid->nx;  ist->ny = grid->ny;  ist->nz = grid->nz;
    ist->ngrid = grid->ngrid = ist->nx * ist->ny * ist->nz;
    grid->dkx = par->dkx = (TWOPI) / ((double) grid->nx * grid->dx);
    grid->dky = par->dky = (TWOPI) / ((double) grid->ny * grid->dy);
    grid->dkz = par->dkz = (TWOPI) / ((double) grid->nz * grid->dz);

    ist->ngrid_1 = 1.0 / (double)ist->ngrid;
    
    grid->nx_1 = 1.0 / (double) grid->nx;
    grid->ny_1 = 1.0 / (double) grid->ny;
    grid->nz_1 = 1.0 / (double) grid->nz;
    
    par->dx = grid->dx; par->dy = grid->dy; par->dz = grid->dz;
    grid->dr = par->dr = sqrt( grid->dx*grid->dx + grid->dy*grid->dy + grid->dz*grid->dz );
    grid->dv = par->dv = grid->dx * grid->dy * grid->dz;

    par->xmin = grid->xmin; par->ymin = grid->ymin; par->zmin = grid->zmin;
    
    ist->nspinngrid = 2 * ist->ngrid;// hardocded for complex spinors, need to make general
    ist->complex_idx = 2; // these are hardocde for complex spinors, need to make general
    stlen = ist->complex_idx * ist->nspinngrid;

    printf("Done reading in unsafe_input.par\n\n");
    printf("xmin = %lg ymin = %lg zmin = %lg\n", grid->xmin, grid->ymin, grid->zmin);
    printf("nx = %ld ny = %ld nz = %ld\n", grid->nx, grid->ny, grid->nz);
    printf("dx = %lg dy = %lg dz = %lg\n", grid->dx, grid->dy, grid->dz);
    printf("grid->dv = %lg \n", grid->dv);
    printf("natoms = %ld\n", ist->natoms);
    printf("mn_states_tot = %ld\n", ist->mn_states_tot);

    // Allocate memory for grid, psi, eigs, sigma_E, and R
    (*gridx)       =   (double*) calloc(grid->nx, sizeof(double));
    (*gridy)       =   (double*) calloc(grid->ny, sizeof(double));
    (*gridz)       =   (double*) calloc(grid->nz, sizeof(double));

    (*psitot)      =   (double*) calloc(ist->mn_states_tot * stlen, sizeof(double));

    (*eig_vals)    =   (double*) calloc(nst_tot, sizeof(double));
    (*sigma_E)     =   (double*) calloc(nst_tot, sizeof(double));

    (*R)           =   (xyz_st*) calloc(ist->natoms, sizeof(xyz_st));

    // Build the grid

    for (j = 0; j < grid->nx; j++){ (*gridx)[j] = grid->xmin + j * grid->dx; }
    for (j = 0; j < grid->ny; j++){ (*gridy)[j] = grid->ymin + j * grid->dy; }
    for (j = 0; j < grid->nz; j++){ (*gridz)[j] = grid->zmin + j * grid->dz; }

    // Define grid maxima
    grid->xmax = (*gridx)[grid->nx - 1];
    grid->ymax = (*gridy)[grid->ny - 1];
    grid->zmax = (*gridz)[grid->nz - 1];
    par->xmax = grid->xmax; par->ymax = grid->ymax; par->zmax = grid->zmax;
    
    // Read in the eig_vals and sigma_E
    if( access( "eval.par", F_OK) != -1 ) {
        pf = fopen("eval.par", "r");
        
        for (j = 0; j < nst_tot; j++){
            fscanf(pf, "%lu %lg %lg", &itmp, &((*eig_vals)[j]),  &((*sigma_E)[j]));
        }

        fclose(pf);
    }  
    else{
        printf("PROGRAM EXITING: eval.par does not exist in directory\n");
        fprintf(stderr, "PROGRAM EXITING: eval.par does not exist in directory\n");
        exit(EXIT_FAILURE);
    }

    // Get the indices of the homo and lumo
    get_fmo_idxs(*eig_vals, *sigma_E, -0.19, 0.001, nst_tot, &homo_idx);

    // Read in psitot from psi.par

    // 1. Get the indices of the first and last states
    vb_mindex = homo_idx - nholes + 1;
    cb_maxdex = vb_mindex + ist->mn_states_tot;

    printf("vb_mindex = %lu cb_maxdex = %lu\n", vb_mindex, cb_maxdex);

    if( access( "psi.par", F_OK) != -1 ) {
        pf = fopen("psi.par", "r");
        fseek(pf, vb_mindex * stlen * sizeof(double), SEEK_SET);
        fread((*psitot), stlen * sizeof(double), ist->mn_states_tot, pf);
        fclose(pf);
    } 
    else{
        printf("PROGRAM EXITING: psi.par does not exist in directory\n");
        fprintf(stderr, "PROGRAM EXITING: psi.par does not exist in directory\n");
        exit(EXIT_FAILURE);
    }

    // Update eig_vals and sigma_E to only have ist->mn_states_tot
    // double *tmpptr;
    // tmpptr = (double*) realloc(*eig_vals, ist->mn_states_tot * sizeof(double));
    // if (tmpptr == NULL) {
    //     fprintf(stderr, "Memory reallocation failed for eig_vals\n");
    //     exit(EXIT_FAILURE);
    // }
    // *eig_vals = tmpptr;

    // tmpptr = (double*) realloc(*sigma_E, ist->mn_states_tot * sizeof(double));
    // if (tmpptr == NULL) {
    //     fprintf(stderr, "Memory reallocation failed for sigma_E\n");
    //     exit(EXIT_FAILURE);
    // }
    // *sigma_E = tmpptr;
    
    // Read in the short eig_vals and sigma_E
    
    pf = fopen("eval.par", "r");
    cntr = 0;
    for (j = 0; j < nst_tot; j++){
        if ( (j >= vb_mindex) && (j < cb_maxdex) ){
            fscanf(pf, "%lu %lg %lg", &itmp, &((*eig_vals)[cntr]),  &((*sigma_E)[cntr]));
            printf("eig_vals[j = %lu] = %lg\n",j, ((*eig_vals)[cntr]) );
            cntr++;
        } else{
            fscanf(pf, "%lu %lg %lg", &itmp, &dbltmp,  &dbltmp);
        }
    }
    fclose(pf);
    
    // Read in atomic configuration
    if( access( "conf.par", F_OK) != -1 ) {
        pf = fopen("conf.par", "r");
        
        for (j = 0; j < ist->natoms; j++){
            fscanf(pf, "%s %lg %lg %lg", endptr, &((*R)[j].x) , &((*R)[j].y), &((*R)[j].z));
        }

        fclose(pf);
    } 
    else{
        printf("PROGRAM EXITING: conf.par does not exist in directory\n");
        fprintf(stderr, "PROGRAM EXITING: conf.par does not exist in directory\n");
        exit(EXIT_FAILURE);
    }

    printf("Done with function read_unsafe_input\n"); fflush(0);

    return;
}

/****************************************************************************/

void read_bseCoeff(int n_xton,int numExcStatesToRead, zomplex* u, FILE* pf){
  int i,j;

  if(numExcStatesToRead>100){numExcStatesToRead=100;}

  for (i = 0; i < n_xton; i++) {
    for (j = 0; j < numExcStatesToRead; j++) {  
    fscanf (pf,"{%lg, %lg}\t", &(u[i*n_xton+j].re), &(u[i*n_xton+j].im));
    } 
  }
}

/****************************************************************************/

void get_fmo_idxs(
    double*          eig_vals,
    double*          sigma_E,
    double           fermiE,
    double           secut,
    unsigned long    n_elems,
    unsigned long*   homo_idx
    ){

    int i;
    int hidx = 0;
    int lidx;
    int cntr = 0;

    // HOMO index 
    for (i = 0; i < n_elems; i++){
        // If the sigma_E val > secut, skip
        if ( (eig_vals[i] < fermiE) && (sigma_E[i] < secut) ){
            hidx = i;
        }
    }
    (*homo_idx) = hidx;

    printf("Found homo_idx = %lu\n", *homo_idx);

    return;
}
