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
  par->KE_max = 10.0; //NOTE: lowered this to 10 for tests
  par->checkpoint_id = 0;
  // Pseudopotential parameters
  flag->LR = 0; // Long range flag. By default, pseudopotentials are short ranged.
  // Spin-orbit and non-local terms
  flag->useSpinors = 0; // default is to use non-spinor wavefunctions
  flag->isComplex = 0;
  ist->nspin = 1; // turning on the useSpinors flag will set nspin to 2.
  flag->SO = 0; // computes the spin-orbit terms in the Hamiltonian
  flag->NL = 0; // computes the non-local terms in the Hamiltonian; automatically on if SO flag on
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

    long stlen;
    long j;
    long itmp;
    
    if( access( "unsafe_input.par", F_OK) != -1 ) {
        pf = fopen("unsafe_input.par", "r");
    
        while (fscanf(pf, "%s = %s", field, tmp) != EOF && i < 100) {
            // ****** ****** ****** ****** ****** ****** 
            // Set parameters&counters for BSE algorithm
            // ****** ****** ****** ****** ****** ****** 
            if (!strcmp(field, "mnStatesTot")) {
                ist->mn_states_tot= strtol(tmp, &endptr, 10);
                if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
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
    grid->dv = grid->dx * grid->dy * grid->dz;
    ist->nx = grid->nx;  ist->ny = grid->ny;  ist->nz = grid->nz;
    ist->ngrid = grid->ngrid = ist->nx * ist->ny * ist->nz;
    ist->nspinngrid = 2 * ist->ngrid;// hardocded for complex spinors, need to make general
    ist->complex_idx = 2; // these are hardocde for complex spinors, need to make general
    stlen = ist->complex_idx * ist->nspinngrid;

    printf("Done reading in unsafe_input.par\n\n");
    printf("xmin = %lg ymin = %lg zmin = %lg\n", grid->xmin, grid->ymin, grid->zmin);
    printf("nx = %ld ny = %ld nz = %ld\n", grid->nx, grid->ny, grid->nz);
    printf("dx = %lg dy = %lg dz = %lg\n", grid->dx, grid->dy, grid->dz);
    printf("natoms = %ld\n", ist->natoms);
    printf("mn_states_tot = %ld\n", ist->mn_states_tot);

    // Allocate memory for grid, psi, eigs, sigma_E, and R
    (*gridx)       =   (double*) calloc(grid->nx, sizeof(double));
    (*gridy)       =   (double*) calloc(grid->ny, sizeof(double));
    (*gridz)       =   (double*) calloc(grid->nz, sizeof(double));

    (*psitot)      =   (double*) calloc(ist->mn_states_tot * stlen, sizeof(double));

    (*eig_vals)    =   (double*) calloc(ist->mn_states_tot, sizeof(double));
    (*sigma_E)     =   (double*) calloc(ist->mn_states_tot, sizeof(double));

    (*R)           =   (xyz_st*) calloc(ist->natoms, sizeof(xyz_st));

    // Build the grid

    for (j = 0; j < grid->nx; j++) (*gridx)[j] = grid->xmin + j * grid->dx;
    for (j = 0; j < grid->ny; j++) (*gridy)[j] = grid->ymin + j * grid->dy;
    for (j = 0; j < grid->nx; j++) (*gridz)[j] = grid->zmin + j * grid->dz;

    // Define grid maxima
    grid->xmax = (*gridx)[grid->nx - 1];
    grid->ymax = (*gridy)[grid->ny - 1];
    grid->zmax = (*gridz)[grid->nz - 1];

    // Read in psitot from psi.par
    if( access( "psi.par", F_OK) != -1 ) {
        pf = fopen("psi.par", "r");
        fread((*psitot), stlen * sizeof(double), ist->mn_states_tot, pf);
        fclose(pf);
    } 
    else{
        printf("PROGRAM EXITING: psi.par does not exist in directory\n");
        fprintf(stderr, "PROGRAM EXITING: psi.par does not exist in directory\n");
        exit(EXIT_FAILURE);
    }

    // Read in the eig_vals and sigma_E
    if( access( "eval.par", F_OK) != -1 ) {
        pf = fopen("eval.par", "r");
        
        for (j = 0; j < ist->mn_states_tot; j++){
            fscanf(pf, "%ld %lg %lg", &itmp, &((*eig_vals)[j]),  &((*sigma_E)[j]));
        }

        fclose(pf);
    } 
    else{
        printf("PROGRAM EXITING: eval.par does not exist in directory\n");
        fprintf(stderr, "PROGRAM EXITING: eval.par does not exist in directory\n");
        exit(EXIT_FAILURE);
    }
    
    // Read in atomic configuration
    if( access( "conf.par", F_OK) != -1 ) {
        pf = fopen("conf.par", "r");
        
        for (j = 0; j < ist->natoms; j++){
            fscanf(pf, "%s %lg %lg %lg", &tmp, &((*R)[j].x) , &((*R)[j].y), &((*R)[j].z));
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
