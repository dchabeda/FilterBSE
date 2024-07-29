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

  // NC configuration parameters
  ist->n_max_atom_types = N_MAX_ATOM_TYPES;
  flag->centerConf = 1; // this should honestly always be 1
  // Basis set configuration
  flag->readGrid = 0; // Read the grid points from an input file
  flag->useGaussianBasis = 0; // Use atom-centered Gaussian basis set
  // Filter algorithm parameters
  par->KE_max = 10.0; //NOTE: lowered this to 10 for tests
  flag->setTargets = 0;
  flag->printPsiFilt = 0;
  flag->printOrtho = 0;
  par->checkpoint_id = 0;
  // Pseudopotential parameters
  ist->max_pot_file_len = 8192;
  flag->useStrain = 0; // By default, do not compute strain dependent terms in pseudopotential
  flag->LR = 0; // Long range flag. By default, pseudopotentials are short ranged.
  strcpy(par->crystal_structure, "unknown"); // set the following parameters for strain dependent potentials to NULL
  strcpy(par->outmost_material, "unknown");
  ist->crystal_structure_int = -1;
  ist->outmost_material_int = -1;
  ist->ngeoms = 1; // number of different psuedopotential geometries (eg. cubic/ortho) for interpolating psuedopotentials
  par->scale_surface_Cs = 1.0; // By default, do not charge balance the surface Cs atoms
  // Spin-orbit and non-local terms
  flag->useSpinors = 0; // default is to use non-spinor wavefunctions
  flag->isComplex = 0;
  ist->nspin = 1; // turning on the useSpinors flag will set nspin to 2.
  flag->SO = 0; // computes the spin-orbit terms in the Hamiltonian
  flag->NL = 0; // computes the non-local terms in the Hamiltonian; automatically on if SO flag on
  ist->nproj = 5; // number of terms to expand projections in. converged by 5
  par->t_rev_factor = 1; // can time rev filt'rd states to get 2X eigst8. mem alloc multiplied by par.t_rev_factor
  // Optional output flags
  flag->saveOutput = 1; // By default, write the formatted output file that can be read by BSE to recreate job state
  flag->calcPotOverlap = 0; // Calculates matrix elements of the potential <i|V|j>
  flag->getAllStates = 1; // which states are written to disk
  flag->calcSpinAngStat = 0; // Are angular momentum statistics computed. Only available with SO coupling
  flag->timeHamiltonian = 0; // print timing info for computing the Hamiltonian
  par->fermi_E = -0.18; // This default value is not good for LR potentials.
  flag->printNorm = 0; // do not print output from normalize_all
  par->sigma_E_cut = 0.01;
  flag->saveCheckpoints = 0; // by default do not save job checkpoints
  flag->restartFromCheckpoint = -1; // by default do not restart from checkpoint
  // Restart job flags
  flag->retryFilter = 0; // When = 1, if no eigstates acquired, then retry the filter job 
  flag->alreadyTried = 0; // gets tripped to 1 after the first time retrying a Filter.
  
  // Parse the input file
  if (access( "input.par", F_OK) != -1 ) {
    pf = fopen("input.par", "r");

    while (fscanf(pf, "%s = %s", field, tmp) != EOF && i < 100) {
      // ****** ****** ****** ****** ****** ****** 
      // Set geometry parameters for grid
      // ****** ****** ****** ****** ****** ****** 
      //printf("field = %s tmp = %s\n", field, tmp);fflush(0);
      if (!strcmp(field, "nx")) {
          grid->nx = strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "ny")) {
          grid->ny = strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "nz")) {
          grid->nz = strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "dGrid")) {
          grid->dx = strtod(tmp, &endptr);
          if (*endptr != '\0') {printf("Error converting string to double.\n"); exit(EXIT_FAILURE);}
          grid->dy = grid->dz = grid->dx; // If the dGrid flag is used, make all grid dimensions equal
      } else if (!strcmp(field, "dx")) {
          grid->dx = strtod(tmp, &endptr);
          if (*endptr != '\0') {printf("Error converting string to double.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "dy")) {
          grid->dy = strtod(tmp, &endptr);
          if (*endptr != '\0') {printf("Error converting string to double.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "dz")) {
          grid->dz = strtod(tmp, &endptr);
          if (*endptr != '\0') {printf("Error converting string to double.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "centerConf")) {
          flag->centerConf = strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "readGrid")) {
          flag->readGrid = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "useGaussianBasis")) {
          flag->useGaussianBasis = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      }
      // ****** ****** ****** ****** ****** ****** 
      // Set options for pseudopotentials
      // ****** ****** ****** ****** ****** ******
      else if (!strcmp(field, "useStrain")) {
          flag->useStrain = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {printf("Error converting string to double.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "interpolatePot")) {
          flag->interpolatePot = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "scaleSurfaceCs")) {
          par->scale_surface_Cs = strtod(tmp, &endptr);
          if (*endptr != '\0') {printf("Error converting string to double.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "longRange")) {
          flag->LR = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {printf("Error converting string to double.\n"); exit(EXIT_FAILURE);}
      } else if (0 == strcmp(field, "crystalStructure")) {
          strcpy(par->crystal_structure, tmp);
      } else if (0 == strcmp(field, "outmostMaterial")) {
          strcpy(par->outmost_material, tmp);
      }
      // ****** ****** ****** ****** ****** ****** 
      // Set parameters&counters for filter algorithm
      // ****** ****** ****** ****** ****** ****** 
      else if (!strcmp(field, "nFilterCycles")) {
          ist->n_filter_cycles = strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "mStatesPerFilter")) {
          ist->m_states_per_filter = strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "nCheby")) {
          ist->ncheby = strtol(tmp, &endptr, 10);
      } else if (!strcmp(field, "VBmin")) {
          par->VBmin = strtod(tmp, &endptr);
          if (*endptr != '\0') {printf("Error converting string to double.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "VBmax")) {
          par->VBmax = strtod(tmp, &endptr);
          if (*endptr != '\0') {printf("Error converting string to double.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "CBmin")) {
          par->CBmin = strtod(tmp, &endptr);
          if (*endptr != '\0') {printf("Error converting string to double.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "CBmax")) {
          par->CBmax = strtod(tmp, &endptr);
          if (*endptr != '\0') {printf("Error converting string to double.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "KEmax")) {
          par->KE_max = strtod(tmp, &endptr);
          if (*endptr != '\0') {printf("Error converting string to double.\n"); exit(EXIT_FAILURE);}
      }else if (!strcmp(field, "fermiEnergy")) {
          par->fermi_E = strtod(tmp, &endptr);
          if (*endptr != '\0') {printf("Error converting string to double.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "setTargets")) {
          flag->setTargets = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
          if (flag->setTargets == 1){
            fscanf(pf, "%ld %ld", &par->n_targets_VB, &par->n_targets_CB);
          }
      } else if (!strcmp(field, "setSeed")) {
          flag->setSeed = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
          if (flag->setSeed == 1){fscanf(pf, "%ld", &par->rand_seed);}
      } else if (!strcmp(field, "printPsiFilt")) {
          flag->printPsiFilt = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "printOrtho")) {
          flag->printOrtho = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
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
      else if (!strcmp(field, "calcPotOverlap")) {
          flag->calcPotOverlap = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "getAllStates")) {
          flag->getAllStates = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
          if (flag->getAllStates == 0){
            fscanf(pf, "%lg", &par->sigma_E_cut);
          }
      } else if (!strcmp(field, "sigmaECut")) {
          par->sigma_E_cut = strtod(tmp, &endptr);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to double.\n"); exit(EXIT_FAILURE);}
          if (0 > par->sigma_E_cut) {fprintf(stderr, "Error: sigma_E_cut is negative.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "timeHamiltonian")) {
          flag->timeHamiltonian = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "printCubes")) {
          flag->printCubes = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
          if (flag->printCubes == 1){fscanf(pf, "%d", &ist->ncubes);}
      } else if (!strcmp(field, "calcSpinAngStat")) {
          flag->calcSpinAngStat = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "printNorm")) {
          flag->printNorm = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "retryFilter")) {
          flag->retryFilter = (int) strtol(tmp, &endptr, 10);
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
          printf("nx = int (number of grid points in x direction)\n");
          printf("ny = int (number of grid points in y direction)\n");
          printf("nz = int (number of grid points in z direction)\n");
          printf("dx = double (grid spacing in x direction, units of Bohr)\n");
          printf("dy = double (grid spacing in y direction, units of Bohr)\n");
          printf("dz = double (grid spacing in z direction, units of Bohr)\n");
          printf("dGrid = double (grid spacing; sets the values of dx, dy, and dz equal)\n");
          printf("centerConf = int (if 1, center the NC atoms at the COM)\n");
          printf("readGrid = int (if 1, read grid from input file grid.par)\n");
          printf("useGaussianBasis = int (if 1, use atom-centered Gaussian basis set)\n");
          printf("mStatesPerFilter = int (number of total energy targets for each filter cycle)\n");
          printf("nFilterCycles = int (number of filter cycles/number random initial wavefunctions)\n");
          printf("nCheby = int (number of terms in the Chebyshev expansion)\n");
          printf("VBmin = double (bottom of valence band energy window)\n");
          printf("VBmax = double (top of valence band energy window)\n");
          printf("CBmin = double (bottom of conduction band energy window)\n");
          printf("CBmax = double (top of conduction band energy window)\n");
          printf("useStrain = int (if 1, calculate strain dependent pseudopots)\n");
          printf("interpolatePot = int (if 1, interpolate between cubic and orthorhombic potentials)\n");
          printf("longRange = int (if 1, the pseudopots have long range terms; no truncate)\n");
          printf("crystalStructure = string (name of crystal structure e.g. wurtzite)\n");
          printf("outmostMaterial = string (name of outmost layer if core-shell e.g. CdS)\n");
          printf("scaleSurfaceCs = double (fractional scaling factor of surface Cs to balance charge)\n");
          printf("nThreads = 1-16 (number of compute threads to parallelize over)\n");
          printf("fermiEnergy = double (fermi_E of the system)\n");
          printf("KEmax = double (maximum kinetic energy value considered)\n");
          printf("spinOrbit = int (0 for no spinOrbit, 1 for spinOrbit)\n");
          printf("NonLocal = int (0 for no non-local, 1 for non-local potential)\n");
          printf("setTargets = int (0 if half/half split of VB/CB targets suffices for your job)\n");
          printf("If setTargets = 1, the next two entries MUST be \'n_targets_VB n_targets_CB\'\n");
          printf("calcPotOverlap = int (calculate matrix overlaps of pseudopotential)\n");
          printf("sigmaECut = double, the cutoff std. dev allowed for printed eigenstates/cubes\n");
          printf("getAllStates = int (0 to only save states with small variance, 1 to save all)\n");
          printf("If getAllStates = 0, the next entry MUST specify the cutoff criteria \'sigma_E_cut\'\n");
          printf("timeHamiltonian = int (1 to print timing information, 0 to not)\n");
          printf("printCubes = int (1 to print wavefunction cube files, 0 to not)\n");
          printf("If printCubes = 1, the next entry MUST specify the number of eigenstates to print\n");
          printf("calcSpinAngStat = int (calculate spin and ang. mom. statistics for eigenstates)\n");
          printf("setSeed = int (if 1, set the random seed in for filter to generate exactly reproducible wavefunctions)\n");
          printf("If setSeed = 1, the next entry MUST specify the random seed as an integer \'rand_seed\'\n");
          printf("printNorm = int, if 1 then norms of wavefunctions are printed every 100 chebyshev iterations\n");
          printf("retryFilter = int, if 1 then if no eigenstates obtained after diag, then filter is restarted.\n");
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

  // ****** ****** ****** ****** ****** ****** 
  // Set dependencies based on input parameters
  // ****** ****** ****** ****** ****** ****** 
  if (flag->SO == 1) {
    flag->useSpinors = 1;
    flag->NL = 1; // SO automatically switches on NL
    par->R_NLcut2 = 1.5 +  6.0 * log(10.0) + 3.0 * grid->dx;
    //par->R_NLcut2 = 0.49 * 6.0 * log(10.0); // radius of grid points around atom for which to compute NL terms
  }
  if (flag->useSpinors == 1) {
    flag->isComplex = 1;
    ist->nspin = 2; // generate spinor wavefunctions
    par->t_rev_factor = 2; // give double the memory allocation to psitot
  }
  if (flag->interpolatePot == 1) {
    ist->ngeoms = 2; // give double memory to pot.r and pot.pseudo vectors
  }
  ist->ngrid = grid->ngrid = grid->nx * grid->ny * grid->nz;
  ist->nspinngrid = ist->nspin * ist->ngrid;
  ist->mn_states_tot = ist->n_filter_cycles * ist->m_states_per_filter;
  ist->nthreads = parallel->nthreads;
  ist->complex_idx = flag->isComplex + 1;

  // Get the number of atoms (needed to initialize the)
  pf = fopen("conf.par" , "r");
  if (pf) {
    fscanf(pf, "%ld", &ist->natoms); 
    fclose(pf);
    } else {
    fprintf(stderr, "ERROR: conf.par must be in current working directory\n");
    exit(EXIT_FAILURE);
  } 

  // Using input parameters, print the current job state
  print_input_state(stdout, flag, grid, par, ist, parallel);

  return;
}



void read_conf(char *file_name, xyz_st *R, atom_info *atom, index_st *ist, par_st *par, flag_st *flag){
  /*******************************************************************
  * This function reads the conf.par file and initializes the NC     *
  * Also, atom-specific parameters (SO/NL, local geom, etc) are set  *
  * Available atom types/IDs that can be successfully read are:      *
  * H  = 1                                                           *
  * P1 = 2                                                           *
  * P2 = 3                                                           *
  * P3 = 4                                                           *
  * P4 = 5                                                           *
  * Si = 14                                                          *
  * P = 15                                                           *
  * S = 16                                                           *
  * Zn = 30                                                          *
  * Ga = 31                                                          *
  * As = 33                                                          *
  * Se = 34                                                          *
  * Cd = 48                                                          *
  * In = 49                                                          *
  * Te = 52                                                          *
  * I = 53                                                           *
  * Cs = 55                                                          *
  * Pb = 82                                                          *
  * inputs:                                                          *
  *  [xyz_st *R] pointer to double arrays containing atom coords     *
  *  [atom_info *atom] ptr to struct containing atom-specific params *
  *  [index_st *ist] pointer to counters, indices, and lengths       *
  *  [par_st *par] ptr to par_st holding VBmin, VBmax... params      *
  *  [flag_st *flag] ptr to flag_st holding job flags                *
  * outputs: void                                                    *
  ********************************************************************/
  
  FILE *pf, *pw;
  long i, j; 
  double xd, yd, zd;

  // initialize the atom identifiers to -1 
  for (j = 0; j < ist->n_max_atom_types; j++){
    ist->atom_types[j] = -1; 
  }

  // initialize counters to 0
  ist->n_atom_types = 0;
  ist->n_NL_atoms = 0;
  xd = yd = zd = 0.0;

  pf = fopen(file_name , "r");
  // This has already been set, but there should be no harm in overwriting it again... famous last words
  fscanf(pf, "%ld", &ist->natoms); // reading first line so that the file pointer moves to the line with coordinates

  // loop over all atoms in conf.par
  for (i = 0; i < ist->natoms; i++){
    // Read each line of the conf.par file
    fscanf (pf, "%s %lf %lf %lf\n", atom[i].atyp, &R[i].x, &R[i].y, &R[i].z);
    printf ("%s %lf %lf %lf\n", atom[i].atyp, R[i].x, R[i].y, R[i].z, atom[i].Zval);
    // Get the atomic number of each atom
    atom[i].Zval = assign_atom_number(atom[i].atyp);
    printf ("%s %lf %lf %lf %d\n", atom[i].atyp, R[i].x, R[i].y, R[i].z, atom[i].Zval);
    // update list of unique atoms in ist->atom_types
    for (j = 0; j <= ist->n_atom_types; j++){
      //if atomtype already in list, add its index to the list and break
      if (atom[i].Zval == ist->atom_types[j]){
        atom[i].idx = j; 
        break;
      }
      // else check if we can write this position
      else if (j == ist->n_atom_types){
        ist->atom_types[j] = atom[i].Zval;
        ist->n_atom_types++;
        atom[i].idx = j;
        break;
      }
      else if (j == ist->n_max_atom_types){
        fprintf(stderr, "Too many distinct atom types! Exiting..."); fflush(0);
        exit(EXIT_FAILURE);
      }
    }

    // increment the n_NL_atoms parameter if NL flag is on
    // (number of atoms for which nonlocal terms will be computed)
    if ((flag->NL == 1) && (atom[i].Zval > 13)) ist->n_NL_atoms++;

    // add the coordinates to the dimension counter for calc'ing the COM
    xd += R[i].x;
    yd += R[i].y;
    zd += R[i].z;
  }
  fclose(pf);

  printf("\tnatoms = %ld\n", ist->natoms);
  printf("\tn_atom_types = %ld\n", ist->n_atom_types);
  printf("\tthe atoms are [ ");
  char atyp[3];
    for (int k = 0; k<ist->n_atom_types; k++){
      assign_atom_type(atyp, ist->atom_types[k]);
      printf("%s ", atyp);
    }
  printf("]\n");
  if (1 == flag->NL) printf("\tn_NL_atoms = %ld\n", ist->n_NL_atoms);
  
  // Center the NC configuration if requested in input file
  // should always happen unless the configuration is already centered at 0.
  // otherwise the configuration will likely be outside of the grid
  if (flag->centerConf){
    xd /= (double)(ist->natoms); // COM coordinate in x
    yd /= (double)(ist->natoms); // COM coordinate in y
    zd /= (double)(ist->natoms); // COM coordinate in z
    for (i = 0; i < ist->natoms; i++){
      // shift the COM to 0
      R[i].x -= xd; 
      R[i].y -= yd;
      R[i].z -= zd;
    }
  }

  // Print out conf.dat, 
  // the file with the atomic configuration used in the calculation
  pw = fopen("conf.dat" , "w");
  fprintf(pw,"%ld\n", ist->natoms);
  for (i = 0; i < ist->natoms; i++) {
    if (flag->SO != 1){
      fprintf(pw, "%s %g %g %g %ld\n", atom[i].atyp, R[i].x, R[i].y, R[i].z, atom[i].idx);
    } 
    if (flag->SO == 1){
      fprintf(pw, "%s %g %g %g %g\n", atom[i].atyp, R[i].x, R[i].y, R[i].z, atom[i].SO_par);
    }
  }
  fclose(pw);
  
  // If strain dependent pseudopotentials requested, determine the crystal structure and outmost material
	if (flag->useStrain){
    ist->crystal_structure_int = assign_crystal_structure(par->crystal_structure);
    
    /*** Assign outmostMaterialInt ***/
    ist->outmost_material_int = assign_outmost_material(par->outmost_material);
    
  }

  return;
}

/*****************************************************************************/

void read_pot(pot_st *pot, xyz_st *R, atom_info *atom, index_st *ist, par_st *par, flag_st *flag){
  /*******************************************************************
  * This function reads the potential files for each atom in conf    *
  * inputs:                                                          *
  *  [pot] struct containing arrays for reading atom pseudopot files *
  *  [R] array of atomic coordinates                                 *
  *  [atom_info *atom] ptr to struct containing atom-specific params *
  *  [index_st *ist] pointer to counters, indices, and lengths       *
  *  [par_st *par] ptr to par_st holding VBmin, VBmax... params      *
  *  [flag_st *flag] ptr to flag_st holding job flags                *
  * outputs: void                                                    *
  ********************************************************************/

  FILE *pf;  
  xyz_st R_min, R_max;
  double dr_check;
  long i, j, atyp_idx; 
  char str[30], tmpstr[30], atype[4], *req; req = malloc(10*sizeof(req[0]));
  // The number of geometries in the pseudopotential for each atom type. usually 1, more if interpolating
  long ngeoms = ist->ngeoms; 
  // pot->file_lens is the number of points in the pseudopotential file (also equivalent to the number of lines in the file)
  // n_atom_types is the number of distinct atoms in the NC configuration. Ex. 2 for CdSe, 3 for CsPbI3
  long ntype = ist->n_atom_types;

  // The ligand potentials are constructed internally by creating a Gaussian on the grid with
  // equation lig = a * exp(- r^2 / b)
  //               P1        P2         P3       P4        PC5    PC6
  double a[6] = { 0.64,    -0.384,     0.050, -0.684,     0.065, 0.080};
  double b[6] = {2.2287033, 2.2287033, 6.000,  2.2287033, 6.000, 6.000};


  // The array pot->r contains the r values in the psuedopotential file
  // The array pot->pseudo contains the Vloc values in the pseudopotential file
  for (atyp_idx = 0; atyp_idx < ngeoms * ntype * ist->max_pot_file_len; atyp_idx++) {
    pot->r[atyp_idx] = pot->pseudo[atyp_idx] = 0;
  }
  // The number of lines in the pseudopotential file can differ for different systems
  // the array pot->file_lens contains the number of lines in the pseudopotential file for each atom type
  for (atyp_idx = 0; atyp_idx < ngeoms * ntype; atyp_idx++) {
    pot->file_lens[atyp_idx] = 0;
  }
  // If the surface Cs atoms will be charge rescaled, then the surface atoms need to be identified
  // Find the min and max value of the configuration to find the edge Cs atoms
  
  if (1.0 != par->scale_surface_Cs){
    printf("\tIdentifying surface Cs atoms\n");
    pf = fopen("surace_Cs.dat", "w");

    R_min.x = R_min.y = R_min.z = 1e4;
    R_max.x = R_max.y = R_max.z = -1e4;

    for (j = 0; j < ist->natoms; j++){
      // Get min of coordinates
      if (R[j].x < R_min.x) R_min.x = R[j].x;
      if (R[j].y < R_min.y) R_min.y = R[j].y;
      if (R[j].z < R_min.z) R_min.z = R[j].z;
      // Get max of coordinates
      if (R[j].x > R_max.x) R_max.x = R[j].x;
      if (R[j].y > R_max.y) R_max.y = R[j].y;
      if (R[j].z > R_max.z) R_max.z = R[j].z;
    }
    
    // Search for surface atoms
    for (j = 0; j < ist->natoms; j++){
      // By default, all atoms should have the full LR term
      atom[j].LR_par = 1.0;
      if (atom[j].Zval == 55){
        // If the Cs is an edge atom, it will be close to one of the extrema and its LR_par should be scaled
        if ((fabs(R[j].x - R_min.x) < PbIBondMax) || (fabs(R[j].x - R_max.x) < PbIBondMax) ||
            (fabs(R[j].y - R_min.y) < PbIBondMax) || (fabs(R[j].y - R_max.y) < PbIBondMax) ||
            (fabs(R[j].z - R_min.z) < PbIBondMax) || (fabs(R[j].z - R_max.z) < PbIBondMax)){
            fprintf(pf, "Cs%ld %lg %lg %lg\n", j, R[j].x, R[j].y, R[j].z);
            atom[j].LR_par = par->scale_surface_Cs;
        }
      }
    }
    fclose(pf);
  }

  printf("\t");
  // ******  ****** ****** ****** ****** ****** ****** ****** ******
  // ******  ****** ****** ****** ****** ****** ****** ****** ******
  // Loop over all of the atom types to read the pseudopot files
  for (atyp_idx = 0; atyp_idx < ntype; atyp_idx++){

    // iatm is the unique numerical index for the atom type (atomic number). 
    long iatm  = ist->atom_types[atyp_idx];
    assign_atom_type(atype, iatm); // atype now contains the atomic symbol. Ex. if iatm = 48, atype = Cd.
    
    /* ****** ****** ****** *****/
    // Handle ligand potentials. 
    if ( (0 == strcmp(atype, "P1")) || (0 == strcmp(atype, "P2")) ||
         (0 == strcmp(atype, "P3")) || (0 == strcmp(atype, "P4")) ||
         (0 == strcmp(atype, "PC5")) || (0 == strcmp(atype, "PC6"))){
      // Get the name of the ligand potential (stored in atype)
      sprintf (str, "pot%c%c%c", atype[0], atype[1], atype[2]);
      strcat(str, ".par");
      // open the ligand potential file if it exists
      pf = fopen(str , "r");
      if (pf != NULL) {
        strcpy(req, "std");
        read_pot_file(pf, pot, atyp_idx, ist->max_pot_file_len, req);
        fclose(pf);
        // print the potential that finished reading
        printf("%s ",str);
      } else {
        // If there is no file for the ligands in the directory, construct the ligand potentials by using hardcoded 
        // Gaussian parameters a and b
        printf("\n\tNo ligand potential file %s\n\t Using default Gaussian parameters a = %lg b = %lg\n", str, a[iatm-2], b[iatm-2]);
        pot->file_lens[atyp_idx] = pot->file_lens[0];
        for (i = 0; i < pot->file_lens[atyp_idx]; i++) {
          pot->r[atyp_idx*ist->max_pot_file_len + i] = pot->r[i];
          pot->pseudo[atyp_idx*ist->max_pot_file_len + i] = (a[iatm-2] * exp(-sqr(pot->r[atyp_idx*ist->max_pot_file_len + i]) / b[iatm-2]));
          
        } 
        pot->dr[atyp_idx] = pot->r[atyp_idx*ist->max_pot_file_len + 1] - pot->r[atyp_idx*ist->max_pot_file_len];
      }
    }
    // ******* ******* ******* ******* ******* ******* ******* *******
    // Handle atomic pseudopotentials here
    // ******* ******* ******* ******* ******* ******* ******* *******
    else {
      // Handle potential reading if potentials will not be interpolated
      if (flag->interpolatePot == 0){
        // If par->scale_surface_Cs is 1.0, then we do not need to scale the surface Cs atoms
        // This is the default path.
        if (1.0 == par->scale_surface_Cs){
          // The standard format for the file names is pot[atype].dat
          sprintf(str, "pot%c%c%c", atype[0], atype[1], atype[2]);
          strcat(str, ".par");

          pf = fopen(str , "r");
          if (pf != NULL) {
            strcpy(req, "std");
            read_pot_file(pf, pot, atyp_idx, ist->max_pot_file_len, req);
            fclose(pf);
            // print the name of the potential file that was successfully read
            printf("%s ",str); 
          } else {
            fprintf(stderr, "Unable to open pot file %s to read! Exiting...\n", str); fflush(0); 
            exit(EXIT_FAILURE);
          }
        }
        // If par->scale_surface_Cs is 1.0, then we do not need to scale the surface Cs atoms
        else if ( (1.0 != par->scale_surface_Cs) && (par->scale_surface_Cs < 1.0) ){
          // We need to read in the short and long range potentials separately
          // in order to scale just the LR portion of the Cs atom potentials
          // The SR potential file name is pot[atype]_SR.dat
          sprintf (str, "pot%c%c%c", atype[0], atype[1], atype[2]);
          strcat(str, "_SR.par");
          pf = fopen(str , "r");
          if (pf != NULL) {
            strcpy(req, "std");
            read_pot_file(pf, pot, atyp_idx, ist->max_pot_file_len, req);
            printf("%s ",str); 
            i = pot->file_lens[atyp_idx];
            dr_check = pot->dr[atyp_idx];
          } else {
            fprintf(stderr, "Unable to open pot file %s to read! Exiting...\n", str); fflush(0); 
            exit(EXIT_FAILURE);
          }

          // ***** LONG RANGE ***** ***** LONG RANGE ***** ***** LONG RANGE *****
          // The LR potential file name is pot[atype]_LR.dat
          sprintf (str, "pot%c%c%c", atype[0], atype[1], atype[2]);
          strcat(str, "_LR.par");
          pf = fopen(str , "r");
          if (pf != NULL) {
            strcpy(req, "LR");
            read_pot_file(pf, pot, atyp_idx, ist->max_pot_file_len, req);
            fclose(pf);
            printf("%s ",str); 
          } else {
            fprintf(stderr, "Unable to open pot file %s to read! Exiting...\n", str); fflush(0); 
            exit(EXIT_FAILURE);
          }

          // Get the r-spacing and file length
          // check that it is the same as for the short range
          if ((dr_check != pot->dr[j]) || (i != pot->file_lens[j]) ){
            fprintf(stderr, "ERROR: short and long range pseudopotential files have different r-spacing or length\n");
            exit(EXIT_FAILURE);
          }
        } else {
          fprintf(stderr, "Invalid value of parameter scaleSurfaceCs = %lg > 1.0\n", par->scale_surface_Cs);
          exit(EXIT_FAILURE);
        }
        // Read in a4 and a5 parameters for strain dependent pseudopotential terms
        if (1 == flag->useStrain){
          // Read parameters for strain dependent pseudopotential terms
          // The standard format for the file names is pot[atype]_a[4 or 5].dat
          sprintf (str, "pot%c%c%c", atype[0], atype[1], atype[2]);
          strcat(str, "_a4.par");
          pf = fopen(str, "r");
          if (pf != NULL) {
              fscanf(pf, "%lg", &pot->a4_params[j]);
              fclose(pf);
          } else {
              pot->a4_params[j] = 0.0;
              fprintf(stderr, "\tWARNING: strain dependent pseudopotential requested; no %s file...\n", str);
              // exit(EXIT_FAILURE);
          }
          
          sprintf (str, "pot%c%c%c", atype[0], atype[1], atype[2]);
          strcat(str, "_a5.par");
          pf = fopen(str, "r");
          if (pf != NULL) {
              fscanf(pf, "%lg", &pot->a5_params[j]);
              fclose(pf);
          } else {
              pot->a5_params[j] = 0.0;
              fprintf(stderr, "\tWARNING: strain dependent pseudopotential requested; no %s file...\n", str);
              // exit(EXIT_FAILURE);
          }
        }
      } 
      // ******* ******* ******* ******* ******* ******* ******* *******
      // Handle potential reading if potentials will be interpolated 
      else if (1 == flag->interpolatePot){
        // The pseudopotentials will be interpolated between geometries
        // so we need to read in two pseudopotential files
        // the format is pot[atype]_[geom].dat where geom is the geometry indicator eg. "cubic/ortho"
        // each atom's entry is spaced apart by ngeoms*j*n, where j is the index of this atom in the array of
        // atoms in conf.par (max is n_atom_types) and n is the max length of a potential file, currently
        // hardcoded to be 8192. The geom1 values are stored at all even indices
        // 2*j*n while geom2 values are stored on all odd indices 2*j*n + 1
        // I assume ngeoms is 2. 
        // If ngeoms is larger than 2, you made a shite pseudopotential.
        if (1.0 == par->scale_surface_Cs){
          // ******* ******* ******* ******* ******* ******* *******
          // CUBIC (Or geom1)
          // ******* ******* ******* ******* ******* ******* *******
          sprintf (str, "pot%c%c%c", atype[0], atype[1], atype[2]);
          strcat(str, "_cubic.par");
          // Get the cubic pseudopotential
          pf = fopen(str , "r");
          if (pf != NULL) {
            strcpy(req, "std");
            read_pot_file(pf, pot, ngeoms*atyp_idx, ist->max_pot_file_len, req);
            fclose(pf);
            printf("%s ",str); // print the potential that finished reading
          } else {
            fprintf(stderr, "Unable to open pot file %s to read! Exiting...\n", str); fflush(0); 
            exit(EXIT_FAILURE);
          }

          // ******* ******* ******* ******* ******* ******* *******
          // ORTHORHOMBIC GEOMETRY
          // ******* ******* ******* ******* ******* ******* *******
          // get ORTHO pseudopotential
          sprintf (str, "pot%c%c%c", atype[0], atype[1], atype[2]);
          strcat(str, "_ortho.par");
          pf = fopen(tmpstr , "r");
          if (pf != NULL) {
            strcpy(req, "ortho");
            read_pot_file(pf, pot, ngeoms*atyp_idx, ist->max_pot_file_len, req);
            fclose(pf);
            printf("%s ",str); // print the potential that finished reading
          }
          else {
            fprintf(stderr, "Unable to open pot file %s to read! Exiting...\n", tmpstr); fflush(0); 
            exit(EXIT_FAILURE);
          }
        }
        // If par->scale_surface_Cs is 1.0, then we do not need to scale the surface Cs atoms
        else if ( (1.0 != par->scale_surface_Cs) && (par->scale_surface_Cs < 1.0) ){
          // ******* ******* ******* ******* ******* ******* *******
          // CUBIC (Or geom1)
          // ******* ******* ******* ******* ******* ******* *******
          // The SR potential file name is pot[atype]_[geom1]_SR.dat
          sprintf (str, "pot%c%c%c", atype[0], atype[1], atype[2]);
          strcat(str, "_cubic_SR.par");
          pf = fopen(str , "r");
          if (pf != NULL) {
            strcpy(req, "std");
            read_pot_file(pf, pot, ngeoms*atyp_idx, ist->max_pot_file_len, req);
            fclose(pf);
            printf("%s ",str); 
            i = pot->file_lens[ngeoms*atyp_idx];
            dr_check = pot->dr[ngeoms*atyp_idx];
          } else {
            fprintf(stderr, "Unable to open pot file %s to read! Exiting...\n", str); fflush(0); 
            exit(EXIT_FAILURE);
          }

          // ***** LONG RANGE ***** ***** LONG RANGE ***** ***** LONG RANGE *****
          // The LR potential file name is pot[atype]_[geom1]_LR.dat
          sprintf (str, "pot%c%c%c", atype[0], atype[1], atype[2]);
          strcat(str, "_cubic_LR.par");
          pf = fopen(str , "r");
          if (pf != NULL) {
            strcpy(req, "LR");
            read_pot_file(pf, pot, ngeoms*atyp_idx, ist->max_pot_file_len, req);
            fclose(pf);
            printf("%s ",str); 
          } else {
            fprintf(stderr, "Unable to open pot file %s to read! Exiting...\n", str); fflush(0); 
            exit(EXIT_FAILURE);
          }

          if ((dr_check != pot->dr[ngeoms*atyp_idx]) || (i != pot->file_lens[ngeoms*atyp_idx]) ){
            fprintf(stderr, "ERROR: short and long range pseudopotential files have different r-spacing or length\n");
            exit(EXIT_FAILURE);
          }
          // ******* ******* ******* ******* ******* ******* *******
          // ORTHORHOMBIC (Or geom2)
          // ******* ******* ******* ******* ******* ******* *******
          // The SR potential file name is pot[atype]_[geom2]_SR.par
          sprintf (str, "pot%c%c%c", atype[0], atype[1], atype[2]);
          strcat(str, "_ortho_SR.par");
          pf = fopen(str , "r");
          if (pf != NULL) {
            strcpy(req, "ortho");
            read_pot_file(pf, pot, ngeoms*atyp_idx, ist->max_pot_file_len, req);
            fclose(pf);
            printf("%s ",str); 
            i = pot->file_lens[ngeoms*atyp_idx+1];
            dr_check = pot->dr[ngeoms*atyp_idx+1];
          } else {
            fprintf(stderr, "Unable to open pot file %s to read! Exiting...\n", str); fflush(0); 
            exit(EXIT_FAILURE);
          }
          // ***** LONG RANGE ***** ***** LONG RANGE ***** ***** LONG RANGE *****
          // The LR potential file name is pot[atype]_[geom2]_LR.par
          sprintf (str, "pot%c%c%c", atype[0], atype[1], atype[2]);
          strcat(str, "_ortho_LR.par");
          pf = fopen(str , "r");
          if (pf != NULL) {
            strcpy(req, "orthoLR");
            read_pot_file(pf, pot, ngeoms*atyp_idx, ist->max_pot_file_len, req);
            fclose(pf);
            printf("%s ",str); 
          } else {
            fprintf(stderr, "Unable to open pot file %s to read! Exiting...\n", str); fflush(0); 
            exit(EXIT_FAILURE);
          }

          if ((dr_check != pot->dr[ngeoms*atyp_idx + 1]) || (i != pot->file_lens[ngeoms*atyp_idx + 1]) ){
            fprintf(stderr, "ERROR: short and long range pseudopotential files have different r-spacing or length\n");
            exit(EXIT_FAILURE);
          }
        }

        // Read strain dependent parameters if necessary
        if (1 == flag->useStrain){
          // CUBIC
          sprintf (str, "pot%c%c%c", atype[0], atype[1], atype[2]);
          strcat(str, "_a4_cubic.par");
          pf = fopen(str, "r");
          if (pf != NULL) {
              fscanf(pf, "%lg", &pot->a4_params[ngeoms*j]);
              fclose(pf);
          } else {
              pot->a4_params[ngeoms*j] = 0.0;
              fprintf(stderr, "\tWARNING: strain dependent pseudopotential requested; no %s file...\n", str);
              // exit(EXIT_FAILURE);
          }
          
          sprintf (str, "pot%c%c%c", atype[0], atype[1], atype[2]);
          strcat(str, "_a5_cubic.par");
          pf = fopen(str, "r");
          if (pf != NULL) {
              fscanf(pf, "%lg", &pot->a5_params[ngeoms*j]);
              fclose(pf);
          } else {
              pot->a5_params[j] = 0.0;
              fprintf(stderr, "\tWARNING: strain dependent pseudopotential requested; no %s file...\n", str);
              // exit(EXIT_FAILURE);
          }
          // ORTHORHOMBIC
          sprintf (str, "pot%c%c%c_a4", atype[0], atype[1], atype[2]);
          strcat(str, "_a4_ortho.par");
          pf = fopen(str, "r");
          if (pf != NULL) {
              fscanf(pf, "%lg", &pot->a4_params[ngeoms*j + 1]);
              fclose(pf);
          } else {
              pot->a4_params[ngeoms*j + 1] = 0.0;
              fprintf(stderr, "\tWARNING: strain dependent pseudopotential requested; no %s file...\n", str);
              // exit(EXIT_FAILURE);
          }
          
          sprintf (str, "pot%c%c%c", atype[0], atype[1], atype[2]);
          strcat(str, "_a5_ortho.par");
          if (pf != NULL) {
              fscanf(pf, "%lg", &pot->a5_params[ngeoms*j + 1]);
              fclose(pf);
          } else {
              pot->a5_params[j] = 0.0;
              fprintf(stderr, "\tWARNING: strain dependent pseudopotential requested; no %s file...\n", str);
              // exit(EXIT_FAILURE);
          }
        }
      }
    }
  }

  printf("\n");
  // ****** ****** ****** ****** ****** ****** 
  // Handle spin-orbit potentials
  // ****** ****** ****** ****** ****** ****** 
  // If the spin-orbit flag is on, then SO parameters need
  // to be read from file SO_[atomtype].par. 
  // If the interpolatePot flag is on, then parameters need to
  // be read in from two files, SO_[atomtype]_[geomtype1].par
  // and SO_[atomtype]_[geomtype2].par. For perovskites, there
  // are cubic and orthorhombic parameters, placed in files
  // SO_[atomtype]_cubic.par and SO_[atomtype]_ortho.par.
  for (i = 0; i < ist->natoms; i++){
    if (flag->SO == 1){
      // Set each atom's spin orbit parameter to zero before assigning
      atom[i].SO_par = 0.0;
      
      // Cd (no spin orbit fits exist yet)
      if (atom[i].Zval == 48) atom[i].SO_par = 0.0;
      // Se (no spin orbit fits exist yet)
      if (atom[i].Zval == 34) atom[i].SO_par = 0.0;
      if (atom[i].Zval == 49) atom[i].SO_par = 3.42974738; // In
      if (atom[i].Zval == 33) atom[i].SO_par = 0.78183193; // As
      if (atom[i].Zval == 14) atom[i].SO_par = 0.0; // Si
      if (atom[i].Zval == 1) atom[i].SO_par = 0.0; // H
      if (atom[i].Zval == 30) atom[i].SO_par = 0.0; // Zn
      if (atom[i].Zval == 16) atom[i].SO_par = 0.0; // S
      if (atom[i].Zval == 2) atom[i].SO_par = 0.0; // P1
      if (atom[i].Zval == 10) atom[i].SO_par = 0.0; // P2
      if (atom[i].Zval == 18) atom[i].SO_par = 0.0; // P3
      if (atom[i].Zval == 36) atom[i].SO_par = 0.0; // P4
      if (atom[i].Zval == 52) atom[i].SO_par = 0.0; // Te
      if (atom[i].Zval == 31) atom[i].SO_par = 0.0; // Ga
      
      // Read input file with spin orbit parameters
      // only if no interpolation requested!
      
      // Generate the first part of the SO file name
      sprintf (str, "SO_");
      for (int strnum = 0; strnum < 3; strnum++){  
        if(atom[i].atyp[strnum] == '\0') break;
        strncat(str, &(atom[i].atyp[strnum]), 1);
      }

      if (1 != flag->interpolatePot){
        // This is a job that uses spin-orbit and NL, but does not interpolate the potentials
        if (i == 0) printf("\tReading SO & non-local pot parameters\n");
        
        strcat(str, ".par");
        pf = fopen(str , "r");
        if (pf != NULL) {
          fscanf(pf, "%lg %*g", &atom[i].SO_par);
          fscanf(pf,"%lg %*g", &atom[i].NL_par[0]);
          fscanf(pf,"%lg %*g", &atom[i].NL_par[1]);
          fclose(pf);
        } else{
          printf("\nNo SO(+NL) pot found in: %s\n", str);
          fprintf(stderr, "No SO(+NL) pot found in: %s\n", str);
          exit(EXIT_FAILURE);}
      }
    }
  }
  // ****** ****** ****** ****** ****** ****** 
  // Handle potential interpolation
  // ****** ****** ****** ****** ****** ****** 
  if (1 == flag->interpolatePot){
    printf("\tInterpolating pseudopotential parameters\n");
    interpolate_pot(R, atom, ist, par);
  }

  free(req);
  return;
}

/*****************************************************************************/

void read_pot_file(FILE *pf, pot_st *pot, long j, long n, char *req){
  /*******************************************************************
  * This function reads the potential files for each atom in conf    *
  * inputs:                                                          *
  *  [FILE *pf] file stream containing the pot file                  *
  *  [pot] struct containing arrays for reading atom pseudopot files *
  *  [long j] index of the atom type (0 to ngeoms * n_atom_types)    *
  *  [long n] max length of a potential file (8192)                  *
  *  [char *req] the requested type of pot file to read              *
  *  [flag_st *flag] ptr to flag_st holding job flags                *
  * outputs: void                                                    *
  ********************************************************************/
  long iscan;
  char *std, *LR, *orthoLR, *ortho;
  std = malloc(4 * sizeof(std[0])); strcpy(std, "std"); 
  LR = malloc(2 * sizeof(std[0])); strcpy(LR, "LR");
  ortho = malloc(6 * sizeof(ortho[0])); strcpy(ortho, "ortho");
  orthoLR = malloc(8 * sizeof(ortho[0])); strcpy(orthoLR, "orthoLR");

  if (0 == strcmp(req, (const char *)std)){
    // Loop until the end of the file pot[atype].par and read the contents 
    // Increment pot->file_lens[j] to count how long the pseudopotential file is
    pot->file_lens[j] = 0;
    for (iscan = 0; iscan != EOF; pot->file_lens[j]++){
      iscan = fscanf(pf, "%lg %lg", &pot->r[j*n + pot->file_lens[j]], &pot->pseudo[j*n + pot->file_lens[j]]);
      // pot->r and pot->pseudo store all the r values and v(r) values for each atomic pseudopotential
      // Each atom's entry is spaced apart by j*n, where j is the index of this atom in the array of
      // atoms in conf.par (max is n_atom_types) and n is the max length of a potential file, currently
      // hardcoded to be 8192.
    }
    pot->file_lens[j]--;
    // The last iteration returns the EOF signal, but still increments pot->file_lens[j]. 
    // Decrement it now to avoid a seg fault by running outside the pot->r/pot->pseudo array

    // Get the r-spacing and add to the pot->dr array 
    // only evenly-spaced pseudopotential files are supported currently
    pot->dr[j] = pot->r[j*n+1] - pot->r[j*n+0]; // r spacing

  } else if ((0 == strcmp(req, (const char *)LR))){
    for (pot->file_lens[j] = iscan = 0; iscan != EOF; pot->file_lens[j]++) {
      iscan = fscanf (pf,"%lg %lg", &pot->r_LR[j*n + pot->file_lens[j]], &pot->pseudo_LR[j*n + pot->file_lens[j]]);
    }
    pot->file_lens[j]--;
    pot->dr[j] = pot->r_LR[j*n+1] - pot->r_LR[j*n+0];

  } else if ((0 == strcmp(req, (const char *)ortho))){
    for (pot->file_lens[j+1] = iscan = 0; iscan != EOF; pot->file_lens[j+1]++) {
      iscan = fscanf (pf,"%lg %lg", &pot->r[j*n + n + pot->file_lens[j+1]], &pot->pseudo[j*n + n + pot->file_lens[j+1]]);
    }
    pot->file_lens[j+1]--;
    
    pot->dr[j+1] = pot->r[j*n+n+1] - pot->r[j*n+n+0];

  } else if ((0 == strcmp(req, (const char *)orthoLR))){
    for (pot->file_lens[j+1] = iscan = 0; iscan != EOF; pot->file_lens[j+1]++) {
      iscan = fscanf (pf,"%lg %lg", &pot->r_LR[j*n + n + pot->file_lens[j+1]], &pot->pseudo_LR[j*n + n + pot->file_lens[j+1]]);
    }
    pot->file_lens[j+1]--;
    pot->dr[j+1] = pot->r_LR[j*n+n+1] - pot->r_LR[j*n+n+0];
  }

  // Free memory for all strings
  free(std); free(LR); free(ortho); free(orthoLR);

  return;
}

void interpolate_pot(xyz_st *R, atom_info *atom, index_st *ist, par_st *par){
  /*******************************************************************
  * If there are multiple geometries for an NC (eg. cubic/ortho),    *
  * then this function calcs the potential felt by an atom at some   *
  * intermediate geometry as an interpolation between the fitt'd pots*
  * at each geometry.                                                *
  * inputs:                                                          *
  *  [xyz_st *R] pointer to double arrays containing atom coords     *
  *  [atom_info *atom] ptr to struct containing atom-specific params *
  *  [index_st *ist] pointer to counters, indices, and lengths       *
  *  [par_st *par] ptr to par_st holding VBmin, VBmax... params      *
  * outputs: void                                                    *
  ********************************************************************/

  FILE *pw;
  coeff_st *coeff;
  int atm_id, i;
  char tmpstr[100], str[100], atyp[3];;

  // Calculate the geometry parameters of the configuration
  calc_geom_par(R, atom, ist);
  
  // Make containers for the cubic and orthorhombic parameters:
  // This array only needs ist->n_atom_types number of elements. 
  // However, the next loop is easier if we can place the coefficients
  // in the element corresponding to the atom.Zval, which goes up to at least 82.
  coeff = (coeff_st*) calloc(100, sizeof(coeff_st)); 

  for (i = 0; i < ist->n_atom_types; i++){
      // the ith element of atom_types a unique Zval
    atm_id = ist->atom_types[i];
    // Make the first part of the file name for this atom type
    assign_atom_type(atyp, atm_id);
    sprintf(str, "SO_");
    strcat(str, atyp);
    strcpy(tmpstr, str); // copy the base name including the atom type
    
    // Get the orthorhombic parameters
    strcat(str, "_ortho.par");
    pw = fopen(str , "r");
    if (pw != NULL) {
      // Place the parameters into the atm_id'th index of coeff
      fscanf(pw, "%lg %*g", &coeff[atm_id].SO[0]); 
      fscanf(pw,"%lg %*g", &coeff[atm_id].NL1[0]);
      fscanf(pw,"%lg %*g", &coeff[atm_id].NL2[0]);
      fclose(pw);
    } else{
      printf("\nNo SO(+NL) pot found in: %s\n", str);
      fprintf(stderr, "No SO(+NL) pot found in: %s\n", str);
      exit(EXIT_FAILURE);
    }
    // Get the cubic parameters
    strcat(tmpstr, "_cubic.par");
    pw = fopen(tmpstr , "r");
    if (pw != NULL) {
      // Place the parameters into the atm_id'th index of coeff
      fscanf(pw, "%lg %*g", &coeff[atm_id].SO[1]); 
      fscanf(pw,"%lg %*g", &coeff[atm_id].NL1[1]);
      fscanf(pw,"%lg %*g", &coeff[atm_id].NL2[1]);
      fclose(pw);
    } else{
      printf("\nNo SO(+NL) pot found in: %s\n", str);
      fprintf(stderr, "No SO(+NL) pot found in: %s\n", str);
      exit(EXIT_FAILURE);
    }
  }

  // Loop through all of the atom_infos and update their parameters based on the geometry

  for (i = 0; i < ist->natoms; i++){
    atm_id = atom[i].Zval;
    atom[i].SO_par = atom[i].geom_par*coeff[atm_id].SO[0] + (1.0-atom[i].geom_par)*coeff[atm_id].SO[1];
    atom[i].NL_par[0] = atom[i].geom_par*coeff[atm_id].NL1[0] + (1.0-atom[i].geom_par)*coeff[atm_id].NL1[1];
    atom[i].NL_par[1] = atom[i].geom_par*coeff[atm_id].NL2[0] + (1.0-atom[i].geom_par)*coeff[atm_id].NL2[1];
  }

  return;
}

/*****************************************************************************/
void calc_geom_par(xyz_st *R, atom_info *atom, index_st *ist){
  /*******************************************************************
  * This function calculates the local geometry between cubic/ortho  *
  * for interpolating perovskite potentials                          *
  * inputs:                                                          *
  *  [xyz_st *R] pointer to double arrays containing atom coords     *
  *  [atom_info *atom] ptr to struct containing atom-specific params *
  *  [index_st *ist] pointer to counters, indices, and lengths       *
  * outputs: void                                                    *
  ********************************************************************/

  long j;
  double dist_sqr, bond_angle, avg_I_par;
  long bonded[10];
  int n_bonded;

  for (long i = 0; i < ist->natoms; i++){
  
  // Check what type of atom this is. (Currently only support Pb, I, Cs)  
    switch (atom[i].Zval){

      // ****** ****** ****** ****** ****** ****** ****** ****** 
      // Handle I atoms
      case 53:
      { 
        n_bonded = 0;

        //get list of bonded lead atoms
        for(j = 0; j < ist->natoms; j++){
          if (atom[j].Zval == 82){
            dist_sqr = sqr(R[i].x - R[j].x) + sqr(R[i].y - R[j].y) + sqr(R[i].z - R[j].z);
            if (dist_sqr < sqr(PbIBondMax)){
              bonded[n_bonded] = j;
              n_bonded++;
            }
          }
        }

        // Calculate the bond angle for I atoms
        if(n_bonded != 2){
          // This is an edge I atom
          atom[i].geom_par = 0.0; 
          //default to cubic if edge atom
          //printf("atom %ld (I) is an edge atom (n_bonded=%d)\n",i, n_bonded);
        } 
        else{
          bond_angle = calc_bond_angle(bonded[0], i, bonded[1], R);
          //printf("atom %ld (I) has bond angle %g\n",i, bond_angle);

          if (bond_angle < orthoBondAngle) {
            atom[i].geom_par=1.0;
          }
          if (bond_angle > minCubicBondAngle){
            atom[i].geom_par = 0.0;
          }
          else{
            atom[i].geom_par = 1.0 - ((bond_angle - orthoBondAngle) / (minCubicBondAngle - orthoBondAngle));
          } 
        }

        break;
      } 
      

      // ****** ****** ****** ****** ****** ****** ****** ****** 
      // Handle Cs atoms
      case 55: 
      { 
        atom[i].geom_par = 0.0; // This might need to change now that Cs potentials are nonzero

        break;
      }
      
      // ****** ****** ****** ****** ****** ****** ****** ****** 
      // Handle Pb atoms later after all geom_pars are set for I
      case 82:
        break;

      default:
      {
        fprintf(stderr, "Unknown atom with Zval %d!\nExiting...", atom[i].Zval);
        fflush(0);
        exit(EXIT_FAILURE);
      } 
    }
  }
    
  // ****** ****** ****** ****** ****** ****** ****** ****** 
  // Handle Pb atoms
  for (long i = 0; i < ist->natoms; i++){
    
    switch (atom[i].Zval){
      case 53: 
        break;
      
      case 55:
        break;

      case 82: 
      {
        n_bonded = 0;

        // Get list of bonded I atoms
        for(j = 0; j < ist->natoms; j++){
          if (atom[j].Zval == 53){
            dist_sqr = sqr(R[i].x - R[j].x) + sqr(R[i].y - R[j].y) + sqr(R[i].z - R[j].z);
            if (dist_sqr < sqr(PbIBondMax)){
              bonded[n_bonded] = j;
              n_bonded++;
            }
          }
        }

        // Define the geom_par for this Pb as the average of the bonded I atoms
        if(n_bonded != 6){
          // default to cubic if edge atom
          atom[i].geom_par = 0.0;  
          //printf("atom %ld (Pb) is an edge atom (n_bonded=%d)\n",i, n_bonded);
        } else{
          avg_I_par = 0.0;
          for(j = 0;j < 6; j++){
            avg_I_par += atom[bonded[j]].geom_par;
          }
        
          avg_I_par /= 6.0;
          atom[i].geom_par = avg_I_par;
        }
      }

      // If an unknown atom is encountered that was not caught before
      default:
        fprintf(stderr, "Unknown atom with Zval %d!\nExiting...", atom[i].Zval);
        fflush(0);
        exit(EXIT_FAILURE);
    
    }
  }

  return;
}

/*****************************************************************************/
double calc_bond_angle(long index1, long index2, long index3, xyz_st *R){
  /*******************************************************************
  * This function calculates the angle between three atoms           *
  * inputs:                                                          *
  *  [index1] the index of atom 1 in the R array                     *
  *  [index2] the index of atom 2 in the R array                     *
  *  [index3] the index of atom 3 in the R array                     *
  *  [xyz_st *R] pointer to double arrays containing atom coords     *
  * outputs: [double] bond angle in degrees                          *
  ********************************************************************/

  double dx1 = R[index1].x-R[index2].x;
  double dy1 = R[index1].y-R[index2].y;
  double dz1 = R[index1].z-R[index2].z; 
  double l1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);

  double dx2 = R[index3].x-R[index2].x;
  double dy2 = R[index3].y-R[index2].y;
  double dz2 = R[index3].z-R[index2].z;
  double l2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2); 

  double dot = dx1*dx2 + dy1*dy2 + dz1*dz2;
  if (-1.0>dot/(l1*l2)||1.0<dot/(l1*l2)){return 180.0;}  
  if (-1.0>dot/(l1*l2)||1.0<dot/(l1*l2)){
    printf("Error in Bond Angle (%lf)!!\n", dot/(l1*l2));
    printf("1: %f %f %f\n", R[index1].x, R[index1].y, R[index1].z );
    printf("2: %f %f %f\n", R[index2].x, R[index2].y, R[index2].z );
    printf("3: %f %f %f\n", R[index3].x, R[index3].y, R[index3].z );
  }

  return acos(dot/(l1*l2))*180.0/PIE;

}

/*****************************************************************************/

long assign_atom_number(char atyp[4]){
  /*******************************************************************
  * This function assigns the atomic number (ID) based on name       *
  * inputs:                                                          *
  *  [char[3]] name of atom as string literal (e.g. "Cd")            *
  * outputs: [long] atom ID (Zval or fictitious for ligands)         *
  ********************************************************************/
  
  if ((atyp[0] == 'H') && (atyp[1] == '\0')  && (atyp[2] == '\0')) return 1;
  else if ((atyp[0] == 'P') && (atyp[1] == '1')  && (atyp[2] == '\0')) return 2;
  else if ((atyp[0] == 'P') && (atyp[1] == '2')  && (atyp[2] == '\0')) return 3;
  else if ((atyp[0] == 'P') && (atyp[1] == '3')  && (atyp[2] == '\0')) return 4;
  else if ((atyp[0] == 'P') && (atyp[1] == '4')  && (atyp[2] == '\0')) return 5;
  else if ((atyp[0] == 'P') && (atyp[1] == 'C')  && (atyp[2] == '5'))  return 6;
  else if ((atyp[0] == 'P') && (atyp[1] == 'C')  && (atyp[2] == '6'))  return 7;
  else if ((atyp[0] == 'S') && (atyp[1] == 'i')  && (atyp[2] == '\0')) return 14;
  else if ((atyp[0] == 'P') && (atyp[1] == '\0')  && (atyp[2] == '\0')) return 15;
  else if ((atyp[0] == 'S') && (atyp[1] == '\0')  && (atyp[2] == '\0')) return 16;
  else if ((atyp[0] == 'Z') && (atyp[1] == 'n')  && (atyp[2] == '\0')) return 30;
  else if ((atyp[0] == 'G') && (atyp[1] == 'a')  && (atyp[2] == '\0')) return 31;
  else if ((atyp[0] == 'A') && (atyp[1] == 's')  && (atyp[2] == '\0')) return 33;
  else if ((atyp[0] == 'S') && (atyp[1] == 'e')  && (atyp[2] == '\0')) return 34;
  else if ((atyp[0] == 'C') && (atyp[1] == 'd')  && (atyp[2] == '\0')) return 48;
  else if ((atyp[0] == 'I') && (atyp[1] == 'n')  && (atyp[2] == '\0')) return 49;
  else if ((atyp[0] == 'T') && (atyp[1] == 'e')  && (atyp[2] == '\0')) return 52;
  else if ((atyp[0] == 'I') && (atyp[1] == '\0')  && (atyp[2] == '\0')) return 53;
  else if ((atyp[0] == 'C') && (atyp[1] == 's')  && (atyp[2] == '\0')) return 55;
  else if ((atyp[0] == 'P') && (atyp[1] == 'b')  && (atyp[2] == '\0')) return 82;
  else {
    fprintf(stderr, "atom type %s not in current list", atyp);
    exit(EXIT_FAILURE);
  }

  return (0);
}

/*****************************************************************************/

void assign_atom_type(char *atyp, long j){
  /*******************************************************************
  * This function assigns the atomic name based on number (ID)       *
  * inputs:                                                          *
  *  [atyp] ptr to string literal which will become atom name        *
  *  [j] atom number (ID)                                            *
  * outputs: void                                                    *
  ********************************************************************/

  if (j == 1) {atyp[0] = 'H'; atyp[1] = '\0'; atyp[2] = '\0';}
  else if (j == 2) {atyp[0] = 'P'; atyp[1] = '1'; atyp[2] = '\0';}
  else if (j == 3) {atyp[0] = 'P'; atyp[1] = '2'; atyp[2] = '\0';}
  else if (j == 4) {atyp[0] = 'P'; atyp[1] = '3'; atyp[2] = '\0';}
  else if (j == 5) {atyp[0] = 'P'; atyp[1] = '4'; atyp[2] = '\0';}
  else if (j == 6) {atyp[0] = 'P'; atyp[1] = 'C'; atyp[2] = '5'; atyp[3] = '\0';}
  else if (j == 7) {atyp[0] = 'P'; atyp[1] = 'C'; atyp[2] = '6'; atyp[3] = '\0';}
  else if (j == 14) {atyp[0] = 'S'; atyp[1] = 'i'; atyp[2] = '\0';}
  else if (j == 15) {atyp[0] = 'P'; atyp[1] = '\0'; atyp[2] = '\0';}
  else if (j == 16) {atyp[0] = 'S'; atyp[1] = '\0'; atyp[2] = '\0';}
  else if (j == 30) {atyp[0] = 'Z'; atyp[1] = 'n'; atyp[2] = '\0';}
  else if (j == 31) {atyp[0] = 'G'; atyp[1] = 'a'; atyp[2] = '\0';}
  else if (j == 33) {atyp[0] = 'A'; atyp[1] = 's'; atyp[2] = '\0';}
  else if (j == 34) {atyp[0] = 'S'; atyp[1] = 'e'; atyp[2] = '\0';}
  else if (j == 48) {atyp[0] = 'C'; atyp[1] = 'd'; atyp[2] = '\0';}
  else if (j == 49) {atyp[0] = 'I'; atyp[1] = 'n'; atyp[2] = '\0';}
  else if (j == 52) {atyp[0] = 'T'; atyp[1] = 'e'; atyp[2] = '\0';}
  else if (j == 53) {atyp[0] = 'I'; atyp[1] = '\0'; atyp[2] = '\0';}
  else if (j == 55) {atyp[0] = 'C'; atyp[1] = 's'; atyp[2] = '\0';}
  else if (j == 82) {atyp[0] = 'P'; atyp[1] = 'b'; atyp[2] = '\0';}
  else {
    fprintf(stderr, "atomic number %ld not in current list\n", j);
    exit(EXIT_FAILURE);
  }
  return;
}

/*****************************************************************************/

long get_number_of_atom_types(atom_info *atom, index_st *ist, long *list){
  /*******************************************************************
  * This function counts the number of atom types (deprecated)       *
  * outputs: [long] number of atom types                             *
  ********************************************************************/

  long jatom, jlist, nlist, flag;

  for (jlist = 0; jlist < 200; jlist++) list[jlist] = -1;
  list[0] = 8; list[1] = 9; list[2] = 10; list[3] = 11;
  for (nlist = 4, jatom = 0; jatom < ist->natoms; jatom++){
    for (flag = 0, jlist = 0; jlist < nlist; jlist++)
      if (atom[jatom].idx == list[jlist]) flag = 1;
    if (flag != 1) {list[nlist] = atom[jatom].idx; nlist++;}
  }
  //for (jlist = 0; jlist < nlist-4; jlist++) printf ("atom in list %ld\n",list[jlist+4]);
  return(nlist-4);
}

/*****************************************************************************/

int assign_crystal_structure(char *crystal_structure){
  /*******************************************************************
  * This function assigns the crystal structure identity             *
  * inputs:                                                          *
  *  [crystal_structure] string of structure name (e.g. "wurtzite")  *
  * outputs: [long] crystal_structure_int                            *
  ********************************************************************/
  char *wurtzite; wurtzite = malloc(10*sizeof(wurtzite[0]));
  strcpy(wurtzite, "wurtzite");
  char *zincblende; zincblende = malloc(11*sizeof(zincblende[0]));
  strcpy(zincblende, "zincblende");

  if ( 0 == strcmp(crystal_structure, (const char *) wurtzite) ) return 0;
  else if ( 0 == strcmp(crystal_structure, (const char *) zincblende)  ) return 1;
  else {
    fprintf(stderr, "Crystal structure %s not recognized", crystal_structure);
    exit(EXIT_FAILURE);
  }

  return (0);
}

/*****************************************************************************/

int assign_outmost_material(char *outmost_material){
  /*******************************************************************
  * This function assigns the material identity of the outer layer   *
  * inputs:                                                          *
  *  [outmost_material] name of material as string (e.g. "CdS")      *
  * outputs: [long] atom ID (Zval or fictitious for ligands)         *
  ********************************************************************/

  if (! strcmp(outmost_material, "CdS")) {  return 0;}
	else if (! strcmp(outmost_material, "CdSe")) {  return 1;}
	else if (! strcmp(outmost_material, "InP")) {  return 2;}
	else if (! strcmp(outmost_material, "InAs")) {  return 3;}
	else if (! strcmp(outmost_material, "alloyInGaP")) {  return 4;} //cation terminated surface only
	else if (! strcmp(outmost_material, "alloyInGaAs")) {  return 5;} // cation terminated surface only
	else if (! strcmp(outmost_material, "GaAs")) {  return 6;}
	else if (! strcmp(outmost_material, "ZnSe")) {   return 7;}
	else if (! strcmp(outmost_material, "ZnS")) {   return 8;}
	else if (! strcmp(outmost_material, "alloyGaAsP")) {  return 9;} // cation terminated surface only
	else if (! strcmp(outmost_material, "alloyGaAsSb")) {  return 10; } // cation terminated surface only
  else if (! strcmp(outmost_material, "GaP")) {  return 11;}
	else {
		fprintf(stderr, "\n\nOutmost material %s not recognized -- the program is exiting!!!\n\n", outmost_material);
		fflush(stdout);
		exit(EXIT_FAILURE);
	}

  return (0);
}

/*****************************************************************************/

double ret_ideal_bond_len(long natyp_1, long natyp_2, int crystal_structure_int){
	
  if (((natyp_1==48) && (natyp_2==34)) && (crystal_structure_int==0)) return(2.6326);  // CdSe, wz
  else if (((natyp_1==34) && (natyp_2==48)) && (crystal_structure_int==0)) return(2.6326); 
  else if (((natyp_1==48) && (natyp_2==34)) && (crystal_structure_int==1)) return(2.6233); // CdSe, zb
  else if (((natyp_1==34) && (natyp_2==48)) && (crystal_structure_int==1)) return(2.6233); 
  else if (((natyp_1==48) && (natyp_2==16)) && (crystal_structure_int==0)) return(2.5292); // CdS, wz
  else if (((natyp_1==16) && (natyp_2==48)) && (crystal_structure_int==0)) return(2.5292); 
  else if (((natyp_1==48) && (natyp_2==16)) && (crystal_structure_int==1)) return(2.5193); // CdS, zb
  else if (((natyp_1==16) && (natyp_2==48)) && (crystal_structure_int==1)) return(2.5193); 
  else if (((natyp_1==49) && (natyp_2==33)) && (crystal_structure_int==1)) return(2.6228); // InAs, zb
  else if (((natyp_1==33) && (natyp_2==49)) && (crystal_structure_int==1)) return(2.6228); 
  else if (((natyp_1==49) && (natyp_2==15)) && (crystal_structure_int==1)) return(2.5228); // InP, zb
  else if (((natyp_1==15) && (natyp_2==49)) && (crystal_structure_int==1)) return(2.5228); 
  else if (((natyp_1==31) && (natyp_2==33)) && (crystal_structure_int==1)) return(2.4480); // GaAs, zb
  else if (((natyp_1==33) && (natyp_2==31)) && (crystal_structure_int==1)) return(2.4480); 
  else if (((natyp_1==31) && (natyp_2==15)) && (crystal_structure_int==1)) return(2.360); // GaP, zb
  else if (((natyp_1==15) && (natyp_2==31)) && (crystal_structure_int==1)) return(2.360); 
  else if (((natyp_1==30) && (natyp_2==34)) && (crystal_structure_int==1)) return(2.45); // ZnSe, zb
  else if (((natyp_1==34) && (natyp_2==30)) && (crystal_structure_int==1)) return(2.45); 
  else if (((natyp_1==30) && (natyp_2==16)) && (crystal_structure_int==1)) return(2.33); // ZnS, zb
  else if (((natyp_1==16) && (natyp_2==30)) && (crystal_structure_int==1)) return(2.33); 
  else {
    fprintf(stderr, "Atom pair type %ld %ld with crystal_structure_int %d not in current list of bond lengths.\n", natyp_1, natyp_2, crystal_structure_int);
    exit(EXIT_FAILURE);
  }
  return(0);
}

/****************************************************************************/

