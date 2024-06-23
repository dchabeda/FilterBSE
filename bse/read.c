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

  // Filter algorithm parameters
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
  flag->calcSpinAngStat = 0; // Are angular momentum statistics computed. Only available with SO coupling
  flag->timingSpecs = 0; // print timing info for computing the Hamiltonian
  par->fermi_E = -0.18; // This default value is not good for LR potentials.
  par->sigma_E_cut = 0.0001;
  flag->saveCheckpoints = 0; // by default do not save job checkpoints
  flag->restartFromCheckpoint = -1; // by default do not restart from checkpoint
  
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
      } else if (!strcmp(field, "setSeed")) {
          flag->setSeed = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
          if (flag->setSeed == 1){fscanf(pf, "%ld", &par->rand_seed);}
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
          flag->timeSpecs = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
      } else if (!strcmp(field, "printFPDensity")) {
          flag->printFPDensity = (int) strtol(tmp, &endptr, 10);
          if (*endptr != '\0') {fprintf(stderr, "Error converting string to long.\n"); exit(EXIT_FAILURE);}
          if (flag->printFPDensity == 1){fscanf(pf, "%d", &ist->n_FP_density);}
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
          printf("calcSpinAngStat = int (calculate spin and ang. mom. statistics for eigenstates)\n");
          printf("setSeed = int (if 1, set the random seed in for filter to generate exactly reproducible wavefunctions)\n");
          printf("If setSeed = 1, the next entry MUST specify the random seed as an integer \'rand_seed\'\n");
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

  // Set grid defaults
  ist->nthreads = parallel->nthreads;
  ist->complex_idx = flag->isComplex + 1;

  // Get the number of atoms (needed to initialize the R array)
  pf = fopen("conf.dat" , "r");
  if (pf) {
    fscanf(pf, "%ld", &ist->natoms); 
    fclose(pf);
    } else {
    fprintf(stderr, "ERROR: conf.dat must be in current working directory\n");
    exit(EXIT_FAILURE);
  } 

  // Using input parameters, print the current job state
  print_input_state(stdout, flag, grid, par, ist, parallel);

  return;
}

/****************************************************************************/

long assign_atom_number(char atyp[3])
{
  char strerror[100];
  
  if ((atyp[0] == 'C') && (atyp[1] == 'd')  && (atyp[2] == '\0')) return(0);
  else if ((atyp[0] == 'S') && (atyp[1] == 'e') && (atyp[2] == '\0')) return(1);
  else if ((atyp[0] == 'I') && (atyp[1] == 'n') && (atyp[2] == '\0')) return(2);
  else if ((atyp[0] == 'A') && (atyp[1] == 's') && (atyp[2] == '\0')) return(3);
  else if ((atyp[0] == 'S') && (atyp[1] == 'i') && (atyp[2] == '\0')) return(4);
  else if ((atyp[0] == 'H') && (atyp[1] == '\0') && (atyp[2] == '\0'))  return(5);
  else if ((atyp[0] == 'Z') && (atyp[1] == 'n') && (atyp[2] == '\0'))  return(6);
  else if ((atyp[0] == 'S') && (atyp[1] == '\0') && (atyp[2] == '\0'))  return(7);
  else if ((atyp[0] == 'P') && (atyp[1] == '1') && (atyp[2] == '\0'))  return(8);
  else if ((atyp[0] == 'P') && (atyp[1] == '2') && (atyp[2] == '\0'))  return(9);
  else if ((atyp[0] == 'P') && (atyp[1] == '3') && (atyp[2] == '\0'))  return(10);
  else if ((atyp[0] == 'P') && (atyp[1] == '4') && (atyp[2] == '\0'))  return(11);
  else if ((atyp[0] == 'T') && (atyp[1] == 'e') && (atyp[2] == '\0')) return(12);
  else if ((atyp[0] == 'C') && (atyp[1] == 'd') && (atyp[2] == 'z')) return(13);
  else if ((atyp[0] == 'S') && (atyp[1] == 'e') && (atyp[2] == 'z')) return(14);
  else if ((atyp[0] == 'G') && (atyp[1] == 'a') && (atyp[2] == '\0')) return(15);
  else if ((atyp[0] == 'I') && (atyp[1] == '\0') && (atyp[2] == '\0')) return(53);
  else if ((atyp[0] == 'C') && (atyp[1] == 's') && (atyp[2] == '\0')) return(55);
  else if ((atyp[0] == 'P') && (atyp[1] == 'b') && (atyp[2] == '\0')) return(82);


  else {
    sprintf (strerror,"atom type %s not in current list",atyp);
    nerror (strerror);
  }
  return(0);
}

/****************************************************************************/

void assign_atom_type(char *atyp,long j)
{
  if (j == 0) {atyp[0] = 'C'; atyp[1] = 'd'; atyp[2] = '\0';}
  else if (j == 1) {atyp[0] = 'S'; atyp[1] = 'e'; atyp[2] = '\0';}
  else if (j == 2) {atyp[0] = 'I'; atyp[1] = 'n'; atyp[2] = '\0';}
  else if (j == 3) {atyp[0] = 'A'; atyp[1] = 's'; atyp[2] = '\0';}
  else if (j == 4) {atyp[0] = 'S'; atyp[1] = 'i'; atyp[2] = '\0';}
  else if (j == 5) {atyp[0] = 'H'; atyp[1] = '\0'; atyp[2] = '\0';}
  else if (j == 6) {atyp[0] = 'Z'; atyp[1] = 'n'; atyp[2] = '\0';}
  else if (j == 7) {atyp[0] = 'S'; atyp[1] = '\0'; atyp[2] = '\0';}
  else if (j == 8) {atyp[0] = 'P'; atyp[1] = '1'; atyp[2] = '\0';}
  else if (j == 9) {atyp[0] = 'P'; atyp[1] = '2'; atyp[2] = '\0';}
  else if (j == 10) {atyp[0] = 'P'; atyp[1] = '3'; atyp[2] = '\0';}
  else if (j == 11) {atyp[0] = 'P'; atyp[1] = '4'; atyp[2] = '\0';}
  else if (j == 12) {atyp[0] = 'T'; atyp[1] = 'e'; atyp[2] = '\0';}
  else if (j == 13) {atyp[0] = 'C'; atyp[1] = 'd'; atyp[2] = 'z';}
  else if (j == 14) {atyp[0] = 'S'; atyp[1] = 'e'; atyp[2] = 'z';}
  else if (j == 15) {atyp[0] = 'G'; atyp[1] = 'a'; atyp[2] = '\0';}
  else if (j == 53) {atyp[0] = 'I'; atyp[1] = '\0'; atyp[2] = '\0';}
  else if (j == 55) {atyp[0] = 'C'; atyp[1] = 's'; atyp[2] = '\0';}
  else if (j == 82) {atyp[0] = 'P'; atyp[1] = 'b'; atyp[2] = '\0';}


  return;
}

/****************************************************************************/

void read_bseCoeff(int ms2,int numExcStatesToRead, zomplex* u, FILE* pf){
  int i,j;

  if(numExcStatesToRead>100){numExcStatesToRead=100;}

  for (i = 0; i < ms2; i++) {
    for (j = 0; j < numExcStatesToRead; j++) {  
    fscanf (pf,"{%lg, %lg}\t", &(u[i*ms2+j].re), &(u[i*ms2+j].im));
    } 
  }
}
