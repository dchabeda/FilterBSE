/****************************************************************************/

#include "read.h"

/*****************************************************************************/

void read_input(
    flag_st *flag,
    grid_st *grid,
    index_st *ist,
    par_st *par,
    parallel_st *parallel)
{

  /*******************************************************************
   * This function reads the input.par file and initializes the job   *
   * inputs:                                                          *
   *  [flag_st *par] ptr to flag_st holding job flags                 *
   *  [grid_st *grid] pointer to grid struct, holds the grid          *
   *  [index_st *ist] pointer to counters, indices, and lengths       *
   *  [par_st *par] ptr to par_st holding VBmin, VBmax... params      *
   *  [parallel] struct holding parallelization options (for norm)    *
   *  outputs: void                                                   *
   ********************************************************************/

  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  FILE *pf;

  const int mpir = parallel->mpi_rank;

  int i = 0;

  char field[1000] = {0};
  char tmp[1000] = {0};

  /************************************************************/
  /*******************     SET DEFAULTS     *******************/
  /************************************************************/

  //                           Quasiparticle basis params
  ist->max_hole_states = -1;
  ist->max_elec_states = -1;
  par->KE_max = 20.0; // NOTE: This is Filter default
  par->checkpoint_id = 0;

  //                           Pseudopotential parameters
  flag->LR = 0; // By default, pseudopotentials are short ranged.

  //                           Spin-orbit and non-local terms
  flag->useSpinors = 0; // default is to use non-spinor wavefunctions
  flag->isComplex = 0;
  ist->nspin = 1; // turning on the useSpinors flag will set nspin to 2.
  flag->SO = 0;   // computes the spin-orbit terms in the Hamiltonian
  flag->NL = 0;   // automatically on if SO flag on
  flag->restartCoulomb = 0;
  flag->coulombDone = 0;
  flag->calcCoulombOnly = 0;
  flag->noCalcExciton = 0;

  //                           Optional output flags
  flag->calcDarkStates = 0;
  flag->calcSpinAngStat = 1; // By default compute angular momentum statistics computed. Only available with SO coupling
  flag->timingSpecs = 0;     // print timing info for computing the Hamiltonian
  par->fermi_E = -0.18;      // This default value is not good for LR potentials.
  par->sigma_E_cut = 0.0001;
  flag->saveCheckpoints = 0; // by default do not save job checkpoints
  flag->restartFromChk = -1; // by default do not restart from checkpoint
  flag->printFPDensity = 0;
  ist->n_FP_density = 0;

  /************************************************************/
  /******************    PARSE INPUT FILE    ******************/
  /************************************************************/

  if (access("input.par", F_OK) != -1)
  {
    int fk; // "found key"

    pf = fopen("input.par", "r");

    while (fscanf(pf, "%s = %s", field, tmp) != EOF && i < 100)
    {
      fk = 0;
      // ****** ****** ****** ****** ****** ******
      // Set parameters&counters for BSE algorithm
      // ****** ****** ****** ****** ****** ******
      read_field(field, "maxHoleStates", &ist->max_hole_states, LONG_TYPE, tmp, &fk);
      read_field(field, "maxElecStates", &ist->max_elec_states, LONG_TYPE, tmp, &fk);
      read_field(field, "sigmaECut", &par->sigma_E_cut, DOUBLE_TYPE, tmp, &fk);
      read_field(field, "deltaEhole", &par->delta_E_hole, DOUBLE_TYPE, tmp, &fk);
      read_field(field, "deltaEelec", &par->delta_E_elec, DOUBLE_TYPE, tmp, &fk);
      read_field(field, "KEmax", &par->KE_max, DOUBLE_TYPE, tmp, &fk);
      read_field(field, "fermiEnergy", &par->fermi_E, DOUBLE_TYPE, tmp, &fk);

      // ****** ****** ****** ****** ****** ******
      // Set options for dielectric properties
      // ****** ****** ****** ****** ****** ******
      read_field(field, "epsX", &par->epsX, DOUBLE_TYPE, tmp, &fk);
      read_field(field, "epsY", &par->epsY, DOUBLE_TYPE, tmp, &fk);
      read_field(field, "epsZ", &par->epsZ, DOUBLE_TYPE, tmp, &fk);
      read_field(field, "longRange", &flag->LR, INT_TYPE, tmp, &fk);

      // ****** ****** ****** ****** ****** ******
      // Set options for parallelization
      // ****** ****** ****** ****** ****** ******
      read_field(field, "nThreads", &parallel->nthreads, LONG_TYPE, tmp, &fk);

      // ****** ****** ****** ****** ****** ******
      // Set options for spin-orbit calculation
      // ****** ****** ****** ****** ****** ******
      read_field(field, "useSpinors", &flag->useSpinors, INT_TYPE, tmp, &fk);
      read_field(field, "spinOrbit", &flag->SO, INT_TYPE, tmp, &fk);
      read_field(field, "NonLocal", &flag->NL, INT_TYPE, tmp, &fk);

      // ****** ****** ****** ****** ****** ******
      // Set options for additional output
      // ****** ****** ****** ****** ****** ******
      read_field(field, "timingSpecs", &flag->timingSpecs, INT_TYPE, tmp, &fk);
      read_field(field, "printFPDensity", &flag->printFPDensity, INT_TYPE, tmp, &fk);
      read_field(field, "calcDarkStates", &flag->calcDarkStates, INT_TYPE, tmp, &fk);
      read_field(field, "calcSpinAngStat", &flag->calcSpinAngStat, INT_TYPE, tmp, &fk);
      read_field(field, "saveCheckpoints", &flag->saveCheckpoints, INT_TYPE, tmp, &fk);
      read_field(field, "restartFromChk", &flag->restartFromChk, INT_TYPE, tmp, &fk);
      read_field(field, "restartCoulomb", &flag->restartCoulomb, INT_TYPE, tmp, &fk);
      read_field(field, "coulombDone", &flag->coulombDone, INT_TYPE, tmp, &fk);
      read_field(field, "calcCoulombOnly", &flag->calcCoulombOnly, INT_TYPE, tmp, &fk);
      read_field(field, "noCalcExciton", &flag->noCalcExciton, INT_TYPE, tmp, &fk);
      // ****** ****** ****** ****** ****** ******
      // Handle exceptions
      // ****** ****** ****** ****** ****** ******
      if (fk != 1)
      {
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
        printf("noCalcExciton = int, (if 1, then job will not compute Coulomb integrals or BSE).\n");

        fflush(stdout);
        exit(EXIT_FAILURE);
      }
      i++;
    }
  }
  else
  {
    printf("PROGRAM EXITING: input.par does not exist in directory\n");
    fprintf(stderr, "PROGRAM EXITING: input.par does not exist in directory\n");
    exit(EXIT_FAILURE);
  }
  fclose(pf);

  /************************************************************/
  /******************     SANITY CHECKS     *******************/
  /************************************************************/

  if ((par->epsX < 0.0) || (par->epsY < 0.0) || (par->epsZ < 0.0))
  {
    fprintf(stderr, "ERROR: one of the dielectric constants is negative!");
    exit(EXIT_FAILURE);
  }

  /************************************************************/
  /******************    SET DEPENDENCIES    ******************/
  /************************************************************/

  if (flag->SO == 1)
  {
    flag->useSpinors = 1; // generate spinor wavefunctions
    flag->NL = 1;         // SO automatically switches on NL
  }
  if (flag->useSpinors == 1)
  {
    flag->isComplex = 1; // psi is complex: doubles x2
    ist->nspin = 2;      // two-spins per state: ngrid x2
  }

  ist->ngrid_1 = 1.0 / (double)(ist->ngrid);
  ist->nthreads = parallel->nthreads;
  ist->complex_idx = flag->isComplex + 1;
  par->dv = grid->dv;
  ist->nx = grid->nx;
  ist->ny = grid->ny;
  ist->nz = grid->nz;

  // Using input parameters, print the current job state
  if (mpir == 0)
  {
    print_input_state(stdout, flag, grid, par, ist, parallel);
  }

  return;
}

/****************************************************************************/

void read_unsafe_input(
    double complex **psitot,
    double **eig_vals,
    double **sigma_E,
    xyz_st **R,
    grid_st *grid,
    double **gridx,
    double **gridy,
    double **gridz,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel)
{

  // Read in psitot, set up the grid, etc.
  // but in an unsafe way without using output.dat
  // Wavefunctions might not be aligned with the grid,
  // leading to uncontrolled errors
  // USE AT YOUR OWN RISK

  FILE *pf;

  const int mpir = parallel->mpi_rank;

  int i = 0;
  char field[1000], tmp[1000];
  char *endptr;
  endptr = malloc(10 * sizeof(endptr[0]));

  unsigned long stlen;
  unsigned long j;
  unsigned long cntr;
  unsigned long itmp;
  unsigned long homo_idx;
  unsigned long vb_mindex;
  unsigned long cb_maxdex;
  unsigned long nst_tot;
  unsigned long nholes = 0;
  unsigned long nelecs = 0;

  double dbltmp;
  double Rx, Ry, Rz;

  if (access("unsafe_input.par", F_OK) != -1)
  {
    int fk; // "found key"

    pf = fopen("unsafe_input.par", "r");

    while (fscanf(pf, "%s = %s", field, tmp) != EOF && i < 100)
    {
      fk = 0;

      // ****** ****** ****** ****** ****** ******
      // Read unsafe input for BSE algorithm
      // ****** ****** ****** ****** ****** ******
      read_field(field, "mnStatesTot", &nst_tot, LONG_TYPE, tmp, &fk);
      read_field(field, "nHoles", &nholes, LONG_TYPE, tmp, &fk);
      read_field(field, "nElecs", &nelecs, LONG_TYPE, tmp, &fk);
      read_field(field, "nAtoms", &ist->natoms, LONG_TYPE, tmp, &fk);
      read_field(field, "nx", &grid->nx, LONG_TYPE, tmp, &fk);
      read_field(field, "ny", &grid->ny, LONG_TYPE, tmp, &fk);
      read_field(field, "nz", &grid->nz, LONG_TYPE, tmp, &fk);
      read_field(field, "xmin", &grid->xmin, DOUBLE_TYPE, tmp, &fk);
      read_field(field, "ymin", &grid->ymin, DOUBLE_TYPE, tmp, &fk);
      read_field(field, "zmin", &grid->zmin, DOUBLE_TYPE, tmp, &fk);
      read_field(field, "dx", &grid->dx, DOUBLE_TYPE, tmp, &fk);
      read_field(field, "dy", &grid->dy, DOUBLE_TYPE, tmp, &fk);
      read_field(field, "dz", &grid->dz, DOUBLE_TYPE, tmp, &fk);
      read_field(field, "fermiE", &par->fermi_E, DOUBLE_TYPE, tmp, &fk);
      read_field(field, "sigmaECut", &par->sigma_E_cut, DOUBLE_TYPE, tmp, &fk);
      // ****** ****** ****** ****** ****** ******
      // Handle exceptions
      // ****** ****** ****** ****** ****** ******
      if (fk != 1)
      {
        printf("\nFIELD NOT RECOGNIZED in unsafe_input.par: %s\n", field);
        printf("\nInvalid input field and/or format -> equal sign required after each field\n");
        printf("Only allowed fields are:\n");
        printf("mnStatesTot = long (total number of states in system)\n");
        printf("nHoles = long (number of hole states)\n");
        printf("nElecs = long (number of electron states)\n");
        printf("nAtoms = long (number of atoms in system)\n");
        printf("nx = double (grid points in x direction)\n");
        printf("ny = double (grid points in y direction)\n");
        printf("nz = double (grid points in z direction)\n");
        printf("xmin = double (minimum x coordinate)\n");
        printf("ymin = double (minimum y coordinate)\n");
        printf("zmin = double (minimum z coordinate)\n");
        printf("dx = double (grid spacing in x direction)\n");
        printf("dy = double (grid spacing in y direction)\n");
        printf("dz = double (grid spacing in z direction)\n");
        printf("fermiE = double (Fermi energy)\n");
        printf("sigmaECut = double (eig val variance cutoff)\n");

        fflush(stdout);
        exit(EXIT_FAILURE);
      }
      i++;
    }
    fclose(pf);
  }
  else
  {
    printf("PROGRAM EXITING: unsafe_input.par does not exist in directory\n");
    fprintf(stderr, "PROGRAM EXITING: unsafe_input.par does not exist in directory\n");
    exit(EXIT_FAILURE);
  }

  // Assign and construct system size variables
  if ((nholes != 0) && (nelecs != 0))
  {
    ist->mn_states_tot = nholes + nelecs;
  }
  else
  {
    printf("ERROR: unsafe_input has 0 nholes or 0 nelecs!\n");
    exit(EXIT_FAILURE);
  }

  ist->nx = grid->nx;
  ist->ny = grid->ny;
  ist->nz = grid->nz;
  ist->ngrid = grid->ngrid = ist->nx * ist->ny * ist->nz;
  grid->dkx = par->dkx = (TWOPI) / ((double)grid->nx * grid->dx);
  grid->dky = par->dky = (TWOPI) / ((double)grid->ny * grid->dy);
  grid->dkz = par->dkz = (TWOPI) / ((double)grid->nz * grid->dz);

  ist->ngrid_1 = 1.0 / (double)ist->ngrid;

  grid->nx_1 = 1.0 / (double)grid->nx;
  grid->ny_1 = 1.0 / (double)grid->ny;
  grid->nz_1 = 1.0 / (double)grid->nz;

  par->dx = grid->dx;
  par->dy = grid->dy;
  par->dz = grid->dz;
  grid->dr = par->dr = sqrt(grid->dx * grid->dx + grid->dy * grid->dy + grid->dz * grid->dz);
  grid->dv = par->dv = grid->dx * grid->dy * grid->dz;

  par->xmin = grid->xmin;
  par->ymin = grid->ymin;
  par->zmin = grid->zmin;

  ist->nspinngrid = 2 * ist->ngrid; // hardocded for complex spinors, need to make general
  ist->complex_idx = 2;             // these are hardocde for complex spinors, need to make general
  stlen = ist->nspinngrid;

  if (mpir == 0)
  {
    printf("Done reading in unsafe_input.par\n\n");
    printf("mnStatesTot = %ld\n", nst_tot);
    printf("xmin = %lg ymin = %lg zmin = %lg\n", grid->xmin, grid->ymin, grid->zmin);
    printf("nx = %ld ny = %ld nz = %ld\n", grid->nx, grid->ny, grid->nz);
    printf("dx = %lg dy = %lg dz = %lg\n", grid->dx, grid->dy, grid->dz);
    printf("grid->dv = %lg \n", grid->dv);
    printf("natoms = %ld\n", ist->natoms);
    printf("mn_states_tot = %ld\n", ist->mn_states_tot);
  }

  // Allocate memory for grid, psi, eigs, sigma_E, and R
  (*gridx) = (double *)calloc(grid->nx, sizeof(double));
  (*gridy) = (double *)calloc(grid->ny, sizeof(double));
  (*gridz) = (double *)calloc(grid->nz, sizeof(double));

  (*psitot) = (double complex *)calloc(ist->mn_states_tot * stlen, sizeof(double complex));

  (*eig_vals) = (double *)calloc(nst_tot, sizeof(double));
  (*sigma_E) = (double *)calloc(nst_tot, sizeof(double));

  ALLOCATE(R, ist->natoms, "R in unsafe_input");
  //(*R)           =   (xyz_st*) calloc(ist->natoms, sizeof(xyz_st));

  // Build the grid

  for (j = 0; j < grid->nx; j++)
  {
    (*gridx)[j] = grid->xmin + j * grid->dx;
  }
  for (j = 0; j < grid->ny; j++)
  {
    (*gridy)[j] = grid->ymin + j * grid->dy;
  }
  for (j = 0; j < grid->nz; j++)
  {
    (*gridz)[j] = grid->zmin + j * grid->dz;
  }

  // Define grid maxima
  grid->xmax = (*gridx)[grid->nx - 1];
  grid->ymax = (*gridy)[grid->ny - 1];
  grid->zmax = (*gridz)[grid->nz - 1];
  par->xmax = grid->xmax;
  par->ymax = grid->ymax;
  par->zmax = grid->zmax;

  // Read in the eig_vals and sigma_E
  if (access("eval.par", F_OK) != -1)
  {
    pf = fopen("eval.par", "r");

    for (j = 0; j < nst_tot; j++)
    {
      fscanf(pf, "%lu %lg %lg", &itmp, &((*eig_vals)[j]), &((*sigma_E)[j]));
    }

    fclose(pf);
  }
  else
  {
    printf("PROGRAM EXITING: eval.par does not exist in directory\n");
    fprintf(stderr, "PROGRAM EXITING: eval.par does not exist in directory\n");
    exit(EXIT_FAILURE);
  }

  // Get the indices of the homo and lumo
  get_fmo_idxs(*eig_vals, *sigma_E, par->fermi_E, par->sigma_E_cut, nst_tot, &homo_idx);

  // Read in psitot from psi.par

  // 1. Get the indices of the first and last states
  vb_mindex = homo_idx - nholes + 1;
  cb_maxdex = vb_mindex + ist->mn_states_tot;

  // printf("vb_mindex = %lu cb_maxdex = %lu\n", vb_mindex, cb_maxdex);

  if (access("psi.par", F_OK) != -1)
  {
    pf = fopen("psi.par", "r");
    fseek(pf, vb_mindex * stlen * sizeof(double complex), SEEK_SET);
    fread((*psitot), stlen * sizeof(double complex), ist->mn_states_tot, pf);
    fclose(pf);
  }
  else
  {
    printf("PROGRAM EXITING: psi.par does not exist in directory\n");
    fprintf(stderr, "PROGRAM EXITING: psi.par does not exist in directory\n");
    exit(EXIT_FAILURE);
  }

  // Read in the short eig_vals and sigma_E

  pf = fopen("eval.par", "r");
  cntr = 0;
  for (j = 0; j < nst_tot; j++)
  {
    if ((j >= vb_mindex) && (j < cb_maxdex))
    {
      fscanf(pf, "%lu %lg %lg", &itmp, &((*eig_vals)[cntr]), &((*sigma_E)[cntr]));
      // printf("eig_vals[j = %lu] = %lg\n",j, ((*eig_vals)[cntr]) );
      cntr++;
    }
    else
    {
      fscanf(pf, "%lu %lg %lg", &itmp, &dbltmp, &dbltmp);
    }
  }
  fclose(pf);

  // Read in atomic configuration
  if (access("conf.par", F_OK) != -1)
  {
    pf = fopen("conf.par", "r");

    for (j = 0; j < ist->natoms; j++)
    {
      fscanf(pf, "%s %lg %lg %lg", endptr, &Rx, &Ry, &Rz);
      (*R)[j].x = Rx;
      (*R)[j].y = Ry;
      (*R)[j].z = Rz;
    }

    fclose(pf);
  }
  else
  {
    printf("PROGRAM EXITING: conf.par does not exist in directory\n");
    fprintf(stderr, "PROGRAM EXITING: conf.par does not exist in directory\n");
    exit(EXIT_FAILURE);
  }

  if (mpir == 0)
  {
    printf("Done with function read_unsafe_input\n");
    fflush(0);
  }

  for (j = 0; j < ist->mn_states_tot; j++)
  {
    printf("%ld %lg %lg\n", j, (*eig_vals)[j], (*sigma_E)[j]);
  }

  return;
}

/****************************************************************************/

// Function to read and set an int, long, or double field
void read_field(
    const char *field,
    const char *key,
    void *var,
    VarType type,
    const char *tmp,
    int *fk)
{

  char *endptr;

  if (!strcmp(field, key))
  {

    if (type == INT_TYPE)
    {
      *(int *)var = (int)strtol(tmp, &endptr, 10);
    }

    else if (type == LONG_TYPE)
    {
      *(long *)var = strtol(tmp, &endptr, 10);
    }

    else if (type == DOUBLE_TYPE)
    {
      *(double *)var = strtod(tmp, &endptr);
    }

    if (*endptr != '\0')
    {
      fprintf(stderr, "Error converting string to %s for %s.\n",
              (type == INT_TYPE) ? "int" : (type == LONG_TYPE) ? "long"
                                                               : "double",
              key);
      exit(EXIT_FAILURE);
    }

    *fk = 1;
  }
}

/****************************************************************************/

void read_bseCoeff(int n_xton, int numExcStatesToRead, zomplex *u, FILE *pf)
{
  int i, j;

  if (numExcStatesToRead > 100)
  {
    numExcStatesToRead = 100;
  }

  for (i = 0; i < n_xton; i++)
  {
    for (j = 0; j < numExcStatesToRead; j++)
    {
      fscanf(pf, "{%lg, %lg}\t", &(u[i * n_xton + j].re), &(u[i * n_xton + j].im));
    }
  }
}

/****************************************************************************/

void get_fmo_idxs(
    double *eig_vals,
    double *sigma_E,
    double fermiE,
    double secut,
    unsigned long n_elems,
    unsigned long *homo_idx)
{

  int i;
  int hidx = 0;

  // HOMO index
  for (i = 0; i < n_elems; i++)
  {
    // If the sigma_E val > secut, skip
    if ((eig_vals[i] < fermiE) && (sigma_E[i] < secut))
    {
      hidx = i;
    }
  }
  (*homo_idx) = hidx;

  // printf("Found homo_idx = %lu\n", *homo_idx);

  return;
}
