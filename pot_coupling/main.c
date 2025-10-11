#include "main.h"

/****************************************************************
* Daniel Chabeda 04.21.2025                                     *
*****************************************************************/

/*****************************************************************/
/*****************************************************************
* This is the main function for pot_coupling.x. It is the driver *
* for computing integrals of <psi_i|U(r;R)|psi_j> used in the    *
* derivation of nonadiabatic coupling expressions (Jasrasaria    *
* thesis)                                                        *
******************************************************************
******************************************************************/

int main(int argc, char *argv[]){

  /************************************************************/
	/*******************  DECLARE VARIABLES   *******************/
	/************************************************************/

  //              zomplex types
  zomplex*        psi;            // Wavefunction on the grid
  zomplex*        phi;            // Wavefunction on the grid (helper)
  zomplex*        an;             // Newton interpolation coefficients
  zomplex*        LS = NULL;      // Spin-orbit matrix
  
  // custom structs 
  index_st        ist;            // Indexes for the program
  par_st          par;            // Parameters for the program
  flag_st         flag;           // Flags for options in the program
  parallel_st     parallel;       // MPI parallelization information
  xyz_st*         R;              // Atomic positions in the x, y, and z
  atom_info*      atom;           // Atom specific information
  pot_st          pot;            // Atomic pseudopotential information
  grid_st         grid;           // Grid params and the real space grid
  grid.x =        NULL;
  grid.y =        NULL;
  grid.z =        NULL;
  nlc_st*         nlc = NULL;     // Non-local pseudopotential info

  // double arrays
  double*         psitot = NULL;  // All filtered states
  double*         pot_local;      // Local pseudopotential on the grid
  double*         SO_projectors;  // Spin-orbit projectors
  double*         ksqr;           // Kinetic energy sqr on the grid
  double*         zn;             // Chebyshev support points
  double*         eig_vals;       // Quasiparticle energies
  double*         ene_targets;    // Energy targets for filtering
  double*         sigma_E;        // Standard dev. of the energies
  
  // long int arrays and counters
  long*           nl = NULL;

  //                  Double Arrays
  double*             psi_qp; 

  //                  Clock/Wall Time Output and Stdout Formatting
  time_t              start_time = time(NULL);  // Get the actual time for total wall runtime
  time_t              start_clock = clock();    // Get the starting CPU clock time for total CPU runtime

  char*               endptr;
  
  ist.atom_types    = malloc(N_MAX_ATOM_TYPES * sizeof(ist.atom_types[0])); 

  /************************************************************/
	/*********************     INIT MPI     *********************/
	/************************************************************/

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &parallel.mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &parallel.mpi_size);

  const int mpir    = parallel.mpi_rank;
  parallel.mpi_root = 0;
    

  /************************************************************/
	/*********************    START JOB     *********************/
	/************************************************************/

  if (mpir == 0){
    fprintf(stdout, "******************************************************************************\n");
    fprintf(stdout, "\nRUNNING PROGRAM: BETHE-SALPETHER\n");
    fprintf(stdout, "This calculation began at: %s", ctime(&start_time)); 
    write_separation(stdout, "B");
    fflush(stdout);
  }

  /************************************************************/
	/**********************    MOD INIT    **********************/
	/************************************************************/

  // Init unsafe is a temporary back-compatability measure to enable
  // reading in eigenvalues and eigenfunctions from eval.par & psi.par
  // In future, all output should be formatted output.dat from filter
  flag.initUnsafe = 0;
  if (argv[1] != NULL) flag.initUnsafe = (int) strtol(argv[1], &endptr, 10);

  mod_init(
    &psitot, &psi_qp, &eig_vals, &sigma_E, &R, &grid, 
    &(grid.x), &(grid.y), &(grid.z), &ist, &par, &flag, &parallel
  );

  /************************************************************/
	/**********************    MOD POT     **********************/
	/************************************************************/

  // Initialize the local, nonlocal, and spin-orbit potentials
  mod_pseudopot(
    pot_local, &pot, R, atom, &grid, LS, nlc, nl, SO_projectors, 
    &ist, &par, &flag, &parallel
  );

}