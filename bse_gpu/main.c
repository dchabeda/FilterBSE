#include "main.h"

/****************************************************************
* Core algorithm written by Eran Rabani, John Philbin, Dipti    *
* Jasrasaria, and Daniel Weinberg.                              *
* Current source code structure, MPI parallelism, GPU offload,  *
* scalar/spinor compatability, and documentation written by     *
* Daniel Chabeda 03.24.2025                                     *
*****************************************************************/

/*****************************************************************/
/*****************************************************************
* This is the main function for BSE.x. It is the driver that     *
* controls memory allocation, structs, and the BSE algorithm.    *
* The algorithm computes correlated excitonic states using the   *
* Bethe-Salpether formalism described by Rohlfing and Louie:     *
* Phys. Rev. B 62, 4927                                          *
* It is applied to nanocrystal systems by using an electron-hole *
* basis constructed from a semi-empirical pseudopot. hamiltonian *
******************************************************************/

int main(int argc, char *argv[]){

  /************************************************************/
	/*******************  DECLARE VARIABLES   *******************/
	/************************************************************/

  //                  Complex Arrays
  double*             pot_bare;
  double*             pot_screened;
  zomplex*            bsmat;
  zomplex*            direct;
  zomplex*            exchange; 
  zomplex*            bs_coeff;
  zomplex*            LdotS;
  zomplex*            l2_mom;

  //                  Custom Structs
  par_st              par = {0}; 
  index_st            ist = {0}; 
  grid_st             grid = {0}; 
  flag_st             flag = {0};  
  parallel_st         parallel = {0}; 
  xyz_st*             R = NULL;
  xyz_st*             s_mom; 
  xyz_st*             l_mom; 
  xyz_st*             elec_dip; 
  xyz_st*             mag_dip; 
  double*             rot_strength; 

  //                  Double Arrays
  double*             psitot = NULL; 
  double*             psi_qp; 
  double*             eig_vals = NULL; 
  double*             sigma_E = NULL; 
  double*             xton_ene; 
  double*             h0mat;

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

  // #pragma omp target
  // {
  //   int num_devices = omp_get_num_devices();
  //   printf("Number of OpenMP target devices available: %d\n", num_devices);
  //   printf("Running on OpenMP target device: %d\n", omp_is_initial_device());

  //   int num_threads = omp_get_num_threads();
  //   int max_threads = omp_get_max_threads();
  //   int num_teams = omp_get_num_teams();

  //   printf("GPU Threads per team: %d\n", num_threads);
  //   printf("Max threads per team: %d\n", max_threads);
  //   printf("Number of teams: %d\n", num_teams);
  // }

  // int num_teams = -1, num_threads = -1;

  // #pragma omp target teams num_teams(32) thread_limit(128) map(from:num_teams, num_threads)
  // {
  //   #pragma omp parallel
  //   {
  //       num_teams = omp_get_num_teams();
  //       num_threads = omp_get_num_threads();
  //   }
  // }

  // printf("Number of teams: %d\n", num_teams);
  // printf("Threads per team: %d\n", num_threads);
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
  
  
  if ( (1 == flag.useSpinors) && (mpir == parallel.mpi_root)){
    qp_spin_frac(psi_qp, eig_vals, &grid, &ist, &par, &flag, &parallel);
  }
  
  
  /************************************************************/
	/*********************   MOD DIPOLE     *********************/
	/************************************************************/

  if (mpir == 0){
    mod_dipole(
      psi_qp, eig_vals, &grid, &elec_dip, &mag_dip, &rot_strength,
      &ist, &par, &flag, &parallel
    );
  }

  /************************************************************/
	/**********************    MOD POT     **********************/
	/************************************************************/

  mod_pot(&pot_bare, &pot_screened, &grid, &ist, &par, &flag, &parallel);
  
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  if (1 != flag.coulombDone){
    mod_kernel(
      psi_qp, &direct, &exchange, pot_bare, pot_screened, &ist, &par, &flag, &parallel
    );
  } else if (1 == flag.coulombDone){
    ALLOCATE(&direct,   ist.n_xton * ist.n_xton, "direct");
    ALLOCATE(&exchange, ist.n_xton * ist.n_xton, "exchange");
    load_coulomb_mat(direct, "direct.dat", &ist);
    load_coulomb_mat(exchange, "exchange.dat", &ist);
  } else{
    fprintf(stderr, "ERROR: Invalid value of flag.coulombDone = %d\n", flag.coulombDone);
    exit(EXIT_FAILURE);
  }

  if (parallel.mpi_rank == 0){
    /*************************************************************************/
    /*************************************************************************/
    // 2. Compute single particle properties
    /* */
    /* */
    /* */
    if (1 == flag.SO){
      if (mpir == 0) {
        write_separation(stdout, "T");
        printf("\n  -  COMPUTING ANG.MOM. PROPERTIES | %s\n", get_time());
        write_separation(stdout, "B"); fflush(stdout);
      }
      
      l_mom = (xyz_st *) malloc( (ist.n_holes*ist.n_holes + ist.n_elecs*ist.n_elecs) * sizeof(l_mom[0])); //<psi_r|Lx|psi_s>
      l2_mom =(zomplex *) malloc((ist.n_holes*ist.n_holes + ist.n_elecs*ist.n_elecs) * sizeof(l2_mom[0])); //<psi_r|Lx|psi_s>
      s_mom = (xyz_st *) malloc( (ist.n_holes*ist.n_holes + ist.n_elecs*ist.n_elecs) * sizeof(s_mom[0])); //<psi_r|Sx|psi_s>
      LdotS = (zomplex *) malloc( (ist.n_holes*ist.n_holes+ist.n_elecs*ist.n_elecs) * sizeof(LdotS[0])); //<psi_r|L.S|psi_s>

      // Compute spin matrix elements, e.g. <j|Sx|i>
      calc_spin_mtrx(s_mom, psi_qp, &grid, &ist, &par);
      // Compute angular momentum matrix elements, e.g. <j|Lx|i>
      calc_ang_mom_mtrx(l_mom, l2_mom, LdotS, psi_qp, &grid, &ist, &par);
    }
    /* */
    /* */
    /* */
    ALLOCATE(&bsmat, ist.n_xton * ist.n_xton, "bsmat");
    ALLOCATE(&h0mat, ist.n_xton * ist.n_xton, "h0mat");

    build_BSE_mat(bsmat, direct, exchange, &ist);
    build_h0_mat(h0mat, eig_vals, &ist);
    
    xton_ene = (double *) calloc(ist.n_xton, sizeof(double));
    bs_coeff = (zomplex *) calloc(ist.n_xton*ist.n_xton, sizeof(zomplex));
    bethe_salpeter(bsmat, direct, exchange, bs_coeff, h0mat, xton_ene,(zomplex*) psi_qp, s_mom,l_mom,l2_mom,LdotS, &grid, &ist, &par);
    
    calc_optical_exc(bs_coeff, xton_ene, elec_dip, mag_dip, &ist, &par);
    
    free(xton_ene); 
    free(bs_coeff);
    free(elec_dip); 
    free(mag_dip); 
    free(rot_strength);
    free(s_mom); 
    free(l_mom); 
    free(l2_mom); 
    free(LdotS);
    free(bsmat); 
    free(h0mat);
  }

  /***********************************************************************/
  free(psi_qp); 
  free(eig_vals); 
  free(grid.x); 
  free(grid.y); 
  free(grid.z);
  free(R);
  free(pot_bare); 
  free(pot_screened);
  free(direct); 
  free(exchange);
  free(ist.eval_elec_idxs);
  free(ist.eval_hole_idxs);
  free(ist.atom_types);

  time_t end_time = time(NULL);
  time_t end_clock = clock();

  if (mpir == 0){ 
    write_separation(stdout, "T");
    printf("\nDONE WITH PROGRAM: BETHE-SALPETHER\n");
    printf("This calculation ended at: %s\n", ctime(&end_time)); 
    printf("Total job CPU time (sec) %.4g, wall run time (sec) %.4g",
      ((double)end_clock - (double)start_clock)/(double)(CLOCKS_PER_SEC), 
      (double)end_time - (double)start_time 
    );
    fflush(0);
    write_separation(stdout, "B");
  }

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();
  // exit(0);

  return 0;
}

/*****************************************************************************/
