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
  double complex*     bsmat;
  double complex*     direct;
  double complex*     exchange; 
  double complex*     bs_coeff;
  double complex*     ldots;
  double complex*     psi_qp;
  double complex*     psitot = NULL;
  double complex*     pot_bare;
  double complex*     pot_screened;

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
  double complex*     l2_mom;
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

  // /************************************************************/
	// /**********************    MOD POT     **********************/
	// /************************************************************/

  mod_pot(&pot_bare, &pot_screened, &grid, &ist, &par, &flag, &parallel);

  MPI_Barrier(MPI_COMM_WORLD);
  
  // /************************************************************/
	// /********************    MOD KERNEL     *********************/
	// /************************************************************/

  mod_kernel(
    psi_qp, &direct, &exchange, pot_bare, pot_screened, &ist, &par, &flag, &parallel
  );
  
  
  if (parallel.mpi_rank == 0){
    /*************************************************************************/
    /*************************************************************************/
    mod_bse(
      psi_qp, direct, exchange, &bsmat, &bs_coeff, &h0mat, &xton_ene, eig_vals,
      &grid, &ist, &par, &flag, &parallel
    );
    
    if (mpir == 0) {
      write_separation(stdout, "T");
      printf("\n  -  COMPUTING XTON OPTICAL PROPERTIES | %s\n", get_time());
      write_separation(stdout, "B"); fflush(stdout);
    }
    calc_optical_exc(bs_coeff, xton_ene, eig_vals, elec_dip, mag_dip, &ist, &par);


    if (1 == flag.calcSpinAngStat){
      if (mpir == 0) {
        write_separation(stdout, "T");
        printf("\n  -  COMPUTING ANG.MOM. PROPERTIES | %s\n", get_time());
        write_separation(stdout, "B"); fflush(stdout);
      }
      
      long mat_size = sqr(ist.n_holes) + sqr(ist.n_elecs);

      ALLOCATE(&l_mom,  mat_size, "l_mom");  //<psi_r|L|psi_s>
      ALLOCATE(&l2_mom, mat_size, "l2_mom"); //<psi_r|L2|psi_s>
      ALLOCATE(&s_mom,  mat_size, "s_mom");  //<psi_r|S|psi_s>
      ALLOCATE(&ldots,  mat_size, "ldots");  //<psi_r|L.S|psi_s>

      if (mpir == 0) {
        printf("\n  -  Quasiparticle angular momenta | %s\n", get_time());
      }
      // Compute spin matrix elements, e.g. <j|Sx|i>
      calc_qp_spin_mtrx(psi_qp, s_mom, &grid, &ist, &par);
      // Compute angular momentum matrix elements, e.g. <j|Lx|i>
      calc_qp_ang_mom_mtrx(psi_qp, l_mom, l2_mom, ldots, &grid, &ist, &par);
    
    
      if (mpir == 0) {
        printf("\n  -  Exciton angular momenta | %s\n", get_time());
      }

      calc_xton_spin_mtrx(bs_coeff, s_mom, &ist, &par, &flag, &parallel);
      calc_xton_ang_mom_mtrx(bs_coeff, s_mom, l_mom, l2_mom, ldots, &ist, &par, &flag, &parallel);
      free(s_mom); 
      free(l_mom); 
      free(l2_mom); 
      free(ldots);
    }

    free(elec_dip); 
    free(mag_dip);
    free(rot_strength);
    free(xton_ene); 
    free(bs_coeff); 
    free(bsmat); 
    free(h0mat);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  // /***********************************************************************/
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
  // free(ist.eval_elec_idxs);
  // free(ist.eval_hole_idxs);
  free(ist.atom_types);

  time_t end_time = time(NULL);
  time_t end_clock = clock();

  if (mpir == 0){ 
    write_separation(stdout, "T");
    printf("\nDONE WITH PROGRAM: BETHE-SALPETHER\n");
    printf("This calculation ended at: %s\n", ctime(&end_time)); 
    printf("Total job CPU time (sec) %.4g, wall run time (sec) %.4g | %s",
      ((double)end_clock - (double)start_clock)/(double)(CLOCKS_PER_SEC), 
      (double)end_time - (double)start_time, format_duration((double)end_time - (double)start_time) 
    );
    fflush(0);
    write_separation(stdout, "B");
  }

  

  MPI_Finalize();

  return 0;
}

/*****************************************************************************/
