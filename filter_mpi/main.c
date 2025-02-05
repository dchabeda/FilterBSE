#include "fd.h"
#include <mpi.h>

/*****************************************************************************/

int main(int argc, char *argv[]){
  /*****************************************************************
  * This is the main function for Filter.x. It is the driver that  *
  * controls memory allocation, structs, and the filter algorithm. *
  * The algorithm computes quasiparticle excited states using      *
  * sparse-matrix techniques. It is applied to nanocrystal         *
  * systems through the use of semiempirical pseudopotentials.     *
  ******************************************************************/  

  // DECLARE VARIABLES AND STRUCTS
  // file pointers
  FILE *ppsi, *peig, *pseed; 
  // zomplex types
  zomplex *psi, *phi, *an; 
  int i;
  // custom structs 
  flag_st flag; 
  index_st ist; 
  par_st par; 
  atom_info *atom; 
  pot_st pot; 
  grid_st grid;
  gauss_st *gauss;
  xyz_st *R; 
  nlc_st *nlc = NULL; 
  parallel_st parallel; 
  // double arrays
  double *psitot, *psi_rank;
  double *ksqr, *zn, *pot_local, *rho; 
  double *eig_vals, *ene_targets, *sigma_E, inital_clock_t, initial_wall_t;
  double *SO_projectors;
  // long int arrays and counters
  long *nl = NULL;
  long jstate, jgrid, jgrid_real, jgrid_imag, jspin, jms, jns, rand_seed, thread_id;
  ist.atom_types = malloc(N_MAX_ATOM_TYPES*sizeof(ist.atom_types[0]));
  // Clock/Wall time output and stdout formatting
  time_t start_time = time(NULL); // Get the actual time for total wall runtime
  time_t start_clock = clock(); // Get the starting CPU clock time for total CPU runtime
  time_t current_time;
  char* c_time_string;
  char *top; top = malloc(2*sizeof(top[0])); 
  char *bottom; bottom = malloc(2*sizeof(bottom[0]));
  strcpy(top, "T\0"); strcpy(bottom, "B\0");

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &parallel.mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &parallel.mpi_size);
  parallel.mpi_root = 0;

  if (parallel.mpi_rank == 0) printf("******************************************************************************\n");
  if (parallel.mpi_rank == 0) printf("\nRUNNING PROGRAM: FILTER DIAGONALIZATION\n");
  if (parallel.mpi_rank == 0) printf("This calculation began at: %s\n", ctime(&start_time)); 
  if (parallel.mpi_rank == 0) printf("Printing from root node %d", parallel.mpi_rank); 
  if (parallel.mpi_rank == 0) write_separation(stdout, bottom);
  fflush(stdout);
  
  /*************************************************************************/
  // Initialize job from input file
  
  if (parallel.mpi_rank == 0) write_separation(stdout, top);
  current_time = time(NULL);
  c_time_string = ctime(&current_time);
  if (parallel.mpi_rank == 0) printf("\n1.\tINITIALIZING JOB | %s\n", c_time_string);
  if (parallel.mpi_rank == 0) write_separation(stdout, bottom); fflush(stdout);

  /*** read initial setup from input.par ***/
  if (parallel.mpi_rank == 0) printf("\nReading job specifications from input.par:\n");
  read_input(&flag, &grid, &ist, &par, &parallel);

  /*** allocating memory ***/
  // the positions of the atoms in the x, y, and z directions 
  if ((R = (xyz_st *) calloc(ist.natoms, sizeof(xyz_st))) == NULL) {
    if (parallel.mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: R array\n\n"); exit(EXIT_FAILURE);
  }
  // the atom specific information 
  if ((atom = (atom_info *) calloc(ist.natoms, sizeof(atom_info))) == NULL){
    if (parallel.mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: atom struct\n\n"); exit(EXIT_FAILURE);
  }
  // the energies of each energy target in filter
  if ((ene_targets = (double *) calloc(ist.m_states_per_filter, sizeof(double))) == NULL){
    if (parallel.mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: ene_targets\n\n"); exit(EXIT_FAILURE);
  }
  
  /*** read the nanocrystal configuration ***/
  if (parallel.mpi_rank == 0) printf("\nReading atomic configuration from conf.par:\n");
  read_conf(R, atom, &ist, &par, &flag, &parallel);

  /*** initialize parameters for the grid ***/
  if (parallel.mpi_rank == 0) printf("\nInitializing the grid parameters:\n");
  init_grid_params(&grid, R, &ist, &par, &flag, &parallel);

  // Allocate memory for the grid in the x, y, and z directions ***/
  if ((grid.x = (double *) calloc(grid.nx, sizeof(double))) == NULL){
    if (parallel.mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: grid.x\n\n"); exit(EXIT_FAILURE);
  }
  if ((grid.y = (double *) calloc(grid.ny, sizeof(double))) == NULL){
    if (parallel.mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: grid.y\n\n"); exit(EXIT_FAILURE);
  }
  if ((grid.z = (double *) calloc(grid.nz, sizeof(double))) == NULL){
    if (parallel.mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: grid.z\n\n"); exit(EXIT_FAILURE);
  }
  // the kinetic energy stored on the grid
  if ((ksqr = (double *) calloc(ist.ngrid, sizeof(double))) == NULL){
    if (parallel.mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: ksqr\n\n"); exit(EXIT_FAILURE);
  }

  /*** build the real- and k-space grids ***/
  if (parallel.mpi_rank == 0) printf("\nBuilding the real-space and k-space grids:\n");
  build_grid_ksqr(ksqr, R, &grid, &ist, &par, &flag, &parallel);
  
  /*** set the energy targets ***/
  if (parallel.mpi_rank == 0) printf("\nSetting the filter energy targets:\n");
  set_ene_targets(ene_targets, &ist, &par, &flag, &parallel);

  
  /*************************************************************************/
  /*** allocating memory for the rest of the program ***/
  if (parallel.mpi_rank == 0) printf("\nAllocating memory for FFT, pot, psi, eig_vals...\n");
  
  // 
  // FFT
  // 
  // FFT
  fftw_init_threads();
  fftw_plan_with_nthreads(parallel.n_inner_threads);
  fftw_plan_loc planfw, planbw; fftw_complex *fftwpsi; 
  long fft_flags=0;
  // Initialize FFT threads for parallel Fourier transform
  // if (fftw_init_threads() == 0) {
  //   fprintf(stderr, "FFTW threading initialization failed.\n");
  //   exit(EXIT_FAILURE);
  // }
  
  // fftw_plan_with_nthreads(parallel.nthreads);

  // Load wisdom if available
  // sprintf(par.fftw_wisdom, "%sfftw_wisdom.dat", par.fft_wisdom_dir);
  // if (fftw_import_wisdom_from_filename(par.fftw_wisdom) == 0) {
  //     if (parallel.mpi_rank == 0) printf("No wisdom file found. Planning from scratch.\n");
  // } else {
  //     if (parallel.mpi_rank == 0) printf("FFT wisdom loaded successfully from %s.\n", par.fftw_wisdom);
  // }
  
  fftwpsi = fftw_malloc(sizeof(fftw_complex)*ist.ngrid);
  /*** initialization for the fast Fourier transform ***/
  planfw = fftw_plan_dft_3d(grid.nz, grid.ny, grid.nx, fftwpsi, fftwpsi, FFTW_FORWARD, fft_flags);
  planbw = fftw_plan_dft_3d(grid.nz, grid.ny, grid.nx, fftwpsi, fftwpsi, FFTW_BACKWARD, fft_flags);
  
  // Save wisdom for future runs
  //if (parallel.mpi_rank == 0) printf("wisdom file: %s\n", par.fftw_wisdom);
  // if (fftw_export_wisdom_to_filename(par.fftw_wisdom) == 0) {
  //   perror("FFT wisdom not saved.\n");
  //   if (parallel.mpi_rank == 0) fprintf(stderr, "Failed to write FFTW wisdom to %s\n", par.fftw_wisdom);
    
  // } else {
  //   if (parallel.mpi_rank == 0) printf("FFTW wisdom saved to %s\n", par.fftw_wisdom);
  // } 
  // fftw_export_wisdom_to_filename(par.fftw_wisdom);
  
  // For reading the atomic potentials ***/
  pot.dr = (double *) calloc(ist.ngeoms * ist.n_atom_types, sizeof(double));
  pot.r = (double *) calloc(ist.ngeoms * ist.max_pot_file_len * ist.n_atom_types, sizeof(double));
  if (1.0 != par.scale_surface_Cs){
    // allocate memory for separately read LR potentials if the surface atoms will be charge balanced
    pot.r_LR = (double *) calloc(ist.ngeoms * ist.max_pot_file_len * ist.n_atom_types, sizeof(double));
  }
  pot.pseudo = (double *) calloc(ist.ngeoms * ist.max_pot_file_len * ist.n_atom_types, sizeof(double));
  if (1.0 != par.scale_surface_Cs){
    pot.pseudo_LR = (double *) calloc(ist.ngeoms * ist.max_pot_file_len * ist.n_atom_types, sizeof(double));
  }
  pot.file_lens = (long *) calloc(ist.ngeoms * ist.n_atom_types, sizeof(long));
  
  // Wavefunction-type objects
  if ((psi = (zomplex *)calloc(ist.nspinngrid, sizeof(zomplex))) == NULL){
    if (parallel.mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: psi\n\n"); exit(EXIT_FAILURE);
  }
  if ((phi = (zomplex *)calloc(ist.nspinngrid, sizeof(zomplex))) == NULL){
    if (parallel.mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: phi\n\n"); exit(EXIT_FAILURE);
  }
  if ((pot_local = (double *) calloc(ist.ngrid, sizeof(double))) == NULL){
    if (parallel.mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: pot_local\n\n"); exit(EXIT_FAILURE);
  }
  
  // For Newton interpolation coefficients
  an = (zomplex *) calloc(ist.ncheby*ist.m_states_per_filter, sizeof(zomplex));
  zn = (double *) calloc(ist.ncheby, sizeof(double));
  
  // The value of isComplex is 1 if wavefunctions are real valued, 2 if functions are complex valued
  // The factor of par.t_rev_factor (2 w spinors, 1 w/o spinors) in the psitot memory allocation
  // is because we time reverse the spinors to get double the orthogonal states
  // Spinor calculations are 8 times more memory intensive than scalar calculations
  
  long long psi_rank_size = ist.n_states_per_rank*ist.nspinngrid*ist.complex_idx;
  if ((psi_rank = (double *) calloc(psi_rank_size, sizeof(psi_rank[0]))) == NULL){
    if (parallel.mpi_rank == 0) fprintf(stderr,"\nOUT OF MEMORY: psi_rank\n\n"); exit(EXIT_FAILURE);
  }
  // the quasiparticle energies and standard deviations
  if ((eig_vals = (double *) calloc(par.t_rev_factor*ist.mn_states_tot, sizeof(double))) == NULL){
    if (parallel.mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: eig_vals\n\n"); exit(EXIT_FAILURE);
  }
  if ((sigma_E = (double *) calloc(par.t_rev_factor*ist.mn_states_tot, sizeof(double))) == NULL){
    if (parallel.mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: sigma_E\n\n"); exit(EXIT_FAILURE);
  }
  // memory allocation for the spin-orbit potential 
  if ( (flag.SO == 1) || (flag.NL == 1) ){
    if ((SO_projectors = (double*) calloc(PROJ_LEN * ist.nproj, sizeof(double)))==NULL){nerror("mem_SO_projector");}  
  }
  // memory allocation for the non-local potential
  if (1 == flag.NL){
    if ((nlc = (nlc_st *) calloc(ist.n_NL_atoms*ist.n_NL_gridpts, sizeof(nlc_st))) == NULL){ 
      if (parallel.mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: nlc\n\n"); exit(EXIT_FAILURE);
    }
    if ((nl = (long *) calloc(ist.natoms, sizeof(nl[0]))) == NULL) {
      if (parallel.mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: nl\n\n"); exit(EXIT_FAILURE);
    }
  }
  
  if (parallel.mpi_rank == 0) printf("\tdone allocating memory.\n"); fflush(stdout);

  // The filter code supports restarting the job from a saved state. See save.c for formatting
  // of checkpoint files. See read_input in read.c for specifying the checkpoint restart
  switch (flag.restartFromCheckpoint){
    // If no restart requested, the job will run from the beginning with no fast-forward 
    case -1:
      /**************************************************************************/
      if (parallel.mpi_rank == 0) printf("\nInitializing potentials...\n");
      
      if (parallel.mpi_rank == 0) printf("\nLocal pseudopotential:\n");
      build_local_pot(pot_local, &pot, R, ksqr, atom, &grid, &ist, &par, &flag, &parallel);
      
      free(pot.r); pot.r = NULL; 
      free(pot.pseudo); pot.pseudo = NULL; 
      free(pot.dr); pot.dr = NULL; 
      free(pot.file_lens); pot.file_lens = NULL;
      if (1.0 != par.scale_surface_Cs){
        free(pot.r_LR); pot.r_LR = NULL;
        free(pot.pseudo_LR); pot.pseudo_LR = NULL; 
      }
      

      write_cube_file(pot_local, &grid, "localPot.cube");
      
      if(flag.SO==1) {
        if (parallel.mpi_rank == 0) printf("\nSpin-orbit pseudopotential:\n");
        init_SO_projectors(SO_projectors, &grid, R, atom, &ist, &par, &flag, &parallel);
        if (parallel.mpi_rank == 0) printf("\tSO projectors generated.\n");
      }
      /*** initialization for the non-local potential ***/
      if (flag.NL == 1){
        if (parallel.mpi_rank == 0) printf("\nNon-local pseudopotential:\n"); fflush(0);
        init_NL_projectors(nlc, nl, SO_projectors, &grid, R, atom, &ist, &par, &flag, &parallel);
        if (parallel.mpi_rank == 0) printf("\tNL projectors generated.\n");
      }
      // free memory allocated to SO_projectors
      if ( (flag.SO == 1) || (flag.NL == 1) ){
        free(SO_projectors); SO_projectors = NULL;
      }

      // If the job uses a Gaussian basis, run the job in the Gaussian basis
      if (1 == flag.useGaussianBasis){
        //
        //
        gauss_driver(gauss, pot_local, eig_vals, R, atom, &grid, &ist, &par, &flag, &parallel);
        //
        //
        exit(0);
      }

      // If the job uses a Gaussian basis, run the job in the Gaussian basis
      if (1 == flag.periodic){
        //
        //
        periodic_driver(pot_local, nlc, nl, &grid, &ist, &par, &flag, &parallel);
        //
        //
        exit(0);
      }
      
       
      /**************************************************************************/
      /*** calculate the energy range of the hamitonian ***/
      if (parallel.mpi_rank == 0) write_separation(stdout, top);
      current_time = time(NULL);
      c_time_string = ctime(&current_time);
      if (parallel.mpi_rank == 0) printf("\n2. CALCULATING HAMILTONIAN ENERGY RANGE | %s\n", c_time_string);
      if (parallel.mpi_rank == 0) write_separation(stdout, bottom); fflush(stdout);

      inital_clock_t = (double)clock(); 
      initial_wall_t = (double)time(NULL);
      
      get_energy_range(psi, phi, pot_local, &grid, nlc, nl, ksqr, &ist, &par, &parallel, &flag, planfw, planbw, fftwpsi);
      
      if (parallel.mpi_rank == 0) printf("\ndone calculate energy range, CPU time (sec) %g, wall run time (sec) %g\n",
                ((double)clock()-inital_clock_t)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-initial_wall_t); 
      fflush(stdout);

      
      /**************************************************************************/
      // FILTER ALGORITHM
      /**************************************************************************/
      if (parallel.mpi_rank == 0) write_separation(stdout, top);
      current_time = time(NULL);
      c_time_string = ctime(&current_time);
      if (parallel.mpi_rank == 0) printf("\n3. GENERATING COEFFICIENTS | %s\n", c_time_string);
      if (parallel.mpi_rank == 0) write_separation(stdout, bottom); fflush(stdout);
      
      /*** set parameters for the newton interpolation ***/
      par.dt = sqr((double)(ist.ncheby) / (2.5*par.dE));

      if (parallel.mpi_rank == 0) printf("\nCalculating Newton interpolation coefficients for filter functions.\n"); fflush(0);
      gen_newton_coeff(an, zn, ene_targets, &ist, &par, &parallel);

      
      if (flag.setSeed == 1){
        rand_seed = - par.rand_seed;} 
      else {
        Randomize();  rand_seed = -random() + parallel.mpi_rank;
      }

      if (parallel.mpi_rank == 0) printf("ncheby = %ld dt = %g dE = %g\n", ist.ncheby, par.dt, par.dE);
      if (parallel.mpi_rank == 0) printf("\n  Energy width, sigma, of filter function = %.6g a.u.\n", sqrt(1 / (2*par.dt)));
      if (parallel.mpi_rank == 0) printf("  Suggested max span of spectrum for filtering = %.6g a.u.\n", ist.m_states_per_filter * sqrt(1 / (2*par.dt)));
      if (parallel.mpi_rank == 0) printf("  Requested span of spectrum to filter = %.6g a.u.\n", (ene_targets[par.n_targets_VB-1] - ene_targets[0]) + (ene_targets[par.n_targets_VB+par.n_targets_CB-1] - ene_targets[par.n_targets_VB])); fflush(stdout);
      
      // Save checkpoint after computing the potential on the grid and Newton interpolation coeffs
      if (1 == flag.saveCheckpoints){
        par.checkpoint_id = 0;
        if (parallel.mpi_rank == 0) write_separation(stdout, top);
        if (parallel.mpi_rank == 0) printf("****    CHECKPOINT %d *** CHECKPOINT %d ** CHECKPOINT %d *** CHECKPOINT %d    ****", par.checkpoint_id, par.checkpoint_id, par.checkpoint_id, par.checkpoint_id);
        if (parallel.mpi_rank == 0) write_separation(stdout, bottom); 
        if (parallel.mpi_rank == 0) printf("\n"); fflush(stdout);

        save_job_state("checkpoint_0.dat",par.checkpoint_id,psitot,pot_local,ksqr,an,zn,ene_targets,nl,nlc,&grid,&ist,&par,&flag,&parallel);
      }
    
    // If checkpoint_id is 0, then restart job right before filtering step
    case 0:
      if (flag.restartFromCheckpoint == 0){
        par.checkpoint_id = flag.restartFromCheckpoint;
        if (parallel.mpi_rank == 0) write_separation(stdout, top);
        if (parallel.mpi_rank == 0) printf("****    CHECKPOINT %d *** CHECKPOINT %d ** CHECKPOINT %d *** CHECKPOINT %d    ****", par.checkpoint_id, par.checkpoint_id, par.checkpoint_id, par.checkpoint_id);
        if (parallel.mpi_rank == 0) write_separation(stdout, bottom); fflush(stdout);
        
        restart_from_save("checkpoint_0.dat",par.checkpoint_id,psitot,pot_local,ksqr,an,zn,ene_targets,nl,nlc,&grid,&ist,&par,&flag,&parallel);
      } 
      par.checkpoint_id++;
      /**************************************************************************/
      /*** start filtering loop.  we run over n_filter_cycles cycles and calculate ***/
      /*** m_states_per_filter filtered states at each cycle ***/
      
      
      if (parallel.mpi_rank == 0) write_separation(stdout, top);
      current_time = time(NULL);
      c_time_string = ctime(&current_time);
      if (parallel.mpi_rank == 0) printf("\n4. STARTING FILTERING | %s\n", c_time_string);
      if (parallel.mpi_rank == 0) write_separation(stdout, bottom); 

      // Create the initial random wavefunctions for each filter state
      // This is an array of n_filter_cycles * m_states_per_filter states
      // where every block of length m_states_per_filter has the same random wavefunction
      
      if (parallel.mpi_rank == 0) printf("\n  4.1 Initializing random filter states\n");
      init_filter_states(psi_rank, psi, &grid, &rand_seed, &ist, &par, &flag, &parallel);
      
      inital_clock_t = (double)clock(); 
      initial_wall_t = (double)time(NULL);
      
      if (1 == flag.timeHamiltonian){
        // Initialize state to perform Hamiltonian on for test
        memcpy(&phi[0], &psi_rank[0], ist.complex_idx*ist.nspinngrid*sizeof(psi_rank[0]));
        current_time = time(NULL);
        c_time_string = ctime(&current_time);
        if (parallel.mpi_rank == 0) printf("\n  Timing Hamiltonian operator... | %s\n", c_time_string); fflush(0);
        
        time_hamiltonian(phi,psi,pot_local,nlc,nl,ksqr,&ist,&par,&flag,&parallel);
        
        current_time = time(NULL);
        c_time_string = ctime(&current_time);
        if (parallel.mpi_rank == 0) printf("  Done timing Hamiltonian | %s\n", c_time_string); fflush(0);
      }

      if ((parallel.jns = (long *) calloc(ist.mn_states_tot, sizeof(parallel.jns[0]))) == NULL){ 
        fprintf(stderr, "\nOUT OF MEMORY: parallel.jns\n\n"); exit(EXIT_FAILURE);
      }
      if ((parallel.jms = (long *) calloc(ist.mn_states_tot, sizeof(parallel.jms[0]))) == NULL){ 
        fprintf(stderr, "\nOUT OF MEMORY: parallel.jms\n\n"); exit(EXIT_FAILURE);
      }

      if (parallel.mpi_rank == 0) printf("\n  4.2 Running filter cycle\n");
      run_filter_cycle(psi_rank,pot_local,nlc,nl,ksqr,an,zn,ene_targets,&grid,&ist,&par,&flag,&parallel);
      MPI_Barrier(MPI_COMM_WORLD); // Ensure all ranks synchronize here
      
      
      if (parallel.mpi_rank == 0) printf("\ndone calculating filter, CPU time (sec) %g, wall run time (sec) %g\n",
                ((double)clock()-inital_clock_t)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-initial_wall_t); 
      fflush(stdout);

      if (1 == flag.calcFilterOnly){
        exit(0);
      }
      /*************************************************************************/
      /*** read all filtered states ***/
      /*ppsi = fopen("psi-filt.dat" , "r");
      fread (psitot,sizeof(zomplex),ist.mn_states_tot*ist.ngrid,ppsi);
      fclose(ppsi);*/
      if (parallel.mpi_rank == 0) printf("Allocating mem for psitot\n"); fflush(0);
      if (parallel.mpi_rank == 0){
        if ((psitot = (double *) calloc(ist.complex_idx * par.t_rev_factor * ist.nspinngrid * ist.mn_states_tot, sizeof(psitot[0]))) == NULL){
          if (parallel.mpi_rank == 0) fprintf(stderr,"\nOUT OF MEMORY: psitot\n\n"); exit(EXIT_FAILURE);
        }
      }
      
      /*************************************************************************/
      if (parallel.mpi_rank == 0) printf("Gathering psitot from all mpi_ranks\n"); fflush(0);
      MPI_Barrier(MPI_COMM_WORLD); // Ensure all ranks synchronize here
      
      // If the size of each wavefunction is very large,
      // then MPI_Send/Recv is the only way to communicate
      // if (psi_rank_size > INT_MAX){
      //   printf("The size of psi_rank, %lld > %d, too large for MPI_Gather\n", psi_rank_size, INT_MAX);
      //   printf("Using MPI_Send/Recv to send individual psi to root node.\n");
      //   // Because the entire psi_rank was too large to send at once,
      //   // Split psi_rank into two blocks and send each one separately
      //   int n_blocks = 2;
      //   int buf_size = psi_rank_size / n_blocks;
        
      //   // If MPI root node, then Recv from all other nodes into psitot
      //   if (parallel.mpi_root == parallel.mpi_rank){
      //     printf("This is the root rank = %d, it should receive\n", parallel.mpi_rank); fflush(0);
      //     // Loop over all ranks (except root) to Recv data
      //     for (long n_rank = 0; n_rank < parallel.mpi_size; n_rank++){
      //       printf("Iterating to receive from rank %ld out of %d\n", n_rank, parallel.mpi_size); fflush(0);
      //       if (n_rank != parallel.mpi_root){
      //         printf("Non-root rank %ld\n", n_rank); fflush(0);
      //         long long offset = psi_rank_size*n_rank;
      //         MPI_Status *recv_status;
      //         int recv_tag = n_blocks*n_rank;
      //         // Recv all blocks (with appropriate tag) from the ranks
      //         for (int block = 0; block < n_blocks; block++){ 
      //           MPI_Recv(&psitot[offset + block*buf_size], buf_size, MPI_DOUBLE, n_rank, recv_tag+block, MPI_COMM_WORLD, recv_status);
      //         }
      //       }
      //       // If root rank, then just copy psi_rank into psitot 
      //       else {
      //         printf("Memcpying root rank %ld\n", n_rank); fflush(0);
      //         for (int block = 0; block < n_blocks; block++){
      //           memcpy(&psitot[buf_size*block], &psi_rank[buf_size*block], buf_size*sizeof(double));
      //         }
      //       }
      //     }
      //   // If not the MPI root node, then Send psi_rank data
      //   } else {
      //     printf("Sending from non-root rank %d\n", parallel.mpi_rank); fflush(0);
      //     int send_tag = n_blocks*parallel.mpi_rank;
      //     // Send all blocks from this mpi_rank (with appropriate tag) to the root node
      //     for (int block = 0; block < n_blocks; block++){
      //       long long offset = buf_size * block;
      //       MPI_Send(&psi_rank[offset], buf_size, MPI_DOUBLE, parallel.mpi_root, send_tag+block, MPI_COMM_WORLD);
      //     }
      //   }
      // }
      // // If the size of each wavefunction is small,
      // // then MPI_Gather is an efficient way to communicate
      // else {
        MPI_Gather(psi_rank,psi_rank_size,MPI_DOUBLE,psitot,psi_rank_size,MPI_DOUBLE,0,MPI_COMM_WORLD);
      // }
      if (parallel.mpi_rank == 0) printf("Succesfully gathered all states\n"); fflush(0);
      MPI_Barrier(MPI_COMM_WORLD); // Ensure all ranks synchronize here

      /*************************************************************************/
      if (parallel.mpi_rank == 0){
        if (1 == flag.useSpinors){
          if (parallel.mpi_rank == 0) printf("\nTime-reversing all filtered states (doubles number of orthogonal states)");
          time_reverse_all(&psitot[0], &psitot[ist.complex_idx*ist.nspinngrid*ist.mn_states_tot], &ist, &parallel);
        }

        if (1 == flag.saveCheckpoints){
          if (parallel.mpi_rank == 0) write_separation(stdout, top);
          if (parallel.mpi_rank == 0) printf("****    CHECKPOINT %d *** CHECKPOINT %d ** CHECKPOINT %d *** CHECKPOINT %d    ****", par.checkpoint_id, par.checkpoint_id, par.checkpoint_id, par.checkpoint_id);
          if (parallel.mpi_rank == 0) write_separation(stdout, bottom);
          if (parallel.mpi_rank == 0) printf("\n"); fflush(stdout);

          save_job_state("checkpoint_1.dat",par.checkpoint_id,psitot,pot_local,ksqr,an,zn,ene_targets,nl,nlc,&grid,&ist,&par,&flag,&parallel);
        }
      // If checkpoint_id is 1, restart job after filtering (+ time reversal) but before orthogonalization
      case 1:
        if (1 == flag.restartFromOrtho){
        write_separation(stdout, top);
        printf("**** START FROM ORTHO *** START FROM ORTHO ** START FROM ORTHO *** START FROM ORTHO ****");
        write_separation(stdout, bottom); fflush(stdout);
        
        printf("\nLocal pseudopotential:\n");
        build_local_pot(pot_local, &pot, R, ksqr, atom, &grid, &ist, &par, &flag, &parallel);
        
        free(pot.r); pot.r = NULL; 
        free(pot.pseudo); pot.pseudo = NULL; 
        free(pot.dr); pot.dr = NULL; 
        free(pot.file_lens); pot.file_lens = NULL;
        if (1.0 != par.scale_surface_Cs){
          free(pot.r_LR); pot.r_LR = NULL;
          free(pot.pseudo_LR); pot.pseudo_LR = NULL; 
        }
        
        write_cube_file(pot_local, &grid, "localPot.cube");
        
        if(flag.SO==1) {
          printf("\nSpin-orbit pseudopotential:\n");
          init_SO_projectors(SO_projectors, &grid, R, atom, &ist, &par, &flag, &parallel);
        }
        /*** initialization for the non-local potential ***/
        if (flag.NL == 1){
          printf("\nNon-local pseudopotential:\n"); fflush(0);
          init_NL_projectors(nlc, nl, SO_projectors, &grid, R, atom, &ist, &par, &flag, &parallel);
        }
        // free memory allocated to SO_projectors
        if ( (flag.SO == 1) || (flag.NL == 1) ){
          free(SO_projectors); SO_projectors = NULL;
        }
        
        long psitot_size = par.t_rev_factor*ist.complex_idx*ist.nspinngrid*ist.mn_states_tot*sizeof(double);
        printf("\nNumber of states included for orthogonalization = %ld\n", ist.mn_states_tot); fflush(stdout);
        printf("Size of psitot array = %.2g GB\n", (double) psitot_size/1024/1024/1024); fflush(stdout);
        
        // Allocate psitot to have space for all the filtered states
        if ((psitot = (double *) calloc(ist.complex_idx * par.t_rev_factor * ist.nspinngrid * ist.mn_states_tot, sizeof(psitot[0]))) == NULL){
          if (parallel.mpi_rank == 0) fprintf(stderr,"\nOUT OF MEMORY: psitot\n\n"); exit(EXIT_FAILURE);
        }
        
        // Read in the states from psi-filt.dat file
        ppsi = fopen("psi-filt.dat", "r");
        if (ppsi != NULL){
          printf("Reading psi-filt.dat\n"); fflush(stdout);
          fread(&psitot[0], sizeof(double), ist.complex_idx*ist.nspinngrid*ist.mn_states_tot, ppsi);
          fclose(ppsi);
        } else{
          fprintf(stderr, "ERROR: psi-filt.dat could not be opened\n");
          exit(EXIT_FAILURE);
        }
        printf("psitot[max] = %lg\n", psitot[ist.complex_idx*ist.nspinngrid*ist.mn_states_tot - 1]);
        
        

        printf("\nNormalizing filtered states (for safety)\n"); fflush(stdout);
        normalize_all(psitot,ist.mn_states_tot,&ist,&par,&flag,&parallel);
        if (2 == par.t_rev_factor){
          printf("\nTime-reversing all filtered states (doubles number of orthogonal states)\n"); fflush(stdout);
          time_reverse_all(&psitot[0], &psitot[ist.complex_idx*ist.nspinngrid*ist.mn_states_tot], &ist, &parallel);
        }
      } else if (flag.restartFromCheckpoint == 1){
          par.checkpoint_id = flag.restartFromCheckpoint;
          if (parallel.mpi_rank == 0) write_separation(stdout, top);
          if (parallel.mpi_rank == 0) printf("****    CHECKPOINT %d *** CHECKPOINT %d ** CHECKPOINT %d *** CHECKPOINT %d    ****", par.checkpoint_id, par.checkpoint_id, par.checkpoint_id, par.checkpoint_id);
          if (parallel.mpi_rank == 0) write_separation(stdout, bottom); fflush(stdout);
          
          restart_from_save("checkpoint_1.dat",par.checkpoint_id,psitot,pot_local,ksqr,an,zn,ene_targets,nl,nlc,&grid,&ist,&par,&flag,&parallel);
        }
        par.checkpoint_id++;
        /*************************************************************************/
        /*** orthogonalize and normalize the filtered states using an svd routine ***/
        omp_set_num_threads(parallel.nthreads);

        if (parallel.mpi_rank == 0) write_separation(stdout, top);
        current_time = time(NULL);
        c_time_string = ctime(&current_time);
        if (parallel.mpi_rank == 0) printf("\n5. ORTHOGONALIZATING FILTERED STATES | %s\n", c_time_string); 
        if (parallel.mpi_rank == 0) write_separation(stdout, bottom); fflush(stdout);

        inital_clock_t = (double)clock(); initial_wall_t = (double)time(NULL);
        if (parallel.mpi_rank == 0) printf("mn_states_tot before ortho = %ld\n", ist.mn_states_tot);
        if (1 == flag.isComplex){
          ist.mn_states_tot = ortho((MKL_Complex16*)psitot, grid.dv, &ist, &par, &flag, &parallel);      
        } else if (0 == flag.isComplex) {
          ist.mn_states_tot = ortho(psitot, grid.dv, &ist, &par, &flag, &parallel);
        }
        if (parallel.mpi_rank == 0) printf("mn_states_tot after ortho = %ld\n", ist.mn_states_tot);
        psitot = (double *) realloc(psitot, ist.mn_states_tot * ist.nspinngrid * ist.complex_idx * sizeof(psitot[0]) );

        normalize_all(&psitot[0], ist.mn_states_tot, &ist, &par, &flag, &parallel);
        
        if (parallel.mpi_rank == 0) printf("\ndone calculating ortho, CPU time (sec) %g, wall run time (sec) %g\n",
                  ((double)clock()-inital_clock_t)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-initial_wall_t); 
        fflush(stdout);

        if (1 == flag.saveCheckpoints){
          if (parallel.mpi_rank == 0) write_separation(stdout, top);
          if (parallel.mpi_rank == 0) printf("****    CHECKPOINT %d *** CHECKPOINT %d ** CHECKPOINT %d *** CHECKPOINT %d    ****", par.checkpoint_id, par.checkpoint_id, par.checkpoint_id, par.checkpoint_id);
          if (parallel.mpi_rank == 0) write_separation(stdout, bottom); 
          if (parallel.mpi_rank == 0) printf("\n"); fflush(stdout);

          save_job_state("checkpoint_2.dat",par.checkpoint_id,psitot,pot_local,ksqr,an,zn,ene_targets,nl,nlc,&grid,&ist,&par,&flag,&parallel);
        }
      // If checkpoint_id is 2, restart job right after orthogonalization, before diagonalization
      case 2:
        if (flag.restartFromCheckpoint == 2){
          par.checkpoint_id = flag.restartFromCheckpoint;
          if (parallel.mpi_rank == 0) write_separation(stdout, top);
          if (parallel.mpi_rank == 0) printf("****    CHECKPOINT %d *** CHECKPOINT %d ** CHECKPOINT %d *** CHECKPOINT %d    ****", par.checkpoint_id, par.checkpoint_id, par.checkpoint_id, par.checkpoint_id);
          if (parallel.mpi_rank == 0) write_separation(stdout, bottom); fflush(stdout);
          
          restart_from_save("checkpoint_2.dat",par.checkpoint_id,psitot,pot_local,ksqr,an,zn,ene_targets,nl,nlc,&grid,&ist,&par,&flag,&parallel);
        }
        par.checkpoint_id++;
        /***********************************************************************/
        /*** diagonalize the hamiltonian in the subspace spanned by the ***/
        /*** orthogonal filtered states, generating the eigenstates of the ***/
        /*** hamiltonian within the desired energy range ***/
        if (parallel.mpi_rank == 0) write_separation(stdout, top);
        current_time = time(NULL);
        c_time_string = ctime(&current_time);
        if (parallel.mpi_rank == 0) printf("\n6. DIAGONALIZING HAMILTONIAN | %s\n", c_time_string); 
        if (parallel.mpi_rank == 0) write_separation(stdout, bottom); fflush(stdout);
        
        inital_clock_t = (double)clock(); initial_wall_t = (double)time(NULL);
        diag_H(psitot,pot_local,nlc,nl,ksqr,eig_vals,&ist,&par,&flag,&parallel,planfw,planbw,fftwpsi);
        normalize_all(&psitot[0],ist.mn_states_tot,&ist,&par,&flag,&parallel);
        jms = ist.mn_states_tot;
        if (parallel.mpi_rank == 0) printf("\ndone calculating Hmat, CPU time (sec) %g, wall run time (sec) %g\n",
                  ((double)clock()-inital_clock_t)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-initial_wall_t);
        fflush(stdout);

        /***********************************************************************/
        /*** calculate the standard deviation of these states ***/
        /*** this is used to check if there are ghost states ***/
        if (parallel.mpi_rank == 0) write_separation(stdout, top);
        current_time = time(NULL);
        c_time_string = ctime(&current_time);
        if (parallel.mpi_rank == 0) printf("\n7. CALCULATING VARIANCE OF EIGENVALUES | %s\n", c_time_string);
        if (parallel.mpi_rank == 0) write_separation(stdout, bottom); fflush(stdout);
        
        calc_sigma_E(psitot, pot_local, nlc, nl, ksqr, sigma_E, &ist, &par, &flag);

        /*** write the eigenstates/energies to a file ***/
        if (flag.getAllStates == 1) {if (parallel.mpi_rank == 0) printf("getAllStates flag on\nWriting all eigenstates to disk\n");}
        else{ if (parallel.mpi_rank == 0) printf("getAllStates flag off\nWriting eigenstates with sigE < 0.1 to disk\n");}
        
        if (flag.getAllStates == 1){
            /*** write all eigenvalues and their standard deviation to a file ***/
            if ((peig = fopen("eval.dat" , "w"))==NULL){
              if (parallel.mpi_rank == 0) printf("Out of disk space!\n Writing energy levels to stdout...\n\n\n\n");
              for (jms = 0; jms < ist.mn_states_tot; jms++) if (parallel.mpi_rank == 0) printf ("%ld %.16g %g\n", jms, eig_vals[jms], sigma_E[jms]); 
            }
            else{
              for (jms = 0; jms < ist.mn_states_tot; jms++) if (parallel.mpi_rank == 0) fprintf (peig,"%ld %.16g %g\n", jms, eig_vals[jms], sigma_E[jms]); 
              fclose(peig);
            }

            // Write all eigenstates to disk!
            if((ppsi = fopen("psi.dat" , "w"))==NULL){if (parallel.mpi_rank == 0) printf("Out of disk space!");}
            else {
              fwrite (&psitot[0],sizeof(psitot[0]),ist.complex_idx*ist.mn_states_tot*ist.nspinngrid,ppsi);
              fclose(ppsi);
            }
          }
          // Write output in the case that only converged eigenstates will be printed
          else{
            long eig = 0;
            // Write only the eigenstates having variance less than sigma_E_cut to disk (remove ghost states)
            if (((peig = fopen("eval.dat" , "w"))==NULL) || ((ppsi = fopen("psi.dat" , "w"))==NULL)) {
              if (parallel.mpi_rank == 0) printf("Out of disk space!\n Writing energy levels to stdout...\n\n\n\n");
              for (jms = 0; jms < ist.mn_states_tot; jms++) if (parallel.mpi_rank == 0) printf ("%ld %.16g %g\n", jms, eig_vals[jms], sigma_E[jms]); 
            }
            else{
              for (jms = 0; jms < ist.mn_states_tot; jms++){
                if (sigma_E[jms] < par.sigma_E_cut){
                  // First, move the file pointer so it writes all the eigenstates contiguously in psi.dat
                  fseek(ppsi, ist.complex_idx * eig * ist.nspinngrid * sizeof(double), SEEK_SET);
                  // Then, write the jms'th eigenstate to psi.dat by writing nspinngrid doubles out of psitot
                  fwrite(&psitot[ist.complex_idx*jms*ist.nspinngrid], sizeof(double), ist.complex_idx*ist.nspinngrid, ppsi);
                  // Modify psitot to only contain the eigenstates with small variance
                  for (jgrid = 0; jgrid < ist.nspinngrid; jgrid++){
                    jgrid_real = ist.complex_idx * jgrid;
                    jgrid_imag = ist.complex_idx * jgrid + 1;

                    psitot[ist.complex_idx*eig*ist.nspinngrid + jgrid_real] = psitot[ist.complex_idx*jms*ist.nspinngrid + jgrid_real];
                    if (1 == flag.isComplex){
                      psitot[ist.complex_idx*eig*ist.nspinngrid + jgrid_imag] = psitot[ist.complex_idx*jms*ist.nspinngrid + jgrid_imag];
                    }
                  }
                  // Write the jms'th eigenenergy to eval.dat
                  if (parallel.mpi_rank == 0) fprintf (peig,"%ld %.16g %g\n", jms, eig_vals[jms], sigma_E[jms]); 
                  // Modify the eig_vals and sigma_E to carry the values for the appropriate small variance eigenstate
                  eig_vals[eig] = eig_vals[jms];
                  sigma_E[eig] = sigma_E[jms];
                  eig++;
                }
                else{ continue;}
              } 
            }
          // Reset the mn_states_tot value if we choose not to get all the states.
          ist.mn_states_tot = eig; 
          if (parallel.mpi_rank == 0) printf("Number of eigenstates with sigE < 0.1 saved to disk: new_mstot = %ld\n", ist.mn_states_tot);
        }

        // If no eigenstates were obtained, then the job has failed and alternate options should be pursued
        if (ist.mn_states_tot == 0){
          if (parallel.mpi_rank == 0) printf("The filter diagonalization yielded %ld eigenstates\n", ist.mn_states_tot);
          if ( (flag.retryFilter == 1) && (flag.alreadyTried == 0) ){
            if (parallel.mpi_rank == 0) printf("The filter diagonalization will be reattempted with increased filters and Chebyshev length\n");
            ist.n_filter_cycles += 16;
            ist.ncheby += 1024;
            if (parallel.mpi_rank == 0) printf("New n_filter_cycles = %ld, new ncheby = %ld\n", ist.n_filter_cycles, ist.ncheby);
            flag.alreadyTried = 1; // prevent the filtering algorithm from repeating indefinitely
          } else{
            if (parallel.mpi_rank == 0) fprintf(stderr, "The filter diagonalization yielded %ld eigenstates (T_T)\n", ist.mn_states_tot);
            if (parallel.mpi_rank == 0) fprintf(stderr, "Job failed :( exiting reluctantly\n");
            exit(EXIT_FAILURE);
          }
        }


        // The standard portion of the filter diagonalization procedure has concluded
        // Print the output.dat file for use in future workflows unless explicitly requested to NOT save output
        if (0 != flag.saveOutput){
          save_output("output.dat",psitot,eig_vals,sigma_E,R,&grid,&ist,&par,&flag,&parallel);
        }

        /**************************************************************************************************************/
        /**************************************************************************************************************/
        /**************************************************************************************************************/

        if (parallel.mpi_rank == 0) write_separation(stdout, top);
        current_time = time(NULL);
        c_time_string = ctime(&current_time);
        if (parallel.mpi_rank == 0) printf("\nCALCULATING OPTIONAL OUTPUT | %s\n", c_time_string); 
        if (parallel.mpi_rank == 0) write_separation(stdout, bottom); fflush(stdout);

        long i, a, ieof, nval;
        char str[50];
        double evalloc, deloc;
        FILE *pf;

        ist.homo_idx = ist.lumo_idx = 0;
        pf = fopen("eval.dat" , "r");
        for (i = ieof = 0; ieof != EOF; i++){
          ieof = fscanf(pf, "%ld %lg %lg", &a, &evalloc, &deloc);
          if (deloc < par.sigma_E_cut && evalloc < par.fermi_E) ist.homo_idx = i;
          if (i > ist.mn_states_tot){
          if (parallel.mpi_rank == 0) printf("No hole states converged to within %lg a.u.\n", par.sigma_E_cut);
          break;
        }
        }
        fclose(pf);

        // nval = i - 1;
        pf = fopen("eval.dat" , "r");
        for (i = 0; i <= ist.homo_idx; i++) {
          fscanf(pf, "%ld %lg %lg", &a, &evalloc, &deloc);
          if (i > ist.mn_states_tot){
            if (parallel.mpi_rank == 0) printf("No electron states converged to within %lg a.u.\n", par.sigma_E_cut);
            break;
          }
        }
        for (i = ist.homo_idx+1, ieof = 0; ieof != EOF; i++) {
          fscanf(pf, "%ld %lg %lg", &a, &evalloc, &deloc);
          if (deloc < par.sigma_E_cut) {
            ist.lumo_idx = i;
            break;
          }
          if (i > ist.mn_states_tot){
            if (parallel.mpi_rank == 0) printf("No electron states converged to within %lg a.u.\n", par.sigma_E_cut);
            break;
          }
        }
        fclose(pf);

        if (parallel.mpi_rank == 0) printf("index of homo, homo_idx = %ld; index of lumo, lumo_idx = %ld\n", ist.homo_idx, ist.lumo_idx); fflush(0);
        // Set the total number of electron and hole states in order to calculate the potential overlap integrals
        ist.total_homo = ist.homo_idx + 1; ist.total_lumo = ist.mn_states_tot - ist.total_homo;
        if (parallel.mpi_rank == 0) printf("total_homo = %ld total_lumo = %ld\n", ist.total_homo, ist.total_lumo); fflush(0);


        if (flag.printCubes == 1){
          /*** Write homo and lumo cube files ***/
          
          if (parallel.mpi_rank == 0) write_separation(stdout, top);
          if (parallel.mpi_rank == 0) printf("\nWRITING CUBE FILES\n"); 
          if (parallel.mpi_rank == 0) write_separation(stdout, bottom); fflush(stdout);

          if ((ist.homo_idx == 0) || (ist.lumo_idx == 0)){
            if (parallel.mpi_rank == 0) printf("\nDid not converge enough electron or hole states to visualize cube files.\n");
          } else{
          if ((rho = (double *) calloc(ist.ngrid, sizeof(double))) == NULL){
            if (parallel.mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: rho\n\n"); exit(EXIT_FAILURE);
          }

          inital_clock_t = (double)clock(); initial_wall_t = (double)time(NULL);

          for (i = 0; (i < ist.total_homo) && (i < ist.ncubes); i++){
            //Spin Up Wavefunction
            sprintf(str,"homo-%ld-Up.cube",i);
            for (jgrid = 0; jgrid < ist.ngrid; jgrid++){
              jgrid_real = ist.complex_idx * jgrid;
              jgrid_imag = ist.complex_idx * jgrid + 1;
              
              rho[jgrid] = sqr(psitot[ist.complex_idx*(ist.homo_idx-i)*ist.nspinngrid + jgrid_real]);
              if (1 == flag.isComplex) rho[jgrid] += sqr(psitot[ist.complex_idx*(ist.homo_idx-i)*ist.nspinngrid + jgrid_imag]);
            }
            write_cube_file(rho, &grid, str);
            //Spin Down Wavefunction
            if (1 == flag.useSpinors){    
              sprintf(str,"homo-%ld-Dn.cube", i);
              for (jgrid = 0; jgrid < ist.ngrid; jgrid++){
                jgrid_real = ist.complex_idx * jgrid;
                jgrid_imag = ist.complex_idx * jgrid + 1;
                
                rho[jgrid] = sqr(psitot[ist.complex_idx*((ist.homo_idx-i)*ist.nspinngrid+ist.ngrid)+jgrid_real]) 
                    + sqr(psitot[ist.complex_idx*((ist.homo_idx-i)*ist.nspinngrid+ist.ngrid)+jgrid_imag]);    
              }
              write_cube_file(rho, &grid, str);
            } 
          }

          for (i = 0;  (i < ist.total_lumo) && (i < ist.ncubes); i++){
            sprintf(str,"lumo+%ld-Up.cube",i);
            for (jgrid = 0; jgrid < ist.ngrid; jgrid++){
              jgrid_real = ist.complex_idx * jgrid;
              jgrid_imag = ist.complex_idx * jgrid + 1;
              
              rho[jgrid] = sqr(psitot[ist.complex_idx*(ist.lumo_idx+i)*ist.nspinngrid + jgrid_real]);
              if (1 == flag.isComplex) rho[jgrid] += sqr(psitot[ist.complex_idx*(ist.lumo_idx+i)*ist.nspinngrid + jgrid_imag]);
            }
            write_cube_file(rho, &grid, str);

            if (1 == flag.useSpinors){
              sprintf(str,"lumo+%ld-Dn.cube",i);
              for (jgrid = 0; jgrid < ist.ngrid; jgrid++){
                jgrid_real = ist.complex_idx * jgrid;
                jgrid_imag = ist.complex_idx * jgrid + 1;
              
                rho[jgrid] = sqr(psitot[ist.complex_idx*((ist.lumo_idx+i)*ist.nspinngrid+ist.ngrid)+jgrid_real]) 
                    + sqr(psitot[ist.complex_idx*((ist.lumo_idx+i)*ist.nspinngrid+ist.ngrid)+jgrid_imag]);
              }
              write_cube_file(rho, &grid, str);
            }
          }
          free(rho);

          if (parallel.mpi_rank == 0) printf("\ndone calculating cubes, CPU time (sec) %g, wall run time (sec) %g\n",
            ((double)clock()-inital_clock_t)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-initial_wall_t);
          }
        }
        if (flag.calcPotOverlap == 1){

          if (parallel.mpi_rank == 0) write_separation(stdout, top);
          current_time = time(NULL);
          c_time_string = ctime(&current_time);
          if (parallel.mpi_rank == 0) printf("\nCALCULATING POTENTIAL MATRIX ELEMENTS | %s\n", c_time_string); fflush(0);
          if (parallel.mpi_rank == 0) write_separation(stdout, bottom); fflush(stdout);

          calc_pot_overlap(&psitot[0], pot_local, nlc, nl, eig_vals, &par, &ist, &flag);
        } 
        if (flag.calcSpinAngStat == 1) {
          if (parallel.mpi_rank == 0) write_separation(stdout, top);
          current_time = time(NULL);
          c_time_string = ctime(&current_time);
          if (parallel.mpi_rank == 0) printf("\nCALCULATING SPIN & ANG. MOM. STATISTICS | %s\n", c_time_string); 
          if (parallel.mpi_rank == 0) write_separation(stdout, bottom); fflush(stdout);

          calc_angular_exp(psitot, &grid,0, ist.mn_states_tot, &ist, &par, &flag, &parallel, planfw, planbw, fftwpsi);
        
        } 
        if ( (flag.printCubes != 1) && (flag.calcPotOverlap != 1) && (flag.calcSpinAngStat != 1) ) {
          if (parallel.mpi_rank == 0) printf("\nNo optional output requested.\n");
        }
      }
      /*************************************************************************/
      /*** free memory ***/
      free(ist.atom_types);
      free(grid.x); free(grid.y); free(grid.z); 
      free(R); free(atom);
      free(psi); free(phi); 
      free(pot_local); free(ksqr);
      free(nlc); free(nl); 
      free(eig_vals); free(ene_targets); free(sigma_E);
      free(an); free(zn);      
      
      fftw_destroy_plan(planfw);
      fftw_destroy_plan(planbw);
      fftw_free(fftwpsi);
      fftw_cleanup_threads();
    
      time_t end_time = time(NULL);
      time_t end_clock = clock();
      
      if (0 == parallel.mpi_rank){
        free(psitot); 
        double elapsed_seconds = difftime(end_time, start_time);
        char* duration = format_duration(elapsed_seconds);
        if (parallel.mpi_rank == 0) write_separation(stdout, top);
        printf("\nDONE WITH PROGRAM: FILTER DIAGONALIZATION\n");
        printf("This calculation ended at: %s\n", ctime(&end_time)); 
        printf("Total job CPU time (sec) %.4g | %s\n",
                  ((double)end_clock - (double)start_clock)/(double)(CLOCKS_PER_SEC), format_duration(((double)end_clock - (double)start_clock)/(double)(CLOCKS_PER_SEC)) );fflush(0);
        printf("Total wall run time (sec) %.4g | %s", (double)end_time - (double)start_time, duration);
        write_separation(stdout, bottom);
        
        free(top); free(bottom);
      }
      
      MPI_Finalize(); // Finalize the MPI tasks and prepare to exit program
      
      exit(0);
      
  } // End of switch statement
} // End of main

/*****************************************************************************/
