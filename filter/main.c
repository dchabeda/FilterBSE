#include "fd.h"

/*****************************************************************************/

int main(int argc, char *argv[]){
  /*****************************************************************
  * This is the main function for filter.x. It is the driver that  *
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
  // FFT
  fftw_plan_loc planfw, planbw; fftw_complex *fftwpsi; 
  long fft_flags=0;
  // custom structs 
  flag_st flag; index_st ist; par_st par; atom_info *atom; 
  pot_st pot; grid_st grid; xyz_st *R; nlc_st *nlc = NULL; 
  parallel_st parallel; 
  // double arrays
  double *psitot;
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
  char *top; top = malloc(2*sizeof(top[0])); 
  char *bottom; bottom = malloc(2*sizeof(bottom[0]));
  strcpy(top, "T\0"); strcpy(bottom, "B\0");
  
  fprintf(stdout, "******************************************************************************\n");
  printf("\nRUNNING PROGRAM: FILTER DIAGONALIZATION\n");
  printf("This calculation began at: %s", ctime(&start_time)); 
  write_separation(stdout, bottom);
  fflush(stdout);

  /*************************************************************************/
  // Initialize job from input file
  
  write_separation(stdout, top);
  printf("\n1.\tINITIALIZING JOB\n");
  write_separation(stdout, bottom); fflush(stdout);

  /*** read initial setup from input.par ***/
  printf("\nReading job specifications from input.par:\n");
  read_input(&flag, &grid, &ist, &par, &parallel);

  /*** allocating memory ***/
  // the positions of the atoms in the x, y, and z directions 
  if ((R = (xyz_st *) calloc(ist.natoms, sizeof(xyz_st))) == NULL) {
    fprintf(stderr, "\nOUT OF MEMORY: R array\n\n"); exit(EXIT_FAILURE);
  }
  // the atom specific information 
  if ((atom = (atom_info *) calloc(ist.natoms, sizeof(atom_info))) == NULL){
    fprintf(stderr, "\nOUT OF MEMORY: atom struct\n\n"); exit(EXIT_FAILURE);
  }
  // the energies of each energy target in filter
  if ((ene_targets = (double *) calloc(ist.m_states_per_filter, sizeof(double))) == NULL){
    fprintf(stderr, "\nOUT OF MEMORY: ene_targets\n\n"); exit(EXIT_FAILURE);
  }
  
  /*** read the nanocrystal configuration ***/
  printf("\nReading atomic configuration from conf.par:\n");
  read_conf(R, atom, &ist, &par, &flag);

  /*** initialize parameters for the grid ***/
  printf("\nInitializing the grid parameters:\n");
  init_grid_params(&grid, R, &ist, &par);

  // Allocate memory for the grid in the x, y, and z directions ***/
  if ((grid.x = (double *) calloc(grid.nx, sizeof(double))) == NULL){
    fprintf(stderr, "\nOUT OF MEMORY: grid.x\n\n"); exit(EXIT_FAILURE);
  }
  if ((grid.y = (double *) calloc(grid.ny, sizeof(double))) == NULL){
    fprintf(stderr, "\nOUT OF MEMORY: grid.y\n\n"); exit(EXIT_FAILURE);
  }
  if ((grid.z = (double *) calloc(grid.nz, sizeof(double))) == NULL){
    fprintf(stderr, "\nOUT OF MEMORY: grid.z\n\n"); exit(EXIT_FAILURE);
  }
  // the kinetic energy stored on the grid
  if ((ksqr = (double *) calloc(ist.ngrid, sizeof(double))) == NULL){
    fprintf(stderr, "\nOUT OF MEMORY: ksqr\n\n"); exit(EXIT_FAILURE);
  }

  /*** build the real- and k-space grids ***/
  printf("\nBuilding the real-space and k-space grids:\n");
  build_grid_ksqr(ksqr, R, &grid, &ist, &par, &flag);
  
  /*** set the energy targets ***/
  printf("\nSetting the filter energy targets:\n");
  set_ene_targets(ene_targets, &ist, &par, &flag);

  
  /*************************************************************************/
  /*** allocating memory for the rest of the program ***/
  printf("\nAllocating memory for FFT, pot, psi, eig_vals...");
  
  // FFT
  fftwpsi = fftw_malloc(sizeof(fftw_complex) * ist.ngrid);
  /*** initialization for the fast Fourier transform ***/
  planfw = fftw_plan_dft_3d(grid.nz, grid.ny, grid.nx, fftwpsi, fftwpsi, FFTW_FORWARD, fft_flags);
  planbw = fftw_plan_dft_3d(grid.nz, grid.ny, grid.nx, fftwpsi, fftwpsi, FFTW_BACKWARD, fft_flags);
  
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
    fprintf(stderr, "\nOUT OF MEMORY: psi\n\n"); exit(EXIT_FAILURE);
  }
  if ((phi = (zomplex *)calloc(ist.nspinngrid, sizeof(zomplex))) == NULL){
    fprintf(stderr, "\nOUT OF MEMORY: phi\n\n"); exit(EXIT_FAILURE);
  }
  if ((pot_local = (double *) calloc(ist.ngrid, sizeof(double))) == NULL){
    fprintf(stderr, "\nOUT OF MEMORY: pot_local\n\n"); exit(EXIT_FAILURE);
  }
  
  // For Newton interpolation coefficients
  an = (zomplex *) calloc(ist.ncheby*ist.m_states_per_filter, sizeof(zomplex));
  zn = (double *) calloc(ist.ncheby, sizeof(double));
  
  // The value of isComplex is 1 if wavefunctions are real valued, 2 if functions are complex valued
  // The factor of par.t_rev_factor (2 w spinors, 1 w/o spinors) in the psitot memory allocation
  // is because we time reverse the spinors to get double the orthogonal states
  // Spinor calculations are 8 times more memory intensive than scalar calculations
  if ((psitot = (double *) calloc(ist.complex_idx * par.t_rev_factor * ist.nspinngrid * ist.mn_states_tot, sizeof(psitot[0]))) == NULL){
    fprintf(stderr,"\nOUT OF MEMORY: psitot\n\n"); exit(EXIT_FAILURE);
  }
  // the quasiparticle energies and standard deviations
  if ((eig_vals = (double *) calloc(par.t_rev_factor*ist.mn_states_tot, sizeof(double))) == NULL){
    fprintf(stderr, "\nOUT OF MEMORY: eig_vals\n\n"); exit(EXIT_FAILURE);
  }
  if ((sigma_E = (double *) calloc(par.t_rev_factor*ist.mn_states_tot, sizeof(double))) == NULL){
    fprintf(stderr, "\nOUT OF MEMORY: sigma_E\n\n"); exit(EXIT_FAILURE);
  }
  // memory allocation for the spin-orbit potential 
  if ( (flag.SO == 1) || (flag.NL == 1) ){
    if ((SO_projectors = (double*) calloc(PROJ_LEN * ist.nproj, sizeof(double)))==NULL){nerror("mem_SO_projector");}  
  }
  // memory allocation for the non-local potential
  if (1 == flag.NL){
    if ((nlc = (nlc_st *) calloc(ist.n_NL_atoms*ist.n_NL_gridpts, sizeof(nlc_st))) == NULL){ 
      fprintf(stderr, "\nOUT OF MEMORY: nlc\n\n"); exit(EXIT_FAILURE);
    }
    if ((nl = (long *) calloc(ist.natoms, sizeof(nl[0]))) == NULL) {
      fprintf(stderr, "\nOUT OF MEMORY: nl\n\n"); exit(EXIT_FAILURE);
    }
  }
  
  printf("\tdone allocating memory.\n"); fflush(stdout);

  // The filter code supports restarting the job from a saved state. See save.c for formatting
  // of checkpoint files. See read_input in read.c for specifying the checkpoint restart
  switch (flag.restartFromCheckpoint){
    // If no restart requested, the job will run from the beginning with no fast-forward 
    case -1:
      /**************************************************************************/
      printf("\nInitializing potentials...\n");
      
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
        init_SO_projectors(SO_projectors, &grid, R, atom, &ist, &par);
      }
      /*** initialization for the non-local potential ***/
      if (flag.NL == 1){
        printf("\nNon-local pseudopotential:\n"); fflush(0);
        init_NL_projectors(nlc, nl, SO_projectors, &grid, R, atom, &ist, &par, &flag);
      }
      // free memory allocated to SO_projectors
      if ( (flag.SO == 1) || (flag.NL == 1) ){
        free(SO_projectors); SO_projectors = NULL;
      }
      
       
      /**************************************************************************/
      /*** calculate the energy range of the hamitonian ***/
      write_separation(stdout, top);
      printf("\n2. CALCULATING HAMILTONIAN ENERGY RANGE...\n");
      write_separation(stdout, bottom); fflush(stdout);

      inital_clock_t = (double)clock(); 
      initial_wall_t = (double)time(NULL);
      
      get_energy_range(psi, phi, pot_local, &grid, nlc, nl, ksqr, &ist, &par, &parallel, &flag, planfw, planbw, fftwpsi);
      
      printf("\ndone calculate energy range, CPU time (sec) %g, wall run time (sec) %g\n",
                ((double)clock()-inital_clock_t)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-initial_wall_t); 
      fflush(stdout);
    
      /**************************************************************************/
      // FILTER ALGORITHM
      /**************************************************************************/
      write_separation(stdout, top);
      printf("\n3. GENERATING COEFFICIENTS\n");
      write_separation(stdout, bottom); fflush(stdout);
      
      /*** set parameters for the newton interpolation ***/
      par.dt = sqr((double)(ist.ncheby) / (2.5*par.dE));

      printf("\nCalculating Newton interpolation coefficients for filter functions.\n"); fflush(0);
      gen_newton_coeff(an, zn, ene_targets, &ist, &par, &parallel);

      pseed = fopen("seed.dat", "w");
      fprintf(pseed, "idx  seed\n");
      if (flag.setSeed == 1){
        rand_seed = - par.rand_seed;} 
      else {
        Randomize();  rand_seed = -random();
      }

      printf("ncheby = %ld dt = %g dE = %g\n", ist.ncheby, par.dt, par.dE);
      
      // Save checkpoint after computing the potential on the grid and Newton interpolation coeffs
      if (1 == flag.saveCheckpoints){
        par.checkpoint_id = 0;
        write_separation(stdout, top);
        printf("****    CHECKPOINT %d *** CHECKPOINT %d ** CHECKPOINT %d *** CHECKPOINT %d    ****", par.checkpoint_id, par.checkpoint_id, par.checkpoint_id, par.checkpoint_id);
        write_separation(stdout, bottom); 
        printf("\n"); fflush(stdout);

        save_job_state("checkpoint_0.dat",par.checkpoint_id,psitot,pot_local,ksqr,an,zn,ene_targets,nl,nlc,&grid,&ist,&par,&flag,&parallel);
      }
    
    // If checkpoint_id is 0, then restart job right before filtering step
    case 0:
      if (flag.restartFromCheckpoint == 0){
        par.checkpoint_id = flag.restartFromCheckpoint;
        write_separation(stdout, top);
        printf("****    CHECKPOINT %d *** CHECKPOINT %d ** CHECKPOINT %d *** CHECKPOINT %d    ****", par.checkpoint_id, par.checkpoint_id, par.checkpoint_id, par.checkpoint_id);
        write_separation(stdout, bottom); fflush(stdout);
        
        restart_from_save("checkpoint_0.dat",par.checkpoint_id,psitot,pot_local,ksqr,an,zn,ene_targets,nl,nlc,&grid,&ist,&par,&flag,&parallel);
      } 
      par.checkpoint_id++;
      /**************************************************************************/
      /*** start filtering loop.  we run over n_filter_cycles cycles and calculate ***/
      /*** m_states_per_filter filtered states at each cycle ***/
      write_separation(stdout, top);
      printf("\n4. STARTING FILTERING\n");
      write_separation(stdout, bottom); 

      printf("\nEnergy width, sigma, of filter function = %.6g a.u.\n", sqrt(1 / (2*par.dt)));
      printf("Suggested max span of spectrum for filtering = %.6g a.u.\n", ist.m_states_per_filter * sqrt(1 / (2*par.dt)));
      printf("Requested span of spectrum to filter = %.6g a.u.\n", (ene_targets[par.n_targets_VB-1] - ene_targets[0]) + (ene_targets[par.n_targets_VB+par.n_targets_CB-1] - ene_targets[par.n_targets_VB]));
      fflush(stdout);

      for (jns = 0; jns < ist.n_filter_cycles; jns++) {
        for (jspin = 0; jspin < ist.nspin; jspin++) {
          fprintf(pseed, "%ld %ld\n", ist.nspin*jns + jspin, rand_seed);
          init_psi(&psi[jspin*ist.ngrid], &rand_seed, flag.isComplex, &grid, &parallel);
        }

        for (jms = 0; jms < ist.m_states_per_filter; jms++) {
          for (jgrid = 0; jgrid < ist.nspinngrid; jgrid++) {
            // handle indexing of real and imaginary components if complex
            // ist.complex_idx = (flag.isComplex + 1) = 2 when complex valued functions are in use, 1 if real
            jgrid_real = ist.complex_idx * jgrid;
            jgrid_imag = ist.complex_idx * jgrid + 1;
            jstate = ist.complex_idx * ( jns*ist.m_states_per_filter*ist.nspinngrid + ist.nspinngrid*jms );
            
            // if the wavefunction is real-valued, then only the real component is stored
            psitot[jstate + jgrid_real] = psi[jgrid].re;
            
            // if the wavefunction is complex-valued, then the imaginary component should also be stored
            if (1 == flag.isComplex){
              psitot[jstate + jgrid_imag] = psi[jgrid].im;
              // the imaginary components will be stored one double away (16 bytes, or IMAG_IDX) in memory from the real component
            }
          }
        }
      }
      fclose(pseed);

      inital_clock_t = (double)clock(); 
      initial_wall_t = (double)time(NULL);
      omp_set_dynamic(0);
      omp_set_num_threads(parallel.nthreads);
    #pragma omp parallel for private(jns, thread_id)
      for (jns = 0; jns < ist.n_filter_cycles; jns++) {
        thread_id = omp_get_thread_num();	
        run_filter_cycle(&psitot[ist.complex_idx*jns*ist.m_states_per_filter*ist.nspinngrid], pot_local, nlc, nl, ksqr, an, zn,\
        ene_targets, thread_id, jns, &grid, &ist, &par, &flag, &parallel);
      } 
      printf("\ndone calculating filter, CPU time (sec) %g, wall run time (sec) %g\n",
                ((double)clock()-inital_clock_t)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-initial_wall_t); 
      fflush(stdout);

      /*************************************************************************/
      /*** read all filtered states ***/
      /*ppsi = fopen("psi-filt.dat" , "r");
      fread (psitot,sizeof(zomplex),ist.mn_states_tot*ist.ngrid,ppsi);
      fclose(ppsi);*/

      if (1 == flag.useSpinors){
        printf("\nTime-reversing all filtered states (doubles number of orthogonal states)");
        time_reverse_all(&psitot[0], &psitot[ist.complex_idx*ist.nspinngrid*ist.mn_states_tot], &ist, &parallel);
      }

      if (1 == flag.saveCheckpoints){
        write_separation(stdout, top);
        printf("****    CHECKPOINT %d *** CHECKPOINT %d ** CHECKPOINT %d *** CHECKPOINT %d    ****", par.checkpoint_id, par.checkpoint_id, par.checkpoint_id, par.checkpoint_id);
        write_separation(stdout, bottom);
        printf("\n"); fflush(stdout);

        save_job_state("checkpoint_1.dat",par.checkpoint_id,psitot,pot_local,ksqr,an,zn,ene_targets,nl,nlc,&grid,&ist,&par,&flag,&parallel);
      }
    // If checkpoint_id is 1, restart job after filtering (+ time reversal) but before orthogonalization
    case 1:
      if (flag.restartFromCheckpoint == 1){
        par.checkpoint_id = flag.restartFromCheckpoint;
        write_separation(stdout, top);
        printf("****    CHECKPOINT %d *** CHECKPOINT %d ** CHECKPOINT %d *** CHECKPOINT %d    ****", par.checkpoint_id, par.checkpoint_id, par.checkpoint_id, par.checkpoint_id);
        write_separation(stdout, bottom); fflush(stdout);
        
        restart_from_save("checkpoint_1.dat",par.checkpoint_id,psitot,pot_local,ksqr,an,zn,ene_targets,nl,nlc,&grid,&ist,&par,&flag,&parallel);
      }
      par.checkpoint_id++;
      /*************************************************************************/
      /*** orthogonalize and normalize the filtered states using an svd routine ***/
      
      write_separation(stdout, top);
      printf("\n5. ORTHOGONALIZATING FILTERED STATES\n"); 
      write_separation(stdout, bottom); fflush(stdout);

      inital_clock_t = (double)clock(); initial_wall_t = (double)time(NULL);
      printf("mn_states_tot before ortho = %ld\n", ist.mn_states_tot);
      if (1 == flag.isComplex){
        ist.mn_states_tot = ortho((MKL_Complex16 *)psitot, grid.dv, &ist, &par, &flag);      
      } else if (0 == flag.isComplex) {
        ist.mn_states_tot = ortho(psitot, grid.dv, &ist, &par, &flag);
      }
      printf("mn_states_tot after ortho = %ld\n", ist.mn_states_tot);
      psitot = (double *) realloc(psitot, ist.mn_states_tot * ist.nspinngrid * ist.complex_idx * sizeof(psitot[0]) );

      normalize_all(&psitot[0],grid.dv,ist.mn_states_tot,ist.nspinngrid,parallel.nthreads,ist.complex_idx,flag.printNorm);
      
      printf("\ndone calculating ortho, CPU time (sec) %g, wall run time (sec) %g\n",
                ((double)clock()-inital_clock_t)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-initial_wall_t); 
      fflush(stdout);

      if (1 == flag.saveCheckpoints){
        write_separation(stdout, top);
        printf("****    CHECKPOINT %d *** CHECKPOINT %d ** CHECKPOINT %d *** CHECKPOINT %d    ****", par.checkpoint_id, par.checkpoint_id, par.checkpoint_id, par.checkpoint_id);
        write_separation(stdout, bottom); 
        printf("\n"); fflush(stdout);

        save_job_state("checkpoint_2.dat",par.checkpoint_id,psitot,pot_local,ksqr,an,zn,ene_targets,nl,nlc,&grid,&ist,&par,&flag,&parallel);
      }
    // If checkpoint_id is 2, restart job right after orthogonalization, before diagonalization
    case 2:
      if (flag.restartFromCheckpoint == 2){
        par.checkpoint_id = flag.restartFromCheckpoint;
        write_separation(stdout, top);
        printf("****    CHECKPOINT %d *** CHECKPOINT %d ** CHECKPOINT %d *** CHECKPOINT %d    ****", par.checkpoint_id, par.checkpoint_id, par.checkpoint_id, par.checkpoint_id);
        write_separation(stdout, bottom); fflush(stdout);
        
        restart_from_save("checkpoint_2.dat",par.checkpoint_id,psitot,pot_local,ksqr,an,zn,ene_targets,nl,nlc,&grid,&ist,&par,&flag,&parallel);
      }
      par.checkpoint_id++;
      /***********************************************************************/
      /*** diagonalize the hamiltonian in the subspace spanned by the ***/
      /*** orthogonal filtered states, generating the eigenstates of the ***/
      /*** hamiltonian within the desired energy range ***/
      write_separation(stdout, top);
      printf("\n6. DIAGONALIZING HAMILTONIAN\n"); 
      write_separation(stdout, bottom); fflush(stdout);
      
      inital_clock_t = (double)clock(); initial_wall_t = (double)time(NULL);
      diag_H(psi,phi,psitot,pot_local,nlc,nl,ksqr,eig_vals,&ist,&par,&flag,planfw,planbw,fftwpsi);
      normalize_all(&psitot[0],grid.dv,ist.mn_states_tot,ist.nspinngrid,parallel.nthreads,ist.complex_idx,flag.printNorm);
      jms = ist.mn_states_tot;
      printf("\ndone calculating Hmat, CPU time (sec) %g, wall run time (sec) %g\n",
                ((double)clock()-inital_clock_t)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-initial_wall_t);
      fflush(stdout);

      if (1 == flag.saveCheckpoints){
        write_separation(stdout, top);
        printf("****    CHECKPOINT %d *** CHECKPOINT %d ** CHECKPOINT %d *** CHECKPOINT %d    ****", par.checkpoint_id, par.checkpoint_id, par.checkpoint_id, par.checkpoint_id);
        write_separation(stdout, bottom); fflush(stdout);
        printf("\n");

        save_job_state("checkpoint_3.dat",par.checkpoint_id,psitot,pot_local,ksqr,an,zn,ene_targets,nl,nlc,&grid,&ist,&par,&flag,&parallel);
      }
    // If checkpoint_id is 3, restart job right after diagonalization, before calculating output
    case 3:
      if (flag.restartFromCheckpoint == 3){
        par.checkpoint_id = flag.restartFromCheckpoint;
        write_separation(stdout, top);
        printf("****    CHECKPOINT %d *** CHECKPOINT %d ** CHECKPOINT %d *** CHECKPOINT %d    ****", par.checkpoint_id, par.checkpoint_id, par.checkpoint_id, par.checkpoint_id);
        write_separation(stdout, bottom); fflush(stdout);
        
        restart_from_save("checkpoint_3.dat",par.checkpoint_id,psitot,pot_local,ksqr,an,zn,ene_targets,nl,nlc,&grid,&ist,&par,&flag,&parallel);
      }
      par.checkpoint_id++;
      /***********************************************************************/
      /*** calculate the standard deviation of these states ***/
      /*** this is used to check if there are ghost states ***/
      write_separation(stdout, top);
      printf("\n7. CALCULATING VARIANCE OF EIGENVALUES\n");
      write_separation(stdout, bottom); fflush(stdout);
      
      calc_sigma_E(psi, phi, psitot, pot_local, nlc, nl, ksqr, sigma_E, &ist, &par, &flag, planfw, planbw, fftwpsi);

      /*** write the eigenstates/energies to a file ***/
      if (flag.getAllStates == 1) {printf("getAllStates flag on\nWriting all eigenstates to disk\n");}
      else{ printf("getAllStates flag off\nWriting eigenstates with sigE < 0.1 to disk\n");}
      
      if (flag.getAllStates == 1){
          /*** write all eigenvalues and their standard deviation to a file ***/
          if ((peig = fopen("eval.dat" , "w"))==NULL){
            printf("Out of disk space!\n Writing energy levels to stdout...\n\n\n\n");
            for (jms = 0; jms < ist.mn_states_tot; jms++) printf ("%ld %.16g %g\n", jms, eig_vals[jms], sigma_E[jms]); 
          }
          else{
            for (jms = 0; jms < ist.mn_states_tot; jms++) fprintf (peig,"%ld %.16g %g\n", jms, eig_vals[jms], sigma_E[jms]); 
            fclose(peig);
          }

          // Write all eigenstates to disk!
          if((ppsi = fopen("psi.dat" , "w"))==NULL){printf("Out of disk space!");}
          else {
            fwrite (&psitot[0],sizeof(psitot[0]),ist.complex_idx*ist.mn_states_tot*ist.nspinngrid,ppsi);
            fclose(ppsi);
          }
        }
        // Write output in the case that
        else{
          long eig = 0;
          // Write only the eigenstates having variance less than 0.1 to disk (remove ghost states)
          if (((peig = fopen("eval.dat" , "w"))==NULL) || ((ppsi = fopen("psi.dat" , "w"))==NULL)) {
            printf("Out of disk space!\n Writing energy levels to stdout...\n\n\n\n");
            for (jms = 0; jms < ist.mn_states_tot; jms++) printf ("%ld %.16g %g\n", jms, eig_vals[jms], sigma_E[jms]); 
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
                fprintf (peig,"%ld %.16g %g\n", jms, eig_vals[jms], sigma_E[jms]); 
                // Modify the eig_vals and sigma_E to carRy the values for the appropriate small variance eigenstate
                eig_vals[eig] = eig_vals[jms];
                sigma_E[eig] = sigma_E[jms];
                eig++;
              }
              else{ continue;}
            } 
          }
        // Reset the mn_states_tot value if we choose not to get all the states.
        ist.mn_states_tot = eig; 
        printf("Number of eigenstates with sigE < 0.1 saved to disk: new_mstot = %ld\n", ist.mn_states_tot);
      }

      // If no eigenstates were obtained, then the job has failed and alternate options should be pursued
      if (ist.mn_states_tot == 0){
        printf("The filter diagonalization yielded %ld eigenstates\n", ist.mn_states_tot);
        if ( (flag.retryFilter == 1) && (flag.alreadyTried == 0) ){
          printf("The filter diagonalization will be reattempted with increased filters and Chebyshev length\n");
          ist.n_filter_cycles += 16;
          ist.ncheby += 1024;
          printf("New n_filter_cycles = %ld, new ncheby = %ld\n", ist.n_filter_cycles, ist.ncheby);
          flag.alreadyTried = 1; // prevent the filtering algorithm from repeating indefinitely
        } else{
          fprintf(stderr, "The filter diagonalization yielded %ld eigenstates (T_T)\n", ist.mn_states_tot);
          fprintf(stderr, "Job failed :( exiting reluctantly\n");
          exit(EXIT_FAILURE);
        }
      }


      // The standard portion of the filter diagonalization procedure has concluded
      // Print the output.dat file for use in future workflows unless explicitly requested to NOT save output
      if (0 != flag->saveOutput){
        save_job_state("output.dat",par.checkpoint_id,psitot,pot_local,ksqr,an,zn,ene_targets,nl,nlc,&grid,&ist,&par,&flag,&parallel);
      }

      write_separation(stdout, top);
      printf("\nCALCULATING OPTIONAL OUTPUT\n"); 
      write_separation(stdout, bottom); fflush(stdout);

      long i, a, ieof, nval;
      char str[50];
      double evalloc, deloc;
      FILE *pf;

      ist.homo_idx = ist.lumo_idx = 0;
      pf = fopen("eval.dat" , "r");
      for (i = ieof = 0; ieof != EOF; i++){
        ieof = fscanf(pf, "%ld %lg %lg", &a, &evalloc, &deloc);
        if (deloc < par.sigma_E_cut && evalloc < par.fermi_E) ist.homo_idx = i;
      }
      fclose(pf);

      nval = i - 1;
      pf = fopen("eval.dat" , "r");
      for (i = 0; i <= ist.homo_idx; i++) fscanf(pf, "%ld %lg %lg", &a, &evalloc, &deloc);
      for (i = ist.homo_idx+1; i < nval; i++) {
        fscanf(pf, "%ld %lg %lg", &a, &evalloc, &deloc);
        if (deloc < par.sigma_E_cut) {
          ist.lumo_idx = i;
          break;
        }
      }
      fclose(pf);

      printf("index of homo, homo_idx = %ld; index of lumo, nlumo = %ld\n", ist.homo_idx, ist.lumo_idx); fflush(0);
      // Set the total number of electron and hole states in order to calculate the potential overlap integrals
      ist.total_homo = ist.homo_idx + 1; ist.total_lumo = ist.mn_states_tot - ist.total_homo;
      printf("total_homo = %ld total_lumo = %ld\n", ist.total_homo, ist.total_lumo); fflush(0);


      if (flag.printCubes == 1){
        /*** Write homo and lumo cube files ***/
        
        write_separation(stdout, top);
        printf("\nWRITING CUBE FILES\n"); 
        write_separation(stdout, bottom); fflush(stdout);

        if ((ist.homo_idx == 0) || (ist.lumo_idx == 0)){
          printf("\nDid not converge enough electron or hole states to visualize cube files.\n");
        } else{
        if ((rho = (double *) calloc(ist.ngrid, sizeof(double))) == NULL){
          fprintf(stderr, "\nOUT OF MEMORY: rho\n\n"); exit(EXIT_FAILURE);
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

        printf("\ndone calculating cubes, CPU time (sec) %g, wall run time (sec) %g\n",
          ((double)clock()-inital_clock_t)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-initial_wall_t);
        }
      }
      if (flag.calcPotOverlap == 1){

        write_separation(stdout, top);
        printf("\nCALCULATING POTENTIAL MATRIX ELEMENTS\n"); fflush(0);
        write_separation(stdout, bottom); fflush(stdout);

        calc_pot_overlap(&psitot[0], pot_local, nlc, nl, eig_vals, &par, &ist, &flag);
      } 
      if (flag.calcSpinAngStat == 1) {
        write_separation(stdout, top);
        printf("\nCALCULATING SPIN & ANG. MOM. STATISTICS\n"); 
        write_separation(stdout, bottom); fflush(stdout);

        calc_angular_exp(psitot, &grid,0, ist.mn_states_tot, &ist, &par, &flag, &parallel, planfw, planbw, fftwpsi);
      
      } 
      if ( (flag.printCubes != 1) && (flag.calcPotOverlap != 1) && (flag.calcSpinAngStat != 1) ) {
        printf("\nNo optional output requested.\n");
      }


      /*************************************************************************/
      /*** free memory ***/
      free(ist.atom_types);
      free(grid.x); free(grid.y); free(grid.z); 
      free(R); free(atom);
      free(psitot); free(psi); free(phi); 
      free(pot_local); free(ksqr);
      free(nlc); free(nl); 
      free(eig_vals); free(ene_targets); free(sigma_E);
      free(an); free(zn);      
      
      fftw_destroy_plan(planfw);
      fftw_destroy_plan(planbw);
      fftw_free(fftwpsi);
    
      time_t end_time = time(NULL);
      time_t end_clock = clock();

      write_separation(stdout, top);
      printf("\nDONE WITH PROGRAM: FILTER DIAGONALIZATION\n");
      printf("This calculation ended at: %s\n", ctime(&end_time)); 
      printf("Total job CPU time (sec) %.4g, wall run time (sec) %.4g",
                ((double)end_clock - (double)start_clock)/(double)(CLOCKS_PER_SEC), (double)end_time - (double)start_time );fflush(0);
      write_separation(stdout, bottom);
      
      free(top); free(bottom);
      exit(0);

  } // End of switch statement
} // End of main

/*****************************************************************************/
