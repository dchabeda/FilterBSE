#include "fd.h"

long rand_interval(long min, long max, long *seed);

int main(int argc, char *argv[]){
    // DECLARE VARIABLES AND STRUCTS
    // file pointers
    FILE *pf, *ppsi; 
    // zomplex types
    zomplex *psi, *phi; 
    // FFT
    fftw_plan_loc planfw, planbw; fftw_complex *fftwpsi; 
    long fft_flags=0;
    // custom structs 
    flag_st flag; index_st ist; par_st par; 
    atom_info *atom_equil; atom_info *atom;
    pot_st pot; 
    grid_st grid_par; 
    xyz_st *grid;
    xyz_st *R; xyz_st *R_equil;
    nlc_st *nlc = NULL; nlc_st *nlc_equil = NULL; 
    parallel_st parallel; 
    // double arrays
    double *ksqr, *pot_local_equil, *pot_local;
    double *SO_projectors_equil; double *SO_projectors; 
    // long int arrays and counters
    long *nl_equil = NULL; long *nl = NULL;
    long i, jms, j, rand_seed;
    int start, end;
    ist.atom_types = malloc(N_MAX_ATOM_TYPES*sizeof(ist.atom_types[0]));
    // Clock/Wall time output and stdout formatting
    time_t start_time = time(NULL); // Get the actual time for total wall runtime
    // time_t start_clock = clock(); // Get the starting CPU clock time for total CPU runtime
    char *top; top = malloc(2*sizeof(top[0])); 
    char *bottom; bottom = malloc(2*sizeof(bottom[0]));
    strcpy(top, "T\0"); strcpy(bottom, "B\0");
    
    fprintf(stdout, "******************************************************************************\n");
    printf("\nRUNNING PROGRAM: POT MATRIX ELEMS\n");
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
    read_input(&flag, &grid_par, &ist, &par, &parallel);

    /*** allocating memory ***/
    // the positions of the atoms in the x, y, and z directions 
    if ((R_equil = (xyz_st *) calloc(ist.natoms, sizeof(xyz_st))) == NULL) {
        fprintf(stderr, "\nOUT OF MEMORY: R array\n\n"); exit(EXIT_FAILURE);
    }
    // the atom specific information 
    if ((atom_equil = (atom_info *) calloc(ist.natoms, sizeof(atom_info))) == NULL){
        fprintf(stderr, "\nOUT OF MEMORY: atom struct\n\n"); exit(EXIT_FAILURE);
    }
    
    /*** read the nanocrystal configuration ***/
    char str[100];
    
    char *file_name_equil; file_name_equil = malloc(16*sizeof(file_name_equil[0]));
    strcpy(file_name_equil, "conf_equil.par");
    sprintf(str, "\nReading atomic configuration from %s:\n", file_name_equil);
    printf("%s", str);
    read_conf(file_name_equil, R_equil, atom_equil, &ist, &par, &flag);

    /*** initialize parameters for the grid ***/
    printf("\nInitializing the grid parameters:\n");
    init_grid_params(&grid_par, R_equil, &ist, &par, &flag);

    // Allocate memory for the grid in the x, y, and z directions ***/
    // Allocate memory for the grid in the x, y, and z directions ***/
    if ((grid = (xyz_st *) calloc(grid_par.ngrid, sizeof(xyz_st))) == NULL){
        fprintf(stderr, "\nOUT OF MEMORY: grid\n\n"); exit(EXIT_FAILURE);
    }
    // the kinetic energy stored on the grid
    if ((ksqr = (double *) calloc(ist.ngrid, sizeof(double))) == NULL){
        fprintf(stderr, "\nOUT OF MEMORY: ksqr\n\n"); exit(EXIT_FAILURE);
    }

    /*** build the real- and k-space grids ***/
    printf("\nBuilding the real-space and k-space grids:\n");
    build_grid_ksqr(ksqr, R_equil, grid, &grid_par, &ist, &par, &flag);
    
    /*************************************************************************/
    /*** allocating memory for the rest of the program ***/
    printf("\nAllocating memory for FFT, pot, psi, eig_vals...");
    
    // FFT
    fftwpsi = fftw_malloc(sizeof(fftw_complex) * ist.ngrid);
    /*** initialization for the fast Fourier transform ***/
    planfw = fftw_plan_dft_3d(grid_par.nz, grid_par.ny, grid_par.nx, fftwpsi, fftwpsi, FFTW_FORWARD, fft_flags);
    planbw = fftw_plan_dft_3d(grid_par.nz, grid_par.ny, grid_par.nx, fftwpsi, fftwpsi, FFTW_BACKWARD, fft_flags);
    
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
    if ((pot_local_equil = (double *) calloc(ist.ngrid, sizeof(double))) == NULL){
        fprintf(stderr, "\nOUT OF MEMORY: pot_local_equil\n\n"); exit(EXIT_FAILURE);
    }
    
    // memory allocation for the spin-orbit potential 
    if ( (flag.SO == 1) || (flag.NL == 1) ){
        if ((SO_projectors_equil = (double*) calloc(PROJ_LEN * ist.nproj, sizeof(double)))==NULL){nerror("mem_SO_projector");}  
        if ((SO_projectors = (double*) calloc(PROJ_LEN * ist.nproj, sizeof(double)))==NULL){nerror("mem_SO_projector");}  
    }
    // memory allocation for the non-local potential
    if (1 == flag.NL){
        if ((nlc_equil = (nlc_st *) calloc(ist.n_NL_atoms*ist.n_NL_gridpts, sizeof(nlc_st))) == NULL){ 
        fprintf(stderr, "\nOUT OF MEMORY: nlc\n\n"); exit(EXIT_FAILURE);
        }
        if ((nl_equil = (long *) calloc(ist.natoms, sizeof(nl[0]))) == NULL) {
        fprintf(stderr, "\nOUT OF MEMORY: nl\n\n"); exit(EXIT_FAILURE);
        }

        if ((nlc = (nlc_st *) calloc(ist.n_NL_atoms*ist.n_NL_gridpts, sizeof(nlc_st))) == NULL){ 
        fprintf(stderr, "\nOUT OF MEMORY: nlc\n\n"); exit(EXIT_FAILURE);
        }
        if ((nl = (long *) calloc(ist.natoms, sizeof(nl[0]))) == NULL) {
        fprintf(stderr, "\nOUT OF MEMORY: nl\n\n"); exit(EXIT_FAILURE);
        }
    }
    
    printf("\tdone allocating memory.\n"); fflush(stdout);

    /**************************************************************************/
    printf("\nInitializing potentials...\n");
    
    printf("\nLocal pseudopotential:\n");
    build_local_pot(pot_local_equil, &pot, R_equil, ksqr, atom_equil, grid, &grid_par, &ist, &par, &flag, &parallel);
    
    free(pot.r); pot.r = NULL; 
    free(pot.pseudo); pot.pseudo = NULL; 
    free(pot.dr); pot.dr = NULL; 
    free(pot.file_lens); pot.file_lens = NULL;
    if (1.0 != par.scale_surface_Cs){
    free(pot.r_LR); pot.r_LR = NULL;
    free(pot.pseudo_LR); pot.pseudo_LR = NULL; 
    }
    
    write_cube_file(pot_local_equil, &grid_par, "local_pot_equil.cube");
    
    if(flag.SO==1) {
    printf("\nSpin-orbit pseudopotential:\n");
    init_SO_projectors(SO_projectors, R_equil, atom_equil, &grid_par, &ist, &par);
    }
    /*** initialization for the non-local potential ***/
    if (flag.NL == 1){
    printf("\nNon-local pseudopotential:\n"); fflush(0);
    init_NL_projectors(nlc_equil, nl_equil, SO_projectors, R_equil, atom_equil, grid, &grid_par, &ist, &par, &flag);
    }
    // free memory allocated to SO_projectors
    if ( (flag.SO == 1) || (flag.NL == 1) ){
    free(SO_projectors); SO_projectors = NULL;
    }

    // ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 
    // ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
    // Potential mat elems
    // This code computes dU = U(r; R) - U(r; R_equil)
    // the difference between the potential energy surfaces of a Hamiltonian
    // that parametrically depends on the nuclear coordinates vs one that
    // only depends on the equilibrium nuclear configuration.
    // We calculate the variation in the energies as a perturbation expansion
    // E_0(R) = E_0(R_equil) + <phi_0|dU|phi_0> + sum_{i neq 0} |<phi_0|dU|phi_i>|^2/dE + ...
    // We verify that the first order term is sufficient to describe the changes
    // by computing its magnitude in comparison to the energy differences between
    // states in the electronic manifold at the equilibrium geometry.
    // ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 
    // ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    // 1. Read in the equilibrium and distorted nuclear configurations
    /*** allocating memory ***/
    // the distorted positions of the atoms in the x, y, and z directions 
    if ((R = (xyz_st *) calloc(ist.natoms, sizeof(xyz_st))) == NULL) {
        fprintf(stderr, "\nOUT OF MEMORY: R array\n\n"); exit(EXIT_FAILURE);
    }
    // the atom specific information 
    if ((atom = (atom_info *) calloc(ist.natoms, sizeof(atom_info))) == NULL){
        fprintf(stderr, "\nOUT OF MEMORY: atom struct\n\n"); exit(EXIT_FAILURE);
    }
    
    /*** read the nanocrystal configuration ***/
    printf("\nReading atomic configuration from conf.par:\n");
    char *file_name; file_name = malloc(9*sizeof(file_name[0]));
    strcpy(file_name, "conf.par");
    read_conf(file_name, R, atom, &ist, &par, &flag);

    // 2. Calculate their local potentials
    if ((pot_local = (double *) calloc(ist.ngrid, sizeof(double))) == NULL){
        fprintf(stderr, "\nOUT OF MEMORY: pot_local_equil\n\n"); exit(EXIT_FAILURE);
    }
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
    
    build_local_pot(pot_local, &pot, R, ksqr, atom, grid, &grid_par, &ist, &par, &flag, &parallel);
    
    write_cube_file(pot_local, &grid_par, "local_pot.cube");
    
    free(pot.r); pot.r = NULL; 
    free(pot.pseudo); pot.pseudo = NULL; 
    free(pot.dr); pot.dr = NULL; 
    free(pot.file_lens); pot.file_lens = NULL;
    if (1.0 != par.scale_surface_Cs){
    free(pot.r_LR); pot.r_LR = NULL;
    free(pot.pseudo_LR); pot.pseudo_LR = NULL; 
    }

    // // 3. Read in the equilibrium wavefunction.
    // //count number of states found
    // jms = countlines("eval.dat");
    // printf("%ld total states in psi.dat\n", jms); fflush(0);
  
    // //allocate memory for psi
    // if ((psi = (zomplex *) calloc(ist.nspinngrid, sizeof(zomplex))) == NULL) nerror("psi");

    // //read psi from file
	// ppsi = fopen("psi.dat" , "r");

	// char filename[20];
    // for (j = start; j <= end; j++){ 
    //     printf("Reading state %d from psi.dat\n", j);

    //     if(fseek(ppsi,j*ist.nspinngrid*sizeof(zomplex),SEEK_SET)!=0){
    //     printf("Error reading from psi.dat!\n"); exit(EXIT_FAILURE);
    //     }
    //     fread (&psi[0],sizeof(zomplex),ist.nspinngrid,ppsi);


    // // 4. Calculate matrix elements of the equilibrium wavefunction with U(r;R_equil)
    
    // // 5. Calculate matrix elements of the equilibrium wavefunction with U(r;R)
    
    // }
    
    return 0;
}

/*****************************************************************************/
int countlines(char *filename){
  FILE* fp = fopen(filename,"r");
  int lines = 0;
  int ch;
  while(1){
    ch = fgetc(fp);
    if (feof(fp)){ break; }
    if(ch == '\n')
    {
      lines++;
    }
  }
  fclose(fp);
  return lines;
}


