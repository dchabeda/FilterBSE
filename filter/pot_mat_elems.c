#include "fd.h"

long rand_interval(long min, long max, long *seed);

int main(int argc, char *argv[]){
    // DECLARE VARIABLES AND STRUCTS
    // file pointers
    FILE *pf, *pseed, *pr, *ppsi; 
    // zomplex types
    zomplex *psi, *phi; 
    // FFT
    fftw_plan_loc planfw, planbw; fftw_complex *fftwpsi; 
    long fft_flags=0;
    // custom structs 
    flag_st flag; index_st ist; par_st par; atom_info *atom; 
    pot_st pot; grid_st grid; xyz_st *R; nlc_st *nlc = NULL; 
    parallel_st parallel; 
    // double arrays
    double *ksqr, *pot_local;
    double *SO_projectors; 
    // long int arrays and counters
    long *nl = NULL;
    long i, rand_seed;
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
    build_grid_ksqr(ksqr, R, &grid, &ist, &par);
    
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

    // ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 
    // ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
    // Potential mat elems
    // ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 
    // ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    // 1. Initialize list of grid indices within Rnlcut of the atom.
    long g_idx, rand_gridpt;
    long r, r_p;
    int n_r_pts = 50;
    long select_atom = atol(argv[1]);
    zomplex up_val, dn_val;
    printf("The coupling will be calculated on grid points around atom %s %ld\n", atom->atyp[select_atom], select_atom);
    
    
    pseed = fopen("seed.dat", "w");
    if (flag.setSeed == 1){
    rand_seed = - par.rand_seed;} 
    else {
    Randomize();  rand_seed = -random();
    }
    
    printf("\nGrid points of coupling elements <r'|V|r>\n");
    pr = fopen("r_grpts.dat", "w");
    pf = fopen("pot_mat_elems.dat", "w");

    
    for (i = 0; i < n_r_pts; i++){
        fprintf(pseed, "%ld\n", rand_seed);
        rand_gridpt = rand_interval(0, ist.n_NL_gridpts, &rand_seed); rand_seed = - random();
        g_idx = select_atom * ist.n_NL_gridpts + rand_gridpt; 

        r_p = nlc[g_idx].jxyz;
        fprintf(pseed, "%ld\n", rand_seed);

        rand_gridpt = rand_interval(0, ist.n_NL_gridpts, &rand_seed); rand_seed = - random();
        g_idx = select_atom * ist.n_NL_gridpts + rand_gridpt; 

        r = nlc[g_idx].jxyz;
        
        // Initialize a wavefunction with only 1/sqrt(2) at gript r
        psi[r].re = 1/sqrt(2); psi[r].im = 1/sqrt(2);
        psi[r + ist.ngrid].re = 1/sqrt(2); psi[r + ist.ngrid].im = 1/sqrt(2);

        // Apply the potential to this wavefunction to obtain phi = V|psi>
        potential(phi, psi, pot_local, nlc, nl, &ist, &par, &flag);

        // Calculate the inner product <r_p|V|r>
        up_val.re = 1/sqrt(2) * phi[r_p].re + 1/sqrt(2) * phi[r_p].im;
        up_val.im = -1/sqrt(2) * phi[r_p].re + 1/sqrt(2) * phi[r_p].im;
        dn_val.re = 1/sqrt(2) * phi[r_p + ist.ngrid].re + 1/sqrt(2) * phi[r_p + ist.ngrid].im;
        dn_val.im = -1/sqrt(2) * phi[r_p + ist.ngrid].re + 1/sqrt(2) * phi[r_p + ist.ngrid].im;

        fprintf(pf, "%ld %ld %lg %lg %lg %lg\n", r_p, r, up_val.re, up_val.im, dn_val.re, dn_val.im);
    }
    
    fclose(pseed);
    fclose(pf);
    fclose(pr);
    
    return 0;
}

long rand_interval(long min, long max, long *seed) { // min and max included 
  return floor(ran_nrc(seed) * (max - min + 1) + min);
}
