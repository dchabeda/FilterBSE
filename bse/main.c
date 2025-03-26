/*****************************************************************************/

#include "fd.h"
#include <mpi.h>

/*****************************************************************************/

int main(int argc, char *argv[]){
    /*****************************************************************
    * This is the main function for BSE.x. It is the driver that     *
    * controls memory allocation, structs, and the BSE algorithm.    *
    * The algorithm computes correlated excitonic states using the   *
    * Bethe-Salpether formalism described by Rohlfing and Louie:     *
    * Phys. Rev. B 62, 4927                                          *
    * It is applied to nanocrystal systems by using an electron-hole *
    * basis constructed from a semi-empirical pseudopot. hamiltonian *
    ******************************************************************/ 

    // DECLARE VARIABLES AND STRUCTS
    // file pointers
    FILE *pf, *pmem, *ppsi;
    pmem = fopen("mem.dat", "w");
    // zomplex types
    zomplex *pot_bare, *pot_screened, *pot_hartree; 
    zomplex *bsmat, *direct, *exchange, *bs_coeff;
    zomplex *LdotS, *l2_mom;
    // custom structs
    par_st par; index_st ist; grid_st grid; flag_st flag; xyz_st *R = NULL; parallel_st parallel; 
    xyz_st *s_mom, *l_mom;
    xyz_st *trans_dipole, *mag_dipole, *rot_strength;
    xyz_st *S_mom = NULL, *L_mom = NULL, *L2_mom = NULL;
    // FFT 
    fftw_plan_loc *planfw, *planbw; fftw_complex *fftwpsi;
    long fft_flags = 0;
    // double arrays
    double *psitot = NULL, *psi_hole = NULL, *psi_elec = NULL, *psi_qp;
    double *eig_vals = NULL, *sigma_E = NULL, *xton_ene;
    double *h0mat;
    double *gridx = NULL, *gridy = NULL, *gridz = NULL;
    double *rho;
    // long int arrays and counters
    long i, j, state;
    long jgrid, jgrid_real, jgrid_imag;
    long mem = 0;
    ist.atom_types = malloc(N_MAX_ATOM_TYPES*sizeof(ist.atom_types[0]));
    fprintf(pmem, "alloc ist.atom_types %ld B\n", N_MAX_ATOM_TYPES*sizeof(ist.atom_types[0])); mem += N_MAX_ATOM_TYPES*sizeof(ist.atom_types[0]);
    // Clock/Wall time output and stdout formatting
    time_t start_time = time(NULL); // Get the actual time for total wall runtime
    time_t start_clock = clock(); // Get the starting CPU clock time for total CPU runtime
    time_t current_time;
    char* c_time_string;
    char *top; top = malloc(2*sizeof(top[0])); 
    char *endptr;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &parallel.mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &parallel.mpi_size);
    parallel.mpi_root = 0;
    const int mpir = parallel.mpi_rank;
    //

    fprintf(pmem, "alloc top %ld B\n", 2*sizeof(top[0])); mem += 2*sizeof(top[0]);
    char *bottom; bottom = malloc(2*sizeof(bottom[0]));
    fprintf(pmem, "alloc bottom %ld B\n", 2*sizeof(bottom[0])); mem += 2*sizeof(bottom[0]);
    strcpy(top, "T\0"); strcpy(bottom, "B\0");

    /*************************************************************************/
    if (mpir == 0) fprintf(stdout, "******************************************************************************\n");
    if (mpir == 0) printf("\nRUNNING PROGRAM: BETHE-SALPETHER\n");
    if (mpir == 0) printf("This calculation began at: %s", ctime(&start_time)); 
    if (mpir == 0) write_separation(stdout, bottom);
    fflush(stdout);

    /*************************************************************************/
    // 1. Initialize job from input file
    
    if (mpir == 0) write_separation(stdout, top);
    
    if (mpir == 0) printf("\n1.\tINITIALIZING JOB | %s\n", get_time());
    if (mpir == 0) write_separation(stdout, bottom); fflush(stdout);

    /*************************************************************************/
    /*** Read output from filter output.par ***/
    int init_unsafe = 0;
    if (argv[1] != NULL) init_unsafe = (int) strtol(argv[1], &endptr, 10);
    
    if (0 == init_unsafe){
        if (mpir == 0) printf("\nReading filter output from output.dat:\n"); fflush(stdout);
    
    // ******
    // ******
        read_filter_output("output.dat", &psitot, &eig_vals, &sigma_E, &R, &grid, &(grid.x), &(grid.y), &(grid.z), &ist, &par, &flag);
    // ******
    // ******
    } else if (1 == init_unsafe){
        if (mpir == 0) printf("\nReading UNSAFE input from unsafe_input.par:\n"); fflush(stdout);
        read_unsafe_input(&psitot, &eig_vals, &sigma_E, &R, &grid, &(grid.x), &(grid.y), &(grid.z), &ist, &par, &flag, &parallel);
    }

    /*************************************************************************/
    /*** Read initial setup from input.par ***/
    if (mpir == 0) printf("\nReading BSE job specifications from input.par:\n"); fflush(stdout);
    // ******
    // ******
    read_input(&flag, &grid, &ist, &par, &parallel);
    // ******
    // ******
    fflush(stdout);

    /*************************************************************************/
    /*** Determine the configuration of the quasiparticle basis ***/
    if (mpir == 0) printf("\nSetting quasiparticle basis indices:\n"); fflush(stdout);
    
    // Allocate memory for the lists of the indices of eigenstates
    // The maximum possible number of hole states is mn_states_tot from filter.
    // We allocate an entire block of that size for both elecs and holes 
    // because we will reallocate after the get_qp_basis_indices function.
    if((ist.eval_hole_idxs = (long*) malloc(ist.mn_states_tot * sizeof(ist.eval_hole_idxs[0]))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for ist.eval_hole_idxs in main.c\n");
        exit(EXIT_FAILURE);
    }
    if((ist.eval_elec_idxs = (long*) malloc(ist.mn_states_tot * sizeof(ist.eval_elec_idxs[0]))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for ist.eval_elec_idxs in main.c\n");
        exit(EXIT_FAILURE);
    }
    
    // ******
    // ******
    get_qp_basis_indices(eig_vals, sigma_E, &ist.eval_hole_idxs, &ist.eval_elec_idxs, &ist, &par, &flag);
    // ******
    // ******
    
    // Reallocate the eig_vals and sigma_E arrays to only contain the n_qp states
    // eig_vals = realloc(eig_vals, ist.n_qp * sizeof(eig_vals[0]));
    // sigma_E = realloc(sigma_E, ist.n_qp * sizeof(sigma_E[0]));
    
    // Resize the index arrays to tightly contain the indices of eigenstates in VB/CB
    // The conditional checks whether n_holes saturated the memory block (if so, no need to resize)
    // or if they are not equal, then n_holes is smaller and we should shrink the array
    // if (ist.n_holes != ist.max_hole_states){
    //     ist.eval_hole_idxs = realloc(ist.eval_hole_idxs, ist.n_holes * sizeof(long));
    // }
    // if (ist.n_elecs != ist.max_elec_states){
    //     ist.eval_elec_idxs = realloc(ist.eval_elec_idxs, ist.n_elecs * sizeof(long));
    // }
    // fprintf(pmem, "alloc ist.eval_hole_idxs %ld B\n", ist.mn_states_tot * sizeof(ist.eval_hole_idxs[0])); mem += ist.mn_states_tot * sizeof(ist.eval_hole_idxs[0]);
    // fprintf(pmem, "alloc ist.eval_elec_idxs %ld B\n", ist.mn_states_tot * sizeof(ist.eval_elec_idxs[0])); mem += ist.mn_states_tot * sizeof(ist.eval_elec_idxs[0]);
    

    /*************************************************************************/
    // Allocate memory for the quasiparticle wavefunctions
    
    if ((psi_qp = (double *) malloc( ist.complex_idx * ist.nspinngrid * (ist.n_qp) * sizeof(double))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for psi_qp in main.c\n");
        exit(EXIT_FAILURE);
    }
    
    fprintf(pmem, "alloc psi_qp %ld B", ist.complex_idx * ist.nspinngrid * (ist.n_holes + ist.n_elecs) * sizeof(double)); mem += ist.complex_idx * ist.nspinngrid * (ist.n_holes + ist.n_elecs) * sizeof(double);
    if (mpir == 0) write_separation(pmem, top);
    fprintf(pmem, "\ntotal mem usage %ld MB\n", mem / 1000000 );
    if (mpir == 0) write_separation(pmem, bottom); fflush(pmem);

    // if (mpir == 0){
    //     printf("\nTHE HOLES:\n");
    //     for (i = 0; i < ist.n_holes; i++){
    //         printf("%ld\n", ist.eval_hole_idxs[i]);
    //     }
        
    //     printf("\nTHE ELECS:\n");
    //     for (i = 0; i < ist.n_elecs; i++){
    //         printf("%ld\n", ist.eval_elec_idxs[i]);
    //     }
    // }
    
    // ******
    // ******
    if (mpir == 0) printf("\nReading quasiparticle basis wavefunctions:\n"); fflush(stdout);
    get_qp_basis(psi_qp, psitot, eig_vals, sigma_E, &ist, &par, &flag);
    ist.n_xton = ist.n_elecs * ist.n_holes;
    // ******
    // ******
    free(psitot); free(sigma_E); //free(psi_hole); free(psi_elec);
    fprintf(pmem, "free psitot %ld B\n", ist.complex_idx * ist.nspinngrid * ist.mn_states_tot * sizeof(double)); mem -= ist.complex_idx * ist.nspinngrid * ist.mn_states_tot * sizeof(double);
    // fprintf(pmem, "free psi_hole %ld B\n", ist.complex_idx * ist.nspinngrid * ist.n_holes * sizeof(psi_hole[0])); mem -= ist.complex_idx*ist.nspinngrid*ist.n_holes*sizeof(psi_hole[0]);
    // fprintf(pmem, "free psi_elec %ld B\n", ist.complex_idx * ist.nspinngrid * ist.n_elecs * sizeof(psi_elec[0])); mem -= ist.complex_idx*ist.nspinngrid*ist.n_elecs*sizeof(psi_elec[0]);
    fprintf(pmem, "free sigma_E %ld B", ist.mn_states_tot * sizeof(double)); mem -= ist.mn_states_tot * sizeof(double);
    if (mpir == 0) write_separation(pmem, top);
    fprintf(pmem, "\ntotal mem usage %ld MB\n", mem / 1000000 );
    if (mpir == 0) write_separation(pmem, bottom); fflush(pmem);

    
    // char fileName[50];
    // rho = (double*) calloc(ist.ngrid, sizeof(double));

    if (1 == flag.useSpinors){
        if (mpir == 0) printf("\nComputing spin composition of quasiparticle spinors | %s\n", get_time()); fflush(stdout);
        FILE *pf = fopen("qp_spins.dat", "w");
        for (int state  = 0; state < ist.n_qp; state++){
            fprintf(pf, "\nStats on state%d: (E=%lg)\n", state, eig_vals[state]);
            double perUp = 0;
            double perDn = 0;
            // sprintf(fileName, "qp_%d_up.cube", state);
            for (long jgrid = 0; jgrid < ist.ngrid; jgrid++){
                jgrid_real = ist.complex_idx * jgrid;
                jgrid_imag = ist.complex_idx * jgrid + 1;
                
                perUp += sqr(psi_qp[state*ist.nspinngrid*ist.complex_idx+jgrid_real])+sqr(psi_qp[state*ist.nspinngrid*ist.complex_idx+jgrid_imag]);
                perDn += sqr(psi_qp[state*ist.nspinngrid*ist.complex_idx+ist.ngrid*ist.complex_idx+jgrid_real])+sqr(psi_qp[state*ist.nspinngrid*ist.complex_idx+ist.ngrid*ist.complex_idx+jgrid_imag]);

                // rho[jgrid] = sqrt(sqr(psi_qp[state*ist.nspinngrid*ist.complex_idx+jgrid_real])+sqr(psi_qp[state*ist.nspinngrid*ist.complex_idx+jgrid_imag]));
            }
            // write_cube_file(rho, &grid, fileName);
            fprintf(pf, " Spin up fraction: %f\n", perUp * grid.dv); fflush(pf);
            fprintf(pf, " Spin dn fraction: %f\n", perDn * grid.dv); fflush(pf);
        }
        fclose(pf);
    }
    // free(rho);
    /*************************************************************************/
    // 2. Compute electron-hole interaction potentials
    
    if (mpir == 0) write_separation(stdout, top);
    
    if (mpir == 0) printf("\n2.\tCOMPUTING ELECTRON-HOLE INTERACTION POTENTIALS | %s\n", get_time());
    if (mpir == 0) write_separation(stdout, bottom); fflush(stdout);

    /*************************************************************************/
    // Allocate memory for the hartree, screened coulomb and bare exchange pots
    if (mpir == 0) printf("Allocating memory for hartree, direct, and exchange potentials..."); fflush(stdout);
    pot_bare  = (zomplex *) malloc(ist.ngrid * sizeof(pot_bare[0]));
    pot_screened  = (zomplex *) malloc(ist.ngrid * sizeof(pot_screened[0]));
    pot_hartree = (zomplex *) malloc(ist.ngrid*parallel.nthreads*sizeof(pot_hartree[0]));
    fprintf(pmem, "alloc pot_bare %ld B\n", ist.ngrid*sizeof(pot_bare[0])); mem += ist.ngrid * sizeof(pot_bare[0]);
    fprintf(pmem, "alloc pot_screened %ld B\n", ist.ngrid*sizeof(pot_screened[0])); mem += ist.ngrid * sizeof(pot_screened[0]);
    fprintf(pmem, "alloc pot_hartree %ld B\n", ist.ngrid*parallel.nthreads*sizeof(pot_bare[0])); mem += ist.ngrid*parallel.nthreads*sizeof(pot_bare[0]);
    if (mpir == 0) printf(" done.\n"); fflush(stdout);

    // Initialize the FFT arrays for parallel Fourier transform
    fftwpsi = fftw_malloc(ist.ngrid*parallel.nthreads*sizeof(fftw_complex));
    fprintf(pmem, "alloc fftwpsi %ld B\n", ist.ngrid*parallel.nthreads*sizeof(fftw_complex)); mem += ist.ngrid*parallel.nthreads*sizeof(fftw_complex);
    
    // Initialize the parallel FFT
    // fftw_plan_with_nthreads(ist.nthreads);
  
    planfw = (fftw_plan_loc *) malloc(parallel.nthreads * sizeof(fftw_plan_loc));
    planbw = (fftw_plan_loc *) malloc(parallel.nthreads * sizeof(fftw_plan_loc));
    for (i = 0; i < ist.nthreads; i++) { 
        planfw[i] = fftw_plan_dft_3d(grid.nz, grid.ny, grid.nx, &fftwpsi[i*ist.ngrid], 
                                     &fftwpsi[i*ist.ngrid], FFTW_FORWARD, fft_flags);
        planbw[i] = fftw_plan_dft_3d(grid.nz, grid.ny, grid.nx, &fftwpsi[i*ist.ngrid],
                                     &fftwpsi[i*ist.ngrid], FFTW_BACKWARD, fft_flags);
    }
    
    if (mpir == 0) printf("Computing interaction potential on grid...\n"); fflush(stdout);
    init_elec_hole_kernel(pot_bare, pot_screened, &grid, &ist, &par, planfw[0], planbw[0], &fftwpsi[0]);

    /*************************************************************************/
    /*************************************************************************/
    // 2. Compute single particle properties
    
    if (mpir == 0) write_separation(stdout, top);
    
    if (mpir == 0) printf("\n3.\tCOMPUTING SINGLE-PARTICLE PROPERTIES | %s\n", get_time());
    if (mpir == 0) write_separation(stdout, bottom); fflush(stdout);

    /*************************************************************************/
    if (mpir == 0) printf("Allocating memory for single particle matrix elements... "); fflush(stdout);
    // trans_dipole def= <psi_i|mu|psi_a>
    if ((trans_dipole = (xyz_st *) malloc(ist.n_elecs*ist.n_holes * sizeof(trans_dipole[0]))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for trans_dipole in main.c\n");
        exit(EXIT_FAILURE);
    }
    fprintf(pmem, "alloc trans_dipole %ld B\n", ist.n_elecs*ist.n_holes * sizeof(trans_dipole[0])); mem += ist.n_elecs*ist.n_holes * sizeof(trans_dipole[0]);
    // mag_dipole def= <psi_i|-1/2L|psi_a>
    if ((mag_dipole = (xyz_st *) malloc(ist.n_elecs*ist.n_holes * sizeof(mag_dipole[0]))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for mag_dipole in main.c\n");
        exit(EXIT_FAILURE);
    }
    fprintf(pmem, "alloc mag_dipole %ld B\n", ist.n_elecs*ist.n_holes * sizeof(mag_dipole[0])); mem += ist.n_elecs*ist.n_holes * sizeof(mag_dipole[0]);
    // rs def= <psi_i|mu|psi_a><a|m|i>
    if ((rot_strength = (xyz_st *) malloc(ist.n_elecs*ist.n_holes * sizeof(rot_strength[0]))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for rot_strength in main.c\n");
        exit(EXIT_FAILURE);
    }
    fprintf(pmem, "alloc rot_strength %ld B\n", ist.n_elecs*ist.n_holes * sizeof(rot_strength[0])); mem += ist.n_elecs*ist.n_holes * sizeof(rot_strength[0]);
    
    if (mpir == 0) write_separation(pmem, top);
    fprintf(pmem, "\ntotal mem usage %ld MB\n", mem / 1000000 );
    if (mpir == 0) write_separation(pmem, bottom); fflush(pmem);
    if (mpir == 0) printf("done\n"); fflush(stdout);

    if (mpir == 0) printf("\nElectric transition dipole moment...\n");
    // calc_elec_dipole(trans_dipole, psi_qp, eig_vals, &grid, &ist, &par, &flag);
    if (mpir == 0) printf("\nMagnetic transition dipole moment...\n");
    // calc_mag_dipole(mag_dipole, psi_qp, eig_vals, &grid, &ist, &par, &flag);
    // calc_rot_strength(rs, mux, muy, muz, mx, my, mz, eig_vals, &ist);
    
    /* */
    /* */
    /* */
    if (1 == flag.SO){
        if (mpir == 0) write_separation(stdout, top);
        current_time = time(NULL);
        c_time_string = ctime(&current_time);
        if (mpir == 0) printf("\n  -  COMPUTING ANG.MOM. PROPERTIES | %s\n", get_time());
        if (mpir == 0) write_separation(stdout, bottom); fflush(stdout);

        l_mom = (xyz_st *) malloc( (ist.n_holes*ist.n_holes + ist.n_elecs*ist.n_elecs) * sizeof(L_mom[0])); //<psi_r|Lx|psi_s>
        l2_mom =(zomplex *) malloc((ist.n_holes*ist.n_holes + ist.n_elecs*ist.n_elecs) * sizeof(L_mom[0])); //<psi_r|Lx|psi_s>
        s_mom = (xyz_st *) malloc( (ist.n_holes*ist.n_holes + ist.n_elecs*ist.n_elecs) * sizeof(s_mom[0])); //<psi_r|Sx|psi_s>
        LdotS = (zomplex *) malloc( (ist.n_holes*ist.n_holes+ist.n_elecs*ist.n_elecs) * sizeof(LdotS[0])); //<psi_r|L.S|psi_s>

        // Compute spin matrix elements, e.g. <j|Sx|i>
        // calc_spin_mtrx(s_mom, psi_qp, &grid, &ist, &par);
        // Compute angular momentum matrix elements, e.g. <j|Lx|i>
        // calc_ang_mom_mtrx(l_mom, l2_mom, LdotS, psi_qp, &grid, &ist, &par);
    }
    /* */
    /* */
    /* */
    
    // /**************************************************************************/
    // /*** this routine computes the coulomb coupling between
    //      single excitons.  On input - it requires the eigenstates stored in psi_qp,
    //      the eigenvalues stored in eval, and pot_hartree computed in init_elec_hole_kernel.
    //      On output it stores the coulomb matrix elements on the disk
    //      in the following format: a, i, b, j, ene_ai, ene_bj, vjbai, vabji.
    //      a - the index of the electron in exciton Sai.
    //      i - the index of the hole in exciton Sai.
    //      b - the index of the electron in exciton Sbj.
    //      j - the index of the hole in exciton Sbj.
    //      ene_ai - the energy of exciton Sai.
    //      ene_bj - the energy of exciton Sbj.
    //      vjbai and vabji are the coulomb matrix elements need to be used to
    //      generate the spin-depedent matrix elements as described by
    //      the last equation in our codument.  ***/
    if (mpir == 0) write_separation(stdout, top);
    if (mpir == 0) printf("\n4.\tCOMPUTING ELEC-HOLE INTERACTION KERNEL | %s\n", get_time());
    if (mpir == 0) write_separation(stdout, bottom); fflush(stdout);
    
    if (mpir == 0) printf("\nThe number of electron-hole pairs in the exciton basis = %ld\n", ist.n_xton);

    direct = (zomplex *) calloc(ist.n_xton * ist.n_xton, sizeof(zomplex)); 
    exchange = (zomplex *) calloc(ist.n_xton * ist.n_xton, sizeof(zomplex)); 
    bsmat = (zomplex *) calloc(ist.n_xton * ist.n_xton, sizeof(zomplex));
    h0mat = (double *) calloc(ist.n_xton * ist.n_xton, sizeof(double)); 
    
    if (1 == flag.isComplex) {
        calc_eh_kernel_cplx((zomplex*)psi_qp, pot_bare, pot_screened, pot_hartree, bsmat, direct, exchange, h0mat, eig_vals, &ist, &par, &flag, planfw, planbw, fftwpsi, &parallel);
    }
    // else if (0 == flag.isComplex){
    //     calc_eh_kernel_real((zomplex*) psi_qp, pot_bare, pot_screened, pot_hartree, bsmat, direct, exchange, h0mat, eig_vals, &ist, &par, &flag, planfw, planbw, fftwpsi);
    // }
    MPI_Barrier(MPI_COMM_WORLD);

    if (parallel.mpi_rank == 0){
    
    xton_ene = (double *) calloc(ist.n_xton, sizeof(double));
    bs_coeff = (zomplex *) calloc(ist.n_xton*ist.n_xton, sizeof(zomplex));
    bethe_salpeter(bsmat, direct, exchange, bs_coeff, h0mat, xton_ene,(zomplex*) psi_qp, s_mom,l_mom,l2_mom,LdotS, &grid, &ist, &par);
    
    calc_optical_exc(bs_coeff, xton_ene, trans_dipole, mag_dipole, &ist, &par);
    
    free(xton_ene); free(bs_coeff);
    }
    /***********************************************************************/
    free(psi_qp); 
    free(eig_vals); 
    free(grid.x); 
    free(grid.y); 
    free(grid.z);
    free(R);
    free(pot_hartree);
    free(pot_bare); 
    free(pot_screened);
    free(direct); 
    free(exchange); 
    free(bsmat); 
    free(h0mat);
    free(trans_dipole); 
    free(mag_dipole); 
    free(rot_strength);
    free(s_mom); 
    free(l_mom); 
    free(l2_mom); 
    free(LdotS);
    for (i = 0; i < ist.nthreads; i++){
        fftw_destroy_plan(planfw[i]); 
        fftw_destroy_plan(planbw[i]);
    }
    free(planfw);
    free(planbw);
    fftw_free(fftwpsi);
    free(ist.eval_elec_idxs);
    free(ist.eval_hole_idxs);
    free(ist.atom_types);

    time_t end_time = time(NULL);
    time_t end_clock = clock();

    if (mpir == 0) write_separation(stdout, top);
    if (mpir == 0) printf("\nDONE WITH PROGRAM: BETHE-SALPETHER\n");
    if (mpir == 0) printf("This calculation ended at: %s\n", ctime(&end_time)); 
    if (mpir == 0) printf("Total job CPU time (sec) %.4g, wall run time (sec) %.4g",
            ((double)end_clock - (double)start_clock)/(double)(CLOCKS_PER_SEC), (double)end_time - (double)start_time );fflush(0);
    if (mpir == 0) write_separation(stdout, bottom);
    
    free(top); free(bottom);
    MPI_Finalize();
    // exit(0);

    return 0;
}

/*****************************************************************************/
