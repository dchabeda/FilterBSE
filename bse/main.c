/*****************************************************************************/

#include "fd.h"

/*****************************************************************************/

int main(int argc, char *argv[]){
    /*****************************************************************
    * This is the main function for BSE.x. It is the driver that     *
    * controls memory allocation, structs, and the BSE algorithm.    *
    * The algorithm computes correlated excitonic states using the   *
    * Bethe-Salpether formalism described by Rohlfing and Louie:     *
    * Phys. Rev. B 62, 4927                                          *
    * It is applied to nanocrystal systems through the use of semi-  *
    * empirical pseudopotentials.                                    *
    ******************************************************************/ 

    // DECLARE VARIABLES AND STRUCTS
    // file pointers
    FILE *pf, *pmem;
    pmem = fopen("mem.dat", "w");
    // zomplex types
    zomplex *pot_direct, *pot_exchange, *pot_hartree; 
    zomplex *bsmat, *direct, *exchange;
    // custom structs
    par_st par; index_st ist; grid_st grid; flag_st flag; xyz_st *R = NULL; parallel_st parallel; 
    // FFT 
    fftw_plan_loc *planfw, *planbw; fftw_complex *fftwpsi;
    long fft_flags=0;
    // double arrays
    double *psitot = NULL, *psi_hole = NULL, *psi_elec = NULL, *psi;
    double *eig_vals = NULL, *sigma_E = NULL;
    double *h0mat;
    double *gridx = NULL, *gridy = NULL, *gridz = NULL;
    double *rho;
    double *mux, *muy, *muz, mx, my, mz, rsx, rsy, rsz;
    double *sx, *sy, *sz;
    double *lx, *ly, *lz, *lsqr, *ls;
    // long int arrays and counters
    long i, j, state;
    long jgrid, jgrid_real, jgrid_imag;
    long mem = 0;
    ist.atom_types = malloc(N_MAX_ATOM_TYPES*sizeof(ist.atom_types[0]));
    fprintf(pmem, "alloc ist.atom_types %ld B\n", N_MAX_ATOM_TYPES*sizeof(ist.atom_types[0])); mem += N_MAX_ATOM_TYPES*sizeof(ist.atom_types[0]);
    // Clock/Wall time output and stdout formatting
    time_t start_time = time(NULL); // Get the actual time for total wall runtime
    time_t start_clock = clock(); // Get the starting CPU clock time for total CPU runtime
    char *top; top = malloc(2*sizeof(top[0])); 
    fprintf(pmem, "alloc top %ld B\n", 2*sizeof(top[0])); mem += 2*sizeof(top[0]);
    char *bottom; bottom = malloc(2*sizeof(bottom[0]));
    fprintf(pmem, "alloc bottom %ld B\n", 2*sizeof(bottom[0])); mem += 2*sizeof(bottom[0]);
    strcpy(top, "T\0"); strcpy(bottom, "B\0");

    /*************************************************************************/
    fprintf(stdout, "******************************************************************************\n");
    printf("\nRUNNING PROGRAM: BETHE-SALPETHER\n");
    printf("This calculation began at: %s", ctime(&start_time)); 
    write_separation(stdout, bottom);
    fflush(stdout);

    /*************************************************************************/
    // 1. Initialize job from input file
    
    write_separation(stdout, top);
    printf("\n1.\tINITIALIZING JOB\n");
    write_separation(stdout, bottom); fflush(stdout);

    /*************************************************************************/
    /*** Read output from filter output.par ***/
    printf("\nReading filter output from output.dat:\n"); fflush(stdout);
    // ******
    // ******
    read_filter_output("output.dat", &psitot, &eig_vals, &sigma_E, &R, &grid, &gridx, &gridy, &gridz, &ist, &par, &flag);
    // ******
    // ******

    // Move the grid values into the grid struct arrays.
    // The memory allocation got too thorny, so we do this simple transfer
    
    if ((grid.x = malloc( grid.nx * sizeof(grid.x[0]))) == NULL ){
        fprintf(stderr, "ERROR: allocating memory for grid.x in main\n");
        exit(EXIT_FAILURE);
    }
    if ((grid.y = malloc(grid.ny * sizeof(grid.y[0]))) == NULL ){
            fprintf(stderr, "ERROR: allocating memory for grid.y in main\n");
            exit(EXIT_FAILURE);
        }
    if ((grid.z = malloc( grid.nz * sizeof(grid.z[0]))) == NULL ){
        fprintf(stderr, "ERROR: allocating memory for grid.z in main\n");
        exit(EXIT_FAILURE);
    }
    fprintf(pmem, "alloc grid.x %ld B\n", grid.nx * sizeof(grid.x[0])); mem += grid.nx * sizeof(grid.x[0]);
    fprintf(pmem, "alloc grid.y %ld B\n", grid.ny * sizeof(grid.y[0])); mem += grid.ny * sizeof(grid.y[0]);
    fprintf(pmem, "alloc grid.z %ld B\n", grid.nz * sizeof(grid.z[0])); mem += grid.nz * sizeof(grid.z[0]);

    for (i = 0; i < grid.nx; i++){
        grid.x[i] = gridx[i];
    }
    for (i = 0; i < grid.ny; i++){
        grid.y[i] = gridy[i];
    }
    for (i = 0; i < grid.nz; i++){
        grid.z[i] = gridz[i];
    }

    /*************************************************************************/
    /*** Read initial setup from input.par ***/
    printf("\nReading BSE job specifications from input.par:\n"); fflush(stdout);
    // ******
    // ******
    read_input(&flag, &grid, &ist, &par, &parallel);
    // ******
    // ******
    fflush(stdout);

    /*************************************************************************/
    /*** Determine the configuration of the quasiparticle basis ***/
    printf("\nSetting quasiparticle basis indices:\n"); fflush(stdout);
    
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
    eig_vals = realloc(eig_vals, ist.n_qp * sizeof(eig_vals[0]));
    sigma_E = realloc(sigma_E, ist.n_qp * sizeof(sigma_E[0]));
    
    // Resize the index arrays to tightly contain the indices of eigenstates in VB/CB
    // The conditional checks whether n_holes saturated the memory block (if so, no need to resize)
    // or if they are not equal, then n_holes is smaller and we should shrink the array
    if (ist.n_holes != ist.max_hole_states){
        ist.eval_hole_idxs = realloc(ist.eval_hole_idxs, ist.n_holes * sizeof(long));
    }
    if (ist.n_elecs != ist.max_elec_states){
        ist.eval_elec_idxs = realloc(ist.eval_elec_idxs, ist.n_elecs * sizeof(long));
    }
    fprintf(pmem, "alloc ist.eval_hole_idxs %ld B\n", ist.mn_states_tot * sizeof(ist.eval_hole_idxs[0])); mem += ist.mn_states_tot * sizeof(ist.eval_hole_idxs[0]);
    fprintf(pmem, "alloc ist.eval_elec_idxs %ld B\n", ist.mn_states_tot * sizeof(ist.eval_elec_idxs[0])); mem += ist.mn_states_tot * sizeof(ist.eval_elec_idxs[0]);
    

    /*************************************************************************/
    // Allocate memory for the electron and hole wavefunctions
    if ((psi_hole = (double *) malloc( ist.complex_idx * ist.nspinngrid * ist.n_holes * sizeof(psi_hole[0]))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for psi_hole in main.c\n");
        exit(EXIT_FAILURE);
    }
    if ((psi_elec = (double *) malloc( ist.complex_idx * ist.nspinngrid * ist.n_elecs * sizeof(psi_elec[0]))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for psi_elec in main.c\n");
        exit(EXIT_FAILURE);
    }
    if ((psi = (double *) malloc( ist.complex_idx * ist.nspinngrid * (ist.n_holes + ist.n_elecs) * sizeof(double))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for psi_elec in main.c\n");
        exit(EXIT_FAILURE);
    }
    fprintf(pmem, "alloc psi_hole %ld B\n", ist.complex_idx * ist.nspinngrid * ist.n_holes * sizeof(psi_hole[0])); mem += ist.complex_idx*ist.nspinngrid*ist.n_holes*sizeof(psi_hole[0]);
    fprintf(pmem, "alloc psi_elec %ld B\n", ist.complex_idx * ist.nspinngrid * ist.n_elecs * sizeof(psi_elec[0])); mem += ist.complex_idx*ist.nspinngrid*ist.n_elecs*sizeof(psi_elec[0]);
    fprintf(pmem, "alloc psi %ld B", ist.complex_idx * ist.nspinngrid * (ist.n_holes + ist.n_elecs) * sizeof(double)); mem += ist.complex_idx * ist.nspinngrid * (ist.n_holes + ist.n_elecs) * sizeof(double);
    write_separation(pmem, top);
    fprintf(pmem, "\ntotal mem usage %ld GB\n", mem * ((long)1e-9) );
    write_separation(pmem, bottom); fflush(pmem);

    printf("\nTHE HOLES:\n");
    for (i = 0; i < ist.n_holes; i++){
        printf("%ld\n", ist.eval_hole_idxs[i]);
    }
    printf("\nTHE ELECS:\n");
    for (i = 0; i < ist.n_elecs; i++){
        printf("%ld\n", ist.eval_elec_idxs[i]);
    }
    
    // ******
    // ******
    get_qp_basis(psi, psitot, psi_hole, psi_elec, eig_vals, sigma_E, &ist, &par, &flag);
    // ******
    // ******
    free(psitot); free(psi_hole); free(psi_elec); free(sigma_E);
    fprintf(pmem, "free psitot %ld B", ist.complex_idx * ist.nspinngrid * ist.mn_states_tot * sizeof(double)); mem -= ist.complex_idx * ist.nspinngrid * ist.mn_states_tot * sizeof(double);
    fprintf(pmem, "free psi_hole %ld B\n", ist.complex_idx * ist.nspinngrid * ist.n_holes * sizeof(psi_hole[0])); mem -= ist.complex_idx*ist.nspinngrid*ist.n_holes*sizeof(psi_hole[0]);
    fprintf(pmem, "free psi_elec %ld B\n", ist.complex_idx * ist.nspinngrid * ist.n_elecs * sizeof(psi_elec[0])); mem -= ist.complex_idx*ist.nspinngrid*ist.n_elecs*sizeof(psi_elec[0]);
    fprintf(pmem, "free sigma_E %ld B", ist.mn_states_tot * sizeof(double)); mem -= ist.mn_states_tot * sizeof(double);
    write_separation(pmem, top);
    fprintf(pmem, "\ntotal mem usage %ld GB\n", mem * ((long)1e-9) );
    write_separation(pmem, bottom); fflush(pmem);

    // print cube files for debug

    // rho = malloc(ist.ngrid * sizeof(rho[0]));
    // char str[100];

    // for (i = 0; i < 2; i++){
    //     //Spin Up Wavefunction
    //     sprintf(str,"hole-%ld-Up.cube",ist.eval_hole_idxs[i]);
    //     for (jgrid = 0; jgrid < ist.ngrid; jgrid++){
    //     jgrid_real = ist.complex_idx * jgrid;
    //     jgrid_imag = ist.complex_idx * jgrid + 1;
        
    //     rho[jgrid] = sqr(psi_hole[ist.complex_idx*i*ist.nspinngrid + jgrid_real]);
    //     if (1 == flag.isComplex) rho[jgrid] += sqr(psitot[ist.complex_idx*i*ist.nspinngrid + jgrid_imag]);
    //     }
    //     write_cube_file(rho, &grid, str);
    //     //Spin Down Wavefunction
    //     if (1 == flag.useSpinors){    
    //     sprintf(str,"hole-%ld-Dn.cube", ist.eval_hole_idxs[i]);
    //     for (jgrid = 0; jgrid < ist.ngrid; jgrid++){
    //         jgrid_real = ist.complex_idx * jgrid;
    //         jgrid_imag = ist.complex_idx * jgrid + 1;
            
    //         rho[jgrid] = sqr(psitot[ist.complex_idx*(i*ist.nspinngrid+ist.ngrid)+jgrid_real]) 
    //             + sqr(psitot[ist.complex_idx*(i*ist.nspinngrid+ist.ngrid)+jgrid_imag]);    
    //     }
    //     write_cube_file(rho, &grid, str);
    //     } 
    // }

    // for (i = 0;  i < 2; i++){
    //     sprintf(str,"elec-%ld-Up.cube", ist.eval_elec_idxs[i]);
    //     for (jgrid = 0; jgrid < ist.ngrid; jgrid++){
    //     jgrid_real = ist.complex_idx * jgrid;
    //     jgrid_imag = ist.complex_idx * jgrid + 1;
        
    //     rho[jgrid] = sqr(psitot[ist.complex_idx*i*ist.nspinngrid + jgrid_real]);
    //     if (1 == flag.isComplex) rho[jgrid] += sqr(psitot[ist.complex_idx*i*ist.nspinngrid + jgrid_imag]);
    //     }
    //     write_cube_file(rho, &grid, str);

    //     if (1 == flag.useSpinors){
    //     sprintf(str,"elec-%ld-Dn.cube", ist.eval_elec_idxs[i]);
    //     for (jgrid = 0; jgrid < ist.ngrid; jgrid++){
    //         jgrid_real = ist.complex_idx * jgrid;
    //         jgrid_imag = ist.complex_idx * jgrid + 1;
        
    //         rho[jgrid] = sqr(psitot[ist.complex_idx*(i*ist.nspinngrid+ist.ngrid)+jgrid_real]) 
    //             + sqr(psitot[ist.complex_idx*(i*ist.nspinngrid+ist.ngrid)+jgrid_imag]);
    //     }
    //     write_cube_file(rho, &grid, str);
    //     }
    // }
    // free(rho);

    /*************************************************************************/
    // 2. Compute electron-hole interaction potentials
    
    write_separation(stdout, top);
    printf("\n2.\tCOMPUTING ELECTRON-HOLE INTERACTION KERNEL\n");
    write_separation(stdout, bottom); fflush(stdout);

    /*************************************************************************/
    // Allocate memory for the hartree, screened coulomb and bare exchange pots
    printf("Allocating memory for hartree, direct, and exchange potentials..."); fflush(stdout);
    pot_direct  = (zomplex *) malloc(ist.ngrid * sizeof(pot_direct[0]));
    pot_exchange  = (zomplex *) malloc(ist.ngrid * sizeof(pot_exchange[0]));
    pot_hartree = (zomplex *) malloc(ist.ngrid*parallel.nthreads*sizeof(pot_hartree[0]));
    fprintf(pmem, "alloc pot_direct %ld B\n", ist.ngrid*sizeof(pot_direct[0])); mem += ist.ngrid * sizeof(pot_direct[0]);
    fprintf(pmem, "alloc pot_exchange %ld B\n", ist.ngrid*sizeof(pot_exchange[0])); mem += ist.ngrid * sizeof(pot_exchange[0]);
    fprintf(pmem, "alloc pot_hartree %ld B\n", ist.ngrid*parallel.nthreads*sizeof(pot_direct[0])); mem += ist.ngrid*parallel.nthreads*sizeof(pot_direct[0]);
    printf(" done.\n"); fflush(stdout);

    // Initialize the FFT arrays for parallel Fourier transform
    fftwpsi = fftw_malloc(ist.ngrid*parallel.nthreads*sizeof(fftw_complex));
    fprintf(pmem, "alloc fftwpsi %ld B\n", ist.ngrid*parallel.nthreads*sizeof(fftw_complex)); mem += ist.ngrid*parallel.nthreads*sizeof(fftw_complex);
    
    // Initialize the parallel FFT
    fftw_plan_with_nthreads(ist.nthreads);
  
    planfw = (fftw_plan_loc *) malloc(parallel.nthreads * sizeof(fftw_plan_loc));
    planbw = (fftw_plan_loc *) malloc(parallel.nthreads * sizeof(fftw_plan_loc));
    for (i = 0; i < ist.nthreads; i++) { 
        planfw[i] = fftw_plan_dft_3d(grid.nz, grid.ny, grid.nx, &fftwpsi[i*ist.ngrid], 
                                     &fftwpsi[i*ist.ngrid], FFTW_FORWARD, fft_flags);
        planbw[i] = fftw_plan_dft_3d(grid.nz, grid.ny, grid.nx, &fftwpsi[i*ist.ngrid],
                                     &fftwpsi[i*ist.ngrid], FFTW_BACKWARD, fft_flags);
    }
    
    init_elec_hole_kernel(pot_direct, pot_exchange, &grid, &par, &ist, planfw[0], planbw[0], &fftwpsi[0]);

    pf = fopen("qp_spins.dat", "w");
    for (int state  = 0; state < ist.n_qp; state++){
        fprintf(pf, "\nStats on state%d: (E=%lg)\n", state, eig_vals[state]);
        double perUp = 0;
        double perDn = 0;
        for (long jgrid = 0; jgrid < ist.ngrid; jgrid++){
            jgrid_real = ist.complex_idx * jgrid;
            jgrid_imag = ist.complex_idx * jgrid + 1;
        
            perUp+=sqr(psi[state*ist.nspinngrid+jgrid_real])+sqr(psi[state*ist.nspinngrid+jgrid_imag]);
            perDn+=sqr(psi[state*ist.nspinngrid+ist.ngrid+jgrid_real])+sqr(psi[state*ist.nspinngrid+ist.ngrid+jgrid_imag]);
        }
        fprintf(pf, " Spin up fraction: %f\n", perUp * par.dv);
        fprintf(pf, " Spin dn fraction: %f\n", perDn * par.dv);
    }
    fclose(pf);
    /*************************************************************************/
    
    // /**************************************************************************/
    // /*** this routine computes the coulomb coupling between
    //      single excitons.  On input - it requires the eigenstates stored in psi,
    //      the eigenvalues stored in eval and pot_hartree computed in init_elec_hole_kernel.
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

    // ist.ms2 = ist.n_elecs * ist.n_holes;
    // bsmat = (zomplex *) calloc(ist.ms2*ist.ms2, sizeof(zomplex));
    // direct = (zomplex *) calloc(ist.ms2*ist.ms2, sizeof(zomplex)); 
    // exchange = (zomplex *) calloc(ist.ms2*ist.ms2, sizeof(zomplex)); 
    // h0mat = (double *) calloc(ist.ms2*ist.ms2, sizeof(double)); 
    
    // single_coulomb_openmp(psi, pot_direct, pot_exchange, pot_hartree, eig_vals, ist, par, planfw, planbw, fftwpsi, bsmat,direct,exchange, h0mat);

    // ppsi = fopen("bsRE.dat", "w");
    // for (i = 0; i < ist.ms2; i++, fprintf(ppsi,"\n"))
    //     for (j = 0; j < ist.ms2; j++)
    //         fprintf(ppsi,"%.*g\t", DBL_DIG, bsmat[i*ist.ms2+j].re);
    // fclose(ppsi);

    // ppsi = fopen("bsIM.dat", "w");
    // for (i = 0; i < ist.ms2; i++, fprintf(ppsi,"\n"))
    //     for (j = 0; j < ist.ms2; j++)
    //         fprintf(ppsi,"%.*g\t", DBL_DIG, bsmat[i*ist.ms2+j].im);
    // fclose(ppsi);

    // ppsi = fopen("h0.dat", "w");
    // for (i = 0; i < ist.ms2; i++, fprintf(ppsi,"\n"))
    //     for (j = 0; j < ist.ms2; j++)
    //          fprintf(ppsi,"%.*g\t", DBL_DIG, h0mat[i*ist.ms2+j]);
    // fclose(ppsi);


    // sx = (zomplex *) calloc(ist.n_holes*ist.n_holes+ist.n_elecs*ist.n_elecs, sizeof(zomplex)); //<psi_r|Sx|psi_s>
    // sy = (zomplex *) calloc(ist.n_holes*ist.n_holes+ist.n_elecs*ist.n_elecs, sizeof(zomplex)); //<psi_r|Sy|psi_s>
    // sz = (zomplex *) calloc(ist.n_holes*ist.n_holes+ist.n_elecs*ist.n_elecs, sizeof(zomplex)); //<psi_r|Sz|psi_s>
    // lx = (zomplex *) calloc(ist.n_holes*ist.n_holes+ist.n_elecs*ist.n_elecs, sizeof(zomplex)); //<psi_r|Lx|psi_s>
    // ly = (zomplex *) calloc(ist.n_holes*ist.n_holes+ist.n_elecs*ist.n_elecs, sizeof(zomplex)); //<psi_r|Ly|psi_s>
    // lz = (zomplex *) calloc(ist.n_holes*ist.n_holes+ist.n_elecs*ist.n_elecs, sizeof(zomplex)); //<psi_r|Lz|psi_s>
    // lsqr = (zomplex *) calloc(ist.n_holes*ist.n_holes+ist.n_elecs*ist.n_elecs, sizeof(zomplex)); //<psi_r|L^2|psi_s>
    // ls  = (zomplex *) calloc(ist.n_holes*ist.n_holes+ist.n_elecs*ist.n_elecs, sizeof(zomplex)); //<psi_r|L*S|psi_s>
    // mux = (zomplex *) calloc(ist.n_elecs*ist.n_holes, sizeof(zomplex)); // <psi_i|ux|psi_a>
    // muy = (zomplex *) calloc(ist.n_elecs*ist.n_holes, sizeof(zomplex)); // <psi_i|uy|psi_a>
    // muz = (zomplex *) calloc(ist.n_elecs*ist.n_holes, sizeof(zomplex)); // <psi_i|uz|psi_a>
    // //mx  = (double *) calloc(ist.n_elecs*ist.n_holes, sizeof(double)); // <psi_i|mx|psi_a>
    // //my  = (double *) calloc(ist.n_elecs*ist.n_holes, sizeof(double)); // <psi_i|my|psi_a>
    // //mz  = (double *) calloc(ist.n_elecs*ist.n_holes, sizeof(double)); // <psi_i|mz|psi_a>
    // //rs  = (double *) calloc(ist.n_elecs*ist.n_holes, sizeof(double)); // <psi_a|u|psi_i>.<psi_a|m|psi_i>

    // spins(sx,sy,sz,psi,ist,par);
    // angular(lx,ly,lz,lsqr,ls,psi,&grid,&ist,&par,planfw[0], planbw[0], &fftwpsi[0]);
    // dipole(psi, mux, muy, muz, eig_vals, &grid, &ist, &par);
    // //TODO? make this work for spinors?//mag_dipole(vx, vy, vz, psi, mx, my, mz, eig_vals, planfw, planbw, fftwpsi, ist, par);
    // //TODO? make this work for spinors?//rotational_strength(rs, mux, muy, muz, mx, my, mz, eig_vals, ist);
    // bethe_salpeter(bsmat, direct, exchange, h0mat, psi, grid.z, mux, muy, muz, mx, my, mz,sx,sy,sz,lx,ly,lz,lsqr,ls, ist, par);
    
    
    
    /***********************************************************************/
    free(psi); 
    free(eig_vals); 
    free(pot_hartree); free(pot_direct); free(pot_exchange);
    // free(bsmat); free(h0mat); free(mux); free(muy); free(muz);
    // //free(mx); free(my); free(mz); free(rs);
    // free(sx); free(sy); free(sz);
    // free(lx); free(ly); free(lz); free(lsqr);
    free(planfw); free(planbw);

    time_t end_time = time(NULL);
    time_t end_clock = clock();

    write_separation(stdout, top);
    printf("\nDONE WITH PROGRAM: BETHE-SALPETHER\n");
    printf("This calculation ended at: %s\n", ctime(&end_time)); 
    printf("Total job CPU time (sec) %.4g, wall run time (sec) %.4g",
            ((double)end_clock - (double)start_clock)/(double)(CLOCKS_PER_SEC), (double)end_time - (double)start_time );fflush(0);
    write_separation(stdout, bottom);
    
    free(top); free(bottom);
    // exit(0);

    return 0;
}

/*****************************************************************************/
