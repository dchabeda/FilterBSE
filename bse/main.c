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
    FILE *ppsi, *peval, *pf;  
    // zomplex types
    zomplex *potq, *potqx, *poth, *psi;
    zomplex *mux, *muy, *muz, mx, my, mz, rsx, rsy, rsz; 
    zomplex *bsmat, *direct, *exchange;
    zomplex *sx, *sy, *sz;
    zomplex *lx, *ly, *lz, *lsqr, *ls;
    // custom structs
    par_st par; index_st ist; grid_st grid; flag_st flag; xyz_st *R = NULL; parallel_st parallel; 
    // FFT 
    fftw_plan_loc *planfw, *planbw; fftw_complex *fftwpsi;
    long fft_flags=0;
    // double arrays
    double *psitot = NULL, *psihomo = NULL, *psilumo = NULL;;
    double *eig_vals = NULL, *sigma_E = NULL;
    double *h0mat;
    double *gridx = NULL, *gridy = NULL, *gridz = NULL;
    double *rho;
    // long int arrays and counters
    long i, a, j, thomo, tlumo, indexfirsthomo;
    long jgrid, jgrid_real, jgrid_imag;
    ist.atom_types = malloc(N_MAX_ATOM_TYPES*sizeof(ist.atom_types[0]));
    // Clock/Wall time output and stdout formatting
    time_t start_time = time(NULL); // Get the actual time for total wall runtime
    time_t start_clock = clock(); // Get the starting CPU clock time for total CPU runtime
    char *top; top = malloc(2*sizeof(top[0])); 
    char *bottom; bottom = malloc(2*sizeof(bottom[0]));
    strcpy(top, "T\0"); strcpy(bottom, "B\0");

    /*************************************************************************/
    fprintf(stdout, "******************************************************************************\n");
    printf("\nRUNNING PROGRAM: BETHE-SALPETHER\n");
    printf("This calculation began at: %s", ctime(&start_time)); 
    write_separation(stdout, bottom);
    fflush(stdout);

    /*************************************************************************/
    // Initialize job from input file
    
    write_separation(stdout, top);
    printf("\n1.\tINITIALIZING JOB\n");
    write_separation(stdout, bottom); fflush(stdout);

    /*** Read output from filter output.par ***/
    printf("\nReading filter output from output.dat:\n");

    read_filter_output("output.dat", &psitot, &eig_vals, &sigma_E, &R, &grid, &gridx, &gridy, &gridz, &ist, &par, &flag);

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

    for (i = 0; i < grid.nx; i++){
        grid.x[i] = gridx[i];
    }
    for (i = 0; i < grid.ny; i++){
        grid.y[i] = gridy[i];
    }
    for (i = 0; i < grid.nz; i++){
        grid.z[i] = gridz[i];
    }

    /*** Read initial setup from input.par ***/
    printf("\nReading BSE job specifications from input.par:\n");
    read_input(&flag, &grid, &ist, &par, &parallel);
    fflush(stdout);

    /*** Determine the configuration of the quasiparticle basis ***/
    printf("\nSetting quasiparticle basis indices:\n");
    get_qp_basis_indices(eig_vals, sigma_E, &ist, &par, &flag);
    // Reallocate the eig_vals and sigma_E arrays to only contain the n_qp states
    eig_vals = realloc(eig_vals, ist.n_qp * sizeof(eig_vals[0]));
    sigma_E = realloc(sigma_E, ist.n_qp * sizeof(sigma_E[0]));

    for (i = 0; i < ist.n_qp; i++){
        printf("% .12f %lg\n", eig_vals[i], sigma_E[i]);
    }

    // // init_size(argc, argv, &par, &ist);
  
    // /*************************************************************************/
    // fftwpsi = fftw_malloc(sizeof (fftw_complex )*ist.ngrid*ist.nthreads);
    // potq  = (zomplex *) calloc(ist.ngrid, sizeof(zomplex));
    // potqx  = (zomplex *) calloc(ist.ngrid, sizeof(zomplex));
    // poth = (zomplex *) calloc(ist.ngrid*ist.nthreads, sizeof(zomplex));
    
    // /**************************************************************************/
    // //Seems to initilize con figuration, grid, and ksqr grid
    // // init(potl, vx, vy, vz, ksqr, rx, ry, rz, &par, &ist);
    // //printf("Done with init"); fflush(0);
    // /*** initialization for the fast Fourier transform ***/
    // fftw_plan_with_nthreads(ist.nthreads);
  
    // planfw = (fftw_plan_loc *) calloc(ist.nthreads, sizeof(fftw_plan_loc));
    // planbw = (fftw_plan_loc *) calloc(ist.nthreads, sizeof(fftw_plan_loc));
    // for (i = 0; i < ist.nthreads; i++) { 
    //     planfw[i] = fftw_plan_dft_3d(ist.nz, ist.ny, ist.nx, &fftwpsi[i*ist.ngrid], 
    //                                  &fftwpsi[i*ist.ngrid], FFTW_FORWARD, fft_flags);
    //     planbw[i] = fftw_plan_dft_3d(ist.nz, ist.ny, ist.nx, &fftwpsi[i*ist.ngrid],
    //                                  &fftwpsi[i*ist.ngrid], FFTW_BACKWARD, fft_flags);
    // }
    
    // init_pot(potq, potqx, &grid, &par, &ist, planfw[0], planbw[0], &fftwpsi[0]);

    // /*************************************************************************/
    // eig_vals = (double *) calloc(ist.nhomo+1, sizeof(double)); 
    // sigma_E = (double *) calloc(ist.nhomo+1, sizeof(double)); 

    // peval = fopen("eval.par" , "r");
    // for (i = 0; i < ist.nhomo+1; i++)
    //     fscanf (peval,"%ld %lg %lg",&a,&eig_vals[i],&sigma_E[i]);
    // fclose (peval);

    // peval = fopen("eval.dat" , "w");
    // for (i = 0; i < ist.nhomo+1; i++){
    //     fprintf (peval,"%ld %g %g\n",i,eig_vals[i],sigma_E[i]);
    // }
    // fclose(peval);
  
    // for (thomo = 0, i = ist.nhomo; i >= 0; i--){
    //     if (sigma_E[i] < par.sigma_E_cut) thomo++;
    //     if (thomo == ist.total_homo) {
    //         indexfirsthomo = i;
    //         break;
    //     }
    // }

    // printf("The index of lowest energy occupied level used = %ld\n", indexfirsthomo); 

    // free(eig_vals);
    // free(sigma_E);

    // psidummy = (zomplex *) calloc(ist.nspinngrid, sizeof(zomplex));
    // eig_vals = (double *)calloc(ist.n_qp, sizeof(double)); 
    // sigma_E = (double *)calloc(ist.n_qp, sizeof(double)); 

    // /**********************************************************************/
    // /*              CHANGES TO THE FOLLOWING CODE START HERE              */ 
    // /**********************************************************************/
    // /* Ideally I would prefer it to be outside of main, as its logic is more 
    //  * specific than "run the BSE Calc", but in the interest of time I'm 
    //  * going to leave it like this and focus back on calculating auger rates
    //  */

    // /**********************************************************************/
    // /*               Read in eval.par and store it in memory              */
    // /**********************************************************************/
    // struct eindex { // Struct I used to store index evals and sigma_E_cut in memory
    //     long index;
    //     double evalue;
    //     double sigma_E_cut;
    // }; 
    
    // struct eindex *evalindex;
    // int size = 1024; // Initial size of array in bytes, is allowed to resize
    // long counter = 0;

    // if ((evalindex = calloc(size, sizeof(struct eindex))) == NULL) {
    //     fprintf(stderr, "Memory error in main\n");
    //     exit(EXIT_FAILURE);
    // }

    // peval = fopen("eval.par" , "r");
    // if (peval) {
    //      while (fscanf(peval, "%ld %lg %lg", &evalindex[counter].index, 
    //             &evalindex[counter].evalue, &evalindex[counter].sigma_E_cut) == 3) {
    //         counter++;
    //         if (counter == size - 1) {
    //             size *= 2;
    //             // TODO: This does not check if realloc succeeded or not
    //             evalindex = realloc(evalindex, size * sizeof(struct eindex));
    //         }
    //     }
    // }
    // fclose(peval);

    // //TODO???: Check if this was a real-scalar filter calc and if so convert to complex spinor rep



    // // Read ead in psi.par as we have already read in psi.par
    // ppsi = fopen("psi.par" , "r");
    // if (ppsi ==NULL) {
    //     printf("no psi.par in cwd\n");
    //     exit(EXIT_FAILURE);
    // }
    // printf("allocating memory for %ld hole and %ld electron states\n",(ist.total_homo),(ist.total_lumo));
    // psihomo = calloc((ist.total_homo) * ist.nspinngrid, sizeof(zomplex));
    // psilumo = calloc((ist.total_lumo) * ist.nspinngrid, sizeof(zomplex));
    // if (!psihomo || !psilumo) {printf("Failed to allocate memory for psihomo/psilumo\n");exit(EXIT_FAILURE);}

    // long foffset = ist.nspinngrid * sizeof(zomplex);  // for random access 
    // char fname[80] = {0};
    // long nstates = 0;  // Total number of states
    
    // fseek(ppsi, foffset * ist.nhomo, SEEK_SET); 
    // counter = ist.nhomo; // Set loop counter to HOMO index
    // thomo = 0;           // Counter for hole states used
    // while (counter >= indexfirsthomo && thomo < ist.max_hole_states) {
    //     fread(psidummy, sizeof(zomplex), ist.nspinngrid, ppsi);
    //     if (evalindex[counter].sigma_E_cut < par.sigma_E_cut) {    
    //         normalize_zomplex(psidummy, par.dv, ist.nspinngrid);
    //         //sprintf(fname, "pzv%ld.dat", counter);
    //         eig_vals[nstates] = evalindex[counter].evalue;
    //         sigma_E[nstates] = evalindex[counter].sigma_E_cut;
    //         for (j = 0; j < ist.nspinngrid; j++) {
    //             psihomo[thomo * ist.nspinngrid + j] = psidummy[j];
    //         }
    //         nstates++;
    //         thomo++;
			
    //     }
    //     counter--;
        
        
    //     int ret =0;
    //     if (counter >=0) {
    //         ret= fseek(ppsi, -2 * foffset, SEEK_CUR);
            
    //     }
    //     if(ret) printf("nonzero exit in fseek!\n"); fflush(0); 
        
    // }

    // fseek(ppsi, foffset * (ist.nhomo + 1), SEEK_SET);
    // counter = ist.nhomo + 1;
    // tlumo = 0;

    // while (tlumo < ist.max_elec_states && tlumo < ist.total_lumo) {
    //     fread(psidummy, sizeof(zomplex), ist.nspinngrid, ppsi);
    //     if (evalindex[counter].sigma_E_cut < par.sigma_E_cut) {
    //         normalize_zomplex(psidummy, par.dv, ist.nspinngrid);
    //         //sprintf(fname, "pzc%ld.dat", counter);
            
	// 	    eig_vals[nstates] = evalindex[counter].evalue;
    //         sigma_E[nstates] = evalindex[counter].sigma_E_cut;
    //         for (j = 0; j < ist.nspinngrid; j++) {      
    //             psilumo[tlumo * ist.nspinngrid + j] = psidummy[j];
    //         }
    //     	nstates++;
    //     	tlumo++;
			
    //     }
    //     counter++;
    //     // printf("tlumo: %d nstates:%d counter:%d\n",tlumo, nstates, counter); fflush(0);

    // }
    // free(psidummy);  printf("freeing psidummy\n");fflush(0);
    // psi = (zomplex *) calloc(nstates, ist.nspinngrid*sizeof(zomplex));
    // if (!psi) {printf("Failed to allocate memory for psi"); fflush(0); exit(EXIT_FAILURE);}
    // fclose(ppsi); printf("closed ppsi");fflush(0);
    
    // /**********************************************************************/

    // for (i = thomo - 1, a = 0; i >= 0 && a < thomo; i--, a++) {
    //     for (j = 0; j < ist.nspinngrid; j++) {
    //         psi[a * ist.nspinngrid + j].re = psihomo[i * ist.nspinngrid + j].re;
    //         psi[a * ist.nspinngrid + j].im = psihomo[i * ist.nspinngrid + j].im;
    //     }
    // }

    // for (i = thomo, a = 0; i < tlumo + thomo && a < tlumo; i++, a++) {
    //     for (j = 0; j < ist.nspinngrid; j++) {
    //         psi[i * ist.nspinngrid + j].re = psilumo[a * ist.nspinngrid + j].re;
    //         psi[i * ist.nspinngrid + j].im = psilumo[a * ist.nspinngrid + j].im;
    //     }
    // }

    // double tmp1, tmp2;
    // /* Rearrange eig_vals */
    // for (i = 0, j = thomo - 1; i < j; i++, j--) {
    //     tmp1 = eig_vals[i];    tmp2 = sigma_E[i];
    //     eig_vals[i] = eig_vals[j]; sigma_E[i] = sigma_E[j];
    //     eig_vals[j] = tmp1;    sigma_E[j] = tmp2;
    // }
    // /* Write smaller psi.par and eval.par for augerBSE */
    // FILE *new_eval = fopen("BSEeval.par", "w");
    // FILE *new_psi = fopen("BSEpsi.par", "w");
    // if (!new_eval || !new_psi) {printf("Failed opening BSEeval.par");fflush(0); exit(EXIT_FAILURE);}

    // fwrite(psi, sizeof(zomplex) * ist.nspinngrid, nstates, new_psi);

    // for (i = 0; i < nstates; i++) {
    //     fprintf(new_eval, "%ld %.*g %.*g\n", i, DBL_DIG, eig_vals[i], DBL_DIG, sigma_E[i]); 
    // }
    // fclose(new_eval);
    // fclose(new_psi);
    // printf("closed eval and psi");fflush(0);
    // /**********************************************************************/
    
    
    // free(evalindex); printf("freeing evalindex\n");fflush(0);
    // free(psihomo); printf("freeing psihomo\n");fflush(0);
    // free(psilumo); printf("freeing psidlumo\n");fflush(0);
    
    // printf("freeing Things");fflush(0);

    // printf("Number of hole eigenstates used in the BSE calculation = %ld\n", thomo);
    // printf("Number of electron eigenstates used in the BSE calculation =  %ld\n", tlumo);
    // printf("Total number of eigenstates used in the BSE calculation = %ld\n", nstates);
    // fflush(0);

    // pf = fopen("qp_spins.dat", "w");
    // for (int state  = 0; state<nstates; state++){
    //     fprintf(pf, "\nStats on state%d: (E=%lg)\n",state, eig_vals[state]);
    //     double perUp = 0;
    //     double perDn = 0;
    //     for (long igrid = 0; igrid < ist.ngrid; igrid++){
    //         perUp+=sqr(psi[state*ist.nspinngrid+igrid].re)+sqr(psi[state*ist.nspinngrid+igrid].im);
    //         perDn+=sqr(psi[state*ist.nspinngrid+ist.ngrid+igrid].re)+sqr(psi[state*ist.nspinngrid+ist.ngrid+igrid].im);
    //     }
    //     fprintf(pf, " Spin up fraction: %f\n", perUp*par.dv);
    //     fprintf(pf, " Spin dn fraction: %f\n", perDn*par.dv);
    // }
    // fclose(pf);
    // /*************************************************************************/
    // ist.n_qp = nstates;
    // ist.nlumo = thomo;
    // ist.nhomo = thomo - 1;
    // ist.total_lumo = tlumo;
    // ist.total_homo = thomo;
  
    // /*                              END CHANGES                           */
    // /*************************************************************************/

    // /**************************************************************************/
    
    // /**************************************************************************/
    // /*** this routine computes the coulomb coupling between
    //      single excitons.  On input - it requires the eigenstates stored in psi,
    //      the eigenvalues stored in eval and poth computed in init_pot.
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

    // ist.ms2 = ist.total_lumo * ist.total_homo;
    // bsmat = (zomplex *) calloc(ist.ms2*ist.ms2, sizeof(zomplex));
    // direct = (zomplex *) calloc(ist.ms2*ist.ms2, sizeof(zomplex)); 
    // exchange = (zomplex *) calloc(ist.ms2*ist.ms2, sizeof(zomplex)); 
    // h0mat = (double *) calloc(ist.ms2*ist.ms2, sizeof(double)); 
    
    // single_coulomb_openmp(psi, potq, potqx, poth, eig_vals, ist, par, planfw, planbw, fftwpsi, bsmat,direct,exchange, h0mat);

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


    // sx = (zomplex *) calloc(ist.total_homo*ist.total_homo+ist.total_lumo*ist.total_lumo, sizeof(zomplex)); //<psi_r|Sx|psi_s>
    // sy = (zomplex *) calloc(ist.total_homo*ist.total_homo+ist.total_lumo*ist.total_lumo, sizeof(zomplex)); //<psi_r|Sy|psi_s>
    // sz = (zomplex *) calloc(ist.total_homo*ist.total_homo+ist.total_lumo*ist.total_lumo, sizeof(zomplex)); //<psi_r|Sz|psi_s>
    // lx = (zomplex *) calloc(ist.total_homo*ist.total_homo+ist.total_lumo*ist.total_lumo, sizeof(zomplex)); //<psi_r|Lx|psi_s>
    // ly = (zomplex *) calloc(ist.total_homo*ist.total_homo+ist.total_lumo*ist.total_lumo, sizeof(zomplex)); //<psi_r|Ly|psi_s>
    // lz = (zomplex *) calloc(ist.total_homo*ist.total_homo+ist.total_lumo*ist.total_lumo, sizeof(zomplex)); //<psi_r|Lz|psi_s>
    // lsqr = (zomplex *) calloc(ist.total_homo*ist.total_homo+ist.total_lumo*ist.total_lumo, sizeof(zomplex)); //<psi_r|L^2|psi_s>
    // ls  = (zomplex *) calloc(ist.total_homo*ist.total_homo+ist.total_lumo*ist.total_lumo, sizeof(zomplex)); //<psi_r|L*S|psi_s>
    // mux = (zomplex *) calloc(ist.total_lumo*ist.total_homo, sizeof(zomplex)); // <psi_i|ux|psi_a>
    // muy = (zomplex *) calloc(ist.total_lumo*ist.total_homo, sizeof(zomplex)); // <psi_i|uy|psi_a>
    // muz = (zomplex *) calloc(ist.total_lumo*ist.total_homo, sizeof(zomplex)); // <psi_i|uz|psi_a>
    // //mx  = (double *) calloc(ist.total_lumo*ist.total_homo, sizeof(double)); // <psi_i|mx|psi_a>
    // //my  = (double *) calloc(ist.total_lumo*ist.total_homo, sizeof(double)); // <psi_i|my|psi_a>
    // //mz  = (double *) calloc(ist.total_lumo*ist.total_homo, sizeof(double)); // <psi_i|mz|psi_a>
    // //rs  = (double *) calloc(ist.total_lumo*ist.total_homo, sizeof(double)); // <psi_a|u|psi_i>.<psi_a|m|psi_i>

    // spins(sx,sy,sz,psi,ist,par);
    // angular(lx,ly,lz,lsqr,ls,psi,&grid,&ist,&par,planfw[0], planbw[0], &fftwpsi[0]);
    // dipole(psi, mux, muy, muz, eig_vals, &grid, &ist, &par);
    // //TODO? make this work for spinors?//mag_dipole(vx, vy, vz, psi, mx, my, mz, eig_vals, planfw, planbw, fftwpsi, ist, par);
    // //TODO? make this work for spinors?//rotational_strength(rs, mux, muy, muz, mx, my, mz, eig_vals, ist);
    // bethe_salpeter(bsmat, direct, exchange, h0mat, psi, grid.z, mux, muy, muz, mx, my, mz,sx,sy,sz,lx,ly,lz,lsqr,ls, ist, par);
    
    
    
    /***********************************************************************/
    // free(psi); 
    free(eig_vals); free(sigma_E); 
    // free(poth); free(potq); free(potqx);
    // free(bsmat); free(h0mat); free(mux); free(muy); free(muz);
    // //free(mx); free(my); free(mz); free(rs);
    // free(sx); free(sy); free(sz);
    // free(lx); free(ly); free(lz); free(lsqr);
    // free(planfw); free(planbw);

    time_t end_time = time(NULL);
    time_t end_clock = clock();

    write_separation(stdout, top);
    printf("\nDONE WITH PROGRAM: BETHE-SALPETHER\n");
    printf("This calculation ended at: %s\n", ctime(&end_time)); 
    printf("Total job CPU time (sec) %.4g, wall run time (sec) %.4g",
            ((double)end_clock - (double)start_clock)/(double)(CLOCKS_PER_SEC), (double)end_time - (double)start_time );fflush(0);
    write_separation(stdout, bottom);
    
    free(top); free(bottom);
    exit(0);

    return 0;
}

/*****************************************************************************/
