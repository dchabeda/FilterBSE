#include "mod_kernel.h"

void mod_kernel(
    double*         psi_qp,
    zomplex**       direct,
    zomplex**       exchange,
    double*         pot_bare,
    double*         pot_screened,
    index_st*       ist,
    par_st*         par,
    flag_st*        flag,
    parallel_st*    parallel
    ){
    
    const int          mpir = parallel->mpi_rank;
    
    if (mpir == 0){
        write_separation(stdout, "T");
        printf("\n3.\tCOMPUTING ELEC-HOLE INTERACTION KERNEL | %s\n", get_time());
        write_separation(stdout, "B"); fflush(stdout);
    } 
    
    if (mpir == 0) printf("\nThe number of electron-hole pairs in the exciton basis = %ld\n", ist->n_xton);

    ALLOCATE(direct,   ist->n_xton * ist->n_xton, "direct");
    ALLOCATE(exchange, ist->n_xton * ist->n_xton, "exchange");
     
    if (flag->coulombDone){
        printf("Loading in Coulomb matrix elements from files\n"); 
        fflush(0);

        load_coulomb_mat(*direct,   "direct.dat", ist);
        load_coulomb_mat(*exchange, "exchange.dat", ist);

    } else{
        if (1 == flag->isComplex) {
            printf("Computing complex e-h kernel\n"); 
            fflush(0);

            calc_eh_kernel_cplx(
                psi_qp, pot_bare, pot_screened, *direct, *exchange,
                ist, par, flag, parallel
            );
        }
        // else if (0 == flag->isComplex){
        //     calc_eh_kernel_real((zomplex*) psi_qp, pot_bare, pot_screened, pot_hartree, bsmat, direct, exchange, h0mat, eig_vals, &ist, &par, &flag, planfw, planbw, fftwpsi);
        // }
    }
    
    
    if (flag->calcCoulombOnly == 1){
        printf("Exitng program after computing Coulomb matrix elements | %s\n", get_time());
        exit(0);
    }

    return;
}