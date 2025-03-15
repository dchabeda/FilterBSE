#include "spin.h"

void qp_spin_frac(
    double*          psi_qp,
    double*          eig_vals,
    index_st*        ist,
    par_st*          par,
    flag_st*         flag,
    parallel_st*     parallel
    ){

    /************************************************************/
	/*******************  DECLARE VARIABLES   *******************/
	/************************************************************/

    FILE*                    pf   = fopen("qp_spins.dat", "w");

    const long               stlen = ist->complex_idx * ist->nspinngrid;

    unsigned long            i;
    unsigned long            i_st;
    unsigned long            jgr; // jgrid real
    unsigned long            jgi; // jgrid imag
    unsigned long            jsgr; // jspingrid real
    unsigned long            jsgi; // jspingrid imag

    double                   psi_iur;
    double                   psi_iui;
    double                   psi_idr;
    double                   psi_idi;
    double                   perUp = 0;
    double                   perDn = 0;

    /************************************************************/
	/*******************  COMPUTE SPIN FRAC   *******************/
	/************************************************************/

    printf("\nComputing spin composition of quasiparticle spinors | %s\n", get_time()); fflush(stdout);
    
    for (i  = 0; i < ist->n_qp; i++){
        i_st = i * stlen;
        
        fprintf(pf, "\nStats on state%d: (E=%lg)\n", i, eig_vals[i]);

        for (jgr = 0; jgr < ist->ngrid; jgr++){
            jgr = ist->complex_idx * jgr;
            jgi = jgr + 1;
            jsgr = jgr + ist->ngrid;
            jsgi = jsgr + 1;

            psi_iur = psi_qp[i_st + jgr];
            psi_iui = psi_qp[i_st + jgi];
            psi_idr = psi_qp[i_st + jsgr];
            psi_idi = psi_qp[i_st + jsgi];
            
            perUp += sqr(psi_iur) + sqr(psi_iui);
            perDn += sqr(psi_idr) + sqr(psi_idi);
        }
        
        fprintf(pf, " Spin up fraction: %f\n", perUp * par->dv); fflush(pf);
        fprintf(pf, " Spin dn fraction: %f\n", perDn * par->dv); fflush(pf);
    }
    fclose(pf);

    return;
}