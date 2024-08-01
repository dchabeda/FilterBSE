#include "pot_coupling.h"

void calc_pot_mat_elems(double *psitot, double *pot_local_equil, nlc_st *nlc_equil, long *nl_equil, double *pot_local, nlc_st *nlc, long *nl, double *eval, par_st *par,index_st *ist, flag_st *flag){
  FILE *pf, *pf1; 
  long i, a; 
  zomplex tmp;
  zomplex g, *phi, *psi;
  
  // Output will be written to this file
  pf = fopen("pot_overlap_equil.dat" , "w"); 
  pf1 = fopen("pot_overlap_dist.dat" , "w"); 

  // Allocate memory for the v|psi> object
  if ((phi = (zomplex *)calloc(ist->nspinngrid, sizeof(zomplex))) == NULL){
    fprintf(stderr, "\nOUT OF MEMORy: phi in calc_pot_overlap\n\n"); exit(EXIT_FAILURE);
  }
  if ((psi = (zomplex *)calloc(ist->nspinngrid, sizeof(zomplex))) == NULL){
    fprintf(stderr, "\nOUT OF MEMORy: psi in calc_pot_overlap\n\n"); exit(EXIT_FAILURE);
  }

  // Main computional work of function performed here - must loop over all electron-electrom (a-a)
  // and hole-hole pairs
  printf("Starting calculation of equilibrium geometry matrix elements"); fflush(0);
  // First the hole-hole couplings
  // compute V_ia = <i|Vloc + Vnonlocal + Vso|j> 
#pragma omp parallel for private(i)
  for (i = 0; i < ist->mn_states_tot; i++){
    long istate, astate, jgridup, jgriddn, jgridup_real, jgridup_imag, jgriddn_real, jgriddn_imag;
    
    istate = ist->complex_idx* i *ist->nspinngrid;
    // Copy the current hole wavefunction into |psi>
    memcpy(&psi[0], &psitot[ist->complex_idx*i*ist->nspinngrid], ist->complex_idx*ist->nspinngrid*sizeof(psitot[0]));
    // Initialize the vector that will become |phi> = Vloc + Vnonlocal + Vso|psi>, as 0.0s
    // This next line is technically not necessaRy because calloc initializes memory to 0
    for (jgridup = 0; jgridup < ist->nspinngrid; jgridup++){
      phi[jgridup].re = phi[jgridup].im = 0.0;
    }
    
    // Compute the potential |phi> = Vloc + Vnonlocal + Vso|psi>
    potential(phi, psi, pot_local_equil, nlc_equil, nl_equil, ist, par, flag);
    
    for (a = 0; a < ist->mn_states_tot; a++){
        astate = ist->complex_idx* a *ist->nspinngrid;
        // Make sure g is zero'd to begin with
        g.re = g.im = 0.0;
        for (jgridup = 0; jgridup < ist->ngrid; jgridup++) {
        jgridup_real = ist->complex_idx * jgridup;
        jgridup_imag = ist->complex_idx * jgridup + 1;
        jgriddn = jgridup+ist->ngrid;
        jgriddn_real = ist->complex_idx *(jgridup+ist->ngrid);
        jgriddn_imag = ist->complex_idx *(jgridup+ist->ngrid) + 1;

        tmp.re =  (psitot[astate+jgridup_real] * phi[jgridup].re + psitot[astate+jgridup_imag] * phi[jgridup].im
                    +psitot[astate+jgriddn_real] * phi[jgriddn].re + psitot[astate+jgriddn_imag] * phi[jgriddn].im );
        if (1 == flag->isComplex){
            tmp.im =  (-psitot[astate+jgridup_imag] * phi[jgridup].re + psitot[astate+jgridup_real] * phi[jgridup].im
                        -psitot[astate+jgriddn_imag] * phi[jgriddn].re + psitot[astate+jgriddn_real] * phi[jgriddn].im) ;
        } else if (0 == flag->isComplex){
            tmp.im =  0.0;
        }
        g.re += tmp.re; g.im += tmp.im;
        }
        g.re *= par->dv;
        g.im *= par->dv;

    fprintf(pf,"%ld %ld %.12f %.12f %.12f\n", i, a, eval[a] - eval[i], g.re, g.im);
    }
  }
  fclose(pf);

  // Now the electron-electron couplings
  printf("Starting calculation of distorted geometry matrix elements"); fflush(0);
  // compute V_ab = <a|Vloc + Vnonlocal + Vso|b> 
#pragma omp parallel for private(i)
  for (i = 0; i < ist->mn_states_tot; i++){
    long istate, astate, jgridup, jgriddn, jgridup_real, jgridup_imag, jgriddn_real, jgriddn_imag;
    
    istate = ist->complex_idx*i*ist->nspinngrid;
    // Copy the current hole wavefunction into |psi>
    memcpy(&psi[0], &psitot[ist->complex_idx*i*ist->nspinngrid], ist->nspinngrid*sizeof(psitot[0]));
    // Initialize the vector that will become |phi> = Vloc + Vnonlocal + Vso|psi>, as 0.0s
    // This next line is technically not necessaRy because calloc initializes memoRy to 0
    for (jgridup = 0; jgridup < ist->nspinngrid; jgridup++){
      phi[jgridup].re = phi[jgridup].im = 0.0;
    }
    
    // Compute the potential |phi> = Vloc + Vnonlocal + Vso|psi>
    potential(phi, psi, pot_local, nlc, nl, ist, par, flag);
    
    for (a = 0; a < ist->mn_states_tot; a++){
        // Make sure g is zero'd to begin with
        g.re = g.im = 0.0;
        for (jgridup = 0; jgridup < ist->ngrid; jgridup++) {
        jgridup_real = ist->complex_idx * jgridup;
        jgridup_imag = ist->complex_idx * jgridup + 1;
        jgriddn = jgridup+ist->ngrid;
        jgriddn_real = ist->complex_idx * (jgridup+ist->ngrid);
        jgriddn_imag = ist->complex_idx * (jgridup+ist->ngrid) + 1;

        tmp.re =  (psitot[astate+jgridup_real] * phi[jgridup].re + psitot[astate+jgridup_imag] * phi[jgridup].im
                    +psitot[astate+jgriddn_real] * phi[jgriddn].re + psitot[astate+jgriddn_imag] * phi[jgriddn].im );
        if (1 == flag->isComplex){
        tmp.im =  (-psitot[astate+jgridup_imag] * phi[jgridup].re + psitot[astate+jgridup_real] * phi[jgridup].im
                    -psitot[astate+jgriddn_imag] * phi[jgriddn].re + psitot[astate+jgriddn_real] * phi[jgriddn].im) ;
        } else if (0 == flag->isComplex){
        tmp.im =  0.0 ;
        }
        g.re += tmp.re; g.im += tmp.im;
        }
        g.re *= par->dv;
        g.im *= par->dv;

        fprintf(pf,"%ld %ld %.12f %.12f %.12f\n", i, a, eval[a] - eval[i], g.re, g.im);
        }
  }
  fclose(pf1);
  
  return;
}
