#include "fd.h"

void calc_pot_overlap(double *psitot, double *pot_local, nlc_st *nlc, long *nl, double *eval, par_st *par,index_st *ist, flag_st *flag){
  FILE *pf, *pf1; 
  long i, a, istate, astate, jgridup, jgriddn, jgridup_real, jgridup_imag, jgriddn_real, jgriddn_imag; 
  zomplex tmp;
  zomplex g, *phi, *psi;
  
  // Output will be written to this file
  pf = fopen("pot_overlap_hole.dat" , "w"); 
  pf1 = fopen("pot_overlap_elec.dat" , "w"); 

  // Allocate memoRy for the v|psi> object
  if ((phi = (zomplex *)calloc(ist->nspinngrid, sizeof(zomplex))) == NULL){
    fprintf(stderr, "\nOUT OF MEMORy: phi in calc_pot_overlap\n\n"); exit(EXIT_FAILURE);
  }
  if ((psi = (zomplex *)calloc(ist->nspinngrid, sizeof(zomplex))) == NULL){
    fprintf(stderr, "\nOUT OF MEMORy: psi in calc_pot_overlap\n\n"); exit(EXIT_FAILURE);
  }

  // Main computional work of function performed here - must loop over all electron-electrom (a-a)
  // and hole-hole pairs
  
  // First the hole-hole couplings
  // compute V_ij = <i|Vloc + Vnonlocal + Vso|j> 
  for (i = 0; i < ist->total_homo; i++){
    istate = ist->complex_idx*i*ist->nspinngrid;
    // Copy the current hole wavefunction into |psi>
    memcpy(&psi[0], &psitot[ist->complex_idx*i*ist->nspinngrid], ist->complex_idx*ist->nspinngrid*sizeof(psitot[0]));
    // Initialize the vector that will become |phi> = Vloc + Vnonlocal + Vso|psi>, as 0.0s
    // This next line is technically not necessaRy because calloc initializes memory to 0
    for (jgridup = 0; jgridup < ist->nspinngrid; jgridup++){
      phi[jgridup].re = phi[jgridup].im = 0.0;
    }
    
    // Compute the potential |phi> = Vloc + Vnonlocal + Vso|psi>
    potential(phi, psi, pot_local, nlc, nl, ist, par, flag);
    
    //for (i = 0; i < ist->total_homo; i++){
    // Make sure g is zero'd to begin with
    g.re = g.im = 0.0;
    for (jgridup = 0; jgridup < ist->ngrid; jgridup++) {
      jgridup_real = ist->complex_idx * jgridup;
      jgridup_imag = ist->complex_idx * jgridup + 1;
      jgriddn = jgridup+ist->ngrid;
      jgriddn_real = ist->complex_idx *(jgridup+ist->ngrid);
      jgriddn_imag = ist->complex_idx *(jgridup+ist->ngrid) + 1;

      tmp.re =  (psitot[istate+jgridup_real] * phi[jgridup].re + psitot[istate+jgridup_imag] * phi[jgridup].im
                  +psitot[istate+jgriddn_real] * phi[jgriddn].re + psitot[istate+jgriddn_imag] * phi[jgriddn].im );
      if (1 == flag->isComplex){
        tmp.im =  (-psitot[istate+jgridup_imag] * phi[jgridup].re + psitot[istate+jgridup_real] * phi[jgridup].im
                    -psitot[istate+jgriddn_imag] * phi[jgriddn].re + psitot[istate+jgriddn_real] * phi[jgriddn].im) ;
      } else if (0 == flag->isComplex){
        tmp.im =  0.0;
      }
      g.re += par->dv*tmp.re; g.im +=  par->dv*tmp.im;
    }

    fprintf(pf,"%ld %.12f %.12f %.12f\n", i, eval[i], g.re, g.im);
    //}
  }
  fclose(pf);

  // Now the electron-electron couplings
  
  // compute V_ab = <a|Vloc + Vnonlocal + Vso|b> 
  for (a = ist->lumo_idx; a < ist->lumo_idx + ist->total_lumo; a++){
    astate = ist->complex_idx*a*ist->nspinngrid;
    // Copy the current hole wavefunction into |psi>
    memcpy(&psi[0], &psitot[ist->complex_idx*a*ist->nspinngrid], ist->nspinngrid*sizeof(psitot[0]));
    // Initialize the vector that will become |phi> = Vloc + Vnonlocal + Vso|psi>, as 0.0s
    // This next line is technically not necessaRy because calloc initializes memoRy to 0
    for (jgridup = 0; jgridup < ist->nspinngrid; jgridup++){
      phi[jgridup].re = phi[jgridup].im = 0.0;
    }
    
    // Compute the potential |phi> = Vloc + Vnonlocal + Vso|psi>
    potential(phi, psi, pot_local, nlc, nl, ist, par, flag);
    
    //for (a = ist->lumo_idx; a < ist->lumo_idx + ist->total_lumo; a++){
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
      g.re += par->dv*tmp.re; g.im += par->dv*tmp.im;
    }
    
    fprintf(pf1,"%ld %.12f %.12f %.12f\n", a, eval[a], g.re, g.im);
    //}
  }
  fclose(pf1);
  
  return;
}
