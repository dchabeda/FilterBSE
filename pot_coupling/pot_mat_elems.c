#include "pot_coupling.h"

void calc_pot_mat_elems(double *psitot, double *pot_local_equil, nlc_st *nlc_equil, long *nl_equil, double *pot_local, nlc_st *nlc, long *nl, double *eval, par_st *par, index_st *ist, flag_st *flag, long n_NL_gridpts_equil, long n_NL_gridpts){
  FILE *pf, *pf1; 
  long i; 
  zomplex tmp;
  zomplex *phi, *psi;
  
  // Output will be written to this file
  pf = fopen("pot_overlap_equil.dat" , "w"); 
  pf1 = fopen("pot_overlap_dist.dat" , "w"); 

  // Allocate memory for the v|psi> object
  if ((phi = (zomplex *)calloc(ist->nspinngrid, sizeof(zomplex))) == NULL){
    fprintf(stderr, "\nOUT OF MEMORY: phi in calc_pot_overlap\n\n"); exit(EXIT_FAILURE);
  }
  if ((psi = (zomplex *)calloc(ist->nspinngrid, sizeof(zomplex))) == NULL){
    fprintf(stderr, "\nOUT OF MEMORY: psi in calc_pot_overlap\n\n"); exit(EXIT_FAILURE);
  }

  // Main computional work of function performed here -
  // Technically, we only need to loop over all electron-electrom (a-b)
  // and hole-hole (i-j) pairs, but for verification purposes all quasiparticles are included


  printf("Starting calculation of equilibrium geometry matrix elements\n\n"); fflush(0);
  // First compute matrix elements of the equilibrium geometry potential
  // compute V_ia = <a|Vloc_equil + Vnonlocal_equil + Vso_equil|i> 
#pragma omp parallel for private(i)
  for (i = 0; i < ist->n_qp; i++){
    long jgridup, jgriddn, jgridup_real, jgridup_imag, jgriddn_real, jgriddn_imag;
    long istate, astate, a;
    zomplex g;

    // printf("i = %ld\n", i); fflush(0);
    istate = ist->complex_idx* i *ist->nspinngrid;
    // Copy the current hole wavefunction into |psi>
    
    //printf("Before memcpy\n"); fflush(0);
    memcpy(&psi[0], &psitot[istate], ist->complex_idx*ist->nspinngrid*sizeof(psitot[0]));
    //printf("memcpy successful\n"); fflush(0);
    // Initialize the vector that will become |phi> = (Vloc + Vnonlocal + Vso)|psi>, as 0.0s
    // This next line is technically not necessary because calloc initializes memory to 0
    for (jgridup = 0; jgridup < ist->nspinngrid; jgridup++){
      phi[jgridup].re = phi[jgridup].im = 0.0;
    }
    //printf("populated phi with zeros\n"); fflush(0);
    
    // Compute the potential |phi> = Vloc + Vnonlocal + Vso|psi>
    potential(phi, psi, pot_local_equil, nlc_equil, nl_equil, ist, par, flag, n_NL_gridpts_equil);
    for (a = 0; a < ist->n_qp; a++){
        //printf("a = %ld\n", a); fflush(0);
        astate = ist->complex_idx* a *ist->nspinngrid;
        // Make sure g is zero'd to begin with
        g.re = g.im = 0.0;
        for (jgridup = 0; jgridup < ist->ngrid; jgridup++) {
            jgridup_real = ist->complex_idx * jgridup;
            jgridup_imag = ist->complex_idx * jgridup + 1;
            jgriddn = jgridup+ist->ngrid;
            jgriddn_real = ist->complex_idx *(jgridup+ist->ngrid);
            jgriddn_imag = ist->complex_idx *(jgridup+ist->ngrid) + 1;
            //printf("calculate tmp\n"); fflush(0);
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
        // Print out the potential matrix elements at the equilibrium geometry
        fprintf(pf,"%2ld % .9f %2ld % .9f % .12g % .12g\n", i, eval[i], a, eval[a], g.re, g.im);
    }
  }
  fclose(pf);

  // Next compute matrix elements of the distorted geometry potential
  // (Note: all wavefunctions were computed at the equilibrium geometry)
  printf("Starting calculation of distorted geometry matrix elements"); fflush(0);
  // compute V_ia = <a|Vloc_dist + Vnonlocal_dist + Vso_dist|i> 
#pragma omp parallel for private(i)
  for (i = 0; i < ist->n_qp; i++){
    long jgridup, jgriddn, jgridup_real, jgridup_imag, jgriddn_real, jgriddn_imag;
    long istate, astate, a;
    zomplex g;

    istate = ist->complex_idx*i*ist->nspinngrid;
    // Copy the current hole wavefunction into |psi>
    memcpy(&psi[0], &psitot[ist->complex_idx*i*ist->nspinngrid], ist->nspinngrid*sizeof(psitot[0]));
    // Initialize the vector that will become |phi> = Vloc + Vnonlocal + Vso|psi>, as 0.0s
    // This next line is technically not necessaRy because calloc initializes memoRy to 0
    for (jgridup = 0; jgridup < ist->nspinngrid; jgridup++){
      phi[jgridup].re = phi[jgridup].im = 0.0;
    }
    
    // Compute the potential |phi> = Vloc + Vnonlocal + Vso|psi>
    potential(phi, psi, pot_local, nlc, nl, ist, par, flag, n_NL_gridpts);
    
    for (a = 0; a < ist->n_qp; a++){
        astate = ist->complex_idx* a *ist->nspinngrid;
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
            g.re += tmp.re; 
            g.im += tmp.im;
        }
        g.re *= par->dv;
        g.im *= par->dv;
        // Print out the potential matrix elements at the distorted geometry
        fprintf(pf1,"%2ld % .9f %2ld % .9f % .12g % .12g\n", i, eval[i], a, eval[a], g.re, g.im);
    }
  }
  fclose(pf1);
  
  return;
}
