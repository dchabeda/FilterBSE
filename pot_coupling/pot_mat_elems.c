#include "pot_coupling.h"

void calc_pot_mat_elems(double *psitot, double *pot_local, nlc_st *nlc, long *nl, double *eval, par_st *par, index_st *ist, flag_st *flag, long n_NL_gridpts, FILE *fij, FILE *fab){
  FILE *pf; 
  long i, a; 
  zomplex tmp;
  
  printf("Starting calculation of distorted geometry matrix elements\n"); fflush(0);
  
  // Main computional work of function performed here -
  // We loop over all hole-hole (i-j) and electron-electron (a-b) pairs

  // compute V_ij = <j|V_loc + V_NL + V_SO|i> 
  #pragma omp parallel for private(i)
  for (i = 0; i < ist->n_holes; i++){
    long jgridup, jgriddn, jgridup_real, jgridup_imag, jgriddn_real, jgriddn_imag;
    long istate, jstate, j;
    zomplex g;
    zomplex *phi, *psi;

    istate = i * ist->complex_idx * ist->nspinngrid;

    // Allocate memory for the v|psi> object
    if ((phi = (zomplex *)calloc(ist->nspinngrid, sizeof(zomplex))) == NULL){
      fprintf(stderr, "\nOUT OF MEMORY: phi in calc_pot_overlap\n\n"); exit(EXIT_FAILURE);
    }
    if ((psi = (zomplex *)calloc(ist->nspinngrid, sizeof(zomplex))) == NULL){
      fprintf(stderr, "\nOUT OF MEMORY: psi in calc_pot_overlap\n\n"); exit(EXIT_FAILURE);
    }

    
    // Copy the current hole wavefunction into |psi>
    memcpy(&phi[0], &psitot[istate], ist->complex_idx * ist->nspinngrid*sizeof(psitot[0]));
    
    // Compute the action of potential operator: |phi> = Vloc + Vnonlocal + Vso|psi>
    potential(phi, psi, pot_local, nlc, nl, ist, par, flag, n_NL_gridpts);

    for (j = 0; j < ist->n_holes; j++){
        jstate = ist->complex_idx* j *ist->nspinngrid;

        g.re = g.im = 0.0;
        if (1 == flag->isComplex){
          for (jgridup = 0; jgridup < ist->ngrid; jgridup++) {
            jgridup_real = ist->complex_idx * jgridup;
            jgridup_imag = ist->complex_idx * jgridup + 1;
            jgriddn      = jgridup + ist->ngrid;
            jgriddn_real = ist->complex_idx * (jgridup+ist->ngrid);
            jgriddn_imag = ist->complex_idx * (jgridup+ist->ngrid) + 1;

            tmp.re = 
                psitot[jstate+jgridup_real] * phi[jgridup].re + psitot[jstate+jgridup_imag] * phi[jgridup].im
              + psitot[jstate+jgriddn_real] * phi[jgriddn].re + psitot[jstate+jgriddn_imag] * phi[jgriddn].im;
            
            tmp.im = 
                psitot[jstate+jgridup_real] * phi[jgridup].im  - psitot[jstate+jgridup_imag] * phi[jgridup].re 
              + psitot[jstate+jgriddn_real] * phi[jgriddn].im - psitot[jstate+jgriddn_imag] * phi[jgriddn].re;
            
            g.re += tmp.re; 
            g.im += tmp.im;
          }
        } else {
          for (jgridup = 0; jgridup < ist->ngrid; jgridup++) {
            tmp.re = psitot[jstate+jgridup] * phi[jgridup].re;
            tmp.im = 0.0;
            g.re += tmp.re; 
            g.im += tmp.im;
          }
        }
        g.re *= par->dv;
        g.im *= par->dv;
      
        // Print out the potential matrix elements at the distorted geometry
        // Print out the indices counting from the homo so that 0 is always homo, 1 is always homo-1, etc.
        fprintf(fij,"%d %.2lg %2ld % .9f %2ld % .9f % .12g % .12g\n", par->mode_idx, par->Q_alpha, labs(i - (ist->n_holes - 1)), eval[i], labs(j - (ist->n_holes - 1)) , eval[j], g.re, g.im);
    }
    free(phi);
    free(psi);
  }
  
  // compute V_ab = <b|V_loc + V_NL + V_SO|a> 
  #pragma omp parallel for private(a)
  for (a = ist->lumo_idx; a < ist->lumo_idx + ist->n_elecs; a++){
    long jgridup, jgriddn, jgridup_real, jgridup_imag, jgriddn_real, jgriddn_imag;
    long astate, bstate, b;
    zomplex g;
    zomplex *phi, *psi;
    
    astate = a * ist->complex_idx * ist->nspinngrid;

    // Allocate memory for the v|psi> object
    if ((phi = (zomplex *)calloc(ist->nspinngrid, sizeof(zomplex))) == NULL){
      fprintf(stderr, "\nOUT OF MEMORY: phi in calc_pot_overlap\n\n"); exit(EXIT_FAILURE);
    }
    if ((psi = (zomplex *)calloc(ist->nspinngrid, sizeof(zomplex))) == NULL){
      fprintf(stderr, "\nOUT OF MEMORY: psi in calc_pot_overlap\n\n"); exit(EXIT_FAILURE);
    }

    
    // Copy the current hole wavefunction into |psi>
    memcpy(&phi[0], &psitot[astate], ist->complex_idx*ist->nspinngrid*sizeof(psitot[0]));
    
    // Compute the action of potential operator: |phi> = V_loc + V_NL + V_SO|psi>
    potential(phi, psi, pot_local, nlc, nl, ist, par, flag, n_NL_gridpts);
    
    for (b = ist->lumo_idx; b < ist->lumo_idx + ist->n_elecs; b++){
        bstate = ist->complex_idx* b *ist->nspinngrid;

        g.re = g.im = 0.0;
        if (1 == flag->isComplex){
          for (jgridup = 0; jgridup < ist->ngrid; jgridup++) {
              jgridup_real = ist->complex_idx * jgridup;
              jgridup_imag = ist->complex_idx * jgridup + 1;
              jgriddn      = jgridup+ist->ngrid;
              jgriddn_real = ist->complex_idx * (jgridup+ist->ngrid);
              jgriddn_imag = ist->complex_idx * (jgridup+ist->ngrid) + 1;

              tmp.re = 
                  psitot[bstate+jgridup_real] * phi[jgridup].re + psitot[bstate+jgridup_imag] * phi[jgridup].im
                + psitot[bstate+jgriddn_real] * phi[jgriddn].re + psitot[bstate+jgriddn_imag] * phi[jgriddn].im;

              
              tmp.im = 
                psitot[bstate+jgridup_real] * phi[jgridup].im - psitot[bstate+jgridup_imag] * phi[jgridup].re + 
                psitot[bstate+jgriddn_real] * phi[jgriddn].im - psitot[bstate+jgriddn_imag] * phi[jgriddn].re;
            
            g.re += tmp.re; 
            g.im += tmp.im;
          }
        } else if (0 == flag->isComplex){
          for (jgridup = 0; jgridup < ist->ngrid; jgridup++) {
            tmp.re = psitot[bstate+jgridup] * phi[jgridup].re;
            tmp.im = 0.0;
            g.re += tmp.re; 
            g.im += tmp.im;
          }
        }
        g.re *= par->dv;
        g.im *= par->dv;
        // Print out the potential matrix elements at the distorted geometry
        fprintf(fab,"%d %.2lg %2ld % .9f %2ld % .9f % .12g % .12g\n", par->mode_idx, par->Q_alpha, a-ist->n_elecs, eval[a], b-ist->n_elecs, eval[b], g.re, g.im);
    }
    free(psi);
    free(phi);
  }
  

  printf("Finished calculation of distorted geometry matrix elements\n"); fflush(0);
  return;
}
