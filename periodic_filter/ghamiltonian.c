#include "fd.h"


/*****************************************************************************/
void hamiltonian_k(
  zomplex       *psi_out,  zomplex      *psi_tmp,    double *pot_local,    
  vector        *G_vecs,   vector        k,          grid_st *grid, 
  nlc_st        *nlc,      long         *nl,         index_st *ist, 
  par_st        *par,      flag_st      *flag,       fftw_plan_loc planfw, 
  fftw_plan_loc planbw,    fftw_complex *fftwpsi){
  /*******************************************************************
  * This function applies the Hamiltonian onto a state               *
  * inputs:                                                          *
  *  [psi_tmp] nspinngrid-long arr of orig. wavefnc                  *
  *  [psi_out] nspinngrid-long arr to hold |psi_out> = H|psi_tmp>    *
  *  [pot_local] ngrid-long arr holding the value of the local pot   *
  *  [nlc] nlc struct holding values for computing SO and NL pots    *
  *  [nl] natom-long arr holding the number of NL gridpts per atom   *
  *  [ksqr] ngrid-long arr holding the values of k^2 for KE calc     *
  *  [ist] ptr to counters, indices, and lengths                     *
  *  [par] ptr to par_st holding VBmin, VBmax... params              *
  *  [flag] ptr to flag_st holding job flags                         *
  *  [planfw] FFTW3 plan for executing 3D forward DFT                *
  *  [planfw] FFTW3 plan for executing 3D backwards DFT              *
  *  [fftwpsi] location to store outcome of Fourier transform        *
  * outputs: void                                                    *
  ********************************************************************/

  int jspin; 
  // Copy psi_out into psi_tmp
  memcpy(&psi_tmp[0], &psi_out[0], ist->nspinngrid*sizeof(psi_tmp[0]));
  

  // Calculate the action of the kinetic energy part of the Hamiltonian on psi_tmp: |psi_out> = T|psi_tmp>
  kinetic_k(psi_out, G_vecs, k, planfw, planbw, fftwpsi, ist); //spin up
  if (1 == flag->useSpinors){
  kinetic_k(&psi_out[ist->ngrid], G_vecs, k, planfw, planbw, fftwpsi, ist); //spin down
  }

  // If we are computing a periodic filter calculation,
  // then we apply the Hamiltonian onto Bloch states
  // |psi> = e^(ik.r)phi(r), 
  // where phi(r) are our usual filter functions
  if (1 == flag->periodic){
    e_ikr(psi_tmp, k, grid, ist, par, flag, parallel);
  }
  
  // write_state_dat(psi_out, ist->nspinngrid, "psi_out_kinetic.dat");
  // Calculate the action of the potential on the wavefunction: |psi_out> = V|psi_tmp>
  potential(psi_out, psi_tmp, pot_local, nlc, nl, ist, par, flag);

  return;
}

/*****************************************************************************/
// Calculates T|psi_tmp> via FFT and stores result in psi_out for a periodic H

/*****************************************************************************/
// Calculates T(k)|psi_tmp> via FFT and stores result in psi_out 

void kinetic_k(zomplex *psi_out, vector *G_vecs, vector k, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi, index_st *ist){
  /*******************************************************************
  * This function applies the KE operator onto a state               *
  * T|psi> = FT^-1[ |(k+G)|^2 * FT[|psi>] ]                          *
  * inputs:                                                          *
  *  [psi_out] nspinngrid-long arr to hold |psi_out> = T|psi_tmp>    *
  *  [G_vecs] ngrid-long arr holding the G vectors                   *
  *  [k_vecs] ngrid-long arr holding the k vectors                   *
  *  [ist] ptr to counters, indices, and lengths                     *
  *  [planfw] FFTW3 plan for executing 3D forward DFT                *
  *  [planfw] FFTW3 plan for executing 3D backwards DFT              *
  *  [fftwpsi] location to store outcome of Fourier transform        *
  * outputs: void                                                    *
  ********************************************************************/

  long j;
  vector kplusG;
  double kpG2;

  // Copy inputted psi to fftwpsi
  memcpy(&fftwpsi[0], &psi_out[0], ist->ngrid*sizeof(fftwpsi[0]));
  

  // FT from r-space to k-space
  fftw_execute(planfw);

  // Kinetic energy is diagonal in k-space, just multiply fftwpsi by |k+G_j|^2
  for (j = 0; j < ist->ngrid; j++) {
    kplusG = retAddedVectors(k, G_vecs[j]);
    kpG2 = sqr(kplusG.mag);
    fftwpsi[j][0] *= kpG2;
    fftwpsi[j][1] *= kpG2;
  }

  // Inverse FT back to r-space
  fftw_execute(planbw);
  
  // Copy fftwpsi to psi_out to store T|psi_tmp> into |psi_out>
  memcpy(&psi_out[0], &fftwpsi[0], ist->ngrid*sizeof(fftwpsi[0]));

  return;
}

/*****************************************************************************/


