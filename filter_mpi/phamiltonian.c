#include "fd.h"
#include <complex.h>
/*****************************************************************************/
// File contains the main functions for calcuating the action of the hamiltonian on 
// a wavefunction including non-local spin-orbit contributions 
// parallelized with OpenMP

/*****************************************************************************/

void p_hamiltonian(zomplex *psi_out, zomplex *psi_tmp, double *pot_local, nlc_st *nlc, long *nl, double *ksqr,
  index_st *ist, par_st *par, flag_st *flag, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi, int ham_threads){
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
  for (jspin = 0; jspin < ist->nspin; jspin++){
      kinetic(&psi_out[jspin*ist->ngrid], ksqr, planfw, planbw, fftwpsi, ist); //spin up/down
  } 
  // write_state_dat(psi_out, ist->nspinngrid, "psi_out_pkinetic.dat");
  // Calculate the action of the potential on the wavefunction: |psi_out> = V|psi_tmp>
  
  p_potential(psi_out, psi_tmp, pot_local, nlc, nl, ist, par, flag, ham_threads);

  return;
}

/*****************************************************************************/

// This calculates the total action of Vloc + Vnonloc + Vso|psi_tmp>
void p_potential(zomplex *psi_out, zomplex *psi_tmp, double *pot_local, nlc_st *nlc, long *nl, index_st *ist,
  par_st *par, flag_st *flag, int ham_threads){
  /*******************************************************************
  * This function calculates |psi_out> = [Vloc+V_SO+V_NL]|psi_tmp>   *
  * inputs:                                                          *
  *  [psi_out] nspinngrid-long arr where V|psi_tmp> is stored        *
  *  [psi_tmp] nspinngrid-long arr holding input wavefunc            *
  *  [pot_local] ngrid-long arr holding the value of the local pot   *
  *  [nlc] nlc struct holding values for computing SO and NL pots    *
  *  [nl] natom-long arr holding the number of NL gridpts per atom   *
  *  [ist] ptr to counters, indices, and lengths                     *
  *  [par] ptr to par_st holding VBmin, VBmax... params              *
  *  [flag] ptr to flag_st holding job flags                         *
  * outputs: void                                                    *
  ********************************************************************/

  long j, jtmp, jspin;
  

  if(flag->SO==1){
    // Calculate |psi_out> = V_SO|psi_tmp>
    p_spin_orbit_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par, ham_threads);
    //write_state_dat(psi_out, ist->nspinngrid, "psi_out_pSO.dat");
  }
  if (flag->NL == 1){
    // Calculate |psi_out> += V_NL|psi_tmp>
    p_nonlocal_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par, ham_threads);
    //write_state_dat(psi_out, ist->nspinngrid, "psi_out_pNL.dat");
  }
  
  // Calculate the action of the local potential energy part of the Hamiltonian on psi_tmp
  // No openMP atomics necessary for this because all grid indices are unique
  // (no two threads reads/writes the same j/jtmp ever)
  if (1 == flag->useSpinors){
    for (j = 0; j < ist->ngrid; j++) {
      psi_out[j].re += (pot_local[j] * psi_tmp[j].re);
      psi_out[j].im += (pot_local[j] * psi_tmp[j].im);
      // handle spin down component
      jtmp = ist->ngrid + j;
      psi_out[jtmp].re += (pot_local[j] * psi_tmp[jtmp].re);
      psi_out[jtmp].im += (pot_local[j] * psi_tmp[jtmp].im);
    }
  }
  else if (0 == flag->useSpinors){
    for (j = 0; j < ist->ngrid; j++) {
      psi_out[j].re += (pot_local[j] * psi_tmp[j].re);
    }
  } else {
    printf("ERROR: unrecognized value for flag->useSpinors: %d\n", flag->useSpinors);
    fprintf(stderr, "ERROR: unrecognized value for flag->useSpinors: %d\n", flag->useSpinors);
    exit(EXIT_FAILURE);
  }
  
  //write_state_dat(psi_out, ist->nspinngrid, "psi_out_ploc.dat");
  return;
}


void p_spin_orbit_proj_pot(zomplex *psi_out, zomplex *psi_tmp, nlc_st *nlc, long *nl, index_st *ist, par_st *par, int ham_threads){
  /*******************************************************************
  * This function calculates the action of the spin-orbit nonlocal   *
  * potential using separable radial projector functions             *
  * inputs:                                                          *
  *  [psi_out] nspinngrid-long arr where V_SO|psi_tmp> will be stored *
  *  [psi_tmp] nspinngrid-long arr holding input wavefunc             *
  *  [nlc] nlc struct holding values for computing SO and NL pots    *
  *  [nl] natom-long arr holding the number of NL gridpts per atom   *
  *  [ist] ptr to counters, indices, and lengths                     *
  *  [par] ptr to par_st holding VBmin, VBmax... params              *
  * outputs: void                                                    *
  ********************************************************************/

  long jatom, jatom_offset, scratch_offset;
  long NL_gridpt, r_idx, r, r_p;
  int iproj, s, s_p, m, m_p;
  int spin_arr[ist->n_j_ang_mom], m_arr[ist->n_j_ang_mom];
  int j, j_p, jtot;
  int t_id;
  zomplex proj, LS_loc, PLS;
  zomplex *LS;

  double psi_re, psi_im;
  double y1_re, y1_im;
  double nlcproj;

  for (s = 0; s < ist->n_s_ang_mom; s++){
    for (m = 0; m < ist->n_l_ang_mom; m++){
      j = s*3 + m;
      spin_arr[j] = s;
      m_arr[j] = m; 
    }
  }

  // Define the LS matrix
  LS = (zomplex*) calloc(ist->n_j_ang_mom * ist->n_j_ang_mom, sizeof(zomplex));
  def_LS(LS, ist, par);


  omp_set_num_threads(ham_threads);
  #pragma omp parallel for private(jatom, jatom_offset, NL_gridpt, r_idx, r, r_p, iproj, proj, s, s_p, m, m_p, j, j_p, jtot, psi_re, psi_im, y1_re, y1_im, nlcproj, LS_loc, PLS)
  for (jatom = 0; jatom < ist->n_NL_atoms; jatom++){
    long jatom_offset = jatom * ist->n_NL_gridpts;
    // Compute equation 2.81 of Daniel Weinberg dissertation
    // Action of spin-orbit operator on a real space wavefunctions
    for ( iproj = 0; iproj < ist->nproj; iproj++){
      for ( j_p = 0; j_p < ist->n_j_ang_mom; j_p++){
        s_p = spin_arr[j_p];
        m_p = m_arr[j_p];

        proj.re = proj.im = 0.00;
        // Compute the projection of the real space wavefunction onto the basis of |lmr\sigma> 
        // where the radial variable, r, is not computed on a grid but actually are smooth radial 
        // functions
        
        for ( NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++){
          r_idx = jatom_offset + NL_gridpt;
          r = nlc[r_idx].jxyz + (ist->ngrid) * s_p;
          psi_re = psi_tmp[r].re;
          psi_im = psi_tmp[r].im;
          
          y1_re = nlc[r_idx].y1[m_p].re;
          y1_im = nlc[r_idx].y1[m_p].im;
          nlcproj = nlc[r_idx].proj[iproj];
          //weird signs b/c of Y_{lm}^*
          // Calculate the integral in eq 2.81
          proj.re += psi_re * y1_re * nlcproj + psi_im * y1_im * nlcproj;
          proj.im += psi_im * y1_re * nlcproj - psi_re * y1_im * nlcproj;
        }
        
        proj.re *= par->dv;
        proj.im *= par->dv;
          
        for (j = 0; j < ist->n_j_ang_mom; j++){
          //get L_{m,m'}\cdot S_{s,s'}*P_{n,m,s} = PLS_{n,m,m',s,s'}
          s = spin_arr[j];
          m = m_arr[j];
          jtot = j_p * ist->n_j_ang_mom + j;

          LS_loc.re = LS[jtot].re;
          LS_loc.im = LS[jtot].im;
          if (LS_loc.re == 0.0 && LS_loc.im == 0.0){
            continue;
          }
          PLS.re = LS_loc.re * proj.re - LS_loc.im * proj.im;
          PLS.im = LS_loc.re * proj.im + LS_loc.im * proj.re;
          
          for (NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++){
            r_idx = jatom_offset + NL_gridpt;
            r_p = nlc[r_idx].jxyz + (ist->ngrid)*s;
            
            nlcproj = nlc[r_idx].proj[iproj];
            y1_re = nlc[r_idx].y1[m].re;
            y1_im = nlc[r_idx].y1[m].im;

            #pragma omp atomic
            psi_out[r_p].re += nlcproj * y1_re * PLS.re - nlcproj * y1_im * PLS.im;
            #pragma omp atomic
            psi_out[r_p].im += nlcproj * y1_re * PLS.im + nlcproj * y1_im * PLS.re;
          }
        }
      }
    }
  }

  return;
}

/*****************************************************************************/

void p_nonlocal_proj_pot(zomplex *psi_out, zomplex *psi_tmp, nlc_st *nlc, long *nl, index_st *ist, par_st *par, int ham_threads){
  /*******************************************************************
  * This function calculates the action of the angular nonlocal      *
  * potential using separable radial projector functions             *
  * inputs:                                                          *
  *  [psi_out] nspinngrid-long arr where V_NL|psi_tmp> will be stored *
  *  [psi_tmp] nspinngrid-long arr holding input wavefunc             *
  *  [nlc] nlc struct holding values for computing SO and NL pots    *
  *  [nl] natom-long arr holding the number of NL gridpts per atom   *
  *  [ist] ptr to counters, indices, and lengths                     *
  *  [par] ptr to par_st holding VBmin, VBmax... params              *
  * outputs: void                                                    *
  ********************************************************************/

  long jatom, jatom_offset;
  long NL_gridpt, r_idx, r, r_p;
  int iproj, s, m, j, t_id;
  int spin_arr[ist->n_j_ang_mom], m_arr[ist->n_j_ang_mom];
  zomplex proj;
  double proj_re;
  double proj_im;

  double psi_re, psi_im;
  double y1_re, y1_im;
  double NL_proj;
  
  for (s = 0; s < ist->n_s_ang_mom; s++){
    for (m = 0; m < ist->n_l_ang_mom; m++){
      j = s*ist->n_l_ang_mom + m;
      spin_arr[j] = s;
      m_arr[j] = m; 
    }
  }

  omp_set_num_threads(ham_threads);
  #pragma omp parallel for private(jatom, jatom_offset, NL_gridpt, r_idx, r, r_p, iproj, s, m, j, psi_re, psi_im, y1_re, y1_im, NL_proj)
  for (jatom = 0; jatom < ist->n_NL_atoms; jatom++) {
    
    jatom_offset = jatom * ist->n_NL_gridpts;

    for (iproj = 0; iproj < ist->nproj; iproj++) {
      for (j = 0; j < ist->n_j_ang_mom; j++) {
        s = spin_arr[j];
        m = m_arr[j];

        proj.re = proj.im = 0.0;

        // First loop over NL_gridpt
        for (NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++) {
          r_idx = jatom_offset + NL_gridpt;
          r = nlc[r_idx].jxyz + ist->ngrid * s;

          psi_re = psi_tmp[r].re;
          psi_im = psi_tmp[r].im;
          y1_re = nlc[r_idx].y1[m].re;
          y1_im = nlc[r_idx].y1[m].im;
          
          NL_proj = nlc[r_idx].NL_proj[iproj];

          proj.re += psi_re * y1_re * NL_proj + psi_im * y1_im * NL_proj;
          proj.im += psi_im * y1_re * NL_proj - psi_re * y1_im * NL_proj;
        }
        
        // Apply scaling factors
        proj_re *= nlc[jatom_offset].NL_proj_sign[iproj] * par->dv;
        proj_im *= nlc[jatom_offset].NL_proj_sign[iproj] * par->dv;

        // Second loop over NL_gridpt
        for (NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++) {
          r_idx = jatom_offset + NL_gridpt;
          r_p = nlc[r_idx].jxyz + ist->ngrid * s;

          y1_re = nlc[r_idx].y1[m].re;
          y1_im = nlc[r_idx].y1[m].im;
          NL_proj = nlc[r_idx].NL_proj[iproj];
          
          // Update thread-local psi_scratch
          #pragma omp atomic
          psi_out[r_p].re += (NL_proj * (y1_re * proj_re - y1_im * proj_im));
          #pragma omp atomic
          psi_out[r_p].im += NL_proj * (y1_re * proj_im + y1_im * proj_re);
        }
      }
    }
  }

  return;
}

/*****************************************************************************/

