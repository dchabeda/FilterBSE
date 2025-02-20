#include "hamiltonian.h"

/*****************************************************************************/
// File contains the main functions for calcuating the action of the hamiltonian on 
// a wavefunction including non-local spin-orbit contributions 
// parallelized with OpenMP

/*****************************************************************************/

void p_hamiltonian(
  zomplex*        psi_out, 
  zomplex*        psi_tmp, 
  double*         pot_local, 
  zomplex*        projs,
  zomplex*        LS, 
  nlc_st*         nlc, 
  long*           nl, 
  double*         ksqr,
  index_st*       ist, 
  par_st*         par, 
  flag_st*        flag, 
  fftw_plan_loc   planfw, 
  fftw_plan_loc   planbw, 
  fftw_complex*   fftwpsi,
  int             ham_threads
  ){
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
  
  p_potential(psi_out, psi_tmp, pot_local, projs, LS, nlc, nl, ist, par, flag, ham_threads);

  return;
}

/*****************************************************************************/

// This calculates the total action of Vloc + Vnonloc + Vso|psi_tmp>
void p_potential(
  zomplex*        psi_out, 
  zomplex*        psi_tmp, 
  double*         pot_local, 
  zomplex*        projs,
  zomplex*        LS, 
  nlc_st*         nlc, 
  long*           nl, 
  index_st*       ist,
  par_st*         par, 
  flag_st*        flag,
  int             ham_threads
  ){
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

  long j, jtmp;
  

  if(flag->SO==1){
    // Calculate |psi_out> = V_SO|psi_tmp>
    p_spin_orbit_proj_pot(psi_out, psi_tmp, projs, LS, nlc, nl, ist, par, ham_threads);
    //write_state_dat(psi_out, ist->nspinngrid, "psi_out_pSO.dat");
  }
  if (flag->NL == 1){
    // Calculate |psi_out> += V_NL|psi_tmp>
    p_nonlocal_proj_pot(psi_out, psi_tmp, projs, nlc, nl, ist, par, ham_threads);
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
    //#pragma omp parallel for private(j)
    for (j = 0; j < ist->ngrid; j++) {
      // #pragma omp atomic
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


void p_spin_orbit_proj_pot(
  zomplex*        psi_out, 
  zomplex*        psi_tmp,
  zomplex*        projs,
  zomplex*        LS,
  nlc_st*         nlc, 
  long*           nl, 
  index_st*       ist, 
  par_st*         par,
  int             ham_threads
  ){
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

  long    jat;    
  long    jat_off;
  long    tid_off;
  long    nNL_at = ist->n_NL_atoms;
  long    NL_gpt;
  long    r_idx; 
  long    r;     
  long    r_p;   
  long    alpha;
  long    ngrid = ist->ngrid;

  int     iproj;
  int     ip_up;
  int     ip_dn;
  int     ip; 
  int     s;     
  int     s_p;   
  int     m;     
  int     m_p;   
  int     tid;
  const int     nm = ist->n_l_ang_mom;
  const int     nj = ist->n_j_ang_mom;    
  const int     np = ist->nproj;    
  const int     njnp = nj * np; 
  const int     nmnp = nm * np; 
  int     spin_arr[ist->n_j_ang_mom]; 
  int     m_arr[ist->n_j_ang_mom];    
  int     j;
  int     j_p;   
  int     jtot;  

  double  proj_re;          
  double  proj_im;          
  double  psi_re_up;
  double  psi_im_up;
  double  psi_re_dn;
  double  psi_im_dn;
  double  y1_re; 
  double  y1_im; 
  double  SOproj;


  zomplex LS_loc, PLS;

  // Make arrays for indexing the flattened ang_mom loop
  for (s = 0; s < ist->n_s_ang_mom; s++){
    for (m = 0; m < ist->n_l_ang_mom; m++){
      j = s*3 + m;
      spin_arr[j] = s;
      m_arr[j] = m; 
    }
  }

  // zero out all elements of projs
  alpha = ham_threads * njnp;
  for (j = 0; j < alpha; j++){
    projs[j].re = projs[j].im = 0.0;
  }
  alpha = ham_threads * nm * np;

  // omp_set_num_threads(ham_threads);
  #pragma omp parallel for private(tid, tid_off, jat, jat_off, iproj, ip, j_p, s_p, m_p, NL_gpt, r_idx, r, r_p, psi_re_up, psi_im_up, psi_re_dn, psi_im_dn, y1_re, y1_im, SOproj, j, s, m, jtot, LS_loc, PLS, proj_re, proj_im, alpha)
  for (jat = 0; jat < nNL_at; jat++){
    tid = omp_get_thread_num();

    tid_off = tid * nmnp;
    jat_off = jat * ist->n_NL_gridpts;

    const int nNL_gpt = nl[jat];
    // Compute equation 2.81 of Daniel Weinberg dissertation
    // Action of spin-orbit operator on a real space wavefunctions
    for ( NL_gpt = 0; NL_gpt < nNL_gpt; NL_gpt++){

      r_idx = jat_off + NL_gpt;
      r = nlc[r_idx].jxyz;

      psi_re_up = psi_tmp[r].re;
      psi_im_up = psi_tmp[r].im;
      psi_re_dn = psi_tmp[r + ngrid].re;
      psi_im_dn = psi_tmp[r + ngrid].im;

      for (ip = 0; ip < 5; ip++){

        SOproj = nlc[r_idx].proj[ip];

        for (m_p = 0; m_p < 3; m_p++){
          
          ip_up = tid_off + (ip * nm + m_p);
          ip_dn = tid_off + (ip * nm + m_p) + (ham_threads * nmnp);

          // proj[alpha].re = proj.im = 0.00;
          // Compute the projection of the real space wavefunction onto the basis of |lmr\sigma> 
          // where the radial variable, r, is not computed on a grid but actually are smooth radial 
          // functions
        
          y1_re = nlc[r_idx].y1[m_p].re;
          y1_im = nlc[r_idx].y1[m_p].im;
          
          //weird signs b/c of Y_{lm}^*
          // Calculate the integral in eq 2.81
          // up spin
          projs[ip_up].re += (psi_re_up * y1_re * SOproj + psi_im_up * y1_im * SOproj) * par->dv;
          projs[ip_up].im += (psi_im_up * y1_re * SOproj - psi_re_up * y1_im * SOproj) * par->dv;
          // dn spin
          projs[ip_dn].re += (psi_re_dn * y1_re * SOproj + psi_im_dn * y1_im * SOproj) * par->dv;
          projs[ip_dn].im += (psi_im_dn * y1_re * SOproj - psi_re_dn * y1_im * SOproj) * par->dv;
        }
      }
    }
    
    for (NL_gpt = 0; NL_gpt < nNL_gpt; NL_gpt++){
      r_idx = jat_off + NL_gpt;

      for ( ip = 0; ip < 5; ip++){

        SOproj = nlc[r_idx].proj[ip];

        for ( j_p = 0; j_p < 6; j_p++){
          s_p = spin_arr[j_p];

          iproj = tid_off + (ip * nm + m_p) + s_p * (ham_threads * nmnp);
          
          proj_re = projs[iproj].re;
          proj_im = projs[iproj].im;

          for (j = 0; j < 6; j++){
            //get L_{m,m'}\cdot S_{s,s'}*P_{n,m,s} = PLS_{n,m,m',s,s'}
            jtot = j_p * nj + j;

            LS_loc.re = LS[jtot].re;
            LS_loc.im = LS[jtot].im;

            if (LS_loc.re == 0.0 && LS_loc.im == 0.0){
              // skip this L.S element because it is 0.0
              continue;
            }

            s = spin_arr[j];
            m = m_arr[j];
            
            r_p = nlc[r_idx].jxyz + ist->ngrid * s;

            PLS.re = LS_loc.re * proj_re - LS_loc.im * proj_im;
            PLS.im = LS_loc.re * proj_im + LS_loc.im * proj_re;
         
            y1_re = nlc[r_idx].y1[m].re;
            y1_im = nlc[r_idx].y1[m].im;

            #pragma omp atomic
            psi_out[r_p].re += SOproj * y1_re * PLS.re - SOproj * y1_im * PLS.im;
            #pragma omp atomic
            psi_out[r_p].im += SOproj * y1_re * PLS.im + SOproj * y1_im * PLS.re;
          }
        }
      }
    }
  }

  return;
}

/*****************************************************************************/

void p_nonlocal_proj_pot(
  zomplex*        psi_out, 
  zomplex*        psi_tmp,
  zomplex*        projs, 
  nlc_st*         nlc, 
  long*           nl, 
  index_st*       ist, 
  par_st*         par,
  int             ham_threads
  ){
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

  long jatom, jatom_offset, tid_offset;
  long NL_gridpt, r_idx, r, r_p;
  
  int iproj, s, m, j, tid;
  int nj = ist->n_j_ang_mom;
  int np = ist->nproj;
  int njnp = nj * np;
  int alpha, sgn;
  int spin_arr[ist->n_j_ang_mom], m_arr[ist->n_j_ang_mom];
  
  double proj_re, proj_im;
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

  
  // zero out all elements of projs
  alpha = ham_threads * njnp;
  for (j = 0; j < alpha; j++){
    projs[j].re = projs[j].im = 0.0;
  }

  // omp_set_num_threads(ham_threads);
  #pragma omp parallel for private(tid, tid_offset, jatom, jatom_offset, iproj, j, s, m, NL_gridpt, r_idx, r, r_p, psi_re, psi_im, y1_re, y1_im, NL_proj, proj_re, proj_im, sgn, alpha)
  for (jatom = 0; jatom < ist->n_NL_atoms; jatom++) {
    tid = omp_get_thread_num();

    tid_offset = tid * njnp;
    jatom_offset = jatom * ist->n_NL_gridpts;

    const int nNL_gpt = nl[jatom];

    // First loop over NL_gridpt
    for (NL_gridpt = 0; NL_gridpt < nNL_gpt; NL_gridpt++) {

      r_idx = jatom_offset + NL_gridpt;
      
      // First loop over projector terms
      for (iproj = 0; iproj < 5; iproj++) {

        NL_proj = nlc[r_idx].NL_proj[iproj];
        sgn = nlc[jatom_offset].NL_proj_sign[iproj];
        
        // First loop over ang_mom j
        for (j = 0; j < 6; j++) {
          s = spin_arr[j];
          m = m_arr[j];
          alpha = tid_offset + (iproj * nj + j);

          r = nlc[r_idx].jxyz + ist->ngrid * s;

          psi_re = psi_tmp[r].re;
          psi_im = psi_tmp[r].im;
          y1_re = nlc[r_idx].y1[m].re;
          y1_im = nlc[r_idx].y1[m].im;
          
          projs[alpha].re += sgn * (psi_re * y1_re * NL_proj + psi_im * y1_im * NL_proj) * par->dv;
          projs[alpha].im += sgn * (psi_im * y1_re * NL_proj - psi_re * y1_im * NL_proj) * par->dv;
        } 
      }
      
      // Second loop over iproj
      for (iproj = 0; iproj < 5; iproj++) {
        // Second loop over ang_mom j
        for (j = 0; j < 6; j++) {
          s = spin_arr[j];
          m = m_arr[j];
          alpha = tid_offset + (iproj * nj + j);

          proj_re = projs[alpha].re;
          proj_im = projs[alpha].im;

          r_p = nlc[r_idx].jxyz + ist->ngrid * s;
  
          y1_re = nlc[r_idx].y1[m].re;
          y1_im = nlc[r_idx].y1[m].im;
          NL_proj = nlc[r_idx].NL_proj[iproj];
          
          // Update thread-local psi_out
          #pragma omp atomic
          psi_out[r_p].re += (NL_proj * (y1_re * proj_re - y1_im * proj_im));
          #pragma omp atomic
          psi_out[r_p].im += NL_proj * (y1_re * proj_im + y1_im * proj_re);

        } // end of j
      } // end of iproj
    } // end of NL_gridpt
  } // end of jatom

  return;
}

/*****************************************************************************/

