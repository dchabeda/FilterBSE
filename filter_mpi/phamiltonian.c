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
  *  [planbw] FFTW3 plan for executing 3D backwards DFT              *
  *  [fftwpsi] location to store outcome of Fourier transform        *
  * outputs: void                                                    *
  ********************************************************************/

  int jspin; 
  
  // Copy psi_out into psi_tmp
  memcpy(&psi_tmp[0], &psi_out[0], ist->nspinngrid*sizeof(psi_tmp[0]));
  // write_state_dat(psi_out, ist->nspinngrid, "psi_out_init.dat");
  
  // Calculate the action of the kinetic energy part of the Hamiltonian on psi_tmp: |psi_out> = T|psi_tmp>
  for (jspin = 0; jspin < ist->nspin; jspin++){
      kinetic(&psi_out[jspin*ist->ngrid], ksqr, planfw, planbw, fftwpsi, ist); //spin up/down
  } 
  // write_state_dat(psi_out, ist->nspinngrid, "psi_out_pkinetic.dat");
  // Calculate the action of the potential on the wavefunction: |psi_out> = V|psi_tmp>
  
  p_potential(psi_out, psi_tmp, pot_local, LS, nlc, nl, ist, par, flag, ham_threads, (vector){0});

  return;
}

/*****************************************************************************/

// This calculates the total action of Vloc + Vnonloc + Vso|psi_tmp>
void p_potential(
  zomplex*        psi_out, 
  zomplex*        psi_tmp, 
  double*         pot_local,
  zomplex*        LS, 
  nlc_st*         nlc, 
  long*           nl, 
  index_st*       ist,
  par_st*         par,
  flag_st*        flag,
  int             ham_threads,
  vector          k
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
    p_spin_orbit_proj_pot(psi_out, psi_tmp, LS, nlc, nl, ist, par, ham_threads, k);
    // write_state_dat(psi_out, ist->nspinngrid, "psi_out_pSO.dat");
  }
  if (flag->NL == 1){
    // Calculate |psi_out> += V_NL|psi_tmp>
    p_nonlocal_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par, ham_threads, k);
    // write_state_dat(psi_out, ist->nspinngrid, "psi_out_pNL.dat");
  }
  
  // Calculate the action of the local potential energy part of the Hamiltonian on psi_tmp
  // No openMP atomics necessary for this because all grid indices are unique
  // (no two threads reads/writes the same j/jtmp ever)
  
  // TODO: Spinor version should store pot_local twice, contiguously in memory
  // so that this loop is flat over nspinngrid. Daniel C. 3.2.2025
  // NOTE: parallelizing over grid points here is not efficient. Not worthwhile
  // The local potential must hit the imaginary channel too whenever the wavefunction
  // is complex (periodic k != 0, or spinor/SO runs). Previously this branched on
  // useSpinors, so the periodic non-SO case (isComplex==1, useSpinors==0) applied V
  // only to psi.re -- dropping V|Im(psi)>. That made V (hence H) non-Hermitian on
  // complex states (Im<psi|V|psi> != 0), so diag_H/the filter built eigenvectors of a
  // different operator than calc_sigma_E_k measured. Branch on isComplex to match the
  // serial potential() in hamiltonian.c. nspinngrid == ngrid when nspin == 1, so the
  // real path is unchanged.
  if (1 == flag->isComplex){
    // #pragma omp parallel for private(j)
    for (j = 0; j < ist->nspinngrid; j++) {
      psi_out[j].re += (pot_local[j] * psi_tmp[j].re);
      psi_out[j].im += (pot_local[j] * psi_tmp[j].im);
    }
  }
  else {
    for (j = 0; j < ist->nspinngrid; j++) {
      psi_out[j].re += (pot_local[j] * psi_tmp[j].re);
    }
  }
  
  // write_state_dat(psi_out, ist->nspinngrid, "psi_out_ploc.dat");
  
  return;
}


void p_spin_orbit_proj_pot(
  zomplex*        psi_out,
  zomplex*        psi_tmp,
  zomplex*        LS,
  nlc_st*         nlc,
  long*           nl,
  index_st*       ist,
  par_st*         par,
  int             ham_threads,
  vector          k
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
  long    nNL_at = ist->n_NL_atoms;
  long    NL_gpt;
  long    r_idx; 
  long    r;     
  long    r_p;

  int     ip; 
  int     s;     
  int     s_p;   
  int     m;     
  int     m_p;
  int     spin_arr[ist->n_j_ang_mom]; 
  int     m_arr[ist->n_j_ang_mom];    
  int     j;
  int     j_p;   
  int     jtot;  
         
  double  psi_re;
  double  psi_im;
  double  y1_re; 
  double  y1_im; 
  double  SOproj;


  zomplex proj, LS_loc, PLS;

  // Make arrays for indexing the flattened ang_mom loop
  for (s = 0; s < ist->n_s_ang_mom; s++){
    for (m = 0; m < ist->n_l_ang_mom; m++){
      j = s*3 + m;
      spin_arr[j] = s;
      m_arr[j] = m; 
    }
  }

  // Bloch phase: periodic (k != 0) operator twists projectors by e^{-i k.delta}
  // (ket) / e^{+i k.delta} (bra). has_phase == 0 reduces to the real-space form.
  const int has_phase = (k.x != 0.0) || (k.y != 0.0) || (k.z != 0.0);

  omp_set_num_threads(ham_threads);
  #pragma omp parallel for private(jat, jat_off, proj, ip, j_p, s_p, m_p, NL_gpt, r_idx, r, r_p, psi_re, psi_im, y1_re, y1_im, SOproj, j, s, m, jtot, LS_loc, PLS)
  for (jat = 0; jat < nNL_at; jat++){

    jat_off = jat * ist->n_NL_gridpts;

    const int nNL_gpt = nl[jat];

    // Per-gridpoint phase factors cos/sin(k.delta) for this atom (thread-private).
    double pc[nNL_gpt > 0 ? nNL_gpt : 1], ps[nNL_gpt > 0 ? nNL_gpt : 1];
    if (has_phase)
      for (NL_gpt = 0; NL_gpt < nNL_gpt; NL_gpt++){
        r_idx = jat_off + NL_gpt;
        double phi = k.x * nlc[r_idx].dx + k.y * nlc[r_idx].dy + k.z * nlc[r_idx].dz;
        pc[NL_gpt] = cos(phi);
        ps[NL_gpt] = sin(phi);
      }

    // Compute equation 2.81 of Daniel Weinberg dissertation
    // Action of spin-orbit operator on a real space wavefunctions
    for (ip = 0; ip < 5; ip++){
      for (j_p = 0; j_p < 6; j_p++){
        // Compute the projection of the real space wavefunction onto the basis of |lmr\sigma>
        // where the radial variable, r, is not computed on a grid but actually are smooth radial
        // functions
        s_p = spin_arr[j_p];
        m_p = m_arr[j_p];
        proj.re = proj.im = 0.00;
        for ( NL_gpt = 0; NL_gpt < nNL_gpt; NL_gpt++){

          r_idx = jat_off + NL_gpt;
          r = nlc[r_idx].jxyz + ist->ngrid * s_p;

          psi_re = psi_tmp[r].re;
          psi_im = psi_tmp[r].im;

          SOproj = nlc[r_idx].proj[ip];

          y1_re = nlc[r_idx].y1[m_p].re;
          y1_im = nlc[r_idx].y1[m_p].im;

          //weird signs b/c of Y_{lm}^*: term = SOproj * (psi . conj(y1))
          double tre = SOproj * (psi_re * y1_re + psi_im * y1_im);
          double tim = SOproj * (psi_im * y1_re - psi_re * y1_im);
          if (has_phase){ // bra twist: x e^{+i k.delta}
            double c = pc[NL_gpt], sph = ps[NL_gpt];
            proj.re += tre * c - tim * sph;
            proj.im += tre * sph + tim * c;
          } else {
            proj.re += tre;
            proj.im += tim;
          }
        }

        proj.re *= par->dv;
        proj.im *= par->dv;

        for ( j = 0; j < 6; j++){
          jtot = j_p * 6 + j;

          LS_loc.re = LS[jtot].re;
          LS_loc.im = LS[jtot].im;

          if (LS_loc.re == 0.0 && LS_loc.im == 0.0){
            // skip this L.S element because it is 0.0
            continue;
          }

          s = spin_arr[j];
          m = m_arr[j];

          PLS.re = LS_loc.re * proj.re - LS_loc.im * proj.im;
          PLS.im = LS_loc.re * proj.im + LS_loc.im * proj.re;

          for (NL_gpt = 0; NL_gpt < nNL_gpt; NL_gpt++){
            r_idx = jat_off + NL_gpt;
            r_p = nlc[r_idx].jxyz + ist->ngrid * s;

            SOproj = nlc[r_idx].proj[ip];

            y1_re = nlc[r_idx].y1[m].re;
            y1_im = nlc[r_idx].y1[m].im;

            // out = SOproj * (y1 . PLS)
            double ore = SOproj * (y1_re * PLS.re - y1_im * PLS.im);
            double oim = SOproj * (y1_re * PLS.im + y1_im * PLS.re);
            if (has_phase){ // ket twist: x e^{-i k.delta}
              double c = pc[NL_gpt], sph = ps[NL_gpt];
              double rre = ore * c + oim * sph;
              double rim = oim * c - ore * sph;
              #pragma omp atomic
              psi_out[r_p].re += rre;
              #pragma omp atomic
              psi_out[r_p].im += rim;
            } else {
              #pragma omp atomic
              psi_out[r_p].re += ore;
              #pragma omp atomic
              psi_out[r_p].im += oim;
            }
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
  nlc_st*         nlc,
  long*           nl,
  index_st*       ist,
  par_st*         par,
  int             ham_threads,
  vector          k
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

  long jatom, jatom_offset;
  long NL_gridpt, r_idx, r, r_p;
  
  int iproj, s, m, j;
  int sgn;
  int spin_arr[ist->n_j_ang_mom], m_arr[ist->n_j_ang_mom];
  
  double psi_re, psi_im;
  double y1_re, y1_im;
  double NL_proj;

  zomplex proj;
  
  for (s = 0; s < ist->n_s_ang_mom; s++){
    for (m = 0; m < ist->n_l_ang_mom; m++){
      j = s*ist->n_l_ang_mom + m;
      spin_arr[j] = s;
      m_arr[j] = m; 
    }
  }

  // Bloch phase twist for k != 0; has_phase == 0 reduces to the real-space form.
  const int has_phase = (k.x != 0.0) || (k.y != 0.0) || (k.z != 0.0);

  omp_set_num_threads(ham_threads);
  #pragma omp parallel for private(jatom, jatom_offset, proj, j, s, m, NL_gridpt, r_idx, r, r_p, psi_re, psi_im, y1_re, y1_im, NL_proj, sgn)
  for (jatom = 0; jatom < ist->n_NL_atoms; jatom++) {

    jatom_offset = jatom * ist->n_NL_gridpts;

    const int nNL_gpt = nl[jatom];

    // Per-gridpoint phase factors cos/sin(k.delta) for this atom (thread-private).
    double pc[nNL_gpt > 0 ? nNL_gpt : 1], ps[nNL_gpt > 0 ? nNL_gpt : 1];
    if (has_phase)
      for (NL_gridpt = 0; NL_gridpt < nNL_gpt; NL_gridpt++) {
        r_idx = jatom_offset + NL_gridpt;
        double phi = k.x * nlc[r_idx].dx + k.y * nlc[r_idx].dy + k.z * nlc[r_idx].dz;
        pc[NL_gridpt] = cos(phi);
        ps[NL_gridpt] = sin(phi);
      }

    for (iproj = 0; iproj < 5; iproj++) {
      for (j = 0; j < 6; j++) {
        s = spin_arr[j];
        m = m_arr[j];

        proj.re = proj.im = 0.0;
        for (NL_gridpt = 0; NL_gridpt < nNL_gpt; NL_gridpt++) {

          r_idx = jatom_offset + NL_gridpt;
          r = nlc[r_idx].jxyz + ist->ngrid * s;

          NL_proj = nlc[r_idx].NL_proj[iproj];

          psi_re = psi_tmp[r].re;
          psi_im = psi_tmp[r].im;

          y1_re = nlc[r_idx].y1[m].re;
          y1_im = nlc[r_idx].y1[m].im;

          double tre = NL_proj * (psi_re * y1_re + psi_im * y1_im);
          double tim = NL_proj * (psi_im * y1_re - psi_re * y1_im);
          if (has_phase){ // bra twist: x e^{+i k.delta}
            double c = pc[NL_gridpt], sph = ps[NL_gridpt];
            proj.re += tre * c - tim * sph;
            proj.im += tre * sph + tim * c;
          } else {
            proj.re += tre;
            proj.im += tim;
          }
        }

        sgn = nlc[jatom_offset].NL_proj_sign[iproj];
        proj.re *= sgn * par->dv;
        proj.im *= sgn * par->dv;

        for (NL_gridpt = 0; NL_gridpt < nNL_gpt; NL_gridpt++) {
          r_idx = jatom_offset + NL_gridpt;
          r_p = nlc[r_idx].jxyz + ist->ngrid * s;

          y1_re = nlc[r_idx].y1[m].re;
          y1_im = nlc[r_idx].y1[m].im;

          NL_proj = nlc[r_idx].NL_proj[iproj];

          double ore = NL_proj * (y1_re * proj.re - y1_im * proj.im);
          double oim = NL_proj * (y1_re * proj.im + y1_im * proj.re);
          if (has_phase){ // ket twist: x e^{-i k.delta}
            double c = pc[NL_gridpt], sph = ps[NL_gridpt];
            double rre = ore * c + oim * sph;
            double rim = oim * c - ore * sph;
            #pragma omp atomic
            psi_out[r_p].re += rre;
            #pragma omp atomic
            psi_out[r_p].im += rim;
          } else {
            #pragma omp atomic
            psi_out[r_p].re += ore;
            #pragma omp atomic
            psi_out[r_p].im += oim;
          }
        }
      }
    } // end of NL_gridpt
  } // end of jatom

  return;
}

/*****************************************************************************/

