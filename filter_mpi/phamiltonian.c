#include "fd.h"
#include <complex.h>
/*****************************************************************************/
// File contains the main functions for calcuating the action of the hamiltonian on 
// a wavefunction including non-local spin-orbit contributions 
// parallelized with OpenMP

/*****************************************************************************/

void p_hamiltonian(zomplex *psi_scratch, zomplex *psi_out, zomplex *psi_tmp, double *pot_local, nlc_st *nlc, long *nl, double *ksqr,
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
  p_potential(psi_scratch, psi_out, psi_tmp, pot_local, nlc, nl, ist, par, flag, ham_threads);

  return;
}

/*****************************************************************************/

// This calculates the total action of Vloc + Vnonloc + Vso|psi_tmp>
void p_potential(zomplex *psi_scratch, zomplex *psi_out, zomplex *psi_tmp, double *pot_local, nlc_st *nlc, long *nl, index_st *ist,
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
    p_spin_orbit_proj_pot(psi_scratch, psi_out, psi_tmp, nlc, nl, ist, par, ham_threads);
    //write_state_dat(psi_out, ist->nspinngrid, "psi_out_pSO.dat");
  }
  if (flag->NL == 1){
    // Calculate |psi_out> += V_NL|psi_tmp>
    double* tmp1, *tmp2;
    p_nonlocal_proj_pot(psi_out, psi_tmp,tmp1, tmp2, nlc, nl, ist, par, ham_threads);
    //write_state_dat(psi_out, ist->nspinngrid, "psi_out_pNL.dat");
  }
  
  // Calculate the action of the local potential energy part of the Hamiltonian on psi_tmp
  for (jspin = 0; jspin < ist->nspin; jspin++){
    for (j = 0; j < ist->ngrid; j++) {
      jtmp = ist->ngrid * jspin + j ; // generalized indexing to handle spinors or single-component wavefuncs

      psi_out[jtmp].re += (pot_local[j] * psi_tmp[jtmp].re);
      psi_out[jtmp].im += (pot_local[j] * psi_tmp[jtmp].im);
    }
  }
  //write_state_dat(psi_out, ist->nspinngrid, "psi_out_ploc.dat");
  return;
}


void p_spin_orbit_proj_pot(zomplex *psi_scratch, zomplex *psi_out, zomplex *psi_tmp, nlc_st *nlc, long *nl, index_st *ist, par_st *par, int ham_threads){
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

  long jatom;
  long NL_gridpt, index1, r, index2, r_p;
  int iproj, s, s_p, m, m_p;
  int spin_arr[ist->n_j_ang_mom], m_arr[ist->n_j_ang_mom];
  int j, j_p, jtot;
  int t_id;
  zomplex proj, LS_loc, PLS;
  zomplex *LS;

  for (s = 0; s < ist->n_s_ang_mom; s++){
    for (m = 0; m < ist->n_l_ang_mom; m++){
      j = s*3 + m;
      spin_arr[j] = s;
      m_arr[j] = m; 
    }
  }
  
  double psi_re, psi_im;
  double y1_re, y1_im;
  double nlcproj;

  LS = (zomplex*) calloc(ist->n_j_ang_mom * ist->n_j_ang_mom, sizeof(zomplex));
  def_LS(LS, ist, par);

  // zomplex *psi_scratch = calloc(ist->nthreads * ist->nspinngrid, sizeof(zomplex));
  // if (!psi_scratch) {
  //     fprintf(stderr, "Failed to allocate memory for thread-local psi_out\n");
  //     exit(1);
  // }

  // FILE *pf;
  // pf = fopen("SO.dat", "w");

  omp_set_num_threads(ist->nthreads);
  #pragma omp parallel for private(jatom, NL_gridpt, index1, r, index2, r_p, iproj, proj, s, s_p, m, m_p, j, j_p, jtot, psi_re, psi_im, y1_re, y1_im, nlcproj, LS_loc, PLS, t_id)
  for (jatom = 0; jatom < ist->n_NL_atoms; jatom++){
    t_id = omp_get_thread_num();
    // printf("Inner thread %d\n", t_id);
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
        // FILE *pz;
        // char pz_name[50];
        // sprintf(pz_name, "n_zeros-%ld.dat", jatom);
        // pz = fopen(pz_name, "w");
        // long fiv_ct = 0;
        // long six_ct = 0;
        // long sev_ct = 0;
        // long eig_ct = 0; 
        // long nin_ct = 0;
        // long ten_ct = 0;
        for ( NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++){
          index1 = jatom * ist->n_NL_gridpts + NL_gridpt;
          r = nlc[index1].jxyz + (ist->ngrid) * s_p;
          psi_re = psi_tmp[r].re;
          psi_im = psi_tmp[r].im;
          // If both psi_re and psi_im are zero, then there is no need to continue the loop
          // if ( (psi_re < par->psi_zero_cut) && (psi_im < par->psi_zero_cut) ){
          //   // zero_ct++;
          //   continue;
          // }
          y1_re = nlc[index1].y1[m_p].re;
          y1_im = nlc[index1].y1[m_p].im;
          nlcproj = nlc[index1].proj[iproj];
          //weird signs b/c of Y_{lm}^*
          // Calculate the integral in eq 2.81
          proj.re += psi_re * y1_re * nlcproj + psi_im * y1_im * nlcproj;
          proj.im += psi_im * y1_re * nlcproj - psi_re * y1_im * nlcproj;
        }
        // fprintf(pz, "Five count = %ld / %ld\n", fiv_ct, nl[jatom]);
        // fprintf(pz, "Six count = %ld / %ld\n", six_ct, nl[jatom]);
        // fprintf(pz, "Sev count = %ld / %ld\n", sev_ct, nl[jatom]);
        // fprintf(pz, "Eight count = %ld / %ld\n", eig_ct, nl[jatom]); 
        // fprintf(pz, "Nine count = %ld / %ld\n", nin_ct, nl[jatom]); 
        // fprintf(pz, "Ten count = %ld / %ld\n", ten_ct, nl[jatom]); //fflush(pz); 
        // fclose(pz);
        proj.re *= par->dv;
        proj.im *= par->dv;
          
        for (j = 0; j < ist->n_j_ang_mom; j++){
          //get L_{m,m'}\cdot S_{s,s'}*P_{n,m,s} = PLS_{n,m,m',s,s'}
          s = spin_arr[j];
          m = m_arr[j];
          jtot = j_p * ist->n_j_ang_mom + j;

          PLS.re = PLS.im = 0.00;
          LS_loc.re = LS[jtot].re;
          LS_loc.im = LS[jtot].im;
          if (LS_loc.re == 0.0 && LS_loc.im == 0.0){
            continue;
          }
          PLS.re = LS_loc.re * proj.re - LS_loc.im * proj.im;
          PLS.im = LS_loc.re * proj.im + LS_loc.im * proj.re;
          
          for (NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++){
            index2 = jatom * ist->n_NL_gridpts + NL_gridpt;
            r_p = nlc[index2].jxyz + (ist->ngrid)*s;
            nlcproj = nlc[index2].proj[iproj];
            y1_re = nlc[index2].y1[m].re;
            y1_im = nlc[index2].y1[m].im;

            psi_scratch[t_id * ist->nspinngrid + r_p].re += nlcproj * y1_re * PLS.re - nlcproj * y1_im * PLS.im;
            psi_scratch[t_id * ist->nspinngrid + r_p].im += nlcproj * y1_re * PLS.im + nlcproj * y1_im * PLS.re;
          }
        }
      }
    }
  }

  for (int t_idx = 0; t_idx < ist->nthreads; t_idx++){
    for (long jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
      psi_out[jgrid].re += psi_scratch[t_idx*ist->nspinngrid + jgrid].re;
      psi_out[jgrid].im += psi_scratch[t_idx*ist->nspinngrid + jgrid].im;
    }
  } 
  

  // free(psi_scratch);

  return;
}

/*****************************************************************************/

void p_nonlocal_proj_pot(zomplex *psi_out, zomplex *psi_tmp, double *psi_out_re, double *psi_out_im, nlc_st *nlc, long *nl, index_st *ist, par_st *par, int ham_threads){
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
  long NL_gridpt, index1, index2, r, r_p;
  int iproj, s, m, j, t_id;
  int spin_arr[ist->n_j_ang_mom], m_arr[ist->n_j_ang_mom];
  // zomplex proj;
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

  long nNL = ist->n_NL_gridpts;
   
  // Thread-local storage for scratch work
  // zomplex *scratch = (zomplex *)calloc(ist->nthreads * ist->n_NL_gridpts, sizeof(zomplex));
  // if (!scratch) {
  //     fprintf(stderr, "Error: Unable to allocate thread-local storage.\n");
  //     exit(EXIT_FAILURE);
  // }
  // num_threads(ham_threads)
  for (jatom = 0; jatom < ist->n_NL_atoms; jatom++) {
    
    jatom_offset = jatom * ist->n_NL_gridpts;

    for (iproj = 0; iproj < ist->nproj; iproj++) {
      for (j = 0; j < ist->n_j_ang_mom; j++) {
        s = spin_arr[j];
        m = m_arr[j];

        proj_re = proj_im = 0.0;

        // First loop over NL_gridpt
        
        #pragma omp parallel for private(jatom, jatom_offset, NL_gridpt, index1, r, iproj, s, m, psi_re, psi_im, y1_re, y1_im, NL_proj) reduction(+:proj_re, proj_im)
        for (NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++) {
          // t_id = omp_get_thread_num();
          index1 = jatom_offset + NL_gridpt;
          r = nlc[index1].jxyz + ist->ngrid * s;

          psi_re = psi_tmp[r].re;
          psi_im = psi_tmp[r].im;
          y1_re = nlc[index1].y1[m].re;
          y1_im = nlc[index1].y1[m].im;
          // if ((psi_re < par->psi_zero_cut) && (psi_im < par->psi_zero_cut)){
          //   continue;
          // }
          NL_proj = nlc[index1].NL_proj[iproj];

          proj_re += psi_re * y1_re * NL_proj + psi_im * y1_im * NL_proj;
          proj_im += psi_im * y1_re * NL_proj - psi_re * y1_im * NL_proj;
        }
        // printf("post proj\n"); fflush(0);
        // Collect proj from threads
        // for (t_id = 0; t_id < ist->nthreads; t_id++){
        //   for (NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++){
        //     proj.re += psi_scratch[t_id*ist->n_NL_gridpts + NL_gridpt].re;
        //     proj.im += psi_scratch[t_id*ist->n_NL_gridpts + NL_gridpt].im;
        //   }
        // }

        // Apply scaling factors
        proj_re *= nlc[jatom_offset].NL_proj_sign[iproj] * par->dv;
        proj_im *= nlc[jatom_offset].NL_proj_sign[iproj] * par->dv;

        // Second loop over NL_gridpt
        #pragma omp parallel for private(jatom, jatom_offset, NL_gridpt, index2, iproj, r_p, s, m, y1_re, y1_im, NL_proj)
        for (NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++) {
          // t_id = omp_get_thread_num();
          index2 = jatom_offset + NL_gridpt;
          r_p = nlc[index2].jxyz + ist->ngrid * s;

          y1_re = nlc[index2].y1[m].re;
          y1_im = nlc[index2].y1[m].im;
          NL_proj = nlc[index2].NL_proj[iproj];
          // printf("pre psi_out\n"); fflush(0);
          // Update thread-local psi_scratch
          #pragma omp atomic
          {
          psi_out[r_p].re += (NL_proj * (y1_re * proj_re - y1_im * proj_im));
          psi_out[r_p].im += NL_proj * (y1_re * proj_im + y1_im * proj_re);
          }
          // printf("post psi_out\n"); fflush(0);
        }

        // for (t_id = 0; t_id < ist->nthreads; t_id++){
          
        // }
      }
    }
  }

  // Reduction step: Combine thread-local psi_scratch into the global psi_out
  // for (int t_idx = 0; t_idx < ist->nthreads; t_idx++){
  //   for (long jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
  //     psi_out[jgrid].re += psi_scratch[t_idx * ist->nspinngrid + jgrid].re;
  //     psi_out[jgrid].im += psi_scratch[t_idx * ist->nspinngrid + jgrid].im;
  //   }
  // }

  // for (jatom = 0; jatom < ist->natoms; jatom++){
  //   jatom_offset = jatom * ist->n_NL_gridpts;
  //   for (NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++){
  //     index2 = jatom_offset + NL_gridpt;
  //     r_p = nlc[index2].jxyz + ist->ngrid * s;
  //     psi_out[r_p].re += psi_out_re[NL_gridpt];
  //     psi_out[r_p].im += psi_out_im[NL_gridpt];
  //   }
  // }

  // Free thread-local storage
  // free(psi_scratch);
  

  return;
}

/*****************************************************************************/

