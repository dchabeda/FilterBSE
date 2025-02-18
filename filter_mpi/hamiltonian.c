/*****************************************************************************/
// File contains the main functions for calcuating the action of the hamiltonian on 
// a wavefunction including non-local spin-orbit contributions 

#include "hamiltonian.h"

/*****************************************************************************/

void hamiltonian(
  zomplex*       psi_out, 
  zomplex*       psi_tmp, 
  double*        pot_local, 
  zomplex*       LS, 
  nlc_st*        nlc, 
  long*          nl, 
  double*        ksqr,
  index_st*      ist, 
  par_st*        par, 
  flag_st*       flag, 
  fftw_plan_loc  planfw, 
  fftw_plan_loc  planbw, 
  fftw_complex*  fftwpsi){
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
  
  // write_state_dat(psi_out, ist->nspinngrid, "psi_out_kinetic.dat");
  // Calculate the action of the potential on the wavefunction: |psi_out> = V|psi_tmp>
  potential(psi_out, psi_tmp, pot_local, LS, nlc, nl, ist, par, flag);

  return;
}

/*****************************************************************************/
// Calculates T|psi_tmp> via FFT and stores result in psi_out 

void kinetic(zomplex *psi_out, double *ksqr, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi, index_st *ist){
  /*******************************************************************
  * This function applies the KE operator onto a state               *
  * T|psi> = FT^-1[ k^2 * FT[|psi>] ]                                *
  * inputs:                                                          *
  *  [psi_out] nspinngrid-long arr to hold |psi_out> = T|psi_tmp>     *
  *  [ksqr] ngrid-long arr holding the values of k^2 for KE calc     *
  *  [ist] ptr to counters, indices, and lengths                     *
  *  [planfw] FFTW3 plan for executing 3D forward DFT                *
  *  [planfw] FFTW3 plan for executing 3D backwards DFT              *
  *  [fftwpsi] location to store outcome of Fourier transform        *
  * outputs: void                                                    *
  ********************************************************************/

  long j;
  

  // Copy inputted psi to fftwpsi
  memcpy(&fftwpsi[0], &psi_out[0], ist->ngrid*sizeof(fftwpsi[0]));
  

  // FT from r-space to k-space
  fftw_execute(planfw);

  // Kinetic energy is diagonal in k-space, just multiply fftwpsi by k^2
  for (j = 0; j < ist->ngrid; j++) {
    fftwpsi[j][0] *= ksqr[j];
    fftwpsi[j][1] *= ksqr[j];
  }

  // Inverse FT back to r-space
  fftw_execute(planbw);
  
  // Copy fftwpsi to psi_out to store T|psi_tmp> into |psi_out>
  memcpy(&psi_out[0], &fftwpsi[0], ist->ngrid*sizeof(fftwpsi[0]));
  
  return;
}


/*****************************************************************************/

// This calculates the total action of Vloc + Vnonloc + Vso|psi_tmp>
void potential(zomplex *psi_out, zomplex *psi_tmp, double *pot_local, zomplex *LS, nlc_st *nlc, long *nl, index_st *ist,
  par_st *par, flag_st *flag){
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
    spin_orbit_proj_pot(psi_out, psi_tmp, LS, nlc, nl, ist, par);
    // write_state_dat(psi_out, ist->nspinngrid, "psi_out_SO.dat");
  }
  if (flag->NL == 1){
    // Calculate |psi_out> += V_NL|psi_tmp>
    nonlocal_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par);
    // write_state_dat(psi_out, ist->nspinngrid, "psi_out_NL.dat");
  }
  
  // Calculate the action of the local potential energy part of the Hamiltonian on psi_tmp
  for (jspin = 0; jspin < ist->nspin; jspin++){
    for (j = 0; j < ist->ngrid; j++) {
      jtmp = ist->ngrid * jspin + j ; // generalized indexing to handle spinors or spinless wavefuncs

      psi_out[jtmp].re += (pot_local[j] * psi_tmp[jtmp].re);
      psi_out[jtmp].im += (pot_local[j] * psi_tmp[jtmp].im);
    }
  }
  // write_state_dat(psi_out, ist->nspinngrid, "psi_out_loc.dat");
  return;
}


void spin_orbit_proj_pot(
  zomplex*     psi_out, 
  zomplex*     psi_tmp,
  zomplex*     LS,
  nlc_st*      nlc, 
  long*        nl, 
  index_st*    ist, 
  par_st*      par){
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
  zomplex proj;

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


  for ( jatom = 0; jatom < ist->n_NL_atoms; jatom++){
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
          index1 = jatom * ist->n_NL_gridpts + NL_gridpt;
          r = nlc[index1].jxyz + (ist->ngrid) * s_p;
          psi_re = psi_tmp[r].re;
          psi_im = psi_tmp[r].im;
          y1_re = nlc[index1].y1[m_p].re;
          y1_im = nlc[index1].y1[m_p].im;
          nlcproj = nlc[index1].proj[iproj];
          //weird signs b/c of Y_{lm}^*
          // Calculate the integral in eq 2.81
          proj.re += psi_re * y1_re * nlcproj + psi_im * y1_im * nlcproj;
          proj.im += psi_im * y1_re * nlcproj - psi_re * y1_im * nlcproj;
          
        }
        
        proj.re *= par->dv;
        proj.im *= par->dv;
          
        for (j = 0; j < ist->n_j_ang_mom; j++){
          //get L_{m,m'}\cdot S_{s,s'}*P_{n,m,s} = PLS_{n,m,m',s,s'}
          jtot = j_p * ist->n_j_ang_mom + j;
          if ( LS[jtot].re == 0.0 ){
            continue;
          }
          s = spin_arr[j];
          m = m_arr[j];
          
          zomplex PLS, LS_loc, LS_test;

          PLS.re = PLS.im = 0.00;
          LS_loc = LS[jtot];
          
          PLS.re = LS_loc.re * proj.re - LS_loc.im * proj.im;
          PLS.im = LS_loc.re * proj.im + LS_loc.im * proj.re;
          // fprintf(pf, "LS.re %lg LS.im %lg PLS.re %lg PLS.im %lg\n", LS.re, LS.im, PLS.re, PLS.im);
          
          for (NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++){
            index2 = jatom * ist->n_NL_gridpts + NL_gridpt;
            r_p = nlc[index2].jxyz + (ist->ngrid)*s;
            nlcproj = nlc[index2].proj[iproj];
            y1_re = nlc[index2].y1[m].re;
            y1_im = nlc[index2].y1[m].im;

            psi_out[r_p].re += nlcproj * y1_re * PLS.re - nlcproj * y1_im * PLS.im;
            psi_out[r_p].im += nlcproj * y1_re * PLS.im + nlcproj * y1_im * PLS.re;
          }
        
        }
      }
    }
  }
  

  return;
}

/*****************************************************************************/

void def_LS(zomplex *LS, index_st *ist, par_st *par){
  int spin_arr[6], m_arr[6];
  int j, j_p, jtot, s, s_p, m, m_p;
  zomplex Lx[3][3], Ly[3][3], Lz[3][3];
  zomplex Sx[2][2], Sy[2][2], Sz[2][2];
  double sq2 = sqrt(0.5);
  
  // Create an indexing array for flattened s, m loops
  for (s = 0; s < ist->n_s_ang_mom; s++){
    for (m = 0; m < ist->n_l_ang_mom; m++){
      j = s*ist->n_l_ang_mom + m;
      spin_arr[j] = s;
      m_arr[j] = m; 
    }
  }

  // Define Lx
  Lx[0][0].re=0.00; Lx[0][0].im=0.00;     Lx[0][1].re=sq2;  Lx[0][1].im=0.00;     Lx[0][2].re=0.00; Lx[0][2].im=0.00;
  Lx[1][0].re=sq2;  Lx[1][0].im=0.00;     Lx[1][1].re=0.00; Lx[1][1].im=0.00;     Lx[1][2].re=sq2;  Lx[1][2].im=0.00;
  Lx[2][0].re=0.00; Lx[2][0].im=0.00;     Lx[2][1].re=sq2;  Lx[2][1].im=0.00;     Lx[2][2].re=0.00; Lx[2][2].im=0.00;
  // Define Ly
  Ly[0][0].re=0.00; Ly[0][0].im=0.00;     Ly[0][1].re=0.00; Ly[0][1].im=1.0*sq2;  Ly[0][2].re=0.00; Ly[0][2].im=0.00;
  Ly[1][0].re=0.00; Ly[1][0].im=-1.0*sq2; Ly[1][1].re=0.00; Ly[1][1].im=0.00;     Ly[1][2].re=0.00; Ly[1][2].im=1.0*sq2;
  Ly[2][0].re=0.00; Ly[2][0].im=0.00;     Ly[2][1].re=0.00; Ly[2][1].im=-1.0*sq2; Ly[2][2].re=0.00; Ly[2][2].im=0.00;
  // Define Lz 
  Lz[0][0].re=-1.0 ;Lz[0][0].im=0.00;     Lz[0][1].re=0.00; Lz[0][1].im=0.00;     Lz[0][2].re=0.00; Lz[0][2].im=0.00;
  Lz[1][0].re=0.00; Lz[1][0].im=0.00;     Lz[1][1].re=0.00; Lz[1][1].im=0.00;     Lz[1][2].re=0.00; Lz[1][2].im=0.00;
  Lz[2][0].re=0.00; Lz[2][0].im=0.00;     Lz[2][1].re=0.00; Lz[2][1].im=0.00;     Lz[2][2].re=1.00; Lz[2][2].im=0.00;

  //Define Sx
  Sx[0][0].re=0.00; Sx[0][0].im=0.00;     Sx[0][1].re=0.50; Sx[0][1].im=0.00;
  Sx[1][0].re=0.50; Sx[1][0].im=0.00;     Sx[1][1].re=0.00; Sx[1][1].im=0.00;
  //Define Sy
  Sy[0][0].re=0.00; Sy[0][0].im=0.00;     Sy[0][1].re=0.00; Sy[0][1].im=-0.50;
  Sy[1][0].re=0.00; Sy[1][0].im=0.50;     Sy[1][1].re=0.00; Sy[1][1].im=0.00;
  //Define Sz
  Sz[0][0].re=0.50; Sz[0][0].im=0.00;     Sz[0][1].re=0.00;  Sz[0][1].im=0.00;
  Sz[1][0].re=0.00; Sz[1][0].im=0.00;     Sz[1][1].re=-0.50; Sz[1][1].im=0.00;

  // Calculate the elements of L.S in the s,m,s_p,m_p basis, but flattened over j
  // These are not int the total ang_mom basis, the loops have just been flattened for
  // computational efficiency
  for (j_p = 0; j_p < ist->n_j_ang_mom; j_p++){
    s_p = spin_arr[j_p];
    m_p = m_arr[j_p];
    for (j = 0; j < ist->n_j_ang_mom; j++){
      s = spin_arr[j];
      m = m_arr[j];
      jtot = j_p*ist->n_j_ang_mom + j;

      LS[jtot].re += (  Lx[m][m_p].re * Sx[s][s_p].re - Lx[m][m_p].im * Sx[s][s_p].im);
      LS[jtot].im += (  Lx[m][m_p].re * Sx[s][s_p].im + Lx[m][m_p].im * Sx[s][s_p].re);

      LS[jtot].re += (  Ly[m][m_p].re * Sy[s][s_p].re - Ly[m][m_p].im * Sy[s][s_p].im);
      LS[jtot].im += (  Ly[m][m_p].re * Sy[s][s_p].im + Ly[m][m_p].im * Sy[s][s_p].re);

      LS[jtot].re += (  Lz[m][m_p].re * Sz[s][s_p].re - Lz[m][m_p].im * Sz[s][s_p].im);
      LS[jtot].im += (  Lz[m][m_p].re * Sz[s][s_p].im + Lz[m][m_p].im * Sz[s][s_p].re);
    }
  }

  
  return;
}

/*****************************************************************************/

void nonlocal_proj_pot(zomplex *psi_out, zomplex *psi_tmp, nlc_st *nlc, long *nl, index_st *ist, par_st *par){
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

  long jatom;
  long NL_gridpt, index1, r, index2, r_p;
  int iproj, spin, m;
  zomplex proj;
  
  for ( jatom =0; jatom < ist->n_NL_atoms; jatom++){
    for ( iproj = 0; iproj < ist->nproj; iproj++){
      for ( spin = 0; spin < 2; spin++){
        for ( m = 0; m < 3; m++){
          proj.re = proj.im = 0.00;
          for ( NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++){
            index1 = jatom * ist->n_NL_gridpts + NL_gridpt;
            r = nlc[index1].jxyz + (ist->ngrid)*spin;
            
            //weird signs b/c of Y_{lm}^*
            //TODO: I've checked these iprojs against the plane wave way to generate iprojs and they match..
            //When I use these iprojs with the LdotS I've checked then i get the same as the RS energy. 
            
            proj.re += psi_tmp[r].re * nlc[index1].y1[m].re * nlc[index1].NL_proj[iproj];
            proj.re += psi_tmp[r].im * nlc[index1].y1[m].im * nlc[index1].NL_proj[iproj];
           
            proj.im += psi_tmp[r].im * nlc[index1].y1[m].re * nlc[index1].NL_proj[iproj];
            proj.im -= psi_tmp[r].re * nlc[index1].y1[m].im * nlc[index1].NL_proj[iproj];
            
          }
          proj.re *= nlc[index1].NL_proj_sign[iproj];
          proj.im *= nlc[index1].NL_proj_sign[iproj];   
          
          proj.re *= par->dv;
          proj.im *= par->dv;
          
          for (NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++){
            index2 = jatom * ist->n_NL_gridpts + NL_gridpt;
            r_p = nlc[index2].jxyz + (ist->ngrid)*spin;

            psi_out[r_p].re += nlc[index2].NL_proj[iproj] * nlc[index2].y1[m].re * proj.re;
            psi_out[r_p].re -= nlc[index2].NL_proj[iproj] * nlc[index2].y1[m].im * proj.im;

            psi_out[r_p].im += nlc[index2].NL_proj[iproj] * nlc[index2].y1[m].re * proj.im;
            psi_out[r_p].im += nlc[index2].NL_proj[iproj] * nlc[index2].y1[m].im * proj.re;
            

          }
        }
      }
    }
  }

  return;
}

/*****************************************************************************/

void time_reverse_all(double *psitot, double *dest, index_st *ist, parallel_st *parallel){
  /*******************************************************************
  * This function applies the real-space time reversal operator to   *
  * all filtered states to generate twice the orthogonal states      *
  * THETA psi(r) alpha =  psi(r)^* beta                              *
  * THETA psi(r) beta = - psi(r)^* alpha                             *
  * inputs:                                                          *
  *  [psitot] m*n*nspinngrid-long arr where all filtered states are  *
  *  [dest] ptr to where the new states will be stored               *
  *  [ist] ptr to struct holding counters, indices, and arr lengths  *
  *  [parallel] struct holding parallelization options               *
  * outputs: void                                                    *
  ********************************************************************/

  long jms, jstate, jgrid, jgrid_real, jgrid_imag;
  long ngrid = ist->ngrid, cngrid = ist->complex_idx * ist->ngrid;
  omp_set_dynamic(0);
  omp_set_num_threads(parallel->nthreads);
#pragma omp parallel for private(jms, jstate, jgrid)
  for (jms = 0; jms < ist->mn_states_tot; jms++) {
    jstate = ist->complex_idx*ngrid*jms;
    for (jgrid = 0; jgrid < ngrid; jgrid++){
      jgrid_real = ist->complex_idx * jgrid;
      jgrid_imag = ist->complex_idx * jgrid + 1;
      //dest dn = (psi_tmp up)^*
      dest[jstate+jgrid_real+cngrid] = psitot[jstate+jgrid_real];
      dest[jstate+jgrid_imag+cngrid] = -1.0*psitot[jstate+jgrid_imag];

      //dest up = -(psi_tmp dn)^*
      dest[jstate+jgrid_real] = -1.0*psitot[jstate+jgrid_real+cngrid];
      dest[jstate+jgrid_imag] = psitot[jstate+jgrid_imag+cngrid];

    }
  } 

  return;
}

/*****************************************************************************/
