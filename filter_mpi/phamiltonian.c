#include "fd.h"

/*****************************************************************************/
// File contains the main functions for calcuating the action of the hamiltonian on 
// a wavefunction including non-local spin-orbit contributions 
// parallelized with OpenMP

/*****************************************************************************/

void p_hamiltonian(zomplex *psi_out, zomplex *psi_tmp, double *pot_local, nlc_st *nlc, long *nl, double *ksqr,
  index_st *ist, par_st *par, flag_st *flag, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi){
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

  // Calculate the action of the potential on the wavefunction: |psi_out> = V|psi_tmp>
  p_potential(psi_out, psi_tmp, pot_local, nlc, nl, ist, par, flag);

  return;
}

/*****************************************************************************/

// This calculates the total action of Vloc + Vnonloc + Vso|psi_tmp>
void p_potential(zomplex *psi_out, zomplex *psi_tmp, double *pot_local, nlc_st *nlc, long *nl, index_st *ist,
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
    p_spin_orbit_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par);
  }
  if (flag->NL == 1){
    // Calculate |psi_out> += V_NL|psi_tmp>
    p_nonlocal_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par);
  }
  
  // Calculate the action of the local potential energy part of the Hamiltonian on psi_tmp
  for (jspin = 0; jspin < ist->nspin; jspin++){
    for (j = 0; j < ist->ngrid; j++) {
      jtmp = ist->ngrid * jspin + j ; // generalized indexing to handle spinors or spinless wavefuncs

      psi_out[jtmp].re += (pot_local[j] * psi_tmp[jtmp].re);
      psi_out[jtmp].im += (pot_local[j] * psi_tmp[jtmp].im);
    }
  }

  return;
}


void p_spin_orbit_proj_pot(zomplex *psi_out, zomplex *psi_tmp, nlc_st *nlc, long *nl, index_st *ist, par_st *par){
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
  //TODO: Find way to get the L and S to be not stack'd for mem reasons???
  zomplex Lx[3][3], Ly[3][3], Lz[3][3];
  zomplex Sx[2][2], Sy[2][2], Sz[2][2];
  double sq2 = sqrt(0.5);
  
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

  // Added by Daniel C to test matrix elements of pseudopotential
  // FILE *pf;
  // char str[40];
  // zomplex *rho;
  // long r_p_idx, idx = 1000;
  //
  long NL_gridpt, index1, r, index2, r_p;
  int iproj, spin_p, m_p, spin, m;
  
  #pragma omp parallel private(jatom, NL_gridpt, index1, r, index2, r_p, iproj, spin, spin_p, m, m_p, Lx, Ly, Lz, Sx, Sy, Sz)
{
    zomplex *psi_out_private = calloc(ist->nspinngrid, sizeof(zomplex));
    if (!psi_out_private) {
        fprintf(stderr, "Failed to allocate memory for thread-local psi_out\n");
        exit(1);
    }

    #pragma omp for
    for (jatom = 0; jatom < ist->n_NL_atoms; jatom++) {
      zomplex proj;
        for (iproj = 0; iproj < ist->nproj; iproj++) {
            for (spin_p = 0; spin_p < 2; spin_p++) {
                for (m_p = 0; m_p < 3; m_p++) {
                    proj.re = proj.im = 0.0;
                    for (NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++) {
                        index1 = jatom * ist->n_NL_gridpts + NL_gridpt;
                        r = nlc[index1].jxyz + (ist->ngrid) * spin_p;

                        proj.re += psi_tmp[r].re * nlc[index1].y1[m_p].re * nlc[index1].proj[iproj];
                        proj.re += psi_tmp[r].im * nlc[index1].y1[m_p].im * nlc[index1].proj[iproj];
                        proj.im += psi_tmp[r].im * nlc[index1].y1[m_p].re * nlc[index1].proj[iproj];
                        proj.im -= psi_tmp[r].re * nlc[index1].y1[m_p].im * nlc[index1].proj[iproj];
                    }

                    proj.re *= par->dv;
                    proj.im *= par->dv;
                    
                    
                    for (spin = 0; spin < 2; spin++) {
                        for (m = 0; m < 3; m++) {
                          zomplex LS, PLS;
                            LS.re = LS.im = 0.0;
                            LS.re += Lx[m][m_p].re * Sx[spin][spin_p].re - Lx[m][m_p].im * Sx[spin][spin_p].im;
                            LS.im += Lx[m][m_p].re * Sx[spin][spin_p].im + Lx[m][m_p].im * Sx[spin][spin_p].re;
                            LS.re += Ly[m][m_p].re * Sy[spin][spin_p].re - Ly[m][m_p].im * Sy[spin][spin_p].im;
                            LS.im += Ly[m][m_p].re * Sy[spin][spin_p].im + Ly[m][m_p].im * Sy[spin][spin_p].re;
                            LS.re += Lz[m][m_p].re * Sz[spin][spin_p].re - Lz[m][m_p].im * Sz[spin][spin_p].im;
                            LS.im += Lz[m][m_p].re * Sz[spin][spin_p].im + Lz[m][m_p].im * Sz[spin][spin_p].re;

                            PLS.re = LS.re * proj.re - LS.im * proj.im;
                            PLS.im = LS.re * proj.im + LS.im * proj.re;

                            for (NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++) {
                                index2 = jatom * ist->n_NL_gridpts + NL_gridpt;
                                r_p = nlc[index2].jxyz + (ist->ngrid) * spin;

                                psi_out_private[r_p].re += nlc[index2].proj[iproj] * nlc[index2].y1[m].re * PLS.re;
                                psi_out_private[r_p].re -= nlc[index2].proj[iproj] * nlc[index2].y1[m].im * PLS.im;
                                psi_out_private[r_p].im += nlc[index2].proj[iproj] * nlc[index2].y1[m].re * PLS.im;
                                psi_out_private[r_p].im += nlc[index2].proj[iproj] * nlc[index2].y1[m].im * PLS.re;
                            }
                        }
                    }
                }
            }
        }
    }

  
    #pragma omp critical
    {
        for (int i = 0; i < ist->nspinngrid; i++) {
            psi_out[i].re += psi_out_private[i].re;
            psi_out[i].im += psi_out_private[i].im;
        }
    }

    free(psi_out_private);
}



// #pragma omp parallel for private(jatom, NL_gridpt, index1, r, index2, r_p, iproj, spin, spin_p, m, m_p, proj, Lx, Ly, Lz, Sx, Sy, Sz)
//   for ( jatom = 0; jatom < ist->n_NL_atoms; jatom++){
    
//     // sprintf(str, "atom_%ld_Vr%ld.dat", jatom, idx);
//     // pf = fopen(str, "w");
//     // rho = (zomplex *) calloc(ist->nspinngrid, sizeof(rho[0])); // allocate memory and initialize all to zero
//     // r_p_idx = nlc[jatom * ist->n_NL_gridpts + idx].jxyz; // get the idx-th gridpt around atom jatom
//     // rho[r_p_idx].re = 1.0 / sqrt(par->dv); // Make wavefunction a delta function at r_p
//     // rho[r_p_idx+ist->nspinngrid].re = 1.0 / sqrt(par->dv);

//     // memcpy(&psi_tmp[0], &rho[0], ist->nspinngrid*sizeof(psi_tmp[0]));

//     for ( iproj = 0; iproj < ist->nproj; iproj++){
//       for ( spin_p = 0; spin_p < 2; spin_p++){
//         for ( m_p = 0; m_p < 3; m_p++){
//           proj.re = proj.im = 0.00;
//           for ( NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++){
//             index1 = jatom * ist->n_NL_gridpts + NL_gridpt;
//             r = nlc[index1].jxyz + (ist->ngrid) * spin_p;
            
//             //weird signs b/c of Y_{lm}^*
//             //TODO: I've checked these iprojs against the plane wave way to generate iprojs and they match..
//             //When I use these iprojs with the LdotS I've checked then i get the same as the RS energy. 
            
//             proj.re += psi_tmp[r].re * nlc[index1].y1[m_p].re * nlc[index1].proj[iproj];
//             proj.re += psi_tmp[r].im * nlc[index1].y1[m_p].im * nlc[index1].proj[iproj];
           
//             proj.im += psi_tmp[r].im * nlc[index1].y1[m_p].re * nlc[index1].proj[iproj];
//             proj.im -= psi_tmp[r].re * nlc[index1].y1[m_p].im * nlc[index1].proj[iproj];
            
//           }
          
//           proj.re *= par->dv;
//           proj.im *= par->dv;
            
//           for (spin = 0; spin<2; spin++){
//             for (m = 0; m < 3; m++){
//               //get L_{m,m'}\cdot S_{s,s'}*P_{n,m,s} = PLS_{n,m,m',s,s'}
//               zomplex LS, PLS;

//               LS.re = LS.im = 0.00;
//               PLS.re = PLS.im = 0.00;
//               //TODO: could already have LS 6x6 (real?) matrix loaded... 
//               //checked ordering of m and m_p, printed LdotS and checked against ref. 
//               LS.re += (  Lx[m][m_p].re * Sx[spin][spin_p].re - Lx[m][m_p].im * Sx[spin][spin_p].im);
//               LS.im += (  Lx[m][m_p].re * Sx[spin][spin_p].im + Lx[m][m_p].im * Sx[spin][spin_p].re);

//               LS.re += (  Ly[m][m_p].re * Sy[spin][spin_p].re - Ly[m][m_p].im * Sy[spin][spin_p].im);
//               LS.im += (  Ly[m][m_p].re * Sy[spin][spin_p].im + Ly[m][m_p].im * Sy[spin][spin_p].re);

//               LS.re += (  Lz[m][m_p].re * Sz[spin][spin_p].re - Lz[m][m_p].im * Sz[spin][spin_p].im);
//               LS.im += (  Lz[m][m_p].re * Sz[spin][spin_p].im + Lz[m][m_p].im * Sz[spin][spin_p].re);              
//               //if (parallel->mpi_rank == 0) printf("m:%i m':%i s:%i s':%i  LS: %f+i*%f\n", m, m_p, spin, spin_p, LS.re, LS.im);

//               PLS.re = LS.re * proj.re - LS.im * proj.im;
//               PLS.im = LS.re * proj.im + LS.im * proj.re;

//               //TODO: check why signs funky in update psi_out 214, 215, 217, 218;
//               for (NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++){
//                 index2 = jatom * ist->n_NL_gridpts + NL_gridpt;
//                 r_p = nlc[index2].jxyz + (ist->ngrid)*spin;

//                 psi_out[r_p].re += nlc[index2].proj[iproj] * nlc[index2].y1[m].re * PLS.re;
//                 psi_out[r_p].re -= nlc[index2].proj[iproj] * nlc[index2].y1[m].im * PLS.im;

//                 psi_out[r_p].im += nlc[index2].proj[iproj] * nlc[index2].y1[m].re * PLS.im;
//                 psi_out[r_p].im += nlc[index2].proj[iproj] * nlc[index2].y1[m].im * PLS.re;
                
                
//                 /*
//                 soPE+= psi_tmp[r].re *((nlc[index1].rsProj[iproj] * nlc[index1].y1[m].re * PLS.re)
//                                               -(nlc[index1].rsProj[iproj] * nlc[index1].y1[m].im * PLS.im))*rsParams.dV;
//                 soPE+= psi_tmp[r].im * ((nlc[index1].rsProj[iproj] * nlc[index1].y1[m].re * PLS.im)
//                                               +(nlc[index1].rsProj[iproj] * nlc[index1].y1[m].im * PLS.re))*rsParams.dV;
//                 */
//               }
//             } 
//           }
//         }
//       }
//     }
    
  // }

  return;
}

/*****************************************************************************/

void p_nonlocal_proj_pot(zomplex *psi_out, zomplex *psi_tmp, nlc_st *nlc, long *nl, index_st *ist, par_st *par){
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

//   long jatom;
//   long NL_gridpt, index1, r, index2, r_p;
//   int iproj, spin, m;
//   zomplex proj;
  long jatom, NL_gridpt, index1, index2, r, r_p;
  int iproj, spin, m;
  zomplex proj;

  #pragma omp parallel
  {
      // Thread-local storage for psi_out
      zomplex *psi_out_private = (zomplex *)calloc(ist->ngrid, sizeof(zomplex));
      if (!psi_out_private) {
          fprintf(stderr, "Error: Unable to allocate thread-local storage.\n");
          exit(EXIT_FAILURE);
      }

      #pragma omp for private(jatom, NL_gridpt, index1, index2, r, r_p, iproj, spin, m, proj)
      for (jatom = 0; jatom < ist->n_NL_atoms; jatom++) {
          long jatom_offset = jatom * ist->n_NL_gridpts;

          for (iproj = 0; iproj < ist->nproj; iproj++) {
              for (spin = 0; spin < 2; spin++) {
                  for (m = 0; m < 3; m++) {
                      proj.re = proj.im = 0.0;

                      // First loop over NL_gridpt
                      for (NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++) {
                          index1 = jatom_offset + NL_gridpt;
                          r = nlc[index1].jxyz + ist->ngrid * spin;

                          double psi_re = psi_tmp[r].re;
                          double psi_im = psi_tmp[r].im;
                          double y1_re = nlc[index1].y1[m].re;
                          double y1_im = nlc[index1].y1[m].im;
                          double NL_proj = nlc[index1].NL_proj[iproj];

                          proj.re += psi_re * y1_re * NL_proj + psi_im * y1_im * NL_proj;
                          proj.im += psi_im * y1_re * NL_proj - psi_re * y1_im * NL_proj;
                      }

                      // Apply scaling factors
                      double proj_sign = nlc[jatom_offset].NL_proj_sign[iproj];
                      proj.re *= proj_sign * par->dv;
                      proj.im *= proj_sign * par->dv;

                      // Second loop over NL_gridpt
                      for (NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++) {
                          index2 = jatom_offset + NL_gridpt;
                          r_p = nlc[index2].jxyz + ist->ngrid * spin;

                          double y1_re = nlc[index2].y1[m].re;
                          double y1_im = nlc[index2].y1[m].im;
                          double NL_proj = nlc[index2].NL_proj[iproj];

                          double update_re = NL_proj * (y1_re * proj.re - y1_im * proj.im);
                          double update_im = NL_proj * (y1_re * proj.im + y1_im * proj.re);

                          // Update thread-local psi_out_private
                          psi_out_private[r_p].re += update_re;
                          psi_out_private[r_p].im += update_im;
                      }
                  }
              }
          }
      }

      // Reduction step: Combine thread-local psi_out_private into the global psi_out
      #pragma omp critical
      {
          for (long i = 0; i < ist->ngrid; i++) {
              psi_out[i].re += psi_out_private[i].re;
              psi_out[i].im += psi_out_private[i].im;
          }
      }

      // Free thread-local storage
      free(psi_out_private);
  }


// #pragma omp parallel for private(jatom, NL_gridpt, index1, r, index2, r_p, iproj, spin, m, proj)
//   for ( jatom =0; jatom < ist->n_NL_atoms; jatom++){
//     for ( iproj = 0; iproj < ist->nproj; iproj++){
//       for ( spin = 0; spin < 2; spin++){
//         for ( m = 0; m < 3; m++){
//           proj.re = proj.im = 0.00;
//           for ( NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++){
//             index1 = jatom * ist->n_NL_gridpts + NL_gridpt;
//             r = nlc[index1].jxyz + (ist->ngrid)*spin;
            
//             //weird signs b/c of Y_{lm}^*
//             //TODO: I've checked these iprojs against the plane wave way to generate iprojs and they match..
//             //When I use these iprojs with the LdotS I've checked then i get the same as the RS energy. 
            
//             proj.re += psi_tmp[r].re * nlc[index1].y1[m].re * nlc[index1].NL_proj[iproj];
//             proj.re += psi_tmp[r].im * nlc[index1].y1[m].im * nlc[index1].NL_proj[iproj];
           
//             proj.im += psi_tmp[r].im * nlc[index1].y1[m].re * nlc[index1].NL_proj[iproj];
//             proj.im -= psi_tmp[r].re * nlc[index1].y1[m].im * nlc[index1].NL_proj[iproj];
            
//           }
//           proj.re *= nlc[index1].NL_proj_sign[iproj];
//           proj.im *= nlc[index1].NL_proj_sign[iproj];   
          
//           proj.re *= par->dv;
//           proj.im *= par->dv;
          
//           for (NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++){
//             index2 = jatom * ist->n_NL_gridpts + NL_gridpt;
//             r_p = nlc[index2].jxyz + (ist->ngrid)*spin;

//             psi_out[r_p].re += nlc[index2].NL_proj[iproj] * nlc[index2].y1[m].re * proj.re;
//             psi_out[r_p].re -= nlc[index2].NL_proj[iproj] * nlc[index2].y1[m].im * proj.im;

//             psi_out[r_p].im += nlc[index2].NL_proj[iproj] * nlc[index2].y1[m].re * proj.im;
//             psi_out[r_p].im += nlc[index2].NL_proj[iproj] * nlc[index2].y1[m].im * proj.re;
            

//           }
//         }
//       }
//     }
//   }

  return;
}

/*****************************************************************************/

