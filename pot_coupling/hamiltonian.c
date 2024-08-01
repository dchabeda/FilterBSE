/*****************************************************************************/
// File contains the main functions for calcuating the action of the hamiltonian on 
// a wavefunction including non-local spin-orbit contributions 

#include "pot_coupling.h"

/*****************************************************************************/

// This calculates the total action of Vloc + Vnonloc + Vso|psi_tmp>
void potential(zomplex *psi_out, zomplex *psi_tmp, double *pot_local, nlc_st *nlc, long *nl, index_st *ist,
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
  memcpy(&psi_tmp[0], &psi_out[0], ist->nspinngrid*sizeof(psi_tmp[0]));

  if(flag->SO==1){
    // Calculate |psi_out> = V_SO|psi_tmp>
    spin_orbit_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par);
  }
  if (flag->NL == 1){
    // Calculate |psi_out> += V_NL|psi_tmp>
    nonlocal_proj_pot(psi_out, psi_tmp, nlc, nl, ist, par);
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


void spin_orbit_proj_pot(zomplex *psi_out, zomplex *psi_tmp, nlc_st *nlc, long *nl, index_st *ist, par_st *par){
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

  zomplex proj;
  long jatom, NL_gridpt, index1, r, index2, r_p;
  int iproj, spin_p, m_p, spin, m;
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

  for ( jatom = 0; jatom < ist->n_NL_atoms; jatom++){

    // sprintf(str, "atom_%ld_Vr%ld.dat", jatom, idx);
    // pf = fopen(str, "w");
    // rho = (zomplex *) calloc(ist->nspinngrid, sizeof(rho[0])); // allocate memory and initialize all to zero
    // r_p_idx = nlc[jatom * ist->n_NL_gridpts + idx].jxyz; // get the idx-th gridpt around atom jatom
    // rho[r_p_idx].re = 1.0 / sqrt(par->dv); // Make wavefunction a delta function at r_p
    // rho[r_p_idx+ist->nspinngrid].re = 1.0 / sqrt(par->dv);

    // memcpy(&psi_tmp[0], &rho[0], ist->nspinngrid*sizeof(psi_tmp[0]));

    for ( iproj = 0; iproj < ist->nproj; iproj++){
      for ( spin_p = 0; spin_p < 2; spin_p++){
        for ( m_p = 0; m_p < 3; m_p++){
          proj.re = proj.im = 0.00;
          for ( NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++){
            index1 = jatom * ist->n_NL_gridpts + NL_gridpt;
            r = nlc[index1].jxyz + (ist->ngrid) * spin_p;
            
            //weird signs b/c of Y_{lm}^*
            //TODO: I've checked these iprojs against the plane wave way to generate iprojs and they match..
            //When I use these iprojs with the LdotS I've checked then i get the same as the RS energy. 
            
            proj.re += psi_tmp[r].re * nlc[index1].y1[m_p].re * nlc[index1].proj[iproj];
            proj.re += psi_tmp[r].im * nlc[index1].y1[m_p].im * nlc[index1].proj[iproj];
           
            proj.im += psi_tmp[r].im * nlc[index1].y1[m_p].re * nlc[index1].proj[iproj];
            proj.im -= psi_tmp[r].re * nlc[index1].y1[m_p].im * nlc[index1].proj[iproj];
            
          }
          proj.re *= par->dv;
          proj.im *= par->dv;
           
          for (spin = 0; spin<2; spin++){
            for (m = 0; m < 3; m++){
              //get L_{m,m'}\cdot S_{s,s'}*P_{n,m,s} = PLS_{n,m,m',s,s'}
              zomplex LS, PLS;

              LS.re = LS.im = 0.00;
              PLS.re = PLS.im = 0.00;
              //TODO: could already have LS 6x6 (real?) matrix loaded... 
              //checked ordering of m and m_p, printed LdotS and checked against ref. 
              LS.re += (  Lx[m][m_p].re * Sx[spin][spin_p].re - Lx[m][m_p].im * Sx[spin][spin_p].im);
              LS.im += (  Lx[m][m_p].re * Sx[spin][spin_p].im + Lx[m][m_p].im * Sx[spin][spin_p].re);

              LS.re += (  Ly[m][m_p].re * Sy[spin][spin_p].re - Ly[m][m_p].im * Sy[spin][spin_p].im);
              LS.im += (  Ly[m][m_p].re * Sy[spin][spin_p].im + Ly[m][m_p].im * Sy[spin][spin_p].re);

              LS.re += (  Lz[m][m_p].re * Sz[spin][spin_p].re - Lz[m][m_p].im * Sz[spin][spin_p].im);
              LS.im += (  Lz[m][m_p].re * Sz[spin][spin_p].im + Lz[m][m_p].im * Sz[spin][spin_p].re);              
              //printf("m:%i m':%i s:%i s':%i  LS: %f+i*%f\n", m, m_p, spin, spin_p, LS.re, LS.im);

              PLS.re = LS.re * proj.re - LS.im * proj.im;
              PLS.im = LS.re * proj.im + LS.im * proj.re;

              //TODO: check why signs funky in update psi_out 214, 215, 217, 218;
              for (NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++){
                index2 = jatom * ist->n_NL_gridpts + NL_gridpt;
                r_p = nlc[index2].jxyz + (ist->ngrid)*spin;

                psi_out[r_p].re += nlc[index2].proj[iproj] * nlc[index2].y1[m].re * PLS.re;
                psi_out[r_p].re -= nlc[index2].proj[iproj] * nlc[index2].y1[m].im * PLS.im;

                psi_out[r_p].im += nlc[index2].proj[iproj] * nlc[index2].y1[m].re * PLS.im;
                psi_out[r_p].im += nlc[index2].proj[iproj] * nlc[index2].y1[m].im * PLS.re;
                
                
                /*
                soPE+= psi_tmp[r].re *((nlc[index1].rsProj[iproj] * nlc[index1].y1[m].re * PLS.re)
                                              -(nlc[index1].rsProj[iproj] * nlc[index1].y1[m].im * PLS.im))*rsParams.dV;
                soPE+= psi_tmp[r].im * ((nlc[index1].rsProj[iproj] * nlc[index1].y1[m].re * PLS.im)
                                              +(nlc[index1].rsProj[iproj] * nlc[index1].y1[m].im * PLS.re))*rsParams.dV;
                */
              }
            } 
          }
        }
      }
    }
    // for (NL_gridpt = 0; NL_gridpt < nl[jatom]; NL_gridpt++){
    //   index2 = jatom * ist->n_NL_gridpts + NL_gridpt;
    //   r_p = nlc[index2].jxyz + (ist->ngrid)*spin;
    //   fprintf(pf, "%ld %ld %lg %lg\n", r_p, idx, psi_out[r_p].re, psi_out[r_p].im);
    // }
    // fclose(pf);
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

  zomplex proj;
  long jatom, NL_gridpt, index1, r, index2, r_p;
  int iproj, spin, m;

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
