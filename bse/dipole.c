/****************************************************************************/

#include "fd.h"

/****************************************************************************/

void calc_elec_dipole(xyz_st *trans_dipole, double *psi_qp, double *eig_vals, grid_st *grid, index_st *ist, par_st *par, flag_st *flag){
  /*******************************************************************
  * This function computes the electric transition dipole matrix     *
  * matrix elements.                                                 *
  * inputs:                                                          *
  *  [trans_dipole] array to hold matrix elems in x, y, z direction  *
  *  [psi_qp] array holding all qp_basis states                      *
  *  [eig_vals] array holding the quasiparticle orbital energies     *
  *  [grid] grid_st instance holding values of all grid points       *
  *  [ist] ptr to counters, indices, and lengths                     *
  *  [par] ptr to par_st holding VBmin, VBmax... params              *
  * outputs: void                                                    *
  ********************************************************************/

  FILE *pf1, *pf2, *pf3, *pf; 
  long a, i, jx, jy, jz, jyz; 
  long jgridup_real, jgridup_imag, jgriddn_real, jgriddn_imag;
  long i_up_real, i_up_imag, i_dn_real, i_dn_imag, a_up_real, a_up_imag, a_dn_real, a_dn_imag;
  zomplex sumX, sumY, sumZ, tmp;
  double z, y, x, bohr_freq, osc_strength;
  
  // Output will be written to these files
  pf = fopen("OS0.dat" , "w"); 
  fprintf(pf, "i  a   sqrt(mu2)     Ea-Ei 	  f_osc       mu_x.re     mu_y.re     mu_z.re");
  if (flag->useSpinors){
    fprintf(pf, "     mu_x.im     mu_y.im     mu_z.im");
  }
  //pf1 = fopen("mu_x.dat", "w"); pf2 = fopen("mu_y.dat", "w"); pf3 = fopen("mu_z.dat", "w");

  // Make sure mu_x, mu_y and mu_z are zero to begin with
  // for (i = 0; i < ist->n_elecs*ist->n_holes; i++) mux[i].re = muy[i].re = muz[i].re = mux[i].im = muy[i].im = muz[i].im = 0.0;

  // Main computional work of function performed here - must loop over all electron-hole (i-a) pairs
  // compute $$\vec{\mu} = <a|\vec{r}|i> = \int dr \sum_\simga \psi_{a}^{*}(r,\sigma) \vec{r} \psi_i(r,\sigma)$$

  for (i = 0; i < ist->n_holes; i++){
    // fprintf(pf1,"\n"); fprintf(pf2,"\n"); fprintf(pf3,"\n");
    for (a = ist->lumo_idx; a < ist->lumo_idx+ist->n_elecs; a++) {
      sumX.re = sumY.re = sumZ.re = sumX.im = sumY.im = sumZ.im = 0.0;
      for (jz = 0; jz < ist->nz; jz++) {
        z = grid->z[jz];
        for (jy = 0; jy < ist->ny; jy++) {
          y = grid->y[jy];
          jyz = ist->nx * (ist->ny * jz + jy);
          for (jx = 0; jx < ist->nx; jx++) {
            x = grid->x[jx];
            // If the wavefunctions are scalar, only one of these indices is necessary (jgridup_real)
            // If psi_qp are spinors, then all four indices are necessary
            jgridup_real = ist->complex_idx * (jyz + jx);
            i_up_real = i*ist->nspinngrid*ist->complex_idx+jgridup_real;
            a_up_real = a*ist->nspinngrid*ist->complex_idx+jgridup_real;
            
            // compute matrix element with real, scalar wavefuntions
      	    tmp.re = psi_qp[a_up_real] * psi_qp[i_up_real];
            // sumX.re += tmp.re * x; sumY.re += tmp.re * y; sumZ.re += tmp.re * z;

            // If using spinors, add the down component of the integral (and the remaining imaginary piece from the up spin)
            if (1 == flag->useSpinors){
              // get indices for computing matrix elem with spinor
              jgridup_imag = ist->complex_idx * (jyz + jx) + 1;
              jgriddn_real = jgridup_real + ist->ngrid*ist->complex_idx;
              jgriddn_imag = jgridup_imag + ist->ngrid*ist->complex_idx;
              i_up_imag = i*ist->nspinngrid*ist->complex_idx+jgridup_imag; a_up_imag = a*ist->nspinngrid*ist->complex_idx+jgridup_imag;
              i_dn_real = i*ist->nspinngrid*ist->complex_idx+jgriddn_real; a_dn_real = a*ist->nspinngrid*ist->complex_idx+jgriddn_real;
              i_dn_imag = i*ist->nspinngrid*ist->complex_idx+jgriddn_imag; a_dn_imag = a*ist->nspinngrid*ist->complex_idx+jgriddn_imag;

              // REAL PART OF MATRIX ELEM CONTINUED
              // imaginary component of up spin contribution
              tmp.re += psi_qp[a_up_imag] * psi_qp[i_up_imag];
              // and down spin contribution
              tmp.re += psi_qp[a_dn_real] * psi_qp[i_dn_real]  + psi_qp[a_dn_imag] * psi_qp[i_dn_imag];

              // IMAG PART OF MATRIX ELEM
              tmp.im = psi_qp[a_up_imag] * psi_qp[i_up_real] + psi_qp[a_dn_imag] * psi_qp[i_dn_real]
                     - psi_qp[a_up_real] * psi_qp[i_up_imag] - psi_qp[a_dn_real] * psi_qp[i_dn_imag];
            
              sumX.im +=  tmp.im * x; sumY.im +=  tmp.im * y; sumZ.im +=  tmp.im * z;
            }
            //
            //
            sumX.re += tmp.re * x; sumY.re += tmp.re * y; sumZ.re += tmp.re * z;
            //
            //
          }
        }
      }
      trans_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].x_re = par->dv * sumX.re;  
      trans_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].y_re = par->dv * sumY.re; 
      trans_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].z_re = par->dv * sumZ.re; 
      osc_strength = sqr(par->dv) * (sqr(sumX.re) + sqr(sumY.re) + sqr(sumZ.re));
       
      if (1 == flag->useSpinors){
        trans_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].x_im = par->dv * sumX.im; 
        trans_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].y_im = par->dv * sumY.im;
        trans_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].z_im = par->dv * sumZ.im;
        osc_strength += sqr(par->dv) * (sqr(sumX.im) + sqr(sumY.im) + sqr(sumZ.im));
      }
      // fprintf (pf1,"%ld %ld %g %g\n",i,a,mux[i*ist->n_elecs+(a-ist->lumo_idx)].re, mux[i*ist->n_elecs+(a-ist->lumo_idx)].im);
      // fprintf (pf2,"%ld %ld %g %g\n",i,a,muy[i*ist->n_elecs+(a-ist->lumo_idx)].re, muy[i*ist->n_elecs+(a-ist->lumo_idx)].im);
      // fprintf (pf3,"%ld %ld %g %g\n",i,a,muz[i*ist->n_elecs+(a-ist->lumo_idx)].re, muz[i*ist->n_elecs+(a-ist->lumo_idx)].im);

      bohr_freq = eig_vals[a] - eig_vals[i];
      fprintf(pf,"\n%ld % ld  %.8f %.12f % .8f % .8f % .8f % .8f", i, a, sqrt(osc_strength), bohr_freq, (2.0/3.0)*bohr_freq*osc_strength, trans_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].x_re,
        trans_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].y_re, trans_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].z_re);

      if (1 == flag->useSpinors){
        fprintf(pf," %.8f %.8f %.8f", trans_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].x_im, trans_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].y_im, trans_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].z_im);
      }
    }
  }
  fclose(pf); // fclose(pf1); fclose(pf2); fclose(pf3);

  return;
}

/****************************************************************************/
// This function calculates the magnetic dipole matrix elements between the
// single-particle electron (a) and hole (i) states: <psi_a|m|psi_i>
// where m = -1/2*L = -1/2*r x p where x is the cross product.


void calc_mag_dipole(xyz_st *mag_dipole, double *psi_qp, double *eig_vals, grid_st *grid, index_st *ist, par_st *par, flag_st *flag)
{
  FILE *pf1, *pf2, *pf3, *pf; 
  long a, i, jx, jy, jz, jyz; 
  long jgrid,jgridup_real, jgridup_imag, jgriddn_real, jgriddn_imag;
  long i_up_real, i_up_imag, i_dn_real, i_dn_imag, a_up_real, a_up_imag, a_dn_real, a_dn_imag;
  long istate_idx, astate_idx;
  zomplex sumX, sumY, sumZ, tmp;
  double ms, bohr_freq, osc_strength;
  zomplex *Lxpsi, *Lypsi, *Lzpsi;
  zomplex *psi_i, *psi_a;

  // Allocate memory
  psi_i = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));;
	psi_a = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));;
	
	Lxpsi = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));
	Lypsi = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));
	Lzpsi = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));
  
  // Output will be written to these files
  pf = fopen("M0.dat" , "w"); 
  fprintf(pf, "  i   a    sqrt(ms)       Ea-Ei         m_x.re      m_x.im      m_y.re      m_y.im      m_z.re      m_z.im\n");
  // if (flag->useSpinors){
  //   fprintf(pf, "     m_x.im     m_y.im     m_z.im\n");
  // }
  
  
  // Output will be written to these files
  pf1 = fopen("mx.dat", "w"); 
  pf2 = fopen("my.dat", "w"); 
  pf3 = fopen("mz.dat", "w");

  // Make sure mx, my and mz are zero to begin with
  for (i = 0; i < ist->n_xton; i++) {
    mag_dipole[i].x_re = mag_dipole[i].y_re = mag_dipole[i].z_re = 0.0;
    mag_dipole[i].x_im = mag_dipole[i].y_im = mag_dipole[i].z_im = 0.0;
  }
  
  // Main computional work of function performed here - must loop over all electron-hole (i-a) pairs
  for (i = 0; i < ist->n_holes; i++) {
    fprintf(pf1,"\n"),fprintf(pf2,"\n"),fprintf(pf3,"\n");
    
    istate_idx = i * ist->nspinngrid * ist->complex_idx;
    // Compute the angular momentum of this qp state
    memcpy(&psi_i[0], &psi_qp[istate_idx], ist->nspinngrid*sizeof(psi_i[0]));

    //spin up part
		l_operator(&Lxpsi[0], &Lypsi[0], &Lzpsi[0], &psi_i[0], grid, ist, par);
		//spin dn part
		l_operator(&Lxpsi[ist->ngrid], &Lypsi[ist->ngrid], &Lzpsi[ist->ngrid], &psi_i[ist->ngrid], grid, ist, par);


    for (a = ist->lumo_idx; a < ist->lumo_idx+ist->n_elecs; a++) {
      astate_idx = a * ist->complex_idx * ist->nspinngrid;
      
      // Compy state a into psi_a
      memcpy(&psi_a[0], &psi_qp[astate_idx], ist->nspinngrid*sizeof(psi_i[0]));

      sumX.re = sumY.re = sumZ.re = 0.0;
      sumX.im = sumY.im = sumZ.im = 0.0;
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
        sumX.re += psi_a[jgrid].re * Lxpsi[jgrid].re - psi_a[jgrid].im * Lxpsi[jgrid].im;
        sumX.im += psi_a[jgrid].re * Lxpsi[jgrid].im - psi_a[jgrid].im * Lxpsi[jgrid].re;

        sumY.re += psi_a[jgrid].re * Lypsi[jgrid].re - psi_a[jgrid].im * Lypsi[jgrid].im;
        sumY.im += psi_a[jgrid].re * Lypsi[jgrid].im - psi_a[jgrid].im * Lypsi[jgrid].re;

        sumZ.re += psi_a[jgrid].re * Lzpsi[jgrid].re - psi_a[jgrid].im * Lzpsi[jgrid].im;
        sumZ.im += psi_a[jgrid].re * Lzpsi[jgrid].im - psi_a[jgrid].im * Lzpsi[jgrid].re;
      }

      mag_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].x_re =  0.5 * sumX.re;
      mag_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].y_re =  0.5 * sumY.re;
      mag_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].z_re =  0.5 * sumZ.re;

      mag_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].x_im =  0.5 * sumX.im;
      mag_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].y_im =  0.5 * sumY.im;
      mag_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].z_im =  0.5 * sumZ.im;

      fprintf(pf1,"%ld %ld %g %g\n",i, a, sumX.re, sumX.im);
      fprintf(pf2,"%ld %ld %g %g\n",i, a, sumY.re, sumY.im);
      fprintf(pf3,"%ld %ld %g %g\n",i, a, sumZ.re, sumZ.im);

      ms = \
      ( sqr(mag_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].x_re) + sqr(mag_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].x_im) 
      + sqr(mag_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].y_re) + sqr(mag_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].y_im)
      + sqr(mag_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].z_re) + sqr(mag_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].z_im) );

      bohr_freq = eig_vals[a] - eig_vals[i];

      fprintf(pf,"%3ld %3ld  %+.8f %+.12f %+.8f %+.8f %+.8f %+.8f %+.8f %+.8f\n", i, a, sqrt(ms), bohr_freq,
        mag_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].x_re,
        mag_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].x_im,
        mag_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].y_re,
        mag_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].y_im,
        mag_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].z_re,
        mag_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].z_im);
    }
  }
  fclose(pf); fclose(pf1); fclose(pf2); fclose(pf3);

  // Free dynamically allocated memory
  free(psi_i); free(psi_a); 
  free(Lxpsi); free(Lypsi); free(Lzpsi);

  return;
}

// /****************************************************************************/

// void rotational_strength(double *rs, double *mux, double *muy, double *muz, double *mx, 
//   double *my, double *mz, double *eig_vals, index_st ist) {
//   FILE *pf;
//   int i, a, index;

//   pf = fopen("rs0.dat", "w");
//   for (i = 0; i < ist->n_holes; i++) {
//     for (a = ist->lumo_idx; a < ist->lumo_idx+ist->n_elecs; a++) {
//       index = i*ist->n_elecs + (a-ist->lumo_idx);
//       rs[index] = mux[index]*mx[index] + muy[index]*my[index] + muz[index]*mz[index];
//       fprintf(pf, "%d %d %.12f  %.16f\n", i, a, eig_vals[a] - eig_vals[i], rs[index]);
//     }
//   }
//   fclose(pf);

//   return;
// }

/****************************************************************************/
