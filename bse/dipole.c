/****************************************************************************/

#include "fd.h"

/****************************************************************************/

void calc_electric_dipole(xyz_st *trans_dipole, double *psi_qp, double *eig_vals, grid_st *grid, index_st *ist, par_st *par, flag_st *flag){
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
  fprintf(pf, " i   a   sqrt(mu2)       Ea-Ei 	  f_osc       mu_x.re     mu_y.re     mu_z.re");
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
            i_up_real = i*ist->nspinngrid+jgridup_real;
            a_up_real = a*ist->nspinngrid+jgridup_real;
            
            // compute matrix element with real, scalar wavefuntions
      	    tmp.re = psi_qp[a_up_real] * psi_qp[i_up_real];
            sumX.re += tmp.re * x; sumY.re += tmp.re * y; sumZ.re += tmp.re * z;

            // If using spinors, add the down component of the integral (and the remaining imaginary piece from the up spin)
            if (1 == flag->useSpinors){
              // get indices for computing matrix elem with spinor
              jgridup_imag = ist->complex_idx * (jyz + jx) + 1;
              jgriddn_real = jgridup_real + ist->ngrid;
              jgriddn_imag = jgridup_imag + ist->ngrid;
              i_up_imag = i*ist->nspinngrid+jgridup_imag; a_up_imag = a*ist->nspinngrid+jgridup_imag;
              i_dn_real = i*ist->nspinngrid+jgriddn_real; a_dn_real = a*ist->nspinngrid+jgriddn_real;
              i_dn_imag = i*ist->nspinngrid+jgriddn_imag; a_dn_imag = a*ist->nspinngrid+jgriddn_imag;

              // REAL PART OF MATRIX ELEM CONTINUED
              // imaginary component of up spin and down spin contribution
              tmp.re += psi_qp[a_up_imag] * psi_qp[i_up_imag];
              // and down spin contribution
              tmp.re += + psi_qp[a_dn_real] * psi_qp[i_dn_real]  + psi_qp[a_dn_imag] * psi_qp[i_dn_imag];

              // IMAG PART OF MATRIX ELEM
              tmp.im = (- psi_qp[a_up_imag]) * psi_qp[i_up_real] + psi_qp[a_up_real] * psi_qp[i_up_imag];
              tmp.im += (- psi_qp[a_dn_imag]) * psi_qp[i_dn_real] + psi_qp[a_dn_real] * psi_qp[i_dn_imag];
            
              sumX.re += tmp.re * x; sumY.re += tmp.re * y; sumZ.re += tmp.re * z;
              sumX.im +=  tmp.im * x; sumY.im +=  tmp.im * y; sumZ.im +=  tmp.im * z;
            }
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
      fprintf(pf,"\n%3ld %3ld %.8f %.12f % .8f % .8f % .8f % .8f", i, a, sqrt(osc_strength), bohr_freq, (2.0/3.0)*bohr_freq*osc_strength, trans_dipole[i*ist->n_elecs+(a-ist->lumo_idx)].x_re,
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


// void mag_dipole(double *grid->x, double *grid->y, double *grid->z, double *psi_qp, double *mx, double *my, double *mz, 
//   double *eig_vals, fftw_plan_loc *planfw,fftw_plan_loc *planbw,fftw_complex *fftwpsi, index_st ist, par_st par)
// {
//   FILE *pf1, *pf2, *pf3, *pf; 
//   long a, i, jx, jy, jz, jgrid, jyz; 
//   double sumX, sumY, sumZ;
//   double z, x, y, bohr_freq, ms, tmp;
//   double *kx, *ky, *kz;
//   double *kindex;
//   double *psidx, *psidy, *psidz;
//   zomplex *cpsi;
  
//   // Output will be written to these files
//   pf = fopen("M0.dat" , "w"); 
//   pf1 = fopen("mx.dat", "w"); 
//   pf2 = fopen("my.dat", "w"); 
//   pf3 = fopen("mz.dat", "w");

//   // Make sure mx, my and mz are zero to begin with
//   for (i = 0; i < ist->n_elecs*ist->n_holes; i++) mx[i] = my[i] = mz[i] = 0.0;

//   // Allocate memory 
//   if ((kx = (double *) calloc(ist->nx, sizeof(double)))==NULL) nerror("kx");
//   if ((ky = (double *) calloc(ist->ny, sizeof(double)))==NULL) nerror("ky");
//   if ((kz = (double *) calloc(ist->nz, sizeof(double)))==NULL) nerror("kz");
//   if ((kindex = (double *) calloc(3*ist->ngrid, sizeof(double)))==NULL) nerror("kindex");
//   if ((psidx = (double *) calloc(ist->ngrid, sizeof(double)))==NULL) nerror("psidx");
//   if ((psidy = (double *) calloc(ist->ngrid, sizeof(double)))==NULL) nerror("psidy");
//   if ((psidz = (double *) calloc(ist->ngrid, sizeof(double)))==NULL) nerror("psidz");

//   // Calculate kx 
//   for (jx = 1; jx <= ist->nx / 2; jx++) {
//     kx[jx] = (double)(jx) * par->dkx * ist->nx_1 * ist->ny_1 * ist->nz_1;
//     kx[ist->nx-jx] = -kx[jx];
//   }
//   // Calculate ky 
//   for (jy = 1; jy <= ist->ny / 2; jy++) {
//     ky[jy] = (double)(jy) * par->dky * ist->nx_1 * ist->ny_1 * ist->nz_1;
//     ky[ist->ny-jy] = -ky[jy];
//   }
//   // Calculate kz 
//   for (jz = 1; jz <= ist->nz / 2; jz++) {
//     kz[jz] = (double)(jz) * par->dkz * ist->nx_1 * ist->ny_1 * ist->nz_1;
//     kz[ist->nz-jz] = -kz[jz];
//   }  
//   for (jz = 0; jz < ist->nz; jz++) {
//     for (jy = 0; jy < ist->ny; jy++) {
//       jyz = ist->nx * (ist->ny * jz + jy);
//       for (jx = 0; jx < ist->nx; jx++) {   
//         jgrid = jyz + jx;
//         kindex[3*jgrid]   = kx[jx]; 
//         kindex[3*jgrid+1] = ky[jy];
//         kindex[3*jgrid+2] = kz[jz];
//       } 
//     }
//   }

//   // Allocate memory for Fourier transforms
//   if ((cpsi = (zomplex *) calloc(ist->ngrid, sizeof(zomplex)))==NULL) nerror("cpsi");

//   // Main computional work of function performed here - must loop over all electron-hole (i-a) pairs
//   for (i = 0; i < ist->n_holes; i++, fprintf(pf1,"\n"),fprintf(pf2,"\n"),fprintf(pf3,"\n")) {
//     // Fourier transform the hole wavefunctions
//     for (jgrid = 0; jgrid < ist->ngrid; jgrid++) {
//       cpsi[jgrid].re = psi_qp[i*ist->ngrid+jgrid]; 
//       cpsi[jgrid].im = 0.0;
//     }
//     memcpy(&fftwpsi[0], &cpsi[0], ist->ngrid*sizeof(fftwpsi[0]));
//     fftw_execute(planfw[0]);
//     memcpy(&cpsi[0], &fftwpsi[0], ist->ngrid*sizeof(cpsi[0])); // store the FT of the hole wavefunction in cpsi
//     for(jgrid = 0; jgrid < ist->ngrid; jgrid++) {
//       fftwpsi[jgrid][0] = cpsi[jgrid].im * kindex[3*jgrid];
//       fftwpsi[jgrid][1] = -cpsi[jgrid].re * kindex[3*jgrid];
//     }
//     fftw_execute(planbw[0]);
//     for(jgrid = 0; jgrid < ist->ngrid; jgrid++) {
//       psidx[jgrid] = fftwpsi[jgrid][0]; // x-derivative of the hole wavefunction
//     }

//     for(jgrid = 0; jgrid < ist->ngrid; jgrid++) {
//       fftwpsi[jgrid][0] = cpsi[jgrid].im * kindex[3*jgrid+1];
//       fftwpsi[jgrid][1] = -cpsi[jgrid].re * kindex[3*jgrid+1];
//     }
//     fftw_execute(planbw[0]);
//     for(jgrid = 0; jgrid < ist->ngrid; jgrid++) {
//       psidy[jgrid] = fftwpsi[jgrid][0]; // y-derivative of the hole wavefunction
//     }

//     for(jgrid = 0; jgrid < ist->ngrid; jgrid++) {
//       fftwpsi[jgrid][0] = cpsi[jgrid].im * kindex[3*jgrid+2];
//       fftwpsi[jgrid][1] = -cpsi[jgrid].re * kindex[3*jgrid+2];
//     }
//     fftw_execute(planbw[0]);
//     for(jgrid = 0; jgrid < ist->ngrid; jgrid++) {
//       psidz[jgrid] = fftwpsi[jgrid][0]; // z-derivative of the hole wavefunction
//     }

//     for (a = ist->lumo_idx; a < ist->lumo_idx+ist->n_elecs; a++) {
//       sumX = sumY = sumZ = 0.0;
//       for (jz = 0; jz < ist->nz; jz++) {
//         z = grid->z[jz];
//         for (jy = 0; jy < ist->ny; jy++) {
//           y = grid->y[jy];
//           jyz = ist->nx * (ist->ny * jz + jy);
//           for (jx = 0; jx < ist->nx; jx++) {
//             x = grid->x[jx];
//             jgrid = jyz + jx;
//             tmp = psi_qp[a*ist->ngrid+jgrid] * par->dv * 0.5;
//             sumX += tmp * ( y * psidz[jgrid] - z * psidy[jgrid] );
//             sumY += tmp * ( z * psidx[jgrid] - x * psidz[jgrid] );
//             sumZ += tmp * ( x * psidy[jgrid] - y * psidx[jgrid] );
//           }
//         }
//       }
//       mx[i*ist->n_elecs+(a-ist->lumo_idx)] =  sumX;
//       my[i*ist->n_elecs+(a-ist->lumo_idx)] =  sumY;
//       mz[i*ist->n_elecs+(a-ist->lumo_idx)] =  sumZ;
//       fprintf(pf1,"%ld %ld %g\n",i,a,sumX);
//       fprintf(pf2,"%ld %ld %g\n",i,a,sumY);
//       fprintf(pf3,"%ld %ld %g\n",i,a,sumZ);

//       ms=(sqr(mx[i*ist->n_elecs+(a-ist->lumo_idx)]) + sqr(my[i*ist->n_elecs+(a-ist->lumo_idx)]) + sqr(mz[i*ist->n_elecs+(a-ist->lumo_idx)]));
//       bohr_freq = eig_vals[a] - eig_vals[i];
//       fprintf(pf,"%ld %ld %.8f %.12f %.8f %.8f %.8f %.8f\n", i, a, sqrt(ms), bohr_freq, (4.0/3.0)*bohr_freq*ms,
//          sqr(mx[i*ist->n_elecs+(a-ist->lumo_idx)]),
//          sqr(my[i*ist->n_elecs+(a-ist->lumo_idx)]),
//          sqr(mz[i*ist->n_elecs+(a-ist->lumo_idx)]));
//     }
//   }
//   fclose(pf); fclose(pf1); fclose(pf2); fclose(pf3);

//   // Free dynamically allocated memory
//   free(kx); free(ky); free(kz);
//   free(kindex);
//   free(cpsi);
//   free(psidx); free(psidy); free(psidz);

//   return;
// }

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
