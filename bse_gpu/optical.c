/*****************************************************************************/
#include "fd.h"

/*****************************************************************************/

void calc_optical_exc(zomplex *bs_coeff, double *eval, xyz_st *mu, xyz_st *m, index_st *ist, par_st *par){
  
  FILE *pf, *pf1, *pf2, *pcoeff;
  long a, i, ibs, jgamma;
  zomplex sumx, sumy, sumz;
  zomplex msumx, msumy, msumz;
  zomplex rs;
  double os, mos;
  char str[50];
  // Calculate and print the electric and magnetic dipole strengths and the rotational strength
  pf = fopen("OS.dat", "w");
  pf1 = fopen("M.dat", "w");
  pf2 = fopen("rs.dat", "w");


  for (jgamma = 0; jgamma < ist->n_xton; jgamma++) {
    sumx.re = sumx.im = 0.0; sumy.re = sumy.im = 0.0; sumz.re = sumz.im = 0.0;
    msumx.re = 0.0; msumy.re = 0.0; msumz.re = 0.0;
    msumx.im = 0.0; msumy.im = 0.0; msumz.im = 0.0;

    if (jgamma < 10) {
      sprintf(str, "bs-coeff-%ld.dat", jgamma);
      pcoeff = fopen(str, "w");
    }

    for (ibs = 0, a = ist->lumo_idx; a < ist->lumo_idx+ist->n_elecs; a++) {
      for (i = 0; i < ist->n_holes; i++, ibs++) {
        // Exciton electric dipole
        sumx.re += bs_coeff[ibs*ist->n_xton + jgamma].re * mu[i*ist->n_elecs + (a - ist->lumo_idx)].x_re
                 - bs_coeff[ibs*ist->n_xton + jgamma].im * mu[i*ist->n_elecs + (a - ist->lumo_idx)].x_im;
        sumx.im += bs_coeff[ibs*ist->n_xton + jgamma].re * mu[i*ist->n_elecs + (a - ist->lumo_idx)].x_im
                 + bs_coeff[ibs*ist->n_xton + jgamma].im * mu[i*ist->n_elecs + (a - ist->lumo_idx)].x_re;
        
        sumy.re += bs_coeff[ibs*ist->n_xton + jgamma].re * mu[i*ist->n_elecs + (a - ist->lumo_idx)].y_re
                 - bs_coeff[ibs*ist->n_xton + jgamma].im * mu[i*ist->n_elecs + (a - ist->lumo_idx)].y_im;
        sumy.im += bs_coeff[ibs*ist->n_xton + jgamma].re * mu[i*ist->n_elecs + (a - ist->lumo_idx)].y_im
                 + bs_coeff[ibs*ist->n_xton + jgamma].im * mu[i*ist->n_elecs + (a - ist->lumo_idx)].y_re;
        
        sumz.re += bs_coeff[ibs*ist->n_xton + jgamma].re * mu[i*ist->n_elecs + (a - ist->lumo_idx)].z_re
                 - bs_coeff[ibs*ist->n_xton + jgamma].im * mu[i*ist->n_elecs + (a - ist->lumo_idx)].z_im;
        sumz.im += bs_coeff[ibs*ist->n_xton + jgamma].re * mu[i*ist->n_elecs + (a - ist->lumo_idx)].z_im
                 + bs_coeff[ibs*ist->n_xton + jgamma].im * mu[i*ist->n_elecs + (a - ist->lumo_idx)].z_re;

        // Exciton magnetic dipole
        msumx.re += bs_coeff[ibs*ist->n_xton + jgamma].re * m[i*ist->n_elecs + (a - ist->lumo_idx)].x_re
                 - bs_coeff[ibs*ist->n_xton + jgamma].im * m[i*ist->n_elecs + (a - ist->lumo_idx)].x_im;
        msumx.im += bs_coeff[ibs*ist->n_xton + jgamma].re * m[i*ist->n_elecs + (a - ist->lumo_idx)].x_im
                 + bs_coeff[ibs*ist->n_xton + jgamma].im * m[i*ist->n_elecs + (a - ist->lumo_idx)].x_re;
        
        msumy.re += bs_coeff[ibs*ist->n_xton + jgamma].re * m[i*ist->n_elecs + (a - ist->lumo_idx)].y_re
                 - bs_coeff[ibs*ist->n_xton + jgamma].im * m[i*ist->n_elecs + (a - ist->lumo_idx)].y_im;
        msumy.im += bs_coeff[ibs*ist->n_xton + jgamma].re * m[i*ist->n_elecs + (a - ist->lumo_idx)].y_im
                 + bs_coeff[ibs*ist->n_xton + jgamma].im * m[i*ist->n_elecs + (a - ist->lumo_idx)].y_re;
        
        msumz.re += bs_coeff[ibs*ist->n_xton + jgamma].re * m[i*ist->n_elecs + (a - ist->lumo_idx)].z_re
                 - bs_coeff[ibs*ist->n_xton + jgamma].im * m[i*ist->n_elecs + (a - ist->lumo_idx)].z_im;
        msumz.im += bs_coeff[ibs*ist->n_xton + jgamma].re * m[i*ist->n_elecs + (a - ist->lumo_idx)].z_im
                 + bs_coeff[ibs*ist->n_xton + jgamma].im * m[i*ist->n_elecs + (a - ist->lumo_idx)].z_re;

        
        if (jgamma < 10){
          fprintf(pcoeff, "%ld %ld %lg\n", i, a, sqr(bs_coeff[ibs*ist->n_xton + jgamma].re) + sqr(bs_coeff[ibs*ist->n_xton + jgamma].im));
        }

        // fprintf(po, "%ld %ld %.6lg %.6lg %.6lg %.6lg %.6lg %.6lg\n", i, a, 
        // mu[i*ist->n_elecs + (a - ist->lumo_idx)].x_re, mu[i*ist->n_elecs + (a - ist->lumo_idx)].x_im,
        // mu[i*ist->n_elecs + (a - ist->lumo_idx)].y_re, mu[i*ist->n_elecs + (a - ist->lumo_idx)].y_im,
        // mu[i*ist->n_elecs + (a - ist->lumo_idx)].z_re, mu[i*ist->n_elecs + (a - ist->lumo_idx)].z_im);
        
      }
    } 
    os  = (sqr(sumx.re)+sqr(sumx.im) + sqr(sumy.re)+sqr(sumy.im) + sqr(sumz.re)+sqr(sumz.im));
    mos = (sqr(msumx.re)+sqr(msumx.im) + sqr(msumy.re)+sqr(msumy.im) + sqr(msumz.re)+sqr(msumz.im));
    fprintf(pf,  "%ld % .8f % .8f % .8f % .12f % .12f % .12f % .12f % .12f % .12f\n", jgamma, sqrt(os), eval[jgamma], (2.0/3.0)*eval[jgamma]*os, 
      				sumx.re, sumx.im, sumy.re, sumy.im, sumz.re, sumz.im);
    
    fprintf(pf1, "%ld % .8f % .8f % .8f % .12f % .12f % .12f % .12f % .12f % .12f\n", jgamma, sqrt(mos), eval[jgamma], (4.0/3.0)*eval[jgamma]*mos,
					    msumx.re, msumx.im, msumy.re, msumy.im, msumz.re, msumz.im);
    
    // <a|L|i>*<i|mu|a> = L_ia^* * mu_ia 
    rs.re = msumx.re*sumx.re + msumy.re*sumy.re + msumz.re*sumz.re;
    // L_ia is complex conjugated, so its imag part is negative (overall positive)
    rs.re += msumx.im*sumx.im + msumy.im*sumy.im + msumz.im*sumz.im;

    rs.im = msumx.re*sumx.im + msumy.re*sumy.im + msumz.re*sumz.re;
    rs.im -= msumx.im*sumx.re + msumy.im*sumy.re + msumz.im*sumz.re;

    fprintf(pf2, "%ld %.8f % .16f % .16f\n", jgamma, eval[jgamma], rs.re, rs.im);
    
    if (jgamma < 10){
      fclose(pcoeff);
    }
  }
  
  fclose(pf); fclose(pf1); fclose(pf2);
}
