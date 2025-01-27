/*****************************************************************************/
#include "fd.h"

/*****************************************************************************/

void calc_optical_exc(zomplex *bs_coeff, double *eval, xyz_st *mu, xyz_st *m, index_st *ist, par_st *par){
  
  FILE *pf, *pcoeff;
  long a, i, ibs, jgamma;
  zomplex sumx, sumy, sumz;
  double os;
  char str[50];
  // Calculate and print the electric and magnetic dipole strengths and the rotational strength
  pf = fopen("OS.dat", "w");
  // pf1 = fopen("M.dat", "w");
  // pf2 = fopen("rs.dat", "w");

  FILE *po;

  for (jgamma = 0; jgamma < ist->n_xton; jgamma++) {
    sumx.re = sumx.im = 0.0; sumy.re = sumy.im = 0.0; sumz.re = sumz.im = 0.0;
    //msumx = 0.0; msumy = 0.0; msumz = 0.0;
    po = fopen("mu_test.dat", "w");
    if (jgamma < 10) {
      sprintf(str, "bs-coeff-%ld.dat", jgamma);
      pcoeff = fopen(str, "w");
    }
    for (ibs = 0, a = ist->lumo_idx; a < ist->lumo_idx+ist->n_elecs; a++) {
      for (i = 0; i < ist->n_holes; i++, ibs++) {
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
        //msumx += bs_coeff[jgamma*ist->n_xton + ibs] * mx[i*ist->n_elecs + (a - ist->lumo_idx)];
		    //msumy += bs_coeff[jgamma*ist->n_xton + ibs] * my[i*ist->n_elecs + (a - ist->lumo_idx)];
		    //msumz += bs_coeff[jgamma*ist->n_xton + ibs] * mz[i*ist->n_elecs + (a - ist->lumo_idx)];
        if (jgamma < 10) fprintf(pcoeff, "%ld %ld %lg\n", i, a, sqr(bs_coeff[ibs*ist->n_xton + jgamma].re) + sqr(bs_coeff[ibs*ist->n_xton + jgamma].im));
        fprintf(po, "%ld %ld %.6lg %.6lg %.6lg %.6lg %.6lg %.6lg\n", i, a, mu[i*ist->n_elecs + (a - ist->lumo_idx)].x_re, mu[i*ist->n_elecs + (a - ist->lumo_idx)].x_im,
        mu[i*ist->n_elecs + (a - ist->lumo_idx)].y_re, mu[i*ist->n_elecs + (a - ist->lumo_idx)].y_im,
        mu[i*ist->n_elecs + (a - ist->lumo_idx)].z_re, mu[i*ist->n_elecs + (a - ist->lumo_idx)].z_im);
        
      }
    } 
    os  = (sqr(sumx.re)+sqr(sumx.im) + sqr(sumy.re)+sqr(sumy.im) + sqr(sumz.re)+sqr(sumz.im));
    //mos = (msumx*msumx + msumy*msumy + msumz*msumz);
    fprintf(pf,  "%ld %.8f %.8f % .8f % .12f % .12f % .12f % .12f % .12f % .12f\n", jgamma, sqrt(os), eval[jgamma], (2.0/3.0)*eval[jgamma]*os, 
      				sumx.re, sumx.im, sumy.re, sumy.im, sumz.re, sumz.im);
    //fprintf(pf1, "%ld %.8f %.8f % .8f % .12f % .12f % .12f\n", jgamma, sqrt(mos), eval[jgamma], (4.0/3.0)*eval[jgamma]*mos,
		//			    msumx*msumx, msumy*msumy, msumz*msumz);
    //fprintf(pf2, "%ld %.8f % .16f\n", jgamma, eval[jgamma], (sumx*msumx + sumy*msumy + sumz*msumz));
    fclose(po);
    if (jgamma < 10){
      fclose(pcoeff);
    }
  }
  
  fclose(pf); // fclose(pf1); fclose(pf2);
}