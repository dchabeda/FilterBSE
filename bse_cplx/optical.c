/*****************************************************************************/
#include "fd.h"

/*****************************************************************************/

void calc_optical_exc(
  double complex *bs_coeff, 
  double *eval, 
  xyz_st *mu, 
  xyz_st *m, 
  index_st *ist, 
  par_st *par){
  
  FILE *pf, *pf1, *pf2, *pcoeff;

  long a, i, ibs, j_xton, idx;

  const long      n_el   = ist->n_elecs;
  const long      n_ho   = ist->n_holes;
  const long      lidx   = ist->lumo_idx;
  const long      n_xton = ist->n_xton;

  double os, mos;
  char str[50];

  // Calculate and print the electric and magnetic dipole strengths and the rotational strength
  pf = fopen("OS.dat", "w");
  pf1 = fopen("M.dat", "w");
  pf2 = fopen("rs.dat", "w");


  for (j_xton = 0; j_xton < n_xton; j_xton++) {
    xyz_st mu_sum;
    xyz_st m_sum;
    double complex rs;
    
    if (j_xton < 10) {
      sprintf(str, "bs-coeff-%ld.dat", j_xton);
      pcoeff = fopen(str, "w");
    }
    
    mu_sum.x = mu_sum.y = mu_sum.z = 0.0 + 0.0*I;
    m_sum.x = m_sum.y = m_sum.z = 0.0 + 0.0*I;

    for (ibs = 0, a = lidx; a < lidx+n_el; a++) {
      for (i = 0; i < n_ho; i++, ibs++) {
        // Exciton electric dipole
        idx = i*n_el + (a - lidx);
        mu_sum.x += bs_coeff[ibs*n_xton + j_xton] * mu[idx].x;
        mu_sum.y += bs_coeff[ibs*n_xton + j_xton] * mu[idx].y;
        mu_sum.z += bs_coeff[ibs*n_xton + j_xton] * mu[idx].z;
        
        // Exciton magnetic dipole
        m_sum.x += bs_coeff[ibs*n_xton + j_xton] * m[idx].x;
        m_sum.y += bs_coeff[ibs*n_xton + j_xton] * m[idx].y;
        m_sum.z += bs_coeff[ibs*n_xton + j_xton] * m[idx].z;
        
        if (j_xton < 10){
          fprintf(pcoeff, "%ld %ld %lg\n", i, a, cnorm(bs_coeff[ibs*n_xton + j_xton]));
        }

        // fprintf(po, "%ld %ld %.6lg %.6lg %.6lg %.6lg %.6lg %.6lg\n", i, a, 
        // mu[idx].x_re, mu[idx].x_im,
        // mu[idx].y_re, mu[idx].y_im,
        // mu[idx].z_re, mu[idx].z_im);
        
      }
    } 
    os  = (cnorm(mu_sum.x) + cnorm(mu_sum.y) + cnorm(mu_sum.z));
    mos = (cnorm(m_sum.x) + cnorm(m_sum.y) + cnorm(m_sum.z));

    fprintf(pf,  "%ld % .8f % .8f % .8f % .12f % .12f % .12f % .12f % .12f % .12f\n", j_xton, sqrt(os), eval[j_xton], (2.0/3.0)*eval[j_xton]*os, 
      				mu_sum.x, mu_sum.y, mu_sum.z);
    
    fprintf(pf1, "%ld % .8f % .8f % .8f % .12f % .12f % .12f % .12f % .12f % .12f\n", j_xton, sqrt(mos), eval[j_xton], (4.0/3.0)*eval[j_xton]*mos,
					    m_sum.x, m_sum.y, m_sum.z);
    
    // <a|L|i>*<i|mu|a> = L_ia^* * mu_ia 
    rs = 0.0 + 0.0*I;
    rs += conjmul(m_sum.x, mu_sum.x);
    rs += conjmul(m_sum.y, mu_sum.y);
    rs += conjmul(m_sum.z, mu_sum.z);
    
    fprintf(pf2, "%ld %.8f % .16f % .16f\n", j_xton, eval[j_xton], rs);
    
    if (j_xton < 10){
      fclose(pcoeff);
    }
  }
  
  fclose(pf); fclose(pf1); fclose(pf2);
}
