/*****************************************************************************/
#include "fd.h"

/*****************************************************************************/

void calc_optical_exc(
    double complex *bs_coeff,
    double *xton_ene,
    double *eig_vals,
    xyz_st *mu,
    xyz_st *m,
    index_st *ist,
    par_st *par)
{

  FILE *pf, *pf1, *pf2, *pcoeff;

  long a, i, ibs, j_xton, idx;

  const long n_el = ist->n_elecs;
  const long n_ho = ist->n_holes;
  const long lidx = ist->lumo_idx;
  const long n_xton = ist->n_xton;

  double os, mos;
  char str[50];

  // Sum-rule accumulators (derivation Eq. 25): sum_n R_n must equal sum_ia r_ia.
  double sum_Rn = 0.0;  // total exciton rotational strength
  double sum_ria = 0.0; // total single-particle rotational strength (const in n)

  // Calculate and print the electric and magnetic dipole strengths and the rotational strength
  pf = fopen("OS.dat", "w");
  pf1 = fopen("M.dat", "w");
  pf2 = fopen("rs.dat", "w");

  for (j_xton = 0; j_xton < n_xton; j_xton++)
  {
    xyz_st mu_sum;
    xyz_st m_sum;
    double rs = 0.0;

    if (j_xton < 10)
    {
      sprintf(str, "bs-coeff-%ld.dat", j_xton);
      pcoeff = fopen(str, "w");
    }

    mu_sum.x = mu_sum.y = mu_sum.z = 0.0 + 0.0 * I;
    m_sum.x = m_sum.y = m_sum.z = 0.0 + 0.0 * I;

    for (ibs = 0, a = lidx; a < lidx + n_el; a++)
    {
      for (i = 0; i < n_ho; i++, ibs++)
      {
        // Exciton electric dipole
        idx = i * n_el + (a - lidx);
        mu_sum.x += bs_coeff[ibs * n_xton + j_xton] * mu[idx].x;
        mu_sum.y += bs_coeff[ibs * n_xton + j_xton] * mu[idx].y;
        mu_sum.z += bs_coeff[ibs * n_xton + j_xton] * mu[idx].z;

        // Exciton magnetic dipole
        m_sum.x += bs_coeff[ibs * n_xton + j_xton] * m[idx].x;
        m_sum.y += bs_coeff[ibs * n_xton + j_xton] * m[idx].y;
        m_sum.z += bs_coeff[ibs * n_xton + j_xton] * m[idx].z;

        // Accumulate the total single-particle rotational strength once, for
        // the sum-rule check sum_n R_n = sum_ia r_ia (derivation Eq. 24-25).
        // r_ia = Im(mu_ia . m_ia^*) is independent of the exciton index n.
        if (j_xton == 0)
          sum_ria += cimag(mu[idx].x * conj(m[idx].x))
                   + cimag(mu[idx].y * conj(m[idx].y))
                   + cimag(mu[idx].z * conj(m[idx].z));

        if (j_xton < 10)
        {
          // Print the complex BSE coefficient C^n_ia that mixes electron-hole
          // pair (i,a) into exciton n = j_xton: hole index, electron index,
          // pair transition energy, Re(C), Im(C), and |C|^2.
          double complex c = bs_coeff[ibs * n_xton + j_xton];
          fprintf(pcoeff, "%ld %ld %.12lg % .12lg % .12lg %lg\n",
                  i, a, eig_vals[a] - eig_vals[i], creal(c), cimag(c), cnorm(c));
        }
      }
    }
    os = (cnorm(mu_sum.x) + cnorm(mu_sum.y) + cnorm(mu_sum.z));
    mos = (cnorm(m_sum.x) + cnorm(m_sum.y) + cnorm(m_sum.z));

    // Exciton rotational strength (derivation Eq. 22, compact form):
    //   R_n = Im[ mu_n . m_n^* ],  mu_n = sum_ia C^n_ia mu_ia,  m_n = sum_ia C^n_ia m_ia
    // This is the FULL bilinear contraction: expanding the product reproduces
    //   R_n = sum_ia sum_jb Im[ C^n_ia (C^n_jb)^* mu_ia . m_jb^* ],
    // i.e. the diagonal (|C|^2-weighted single-particle) piece PLUS the
    // off-diagonal (ia != jb) interference terms that dominate excitonic CD.
    // mu_sum = mu_n and m_sum = m_n are the same C-contractions used for OS/M;
    // mag_dip / elec_dip are both stored as <i|O|a> so no extra conjugation is
    // needed beyond the one m_n^* from the Rosenfeld formula.
    rs = cimag(mu_sum.x * conj(m_sum.x)
             + mu_sum.y * conj(m_sum.y)
             + mu_sum.z * conj(m_sum.z));
    sum_Rn += rs;

    fprintf(pf, "%ld % .8f % .8f % .8f % .12f % .12f % .12f % .12f % .12f % .12f\n", j_xton, sqrt(os), xton_ene[j_xton], (2.0 / 3.0) * xton_ene[j_xton] * os,
            creal(mu_sum.x), cimag(mu_sum.x),
            creal(mu_sum.y), cimag(mu_sum.y),
            creal(mu_sum.z), cimag(mu_sum.z));

    fprintf(pf1, "%ld % .8f % .8f % .8f % .12f % .12f % .12f % .12f % .12f % .12f\n", j_xton, sqrt(mos), xton_ene[j_xton], (4.0 / 3.0) * xton_ene[j_xton] * mos,
            creal(m_sum.x), cimag(m_sum.x),
            creal(m_sum.y), cimag(m_sum.y),
            creal(m_sum.z), cimag(m_sum.z));

    // R_n = Im[mu_n . m_n^*], the full bilinear exciton rotational strength.
    fprintf(pf2, "%ld %.8f % .16f\n", j_xton, xton_ene[j_xton], rs);

    if (j_xton < 10)
    {
      fclose(pcoeff);
    }
  }

  fclose(pf);
  fclose(pf1);
  fclose(pf2);

  // Mandatory numerical check (derivation Eq. 25): the BSE only redistributes
  // rotational strength among excitons, so sum_n R_n must equal sum_ia r_ia.
  // (Both the full and diagonal-only formulas satisfy this, so it verifies the
  // implementation is internally consistent rather than distinguishing them.)
  printf("  Exciton CD sum rule: sum_n R_n = % .10g | sum_ia r_ia = % .10g | diff = % .2e\n",
         sum_Rn, sum_ria, sum_Rn - sum_ria);
  fflush(stdout);
}
