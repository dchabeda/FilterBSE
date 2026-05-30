#include "Hmat.h"

/*****************************************************************************/

void diag_H(
    double *psitot,
    double *pot_local,
    vector *G_vecs,
    vector *k_vecs,
    zomplex *LS,
    nlc_st *nlc,
    long *nl,
    double *ksqr,
    double *eval,
    long n_states,
    int k_idx,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel)
{
  /*******************************************************************
   * This function calculates eigenvalues and vectors of the real or  *
   * complex valued matrix, H, where H_ij = <psi_i|H|psi_j>           *
   * It utilizes the MKL LAPACK routines to do the diagonalization    *
   * REAL VALUED                                                      *
   *  https://www.netlib.org/lapack/explore-html-3.6.1/d2/d8a/group__double_s_yeigen_ga442c43fca5493590f8f26cf42fed4044.html#ga442c43fca5493590f8f26cf42fed4044
   * COMPLEX VALUED                                                   *
   *   https://www.netlib.org/lapack/explore-html-3.6.1/d6/dee/zheev_8f_a70c041fd19635ff621cfd5d804bd7a30.html
   * inputs:                                                          *
   *  [psi/phi] container for zomplex wavefunction (for evaluating H) *
   *  [psitot] arr holds all mn_states_tot orthogonal wavefunctions   *
   *  [psims] arr that will hold all ms filtered wavefunctions        *
   *  [pot_local] ngrid-long arr holding the value of the local pot   *
   *  [nlc] nlc struct holding values for computing SO and NL pots    *
   *  [nl] natom-long arr holding the number of NL gridpts per atom   *
   *  [ksqr] ngrid-long arr holding the values of k^2 for KE calc     *
   *  [eval] array that will hold the eigenvalues                     *
   *  [ist] ptr to counters, indices, and lengths                     *
   *  [par] ptr to par_st holding VBmin, VBmax... params              *
   *  [flag] ptr to flag_st holding job flags                         *
   *  [planfw] FFTW3 plan for executing 3D forward DFT                *
   *  [planbw] FFTW3 plan for executing 3D backwards DFT              *
   *  [fftwpsi] location to store outcome of Fourier transform        *
   * outputs: void                                                    *
   ********************************************************************/

  FILE *pg;
  long ims, jms, jgrid, jgrid_real, jgrid_imag;
  const long long mn_states_tot = (int)n_states, lwk = 3 * n_states;
  long long info;
  double *rwork, sumre, sumim;
  zomplex *tpsi;
  double *H, *work;
  MKL_Complex16 *H_z, *work_z;
  time_t current_time;
  char *c_time_string;
  zomplex *psi, *phi;

  const long stlen = ist->nspinngrid * ist->complex_idx;

  // In the periodic path diag_H runs once per k-point, concurrently on each k-group
  // master, so tag the status lines with the k index to keep the interleaved stdout
  // readable. The non-periodic path runs a single diagonalization, so no tag.
  char kpfx[24];
  if (1 == flag->periodic)
    sprintf(kpfx, "ik=%d: ", k_idx);
  else
    kpfx[0] = '\0';

  // FFT
  long fft_flags = FFTW_MEASURE;
  fftw_plan_loc planfw, planbw;
  fftw_complex *fftwpsi;
  // Create FFT structs and plans for Fourier transform
  fftw_init_threads();
  fftw_plan_with_nthreads(par->ham_threads);
  fftwpsi = fftw_malloc(sizeof(fftw_complex) * ist->ngrid);
  planfw = fftw_plan_dft_3d(ist->nz, ist->ny, ist->nx, fftwpsi, fftwpsi, FFTW_FORWARD, fft_flags);
  planbw = fftw_plan_dft_3d(ist->nz, ist->ny, ist->nx, fftwpsi, fftwpsi, FFTW_BACKWARD, fft_flags);

  // Allocate memory specific to scalar or complex arrays
  if (0 == flag->isComplex)
  {
    H = (double *)calloc(n_states * n_states, sizeof(double));
    work = (double *)calloc(lwk, sizeof(double));
    tpsi = (zomplex *)calloc(n_states, sizeof(zomplex));
  }
  if (1 == flag->isComplex)
  {
    H_z = (MKL_Complex16 *)calloc(n_states * n_states, sizeof(MKL_Complex16));
    work_z = (MKL_Complex16 *)calloc(lwk, sizeof(MKL_Complex16));
    rwork = (double *)calloc(3 * n_states, sizeof(double));
    tpsi = (zomplex *)calloc(n_states, sizeof(zomplex));
  }

  if ((psi = (zomplex *)calloc(ist->nspinngrid, sizeof(psi[0]))) == NULL)
  {
    fprintf(stderr, "\nOUT OF MEMORY: psi in diag_H\n\n");
    exit(EXIT_FAILURE);
  }
  if ((phi = (zomplex *)calloc(ist->nspinngrid, sizeof(phi[0]))) == NULL)
  {
    fprintf(stderr, "\nOUT OF MEMORY: phi in diag_H\n\n");
    exit(EXIT_FAILURE);
  }

  printf("%sConstructing Hamiltonian matrix\n", kpfx);
  fflush(stdout);

  omp_set_dynamic(0);
  omp_set_num_threads(parallel->nthreads);
  /*** calculate H|psi_i> ***/

  for (ims = 0; ims < n_states; ims++)
  {

    if (1 == flag->isComplex)
    {
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++)
      {
        jgrid_real = ist->complex_idx * jgrid;
        jgrid_imag = jgrid_real + 1;

        psi[jgrid].re = psitot[ims * stlen + jgrid_real];
        psi[jgrid].im = psitot[ims * stlen + jgrid_imag];
      }
    }
    else if (0 == flag->isComplex)
    {
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++)
      {
        psi[jgrid].re = psitot[ims * stlen + jgrid];
        psi[jgrid].im = 0.0;
      }
    }

    memcpy(&phi[0], &psi[0], ist->nspinngrid * sizeof(phi[0]));
    if (flag->periodic == 0)
    {
      p_hamiltonian(phi, psi, pot_local, LS, nlc, nl, ksqr, ist, par, flag, planfw, planbw, fftwpsi, parallel->nthreads);
    }
    else if (flag->periodic == 1)
    {
      // k_idx is the global k-point index this Hamiltonian block is built at
      // printf("\nk index to diagonalize %d: k = %.6lg\n", k_idx, k_vecs[k_idx].z);
      p_hamiltonian_k(phi, psi, pot_local, G_vecs, k_vecs[k_idx], LS, nlc, nl, ist, par, flag, planfw, planbw, fftwpsi, parallel->nthreads);
    }

/*** calculate <psi_j|H|psi_i> ***/
#pragma omp parallel for private(jms, jgrid, jgrid_real, jgrid_imag) shared(H, H_z)
    for (jms = 0; jms <= ims; jms++)
    {
      // Utilize Hermitian property of Hamiltonian matrix elements to reduce computation to only lower triangle
      if (0 == flag->isComplex)
      {
        H[ims * n_states + jms] = H[jms * n_states + ims] = dotpreal(phi, psitot, jms, ist->nspinngrid, par->dv);
      }
      if (1 == flag->isComplex)
      {
        H_z[ims * n_states + jms] = H_z[jms * n_states + ims] = dotp(phi, psitot, jms, ist->nspinngrid, par->dv);
        H_z[jms * n_states + ims].imag *= -1;
      }
    }

    if ((ims == 0) || (0 == (ims % (n_states / 4 + 1))) || (ims == (n_states - 1)))
    {
      print_progress_bar(ims, n_states, (1 == flag->periodic) ? k_idx : -1);
    }
  }
  // free dynamically allocated memory for psi and phi
  free(psi);
  free(phi);

  // print out Hmat for debugging purposes (per-k so masters do not clobber)
  char hmat_name[64];
  sprintf(hmat_name, "hmat-%d.dat", k_idx);
  pg = fopen(hmat_name, "w");
  for (ims = 0; ims < n_states; ims++)
  {
    for (jms = 0; jms < n_states; jms++)
    {
      if (0 == flag->isComplex)
      {
        fprintf(pg, "%lg ", H[ims * n_states + jms]);
      }
      if (1 == flag->isComplex)
      {
        fprintf(pg, "%lg+i%lg ", H_z[ims * n_states + jms].real, H_z[ims * n_states + jms].imag);
      }
    }
    fprintf(pg, "\n");
  }
  fclose(pg);

  /*** diagonalize the Hamiltonian H ***/
  printf("%sDiagonalizing Hamiltonian | %s\n", kpfx, get_time());
  // Use real, symmetric diagonalization routine for real wavefunctions
  if (0 == flag->isComplex)
  {
    dsyev_("V", "U", &mn_states_tot, &H[0], &mn_states_tot, &eval[0], &work[0], &lwk, &info);
    if (info)
    {
      fprintf(stderr, "error in dsyev_ H\n");
      exit(EXIT_FAILURE);
    }
  }
  // Use Hermitian diagonalization routine for complex wavefunctions
  if (1 == flag->isComplex)
  {
    zheev_("V", "U", &mn_states_tot, &H_z[0], &mn_states_tot, &eval[0], &work_z[0], &lwk, &rwork[0], &info);
    if (info)
    {
      fprintf(stderr, "error in zheev_ H\n");
      exit(EXIT_FAILURE);
    }
  }

  /*** copy the new function into psitot ***/
  printf("%sDiagonalization complete! | %s\n", kpfx, get_time());
  fflush(stdout);

  // The eigenvectors have been computed in the basis of orthogonalized
  // filtered functions (Phi_filter). In order to obtain them in the grid basis
  // as Psi_grid, we must perform a change of basis.
  // The matrix H is holding the eigenvectors, so we perform
  // (Psi_grid)_a = SUM_i H_ai * (Phi_filter)_i
  // for each grid point
  printf("%sWriting out eigenvectors in the grid basis\n", kpfx);
  for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++)
  {
    jgrid_real = ist->complex_idx * jgrid;
    jgrid_imag = ist->complex_idx * jgrid + 1;

#pragma omp parallel for private(jms)
    for (jms = 0; jms < n_states; jms++)
    {
      tpsi[jms].re = psitot[ist->complex_idx * jms * ist->nspinngrid + jgrid_real];
      psitot[ist->complex_idx * jms * ist->nspinngrid + jgrid_real] = 0.0;

      if (1 == flag->isComplex)
      {
        tpsi[jms].im = psitot[ist->complex_idx * jms * ist->nspinngrid + jgrid_imag];
        psitot[ist->complex_idx * jms * ist->nspinngrid + jgrid_imag] = 0.0;
      }
    }

#pragma omp parallel for private(jms, sumre, sumim)
    for (jms = 0; jms < n_states; jms++)
    {

      sumre = sumim = 0.0;
      for (ims = 0; ims < n_states; ims++)
      {
        if (0 == flag->isComplex)
        {
          sumre += H[jms * n_states + ims] * tpsi[ims].re;
        }
        if (1 == flag->isComplex)
        {
          sumre += (H_z[jms * n_states + ims].real * tpsi[ims].re + H_z[jms * n_states + ims].imag * tpsi[ims].im);
          sumim += (H_z[jms * n_states + ims].real * tpsi[ims].im - H_z[jms * n_states + ims].imag * tpsi[ims].re);
        }
      }

      // #pragma omp critical
      psitot[ist->complex_idx * jms * ist->nspinngrid + jgrid_real] = sumre;
      if (1 == flag->isComplex)
      {
        psitot[ist->complex_idx * jms * ist->nspinngrid + jgrid_imag] = sumim;
      }
    }

    if (!(jgrid % (ist->nspinngrid / 4)))
    {
      current_time = time(NULL);
      c_time_string = ctime(&current_time);
      printf("\t%sFinished grid point no. %ld | %s\n", kpfx, jgrid, c_time_string);
      fflush(stdout);
    }
  }

  free(tpsi);
  if (0 == flag->isComplex)
  {
    free(work);
    free(H);
  }
  if (1 == flag->isComplex)
  {
    free(H_z);
    free(work_z);
    free(rwork);
  }

  return;
}

/****************************************************************************/

double dotpreal(zomplex *psi, double *phi, long m, long ngrid, double dv)
{

  long jgrid;
  double sum;

  sum = 0.0;
  for (jgrid = 0; jgrid < ngrid; jgrid++)
  {
    sum += (psi[jgrid].re * phi[m * ngrid + jgrid]);
  }
  sum *= dv;

  return (sum);
}

/****************************************************************************/

MKL_Complex16 dotp(zomplex *psi, double *psitot, long m, long ngrid, double dv)
{

  long jgrid, jgrid_real, jgrid_imag;
  MKL_Complex16 sum;

  sum.real = sum.imag = 0.0;
  for (jgrid = 0; jgrid < ngrid; jgrid++)
  {
    jgrid_real = 2 * jgrid;
    jgrid_imag = 2 * jgrid + 1;

    sum.real += (psi[jgrid].re * psitot[2 * m * ngrid + jgrid_real] + psi[jgrid].im * psitot[2 * m * ngrid + jgrid_imag]);
    sum.imag += (psi[jgrid].re * psitot[2 * m * ngrid + jgrid_imag] - psi[jgrid].im * psitot[2 * m * ngrid + jgrid_real]);
  }
  sum.real *= dv;
  sum.imag *= dv;

  return (sum);
}
