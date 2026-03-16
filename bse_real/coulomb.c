#include "coulomb.h"

/**************************************************************************/
//      this routine computes the coulomb coupling between
//      single excitons.  On input - it requires the eigenstates stored in psi_qp,
//      the eigenvalues stored in eval, and pot_hartree computed in init_elec_hole_kernel.
//      On output it stores the coulomb matrix elements on the disk
//      in the following format: a, i, b, j, ene_ai, ene_bj, vjbai, vabji.
//      a - the index of the electron in exciton Sai.
//      i - the index of the hole in exciton Sai.
//      b - the index of the electron in exciton Sbj.
//      j - the index of the hole in exciton Sbj.
//      ene_ai - the energy of exciton Sai.
//      ene_bj - the energy of exciton Sbj.
//      vjbai and vabji are the coulomb matrix elements need to be used to
//      generate the spin-depedent matrix elements as described by
//      the last equation in our document.  ***/
/**************************************************************************/

void calc_eh_kernel_real(
    double *restrict psi_qp,
    double complex *restrict pot_bare,
    double complex *restrict pot_screened,
    double *restrict direct,
    double *restrict exchange,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel)
{

  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  FILE *pf;

  const int mpir = parallel->mpi_rank;

  //                Indices
  long cntr = 0;

  long i;
  long a;
  long j;
  long b;
  long i_st;
  long a_st;
  long b_st;
  long ibs;
  long jbs;
  long start;
  long loop_idx;
  long ncycles;
  long jg;
  long jsg;

  const long nspngr = ist->nspinngrid;
  const long ngrid = ist->ngrid;
  const long lidx = ist->lumo_idx;
  const long n_el = ist->n_elecs;
  const long n_ho = ist->n_holes;
  const long n_xton = ist->n_xton;

  long *listibs;

  const double dv = par->dv;

  char *fileName;
  fileName = (char *)malloc(30 * sizeof(char) + 1);
  fileName[30] = '\0';

  double complex *restrict rho;
  double *restrict pot_htree;

  ALLOCATE(&rho, parallel->nthreads * ngrid, "rho in coulomb");
  ALLOCATE(&listibs, ist->n_xton, "listibs in coulomb");
  ALLOCATE(&pot_htree, parallel->nthreads * nspngr, "pot_htree");

  /************************************************************/
  /*******************    INIITIALIZE FFT   *******************/
  /************************************************************/

  // Parallel FFT
  fftw_plan_loc *planfw;
  fftw_plan_loc *planbw;
  fftw_complex *fftwpsi;
  long fft_flags = FFTW_MEASURE;

  // Create FFT structs and plans for Fourier transform
  // fftw_init_threads();
  // fftw_plan_with_nthreads(ist->nthreads);

  fftwpsi = fftw_malloc(parallel->nthreads * sizeof(fftw_complex) * ngrid);

  for (j = 0; j < parallel->nthreads; j++)
  {

    planfw[j] = fftw_plan_dft_3d(ist->nz, ist->ny, ist->nx,
                                 &fftwpsi[j * ngrid], fftwpsi[j * ngrid], FFTW_FORWARD, fft_flags);

    planbw[j] = fftw_plan_dft_3d(ist->nz, ist->ny, ist->nx,
                                 fftwpsi[j * ngrid], fftwpsi[j * ngrid], FFTW_BACKWARD, fft_flags);
  }

  /************************************************************/
  /*******************    HANDLE INDEXING   *******************/
  /************************************************************/

  for (ibs = 0, a = lidx; a < lidx + n_el; a++)
  {
    for (i = 0; i < n_ho; i++, ibs++)
    {
      listibs[(a - lidx) * n_ho + i] = ibs;
    }
  }

  // Set the OpenMP parallelization to use nthreads
  omp_set_max_active_levels(1);
  omp_set_num_threads(ist->nthreads);

  /************************************************************/
  /*******************    CALC DIRECT K^D   *******************/
  /************************************************************/

  /*** vabji direct ***/
  // We avoid performing additional computational work
  // 2e integrals have a 4-fold permutation symmetry if real
  // [ij|ab] = [ji|ab] = [ji|ba] = [ij|ba]
  // using Chemist's notation from Szabo and Ostlund

  long ab, ij;
  long ab_tot = n_el * n_el;
  long ij_tot = n_ho * n_ho;
  long last_ab = 0;
  long *lista, *listb, *listi, *listj;

  lista = (long *)calloc(ab_tot, sizeof(long));
  listb = (long *)calloc(ab_tot, sizeof(long));
  listi = (long *)calloc(ij_tot, sizeof(long));
  listj = (long *)calloc(ij_tot, sizeof(long));

  // Generate the indices for flattened loops when a < b
  ab_tot = 0; // Calc ab_tot from the loop trip count
  for (a = lidx; a < lidx + n_el; a++)
  {
    for (b = lidx; b < lidx + n_el; b++)
    {
      lista[ab_tot] = a;
      listb[ab_tot] = b;
      ab_tot++;
    }
  }

  // i, j inner loop indices
  for (i = 0; i < n_ho; i++)
  {
    for (j = 0; j < n_ho; j++)
    {
      listi[i * n_ho + j] = i;
      listj[i * n_ho + j] = j;
    }
  }

  /************************************************************/
  /******************    DO K^D INTEGRAL    *******************/
  /************************************************************/

  printf("Computing direct mat | %s\n", get_time());
  sprintf(fileName, "direct.dat", mpir);
  pf = fopen(fileName, "w");

  cntr = 0;
  ncycles = ab_tot;
  double ab_start_t = omp_get_wtime();
#pragma omp parallel for private(a, b, a_st, b_st, jg, i, j, ij, i_st, j_st, ibs, jbs) // schedule(dynamic, chunk_size)
  for (ab = 0; ab < ab_tot; ab++)
  {
    long thread_id = omp_get_thread_num();
    long thread_st = thread_id * ngrid;

    a = lista[ab];
    b = listb[ab];
    // printf("\n even rank %d ab = %ld a = %ld b = %ld\n", even_rank, ab, a, b); fflush(0);
    // Grab indices of electron-electron states a, b
    a_st = a * nspngr;
    b_st = b * nspngr;

    // Compute hartree potential for a, b density
    // 1) Compute joint density and store in rho

    for (jg = 0; jg < ngrid; jg++)
    {
      rho[thread_st + jg] = psi_qp[a_st + jg] * psi_qp[b_st + jg] + 0.0 * I;
    }

    // Compute the hartree potential and store in pot_htree
    // h_d(r) = \int W(r,r') \rho_{ab}(r') d^3r' via fourier transform

    hartree(&rho[thread_st], pot_screened, &pot_htree[thread_st], ist, planfw[thread_id], planbw[thread_id], &fftwpsi[thread_st]);

    // loop over hole states i, j

    for (ij = 0; ij < ij_tot; ij++)
    {
      // get the matrix indicies for {ai,bj}
      i = listi[ij];
      j = listj[ij];

      long i_st = i * nspngr;
      long j_st = j * nspngr;
      long ibs = listibs[(a - lidx) * n_ho + i];
      long jbs = listibs[(b - lidx) * n_ho + j];

      if (ibs < jbs)
      {
        continue;
      }

      long jg;
      double sum;
      sum = 0.0;

      // K^d_{ai,bj}=\int h_d(r) \sum_\sigma psi_{i}(r,\sigma) psi_{j}^{*}(r,\sigma) d^3r
      for (jg = 0; jg < ngrid; jg++)
      {
        sum += pot_htree[thread_st + jg] * psi_qp[j_st + jg] * psi_qp[i_st + jg];
      }
      sum *= dv;

      direct[ibs * n_xton + jbs] = sum;
    } // end of ij

    // fflush(0);
    for (ij = 0; ij < ij_tot; ij++)
    {
      i = listi[ij];
      j = listj[ij];
      ibs = listibs[(a - lidx) * n_ho + i];
      jbs = listibs[(b - lidx) * n_ho + j];

      fprintf(pf, "%03ld %03ld %03ld %03ld %ld %ld %.12f\n", a, b, i, j, ibs, jbs,
              direct[ibs * n_xton + jbs]);
    }
    // Every 25% of iterations, print output
    if ((cntr == 0) || (0 == cntr % (ncycles / 4 + 1)) || (cntr == (ncycles - 1)))
    {
      // Print out progress
      printf("Direct: ");
      print_progress_bar(cntr, ncycles);

      fflush(0);
    }
    cntr++;
  } // end of ab
  double ab_end_t = omp_get_wtime();

  printf("Done with direct; duration = %lg s (%lg m)\n", (ab_end_t - ab_start_t), (ab_end_t - ab_start_t) / 60.0);

  // Free
  free(lista);
  free(listb);
  free(listi);
  free(listj);
  // end of even MPI ranks

  /*******************************************************************/
  /*******************************************************************/
  /****** BREAK ***** BREAK ***** BREAK ***** BREAK ***** BREAK ******/
  /*******************************************************************/
  /*******************************************************************/

  if (!flag->calcDarkStates)
  {
    for (jg = 0; jg < parallel->nthreads * ngrid; jg++)
    {
      rho[jg] = 0.0 + I * 0.0;
    }

    long ai;
    long bj;
    long ai_tot = n_el * n_ho;
    long bj_tot = n_el * n_ho;
    long last_ai = 0;
    long *lista;
    long *listi;
    long *listb;
    long *listj;
    long jstart;

    ALLOCATE(&lista, ai_tot, "lista");
    ALLOCATE(&listi, ai_tot, "listi");
    ALLOCATE(&listb, ai_tot, "listb");
    ALLOCATE(&listj, ai_tot, "listj");

    // Generate the indices for flattened loop a, i
    for (a = lidx; a < lidx + n_el; a++)
    {
      for (i = 0; i < n_ho; i++)
      {
        lista[(a - lidx) * n_ho + i] = a;
        listi[(a - lidx) * n_ho + i] = i;
      }
    }

    // Generate the indices for flattened loop b, j
    for (b = lidx; b < lidx + n_el; b++)
    {
      for (j = 0; j < n_ho; j++)
      {
        listb[(b - lidx) * n_ho + j] = b;
        listj[(b - lidx) * n_ho + j] = j;
      }
    }

    ncycles = ai_tot;

    printf("\tNumber of exchange ncycles = %lu\n", ncycles);
    fflush(0);

    /************************************************************/
    /******************    DO K^X INTEGRAL   ********************/
    /************************************************************/

    sprintf(fileName, "exchange.dat");
    pf = fopen(fileName, "w");

    // loop over electron states a, i
    cntr = 0;
    double ai_start_t = omp_get_wtime();
#pragma omp parallel for private(a, i, b, j, a_st, i_st, b_st, j_st, jg, bj, ibs, jbs) // schedule(dynamic, chunk_size)
    for (ai = 0; ai < ai_tot; ai++)
    {
      long thread_id = omp_get_thread_num();
      long thread_st = thread_id * ngrid;
      a = lista[ai];
      i = listi[ai];

      a_st = a * nspngr;
      i_st = i * nspngr;

      // 1) Compute joint density and store in rho
      for (jg = 0; jg < ngrid; jg++)
      {
        rho[thread_st + jg] = psi_qp[a_st + jg] * psi_qp[i_st + jg] + 0.0 * I;
      }

      // Compute the hartree potential and store in pot_htree
      // h_d(r) = \int W(r,r') \rho_{ab}(r') d^3r' via fourier transform
      hartree(&rho[thread_st], pot_bare, &pot_htree[thread_st], ist, planfw[thread_id], planbw[thread_id], &fftwpsi[thread_st]);

      // loop over electron-hole pairs b, j
      for (bj = 0; bj < bj_tot; bj++)
      {
        b = listb[bj];
        j = listj[bj];

        long b_st = b * nspngr;
        long j_st = j * nspngr;
        long ibs = listibs[(a - lidx) * n_ho + i];
        long jbs = listibs[(b - lidx) * n_ho + j];

        if (ibs < jbs)
        {
          continue;
        }

        long jg;
        double sum;
        sum = 0.0;

        for (jg = 0; jg < ngrid; jg++)
        {
          sum += pot_htree[jg] * psi_qp[j_st + jg] * psi_qp[b_st + jg];
        }
        sum *= dv;

        exchange[ibs * n_xton + jbs] = -2.0 * sum;
      } // end of bj

      for (bj = 0; bj < bj_tot; bj++)
      {
        b = listb[bj];
        j = listj[bj];
        ibs = listibs[(a - lidx) * n_ho + i];
        jbs = listibs[(b - lidx) * n_ho + j];

        if (ibs < jbs)
        {
          continue;
        }

        fprintf(pf, "%03ld %03ld %03ld %03ld %ld %ld %.12f\n", a, b, i, j, ibs, jbs,
                exchange[ibs * ist->n_xton + jbs]);
      }

      // Every 25% of iterations, print the job progress
      if ((cntr == 0) || (0 == cntr % (ncycles / 4 + 1)) || (cntr == (ncycles - 1)))
      {
        // Print progress
        printf("Exchange: ");
        print_progress_bar(cntr, ncycles);
        fflush(0);
      }

      cntr++;
    } // end of ai
    double ai_end_t = omp_get_wtime();
    fclose(pf);

    printf("Done with exchange integrals on rank %d; duration = %lg s (%lg m)\n", odd_rank, (ai_end_t - ai_start_t), (ai_end_t - ai_start_t) / 60.0);
    fflush(0);

    free(lista);
    free(listi);
    free(listb);
    free(listj);
  }

  free(rho);
  free(listibs);
  free(fileName);
  free(pot_htree);

  fftw_free(fftwpsi);
  for (j = 0; j < parallel->nthreads; j++)
  {
    fftw_destroy_plan(planfw[j]);
    fftw_destroy_plan(planbw[j]);
  }

  return;
}

/***************************************************************************************/

int load_coulomb_mat(double *mat, char *fileName, long *a_max, long *b_max, long *i_max, long *j_max, index_st *ist)
{

  FILE *pf;

  long start;
  long cntr = 0;
  int ieof = 0;
  int done_flag;
  long a, b, i, j, ibs, jbs;
  long tmp_a, tmp_b, tmp_i, tmp_j;
  long max_st_num = ist->n_holes * ist->n_holes * ist->n_elecs * ist->n_elecs;
  double tmp_re, tmp_im;

  // Open the direct/exchange.dat file

  pf = fopen(fileName, "r");

  if (pf == NULL)
  {
    printf("ERROR: could not open file %s\n", fileName);
    fprintf(stderr, "ERROR: could not open file %s\n", fileName);
    exit(EXIT_FAILURE);
  }

  // Scan all the lines and load the values into mat
  // Note, order of a,b,i,j doesn't matter because ibs and jbs
  // are the actual indices of the matrices
  tmp_a = tmp_b = tmp_i = tmp_j = 0;
  while ((ieof != EOF) && (cntr < max_st_num))
  {
    // Scan the file and grab matrix elements
    ieof = fscanf(pf, "%ld %ld %ld %ld %ld %ld %lg %lg", &a, &b, &i, &j, &ibs, &jbs, &tmp_re, &tmp_im);
    // printf("%ld %ld %ld %ld %ld %ld %ld %lg %lg\n", cntr, a, b, i, j, ibs, jbs, tmp_re, tmp_im); fflush(0);
    // Load the matrix elements
    mat[ibs * ist->n_xton + jbs] = CMPLX(tmp_re, tmp_im);

    if (tmp_a < a)
      tmp_a = a;
    if (tmp_b < b)
      tmp_b = b;
    if (tmp_i < i)
      tmp_i = i;
    if (tmp_j < j)
      tmp_j = j;
    cntr++;
  }

  *a_max = tmp_a;
  *b_max = tmp_b;
  *i_max = tmp_i;
  *j_max = tmp_j;

  // *a_max = a;
  // *b_max = b;
  // *i_max = i;
  // *j_max = j;

  printf("Max value of a = %ld\n", *a_max);
  printf("Max value of b = %ld\n", *b_max);
  printf("Max value of i = %ld\n", *i_max);
  printf("Max value of j = %ld\n", *j_max);

  long el_max = ist->lumo_idx + ist->n_elecs - 1;
  long ho_max = ist->n_holes - 1;
  printf("Max value of el = %ld\n", el_max);
  printf("Max value of ho = %ld\n", ho_max);

  done_flag = 0;
  if ((*a_max == el_max) && (*b_max == el_max) && (*i_max == ho_max) && (*j_max == ho_max))
  {
    done_flag = 1;
  }

  return done_flag;
}

/*****************************************************************************/
