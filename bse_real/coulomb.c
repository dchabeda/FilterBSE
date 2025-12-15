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

void calc_eh_kernel_cplx(
    double complex *restrict psi_qp,
    double complex *restrict pot_bare,
    double complex *restrict pot_screened,
    double complex *restrict direct,
    double complex *restrict exchange,
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
  double complex *restrict pot_htree;

  ALLOCATE(&rho, ngrid, "rho in coulomb");
  ALLOCATE(&listibs, ist->n_xton, "listibs in coulomb");
  ALLOCATE(&pot_htree, nspngr, "pot_htree");
  // Note: strictly, only ngrid elements are needed for pot_htree
  // but in order to improve loop indexing and vectorization of
  // Coulomb integrals, the first ngrid elements will be duplicated
  // in the remaining ngrid memory block. This enables a seamless
  // integral over up and down spin parts of the wavefunction

  /************************************************************/
  /*******************    INIITIALIZE FFT   *******************/
  /************************************************************/

  // Parallel FFT
  fftw_plan_loc planfw;
  fftw_plan_loc planbw;
  fftw_complex *fftwpsi;
  long fft_flags = FFTW_MEASURE;

  // Create FFT structs and plans for Fourier transform
  fftw_init_threads();
  fftw_plan_with_nthreads(ist->nthreads);

  fftwpsi = fftw_malloc(sizeof(fftw_complex) * ngrid);

  planfw = fftw_plan_dft_3d(ist->nz, ist->ny, ist->nx,
                            fftwpsi, fftwpsi, FFTW_FORWARD, fft_flags);

  planbw = fftw_plan_dft_3d(ist->nz, ist->ny, ist->nx,
                            fftwpsi, fftwpsi, FFTW_BACKWARD, fft_flags);

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

  /************************************************************/
  /*******************   CONFIG PARALLEL    *******************/
  /************************************************************/

  // Split MPI communicators into even and odd ranks to simultaneously
  // compute direct and exchange integrals

  int rank_parity = mpir % 2;
  MPI_Comm even_comm = MPI_COMM_NULL;
  MPI_Comm odd_comm = MPI_COMM_NULL;

  MPI_Comm_split(MPI_COMM_WORLD, rank_parity, mpir, (rank_parity == 0) ? &even_comm : &odd_comm); // these labels are switched to give exchange int xtra rank

  if ((rank_parity == 0) && (even_comm == MPI_COMM_NULL))
  {
    fprintf(stderr, "ERROR: even_comm is MPI_COMM_NULL\n");
  }
  if ((rank_parity == 1) && (odd_comm == MPI_COMM_NULL))
  {
    fprintf(stderr, "ERROR: odd_comm is MPI_COMM_NULL\n");
  }

  // enable odd_comm if only using 1 MPI rank
  // NOTE: broken, still does not work for just 1 rank
  if (parallel->mpi_size == 1)
  {
    MPI_Comm_split(MPI_COMM_WORLD, rank_parity, mpir, &odd_comm);
  }

  // Set the OpenMP parallelization to use nthreads
  omp_set_max_active_levels(1);
  omp_set_num_threads(ist->nthreads);

  /************************************************************/
  /*******************    CALC DIRECT K^D   *******************/
  /************************************************************/

  /*** vabji direct ***/
  // We avoid performing additional computational work
  // 2e integrals have a 2-fold permutation symmetry if cplx
  // [ij|ab] = [ji|ba]^*
  // 4-fold if real!
  // [ij|ab] = [ji|ab] = [ji|ba] = [ij|ba]
  // using Chemist's notation from Szabo and Ostlund

  if (rank_parity == 0)
  {

    /************************************************************/
    /*****************  RECRUIT 1/2 MPI RANKS  ******************/
    /************************************************************/

    int even_rank, even_size;
    MPI_Comm_rank(even_comm, &even_rank);
    MPI_Comm_size(even_comm, &even_size);
    if (parallel->mpi_size > 1)
    {
      odd_comm = MPI_COMM_NULL;
    }
    if (even_rank == 0)
    {
      printf("Screened direct matrix K^d | %s\n", get_time());
      fflush(0);
    }
    int chunk_size;

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
    /******************    ASSIGN WORKLOADS   *******************/
    /************************************************************/
    // Adjust chunk size based on problem size and thread count
    chunk_size = get_dynamic_process_workload(ij_tot);

    /*** vabji direct ***/
    // loop over electron density (ab) states in a strided manner
    // so that all MPI ranks have similar work loads.

    start = even_rank;
    ncycles = ab_tot / even_size;

    if (even_rank == 0)
    {
      printf("Determining workload when a < b, abtot = %lu | MPI even rank %d\n", ab_tot, even_rank);
      printf("\tncycles = %lu\n", ncycles);
      printf("\tinner loop %ld; %d threads; chunk size %d\n", ij_tot, omp_get_max_threads(), chunk_size);
    }

    /************************************************************/
    /******************     ENABLE RESTART    *******************/
    /************************************************************/

    sprintf(fileName, "direct-%d.dat", even_rank);
    if (flag->restartCoulomb)
    {
      int done_flag;
      long a_tmp, b_tmp, i_tmp, j_tmp;

      // Find start value for continuing computation
      done_flag = load_coulomb_mat(direct, fileName, &a_tmp, &b_tmp, &i_tmp, &j_tmp, ist);

      if (done_flag == 0)
      {
        // Calc the new starting value from the loop trip count
        start = 0;
        for (ab = 0; ab < ab_tot; ab++)
        {
          if ((lista[ab] == a_tmp) && (listb[ab] == b_tmp))
          {
            break;
          }
          start++;
        }
      }
      else
      {
        start = ab_tot;
      }

      // Print out the new matrix elems to this auxiliary file
      sprintf(fileName, "direct-%d_aux.dat", even_rank);
      printf("Even rank %d: continuing direct matrix from ab = %ld a = %ld b = %ld | %s\n",
             even_rank, start, a_tmp, b_tmp, get_time());

      fflush(0);
    }
    else
    {
      if (even_rank == 0)
      {
        printf("Starting at ab = %ld on even rank %d | %s\n", start, even_rank, get_time());
        fflush(0);
      }
    }

    /************************************************************/
    /******************    DO K^D INTEGRAL    *******************/
    /************************************************************/

    pf = fopen(fileName, "w");

    if (even_rank == 0)
    {
      printf("Computing direct mat | %s\n", get_time());
    }
    // Profile expensive loop
    // nvtxRangePushA("Computing BSE direct matrix elements");
    cntr = 0;
    double ab_start_t = omp_get_wtime();
    for (ab = start; ab < ab_tot; ab += even_size)
    {
      a = lista[ab];
      b = listb[ab];
      // printf("\n even rank %d ab = %ld a = %ld b = %ld\n", even_rank, ab, a, b); fflush(0);
      // Grab indices of electron-electron states a, b
      a_st = a * nspngr;
      b_st = b * nspngr;

      // Compute hartree potential for a, b density
      // 1) Compute joint density and store in rho
      // // nvtxRangePushA("Computing ab joint density");

      for (jg = 0; jg < ngrid; jg++)
      {
        rho[jg] = conjmul(psi_qp[a_st + jg], psi_qp[b_st + jg]) + conjmul(psi_qp[a_st + jg + ngrid], psi_qp[b_st + jg + ngrid]);
      }

      // // nvtxRangePop();

      // Compute the hartree potential and store in pot_htree
      // h_d(r) = \int W(r,r') \rho_{ab}(r') d^3r' via fourier transform
      // // nvtxRangePushA("Computing hartree pot");
      hartree(rho, pot_screened, pot_htree, ist, planfw, planbw, fftwpsi);
// // nvtxRangePop();

// loop over hole states i, j
// // nvtxRangePushA("i,j loop of direct");
#pragma omp parallel for private(i, j) // schedule(dynamic, chunk_size)
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
        double complex sum;
        sum = 0.0 + 0.0 * I;

        // K^d_{ai,bj}=\int h_d(r) \sum_\sigma psi_{i}(r,\sigma) psi_{j}^{*}(r,\sigma) d^3r
        for (jg = 0; jg < ngrid; jg++)
        {
          sum += pot_htree[jg] *
                 (conjmul(psi_qp[j_st + jg], psi_qp[i_st + jg]) +
                  conjmul(psi_qp[j_st + jg + ngrid], psi_qp[i_st + jg + ngrid]));
        }
        sum *= dv;

        direct[ibs * n_xton + jbs] = sum;
      } // end of ij
      // // nvtxRangePop();

      // fflush(0);
      for (ij = 0; ij < ij_tot; ij++)
      {
        i = listi[ij];
        j = listj[ij];
        ibs = listibs[(a - lidx) * n_ho + i];
        jbs = listibs[(b - lidx) * n_ho + j];

        fprintf(pf, "%03ld %03ld %03ld %03ld %ld %ld %.12f %.12f\n", a, b, i, j, ibs, jbs,
                creal(direct[ibs * n_xton + jbs]), cimag(direct[ibs * n_xton + jbs]));
      }
      // Every 25% of iterations, print output
      if ((cntr == 0) || (0 == cntr % (ncycles / 4 + 1)) || (cntr == (ncycles - 1)))
      {
        // Print out progress
        if (even_rank == 0)
        {
          printf("Direct: ");
          print_progress_bar(cntr, ncycles);
        }
        fflush(0);
      }
      cntr++;
    } // end of ab
    double ab_end_t = omp_get_wtime();

    if (even_rank == 0)
    {
      printf("Done with direct; duration = %lg s (%lg m)\n", (ab_end_t - ab_start_t), (ab_end_t - ab_start_t) / 60.0);
    }

    // Free
    free(lista);
    free(listb);
    free(listi);
    free(listj);
  } // end of even MPI ranks

  /*******************************************************************/
  /*******************************************************************/
  /****** BREAK ***** BREAK ***** BREAK ***** BREAK ***** BREAK ******/
  /*******************************************************************/
  /*******************************************************************/

  if (rank_parity == 1 || parallel->mpi_size == 1)
  {

    if (parallel->mpi_size == 1)
    {
      for (jg = 0; jg < ngrid; jg++)
        rho[jg] = 0.0 + I * 0.0;
    }
    /************************************************************/
    /*****************  RECRUIT 1/2 MPI RANKS  ******************/
    /************************************************************/

    int odd_rank, odd_size;
    MPI_Comm_rank(odd_comm, &odd_rank);
    MPI_Comm_size(odd_comm, &odd_size);
    if (parallel->mpi_size > 1)
    {
      even_comm = MPI_COMM_NULL;
    }
    if (odd_rank == 0)
    {
      printf("\nBare exchange matrix K^x | %s\n", get_time());
      fflush(0);
    }

    omp_set_max_active_levels(1);
    omp_set_num_threads(ist->nthreads);

    int chunk_size;

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

    /************************************************************/
    /******************    ASSIGN WORKLOADS  ********************/
    /************************************************************/

    // Adjust chunk size based on problem size and thread count
    chunk_size = get_dynamic_process_workload(bj_tot);

    // loop over electron states from start to end
    start = odd_rank;
    ncycles = ai_tot / odd_size;

    if (odd_rank == 0)
    {
      printf("\nDetermining the workload for aitot = %lu on MPI odd rank %d\n", ai_tot, odd_rank);
      printf("\tExchange ncycles = %lu\n", ncycles);
      printf("\t%ld inner loop; %d threads; chunk size %d\n", bj_tot, omp_get_max_threads(), chunk_size);
      fflush(0);
    }

    /************************************************************/
    /******************     ENABLE RESTART   ********************/
    /************************************************************/

    sprintf(fileName, "exchange-%d.dat", odd_rank);

    if (flag->restartCoulomb)
    {
      int done_flag;
      long a_tmp, b_tmp, i_tmp, j_tmp;

      // Find start value for continuing computation
      done_flag = load_coulomb_mat(exchange, fileName, &a_tmp, &b_tmp, &i_tmp, &j_tmp, ist);

      start = 0; // Calc the new starting value from the loop trip count
      for (ai = 0; ai < ai_tot; ai++)
      {
        if ((lista[ai] == a_tmp) && (listi[ai] == i_tmp))
        {
          break;
        }
        start++;
      }

      // Print out the new matrix elems to this auxiliary file
      sprintf(fileName, "exchange-%d_aux.dat", odd_rank);
      printf("Odd rank %d: continuing exchange matrix from ai = %ld a = %ld i = %ld | %s\n",
             odd_rank, start, a_tmp, i_tmp, get_time());

      fflush(0);
    }
    else
    {
      if (odd_rank == 0)
      {
        printf("Starting at ai = %ld on odd rank %d | %s\n", start, odd_rank, get_time());
        fflush(0);
      }
    }

    /************************************************************/
    /******************    DO K^X INTEGRAL   ********************/
    /************************************************************/

    pf = fopen(fileName, "w");

    // loop over electron states a, i
    cntr = 0;
    double ai_start_t = omp_get_wtime();
    for (ai = start; ai < ai_tot; ai += odd_size)
    {

      a = lista[ai];
      i = listi[ai];
      // printf("\n odd rank %d ai = %ld a = %ld i = %ld\n", odd_rank, ai, a, i); fflush(0);
      a_st = a * nspngr;
      i_st = i * nspngr;

      // 1) Compute joint density and store in rho
      for (jg = 0; jg < ngrid; jg++)
      {
        rho[jg] = conjmul(psi_qp[a_st + jg], psi_qp[i_st + jg]) + conjmul(psi_qp[a_st + jg + ngrid], psi_qp[i_st + jg + ngrid]);
      }

      // Compute the hartree potential and store in pot_htree
      // h_d(r) = \int W(r,r') \rho_{ab}(r') d^3r' via fourier transform
      hartree(rho, pot_bare, pot_htree, ist, planfw, planbw, fftwpsi);

// loop over electron-hole pairs b, j
#pragma omp parallel for private(b, j) // schedule(dynamic, chunk_size)
      for (bj = 0; bj < bj_tot; bj++)
      {
        b = listb[bj];
        j = listj[bj];

        long b_st = b * nspngr;
        long j_st = j * nspngr;
        long ibs = listibs[(a - lidx) * n_ho + i];
        long jbs = listibs[(b - lidx) * n_ho + j];

        if (ibs < jbs)
          continue;

        long jg;
        double complex sum;
        sum = 0.0 + 0.0 * I;

        for (jg = 0; jg < ngrid; jg++)
        {
          sum += pot_htree[jg] *
                 (conjmul(psi_qp[j_st + jg], psi_qp[b_st + jg]) +
                  conjmul(psi_qp[j_st + jg + ngrid], psi_qp[b_st + jg + ngrid]));
        }
        sum *= dv;

        exchange[ibs * n_xton + jbs] = -sum;
      } // end of bj

      for (bj = 0; bj < bj_tot; bj++)
      {
        b = listb[bj];
        j = listj[bj];
        ibs = listibs[(a - lidx) * n_ho + i];
        jbs = listibs[(b - lidx) * n_ho + j];

        if (ibs < jbs)
          continue;
        // printf("printing exchange %lu %lu %lu %lu\n", loop_idx, b, i, j);
        fprintf(pf, "%03ld %03ld %03ld %03ld %ld %ld %.12f %.12f\n", a, b, i, j, ibs, jbs,
                creal(exchange[ibs * ist->n_xton + jbs]), cimag(exchange[ibs * ist->n_xton + jbs]));
        // fprintf(pf, "%ld %ld %ld %ld %ld %ld %.16g %.16g\n", a, b, i, j, ibs, jbs,
        //         creal(exchange[ibs * ist->n_xton + jbs]), cimag(exchange[ibs * ist->n_xton + jbs]));
      }

      // Every 25% of iterations, print the job progress
      if ((cntr == 0) || (0 == cntr % (ncycles / 4 + 1)) || (cntr == (ncycles - 1)))
      {
        // Print exchange mat elements
        // for (loop_idx = last_ai; loop_idx <= ai; loop_idx++){
        //   if ((loop_idx - start) % odd_size != 0) continue;
        //   a = lista[loop_idx];
        //   i = listi[loop_idx];
        //   for (bj = 0; bj < bj_tot; bj++){
        //     b = listb[bj];
        //     j = listj[bj];
        //     ibs = listibs[(a - lidx) * n_ho + i];
        //     jbs = listibs[(b - lidx) * n_ho + j];

        //     if (ibs < jbs) continue;
        //     // printf("printing exchange %lu %lu %lu %lu\n", loop_idx, b, i, j);
        //     fprintf(pf,"%ld %ld %ld %ld %ld %ld %.16g %.16g\n", a, b, i, j, ibs, jbs, \
        //         creal(exchange[ibs * ist->n_xton + jbs]), cimag(exchange[ibs * ist->n_xton + jbs])
        //     );
        //   }
        // }
        // last_ai = ai + 1;
        // Print progress
        if (odd_rank == 0)
        {
          printf("Exchange: ");
          print_progress_bar(cntr, ncycles);
        }
        fflush(0);
      }

      cntr++;
    } // end of ai
    double ai_end_t = omp_get_wtime();
    fclose(pf);

    if (odd_rank == 0)
    {
      printf("Done with exchange integrals on rank %d; duration = %lg s (%lg m)\n", odd_rank, (ai_end_t - ai_start_t), (ai_end_t - ai_start_t) / 60.0);
      fflush(0);
    }

    free(lista);
    free(listi);
    free(listb);
    free(listj);

  } // close mpi rank 2

  /************************************************************/
  /****************    SHARE DATA W/ RANKS     ****************/
  /************************************************************/
  MPI_Barrier(MPI_COMM_WORLD);

  // Reduction for even ranks (to rank 0 in even_comm)
  if (rank_parity == 0)
  {
    // Use MPI_Reduce to sum data from all even ranks into rank 0
    if (mpir == 0)
    {
      MPI_Reduce(MPI_IN_PLACE, direct, 2 * sqr(ist->n_xton), MPI_DOUBLE, MPI_SUM, 0, even_comm);
    }
    else
    {
      MPI_Reduce(direct, direct, 2 * sqr(ist->n_xton), MPI_DOUBLE, MPI_SUM, 0, even_comm);
    }
  }
  if (mpir == 0)
    printf("\n\nSuccessfully reduced direct mat from even ranks | %s\n", get_time());
  fflush(0);

  // Reduction for odd ranks (to rank 1 in odd_comm)
  if (rank_parity == 1)
  {
    // Use MPI_Reduce to sum data from all odd ranks into rank 1
    if (mpir == 1)
    {
      MPI_Reduce(MPI_IN_PLACE, exchange, 2 * sqr(ist->n_xton), MPI_DOUBLE, MPI_SUM, 0, odd_comm);
    }
    else
    {
      MPI_Reduce(exchange, exchange, 2 * sqr(ist->n_xton), MPI_DOUBLE, MPI_SUM, 0, odd_comm);
    }
  }
  if (mpir == 1)
    printf("Successfully reduced exchange mat from odd ranks | %s\n", get_time());
  fflush(0);

  MPI_Barrier(MPI_COMM_WORLD);

  // If multiple ranks were used to compute the kernel
  // Send exchange data from rank 1 to rank 0
  if (parallel->mpi_size > 1)
  {
    if (mpir == 1)
    {
      MPI_Send(exchange, 2 * sqr(ist->n_xton), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    if (mpir == 0)
    {
      MPI_Recv(exchange, 2 * sqr(ist->n_xton), MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }
  // else{
  //     printf("ERROR: mpi_size < 1 (how are you using MPI code? u gon seg bruh)\n");
  //     exit(EXIT_FAILURE);
  // }

  if (mpir == 0)
    printf("Successfully sent exchange mat to mpi_rank 0 | %s\n", get_time());
  fflush(0);

  MPI_Barrier(MPI_COMM_WORLD);

  // Clean up communicators
  if (even_comm != MPI_COMM_NULL && even_comm != MPI_COMM_WORLD)
  {
    MPI_Comm_free(&even_comm);
  }
  if (odd_comm != MPI_COMM_NULL && odd_comm != MPI_COMM_WORLD)
  {
    MPI_Comm_free(&odd_comm);
  }

  free(rho);
  free(listibs);
  free(fileName);
  free(pot_htree);

  fftw_free(fftwpsi);
  fftw_destroy_plan(planfw);
  fftw_destroy_plan(planbw);

  return;
}

/***************************************************************************************/

int load_coulomb_mat(double complex *mat, char *fileName, long *a_max, long *b_max, long *i_max, long *j_max, index_st *ist)
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
