#include "coulomb.h"

/**************************************************************************/
// How many quasiparticle states (each `state_bytes` long) may be resident on
// the GPU at once, given a per-rank memory budget. The budget is taken from the
// environment variable BSE_GPU_MEM_GB (gigabytes of usable device memory for
// this MPI rank; default 32 GB, minus `fixed_bytes` of scratch that is always
// resident). Set BSE_GPU_MEM_GB to (device_GB * fraction / ranks_per_gpu) so
// that ranks sharing a device do not oversubscribe it. Returns at least 1.
static long bse_gpu_states_resident(long state_bytes, long fixed_bytes){
  double gb = 32.0;
  const char* e = getenv("BSE_GPU_MEM_GB");
  if (e && *e) gb = atof(e);
  double avail = gb * 1.0e9 - (double) fixed_bytes;
  if (avail < (double) state_bytes) return 1L;
  long n = (long) (avail / (double) state_bytes);
  return (n < 1L) ? 1L : n;
}

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
	double        *psi_qp, 
	double        *pot_bare,
	double        *pot_screened,
	zomplex       *direct,
	zomplex       *exchange,
	index_st*     ist,
	par_st        *par,
	flag_st       *flag,
	parallel_st   *parallel
	){

  // Warm up the GPU: force a compute target region so the CUDA context is
  // established on the device before any target-data mappings. Without this,
  // in the full BSE process the runtime resolves data regions to the host
  // while compute kernels run on the GPU, leaving mapped data absent.
  {
    volatile int warmup = 0;
    #pragma omp target map(tofrom: warmup)
    { warmup = warmup + 1; }
  }

	/************************************************************/
	/*******************  DECLARE VARIABLES   *******************/
	/************************************************************/

  FILE*             pf;

  //                Indices
  int               cntr = 0;

  long              i;
  long              a;
  long              j;
  long              b;
  long              i_st;
  long              a_st;
  long              b_st;
  long              ibs;
  long              jbs;
  long              start;
  long              loop_idx;
  long              ncycles;
  long              jgr;    // jgrid
  long              jgur;   // jgrid_up_real  
  long              jgui;   // jgrid up imag
  long              jgdr;   // jgrid_dn_real 
  long              jgdi;   // jgrid_dn_imag

  const long        nspngr     = ist->nspinngrid;
  const long        ngrid      = ist->ngrid;
  const long        cplx_idx   = ist->complex_idx;
  const long        lidx       = ist->lumo_idx;
  const long        n_el       = ist->n_elecs;
  const long        n_ho       = ist->n_holes;
  const long        n_xton     = ist->n_xton;
  const long        stlen      = nspngr * cplx_idx;
  const long        cngrid     = ngrid * cplx_idx;
  
  long*             listibs;

  const int         mpir       = parallel->mpi_rank;

  const double      dv         = par->dv;

  double*           rho;
  double*           pot_htree;

  char*             fileName;
  fileName          = (char*) malloc(30*sizeof(char)+1);
  fileName[30]      = '\0';

  rho     = (double *)  calloc(cngrid, sizeof(double));
	listibs = (long *)    calloc(ist->n_xton, sizeof(long));
	
  ALLOCATE(&pot_htree, cngrid, "pot_htree");

	/************************************************************/
  /*******************    INIITIALIZE FFT   *******************/
  /************************************************************/

	// Parallel FFT
  fftw_plan_loc     planfw;
  fftw_plan_loc     planbw;
  fftw_complex*     fftwpsi;
  long fft_flags    = FFTW_MEASURE;

  // Create FFT structs and plans for Fourier transform
  fftw_init_threads();
  fftw_plan_with_nthreads(ist->nthreads);
  
  fftwpsi = fftw_malloc(sizeof (fftw_complex) * ist->ngrid);
  
  planfw = fftw_plan_dft_3d(ist->nz, ist->ny ,ist->nx, 
    fftwpsi, fftwpsi, FFTW_FORWARD,fft_flags
  );
  
  planbw = fftw_plan_dft_3d(ist->nz, ist->ny, ist->nx, 
    fftwpsi, fftwpsi, FFTW_BACKWARD, fft_flags
  );
  
  /************************************************************/
  /*******************    HANDLE INDEXING   *******************/
  /************************************************************/

	for (ibs = 0, a = lidx; a < lidx + n_el; a++) {
		for (i = 0; i < n_ho; i++, ibs++) {
			listibs[(a - lidx) * n_ho + i] = ibs;
		}
	}

  /************************************************************/
  /*******************   CONFIG PARALLEL    *******************/
  /************************************************************/

  // Split MPI communicators into even and odd ranks to simultaneously
  // compute direct and exchange integrals

  int rank_parity   = mpir % 2;  
  MPI_Comm          even_comm = MPI_COMM_NULL;
  MPI_Comm          odd_comm = MPI_COMM_NULL;
  
  MPI_Comm_split(MPI_COMM_WORLD, rank_parity, mpir, (rank_parity == 0) ? &even_comm : &odd_comm);

  if ((rank_parity == 0) && (even_comm == MPI_COMM_NULL)) {
      fprintf(stderr, "ERROR: even_comm is MPI_COMM_NULL\n");
  }
  if ((rank_parity == 1) && (odd_comm == MPI_COMM_NULL)) {
      fprintf(stderr, "ERROR: odd_comm is MPI_COMM_NULL\n");
  }

  // enable odd_comm if only using 1 MPI rank
  // NOTE: broken, still does not work for just 1 rank
  if (parallel->mpi_size == 1){
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
	
	if (rank_parity == 0){
    printf("\nComputing screened direct matrix, K^d_(ab,ij) on rank %d\n", mpir); fflush(0);
    
    /************************************************************/
    /*****************  RECRUIT 1/2 MPI RANKS  ******************/
    /************************************************************/
    
    int even_rank, even_size;
    MPI_Comm_rank(even_comm, &even_rank);
    MPI_Comm_size(even_comm, &even_size);
    if (parallel->mpi_size > 1){
      odd_comm = MPI_COMM_NULL;
    }
    

    long        ab;
    long        ab_tot    = n_el * n_el;
    long       *lista, *listb;

    lista   =  (long *) calloc(ab_tot, sizeof(long));
    listb   =  (long *) calloc(ab_tot, sizeof(long));

    for (a = lidx; a < lidx + n_el; a++){
      for (b = lidx; b < lidx + n_el; b++){
        lista[ (a-lidx)*n_el + (b-lidx) ] = a;
        listb[ (a-lidx)*n_el + (b-lidx)] = b;
      }
    }

    /************************************************************/
    /******************    ASSIGN WORKLOADS   *******************/
    /************************************************************/

    /*** vabji direct ***/
    // loop over electron density (ab) states in a strided manner
    // so that all MPI ranks have similar work loads.
    
    printf("Determining the workload for abtot = %lu on MPI even rank %d\n", ab_tot, even_rank);
    start = even_rank;
    cntr = 0;
    ncycles = ab_tot / even_size;
    printf("Even %lu Number of cycles per rank = %lu\n", even_rank, ncycles);
    
    /************************************************************/
    /******************     ENABLE RESTART    *******************/
    /************************************************************/
    
    sprintf(fileName, "direct-%d.dat", even_rank);
    if (flag->restartCoulomb){
      // Find start value for continuing computation
      start = load_coulomb_mat(direct, fileName, ist);

      // Print out the new matrix elems to this auxiliary file
      sprintf(fileName, "direct-%d_aux.dat", even_rank);

      printf("Even rank %d: continuing direct matrix from a = %lu | %s\n", even_rank, start, get_time()); fflush(0);
    }

    /************************************************************/
    /******************    DO K^D INTEGRAL    *******************/
    /************************************************************/


    // Offload the persistent inputs to the GPU exactly once, before the ab
    // loop:
    //   - the hole states of psi_qp  (states 0 .. n_ho-1, each stlen doubles)
    //   - the exciton index map listibs (n_el*n_ho entries)
    // and allocate device scratch that is refilled every ab iteration:
    //   - pot_htree (complex, cngrid = cplx_idx*ngrid doubles)
    //   - blk, the n_ho x n_ho block of results for the current (a,b)
    //     stored interleaved re,im so only a tiny buffer crosses the PCIe
    //     bus each iteration (the full `direct` matrix stays on the host).
    nvtxRangePushA("Offloading direct-kernel data to device");
    int dev = omp_get_default_device();
    int host_dev = omp_get_initial_device();

    // ---- GPU state batching over the hole dimension ----
    // The device kernel needs a "row" hole tile (ii) and a "col" hole tile (jj)
    // resident at once. When all n_ho holes do not fit we tile the hole
    // dimension. Rather than mapping slices of psi_qp (which trips NVHPC's
    // present check when a tile pair straddles an unmapped gap), we allocate
    // fixed device buffers and stream the state tiles into them with
    // omp_target_memcpy; the kernel indexes the buffers via is_device_ptr.
    const long state_bytes = stlen * (long) sizeof(double);
    const long fixed_bytes = (cngrid + 2L * n_ho * n_ho) * (long) sizeof(double);
    long slots = bse_gpu_states_resident(state_bytes, fixed_bytes);
    long Th;
    if (2*n_ho <= slots) { Th = n_ho; }                        // row+col both fit
    else                 { Th = slots / 2; if (Th < 1) Th = 1; } // split row+col
    long n_tiles = (n_ho + Th - 1) / Th;
    printf("Direct kernel on GPU device %d: n_ho=%ld, %ld resident slots, "
           "tile=%ld -> %ld tile(s), %ld tile-pair(s)\n",
           dev, n_ho, slots, Th, n_tiles, n_tiles*(n_tiles+1)/2); fflush(0);

    // Host result block and persistent device buffers.
    double *blk    = (double *) calloc(2 * Th * Th, sizeof(double));
    double *d_row  = (double *) omp_target_alloc(Th * state_bytes, dev);
    double *d_col  = (double *) omp_target_alloc(Th * state_bytes, dev);
    double *d_pot  = (double *) omp_target_alloc(cngrid * sizeof(double), dev);
    double *d_blk  = (double *) omp_target_alloc(2 * Th * Th * sizeof(double), dev);
    long   *d_lidx = (long   *) omp_target_alloc(n_el * n_ho * sizeof(long), dev);
    if (!d_row || !d_col || !d_pot || !d_blk || !d_lidx){
      fprintf(stderr, "ERROR: omp_target_alloc failed in direct kernel\n"); exit(EXIT_FAILURE);
    }
    omp_target_memcpy(d_lidx, listibs, n_el*n_ho*sizeof(long), 0, 0, dev, host_dev);
    nvtxRangePop();

    pf = fopen(fileName , "w");
    printf("Starting at ab = %lu on even rank %d\n", start, even_rank); fflush(0);

    // Profile expensive loop
    nvtxRangePushA("Computing BSE direct matrix elements");
    for (long ti = 0; ti < n_tiles; ti++) {
      long i0 = ti * Th;
      long iN = i0 + Th; if (iN > n_ho) iN = n_ho;
      long ilen = iN - i0;

      // Copy the row hole-tile to the device once; reused for every col tile.
      omp_target_memcpy(d_row, psi_qp, ilen*state_bytes, 0, i0*state_bytes, dev, host_dev);

      for (long tj = 0; tj <= ti; tj++) {
        long j0 = tj * Th;
        long jN = j0 + Th; if (jN > n_ho) jN = n_ho;
        long jlen = jN - j0;

        // Copy the col hole-tile into its own device buffer.
        omp_target_memcpy(d_col, psi_qp, jlen*state_bytes, 0, j0*state_bytes, dev, host_dev);

        cntr = 0;
        for (ab = start; ab < ab_tot; ab += even_size) {
      a = lista[ab];
      b = listb[ab];

      // Grab indices of electron-electron states a, b
      a_st = a * stlen;
      b_st = b * stlen;

      double psi_aur, psi_aui;
      double psi_adr, psi_adi;
      double psi_bur, psi_bui;
      double psi_bdr, psi_bdi;

      // Compute hartree potential for a, b density
      // 1) Compute joint density and store in rho
      nvtxRangePushA("Computing ab joint density");
      #pragma omp simd  
      for (jgr = 0; jgr < ngrid; jgr++){
        // Index this iteration of the loop
        jgur = cplx_idx * jgr;           // jgrid up real
        jgui = jgur + 1;                 // jgrid up imag
        jgdr = jgur + cplx_idx * ngrid;  // jgrid dn real
        jgdi = jgdr + 1;                 // jgrid dn imag

        // Load in the values for this iteration of the loop
        psi_aur = psi_qp[a_st + jgur];
        psi_aui = psi_qp[a_st + jgui];
        psi_adr = psi_qp[a_st + jgdr];
        psi_adi = psi_qp[a_st + jgdi];
        psi_bur = psi_qp[b_st + jgur];
        psi_bui = psi_qp[b_st + jgui];
        psi_bdr = psi_qp[b_st + jgdr];
        psi_bdi = psi_qp[b_st + jgdi];

        // Handle up spin
        rho[jgur] = psi_aur * psi_bur + psi_aui * psi_bui;
        rho[jgui] = psi_aur * psi_bui - psi_aui * psi_bur;
      
        // Handle down spin
        rho[jgur] += psi_adr * psi_bdr + psi_adi * psi_bdi;
        rho[jgui] += psi_adr * psi_bdi - psi_adi * psi_bdr;
      }
      nvtxRangePop();

      // Compute the hartree potential and store in pot_htree 
      // h_d(r) = \int W(r,r') \rho_{ab}(r') d^3r' via fourier transform
      nvtxRangePushA("Computing hartree pot");
      hartree(rho, pot_screened, pot_htree, ist, planfw, planbw, fftwpsi);            
      nvtxRangePop();

      // Push the freshly computed Hartree potential to the device
      nvtxRangePushA("Updating pot Hartree on GPU");
      omp_target_memcpy(d_pot, pot_htree, cngrid*sizeof(double), 0, 0, dev, host_dev);
      nvtxRangePop();

      // K^d_{ai,bj} = dv * \int h_d(r) \sum_sigma psi_i(r,sigma) psi_j^*(r,sigma) d^3r
      //
      // Parallelization: each (i,j) hole pair is handled by one GPU team, and
      // the grid integral is reduced across the threads within that team. This
      // keeps the grid loop (the innermost, memory-bound work) coalesced.
      // States are read from the resident device tile buffers d_row (index ii)
      // and d_col (index jj), addressed by their local offset within the tile.
      nvtxRangePushA("i,j loop of direct");
      #pragma omp target teams distribute collapse(2) \
          is_device_ptr(d_row, d_col, d_pot, d_blk, d_lidx)
      for (long ii = i0; ii < iN; ii++) {
        for (long jj = j0; jj < jN; jj++) {
          long i_loc = (ii - i0) * stlen;   // offset of state ii within d_row
          long j_loc = (jj - j0) * stlen;   // offset of state jj within d_col
          long ibs  = d_lidx[(a - lidx) * n_ho + ii];
          long jbs  = d_lidx[(b - lidx) * n_ho + jj];

          double sum_re = 0.0, sum_im = 0.0;

          // Compute only the upper triangle (symmetry); the rest is left zero
          // and is never read back on the host.
          if (ibs >= jbs) {
            #pragma omp parallel for reduction(+:sum_re, sum_im)
            for (long jgr = 0; jgr < ngrid; jgr++) {
              long jgur = cplx_idx * jgr;
              long jgui = jgur + 1;
              long jgdr = jgur + cplx_idx * ngrid;
              long jgdi = jgdr + 1;

              // Hartree potential at this grid point
              double pot_h_re = d_pot[jgur];
              double pot_h_im = d_pot[jgui];

              // up spin
              double psi_iur = d_row[i_loc + jgur];
              double psi_iui = d_row[i_loc + jgui];
              double psi_jur = d_col[j_loc + jgur];
              double psi_jui = d_col[j_loc + jgui];
              // dn spin
              double psi_idr = d_row[i_loc + jgdr];
              double psi_idi = d_row[i_loc + jgdi];
              double psi_jdr = d_col[j_loc + jgdr];
              double psi_jdi = d_col[j_loc + jgdi];

              double tmp_re, tmp_im;
              // up spin
              tmp_re = (psi_jur * psi_iur + psi_jui * psi_iui);
              tmp_im = (psi_jur * psi_iui - psi_jui * psi_iur);
              sum_re += (pot_h_re * tmp_re - pot_h_im * tmp_im);
              sum_im += (pot_h_re * tmp_im + pot_h_im * tmp_re);
              // dn spin
              tmp_re = (psi_jdr * psi_idr + psi_jdi * psi_idi);
              tmp_im = (psi_jdr * psi_idi - psi_jdi * psi_idr);
              sum_re += (pot_h_re * tmp_re - pot_h_im * tmp_im);
              sum_im += (pot_h_re * tmp_im + pot_h_im * tmp_re);
            }
            sum_re *= dv;
            sum_im *= dv;
          }

          d_blk[2 * ((ii - i0) * jlen + (jj - j0))]     = sum_re;
          d_blk[2 * ((ii - i0) * jlen + (jj - j0)) + 1] = sum_im;
        } // end of jj
      } // end of ii
      nvtxRangePop();

      // Pull the small block back and scatter it into the host direct matrix.
      omp_target_memcpy(blk, d_blk, 2*ilen*jlen*sizeof(double), 0, 0, host_dev, dev);

      for (i = i0; i < iN; i++){
        for (j = j0; j < jN; j++){
          ibs = listibs[(a - lidx) * n_ho + i];
          jbs = listibs[(b - lidx) * n_ho + j];
          if (ibs < jbs){
            continue;
          }

          direct[ibs * n_xton + jbs].re = blk[2 * ((i - i0) * jlen + (j - j0))];
          direct[ibs * n_xton + jbs].im = blk[2 * ((i - i0) * jlen + (j - j0)) + 1];
        }
      }

      // Print progress (only during the first tile-pair to avoid duplicates)
      if (even_rank == 0 && ti == 0 && tj == 0){
        if ( (cntr == 0) || (0 == cntr % (ncycles/8+1)) || (cntr == (ncycles - 1)) ){
          print_progress_bar(cntr, ncycles);
          fflush(0);
        }
        cntr++;
      }
		} // end of ab

      } // end of tj
    } // end of ti
    nvtxRangePop();

    // Write the completed direct matrix for this rank. Rows are complete once
    // all tile-pairs are done; format matches the original per-(a,b) output so
    // load_coulomb_mat / restart can still read it.
    for (ab = start; ab < ab_tot; ab += even_size) {
      a = lista[ab];
      b = listb[ab];
      for (i = 0; i < n_ho; i++){
        for (j = 0; j < n_ho; j++){
          ibs = listibs[(a - lidx) * n_ho + i];
          jbs = listibs[(b - lidx) * n_ho + j];
          if (ibs < jbs) continue;
          fprintf(pf,"%lu %lu %lu %lu %lu %lu %.16g %.16g\n", a, b, i, j, ibs, jbs, \
              direct[ibs * n_xton + jbs].re, direct[ibs * n_xton + jbs].im
          );
        }
      }
    }
    fflush(0);

		fclose(pf);
    printf("  Done computing direct mat\n");
    fflush(0);

    // Free
    free(lista);
    free(listb);

    // Cleanup persistent GPU memory
    omp_target_free(d_row,  dev);
    omp_target_free(d_col,  dev);
    omp_target_free(d_pot,  dev);
    omp_target_free(d_blk,  dev);
    omp_target_free(d_lidx, dev);

    free(blk);

  } // end of even MPI ranks



	/*******************************************************************/
	/*******************************************************************/
  /****** BREAK ***** BREAK ***** BREAK ***** BREAK ***** BREAK ******/
	/*******************************************************************/
  /*******************************************************************/



	if ( rank_parity == 1 || parallel->mpi_size == 1 ){
    printf("\nComputing bare exchange matrix, K^x_(ai,bj) on rank %d\n", mpir); fflush(0);
	
    /************************************************************/
    /*****************  RECRUIT 1/2 MPI RANKS  ******************/
    /************************************************************/

    int odd_rank, odd_size;
    MPI_Comm_rank(odd_comm, &odd_rank);
    MPI_Comm_size(odd_comm, &odd_size);
    if (parallel->mpi_size > 1){
      even_comm = MPI_COMM_NULL;
    }
    

    omp_set_max_active_levels(1);
    omp_set_num_threads(ist->nthreads);

    long  ai;
    long  ai_tot    = n_el * n_ho;
    // long  ns_p_rank = ai_tot / odd_size;
    long  *lista, *listi;
    
    lista   =  (long *) calloc(ai_tot, sizeof(long));
    listi   =  (long *) calloc(ai_tot, sizeof(long));

    for (a = lidx; a < lidx + n_el; a++){
        for (i = 0; i < n_ho; i++){
            lista[ (a-lidx) * n_ho + i] = a;
            listi[ (a-lidx) * n_ho + i] = i;
        }
    }

    /************************************************************/
    /******************    ASSIGN WORKLOADS  ********************/
    /************************************************************/

    //loop over electron states from start to end
    // handle remainder
    printf("Determining the workload for aitot = %lu on MPI odd rank %d\n", ai_tot, odd_rank);
    
    start = odd_rank;
    cntr = 0;
    ncycles = ai_tot / odd_size;
    printf("ncycles = %lu\n", ncycles);

    
    /************************************************************/
    /******************     ENABLE RESTART   ********************/
    /************************************************************/

    sprintf(fileName, "exchange-%d.dat", odd_rank);

    if (flag->restartCoulomb){
      sprintf(fileName, "exchange-%d.dat", odd_rank);
      
      // Find start value for continuing computation
      start = load_coulomb_mat(exchange, fileName, ist);

      strcpy(fileName, "exchange_aux.dat");
      printf("Odd rank %d: continuing exchange matrix from a = %lu | %s\n", odd_rank, start, get_time()); fflush(0);
    }

    /************************************************************/
    /******************    DO K^X INTEGRAL   ********************/
    /************************************************************/

    // Offload persistent inputs once. The exchange kernel reads both the
    // electron states (index b) and hole states (index j) of psi_qp, so the
    // whole basis up to the last electron (lidx + n_el) must be mapped.
    // pot_htree is refilled every ai iteration; blk_x holds the n_el x n_ho
    // block of results for the current (a,i), interleaved re,im.
    nvtxRangePushA("Offloading exchange-kernel data to device");
    int dev = omp_get_default_device();
    int host_dev = omp_get_initial_device();

    // ---- GPU state batching over electron (bb) and hole (jj) dimensions ----
    // The exchange kernel needs one electron tile and one hole tile resident at
    // once. As in the direct kernel we stream state tiles into fixed device
    // buffers (d_elec, d_hole) with omp_target_memcpy and index them via
    // is_device_ptr; electron tiles persist across the hole tiles.
    const long state_bytes = stlen * (long) sizeof(double);
    const long fixed_bytes = (cngrid + 2L * n_el * n_ho) * (long) sizeof(double);
    long slots = bse_gpu_states_resident(state_bytes, fixed_bytes);
    long Te, Th;
    if (n_el + n_ho <= slots) { Te = n_el; Th = n_ho; }
    else {
      Te = slots / 2; Th = slots / 2;
      if (Te > n_el) Te = n_el;  if (Th > n_ho) Th = n_ho;
      if (Te < 1)    Te = 1;     if (Th < 1)    Th = 1;
    }
    long ne_tiles = (n_el + Te - 1) / Te;
    long nh_tiles = (n_ho + Th - 1) / Th;
    printf("Exchange kernel on GPU device %d: n_el=%ld n_ho=%ld, %ld resident slots, "
           "elec-tile=%ld hole-tile=%ld -> %ld x %ld tiles\n",
           dev, n_el, n_ho, slots, Te, Th, ne_tiles, nh_tiles); fflush(0);

    // Host result block and persistent device buffers.
    double *blk_x   = (double *) calloc(2 * Te * Th, sizeof(double));
    double *d_elec  = (double *) omp_target_alloc(Te * state_bytes, dev);
    double *d_hole  = (double *) omp_target_alloc(Th * state_bytes, dev);
    double *d_pot   = (double *) omp_target_alloc(cngrid * sizeof(double), dev);
    double *d_blkx  = (double *) omp_target_alloc(2 * Te * Th * sizeof(double), dev);
    long   *d_lidx  = (long   *) omp_target_alloc(n_el * n_ho * sizeof(long), dev);
    if (!d_elec || !d_hole || !d_pot || !d_blkx || !d_lidx){
      fprintf(stderr, "ERROR: omp_target_alloc failed in exchange kernel\n"); exit(EXIT_FAILURE);
    }
    omp_target_memcpy(d_lidx, listibs, n_el*n_ho*sizeof(long), 0, 0, dev, host_dev);
    nvtxRangePop();

    pf = fopen(fileName , "w");
    printf("Starting at ai = %lu on odd rank %d\n", start, odd_rank); fflush(0);

    for (long te = 0; te < ne_tiles; te++) {
      long e0 = te * Te;
      long eN = e0 + Te; if (eN > n_el) eN = n_el;
      long elen = eN - e0;

      // Copy this electron tile (physical slots lidx+e0 ..) into d_elec.
      omp_target_memcpy(d_elec, psi_qp, elen*state_bytes, 0, (lidx+e0)*state_bytes, dev, host_dev);

      for (long th = 0; th < nh_tiles; th++) {
        long j0 = th * Th;
        long jN = j0 + Th; if (jN > n_ho) jN = n_ho;
        long jlen = jN - j0;

        // Copy this hole tile into d_hole.
        omp_target_memcpy(d_hole, psi_qp, jlen*state_bytes, 0, j0*state_bytes, dev, host_dev);

        cntr = 0;
        //loop over electron states a, i
        for (ai = start; ai < ai_tot; ai += odd_size) {
      // printf("\n\nai = %lu\n", ai);
      a = lista[ai];
      i = listi[ai];

      a_st = a * stlen;
      i_st = i * stlen;

      double psi_aur, psi_aui;
      double psi_adr, psi_adi;
      double psi_iur, psi_iui;
      double psi_idr, psi_idi;
				
      // 1) Compute joint density and store in rho
      for (jgr = 0; jgr < ngrid; jgr++){
        jgur = cplx_idx * jgr;
        jgui = jgur + 1;
        jgdr = jgur + cplx_idx * ngrid;
        jgdi = jgdr + 1;

        // Load values for this iteration
        psi_aur = psi_qp[a_st + jgur]; // up spin
        psi_aui = psi_qp[a_st + jgui];
        psi_iur = psi_qp[i_st + jgur];
        psi_iui = psi_qp[i_st + jgui];
        psi_adr = psi_qp[a_st + jgdr]; // down spin
        psi_adi = psi_qp[a_st + jgdi];
        psi_idr = psi_qp[i_st + jgdr];
        psi_idi = psi_qp[i_st + jgdi];

        // Handle up spin
        rho[jgur] = psi_aur * psi_iur + psi_aui * psi_iui;
        rho[jgui] = psi_aur * psi_iui - psi_aui * psi_iur;
      
        // Handle down spin
        rho[jgur] += psi_adr * psi_idr + psi_adi * psi_idi;
        rho[jgui] += psi_adr * psi_idi - psi_adi * psi_idr;
      }

      // Compute the hartree potential and store in pot_htree
      // h_d(r) = \int W(r,r') \rho_{ab}(r') d^3r' via fourier transform
      hartree(rho, pot_bare, pot_htree, ist, planfw, planbw, fftwpsi);

      // Push the freshly computed potential to the device
      omp_target_memcpy(d_pot, pot_htree, cngrid*sizeof(double), 0, 0, dev, host_dev);

      // K^x_{ai,bj} = -dv * \int h_x(r) \sum_sigma psi_b(r,sigma) psi_j^*(r,sigma) d^3r
      //
      // One team per electron-hole pair (b,j); the grid integral is reduced
      // across the threads of the team. Electron states are read from the
      // resident device tile d_elec (local index bb-e0, physical b=lidx+bb) and
      // hole states from d_hole (local index jj-j0).
      nvtxRangePushA("b,j loop of exchange");
      #pragma omp target teams distribute collapse(2) \
          is_device_ptr(d_elec, d_hole, d_pot, d_blkx, d_lidx)
      for (long bb = e0; bb < eN; bb++) {
        for (long jj = j0; jj < jN; jj++) {
          long b_loc = (bb - e0) * stlen;   // offset of elec bb within d_elec
          long j_loc = (jj - j0) * stlen;   // offset of hole jj within d_hole
          long ibs  = d_lidx[(a - lidx) * n_ho + i];
          long jbs  = d_lidx[bb * n_ho + jj];

          double sum_re = 0.0, sum_im = 0.0;

          // Compute only the upper triangle (symmetry); the rest stays zero
          if (ibs >= jbs) {
            #pragma omp parallel for reduction(+:sum_re, sum_im)
            for (long jgr = 0; jgr < ngrid; jgr++) {
              long jgur = cplx_idx * jgr;
              long jgui = jgur + 1;
              long jgdr = jgur + cplx_idx * ngrid;
              long jgdi = jgdr + 1;

              // Hartree potential at this grid point
              double pot_h_re = d_pot[jgur];
              double pot_h_im = d_pot[jgui];

              // up spin
              double psi_bur = d_elec[b_loc + jgur];
              double psi_bui = d_elec[b_loc + jgui];
              double psi_jur = d_hole[j_loc + jgur];
              double psi_jui = d_hole[j_loc + jgui];
              // dn spin
              double psi_bdr = d_elec[b_loc + jgdr];
              double psi_bdi = d_elec[b_loc + jgdi];
              double psi_jdr = d_hole[j_loc + jgdr];
              double psi_jdi = d_hole[j_loc + jgdi];

              double tmp_re, tmp_im;
              // up spin
              tmp_re = (psi_jur * psi_bur + psi_jui * psi_bui);
              tmp_im = (psi_jur * psi_bui - psi_jui * psi_bur);
              sum_re += (pot_h_re * tmp_re - pot_h_im * tmp_im);
              sum_im += (pot_h_re * tmp_im + pot_h_im * tmp_re);
              // dn spin
              tmp_re = (psi_jdr * psi_bdr + psi_jdi * psi_bdi);
              tmp_im = (psi_jdr * psi_bdi - psi_jdi * psi_bdr);
              sum_re += (pot_h_re * tmp_re - pot_h_im * tmp_im);
              sum_im += (pot_h_re * tmp_im + pot_h_im * tmp_re);
            }
            sum_re *= dv;
            sum_im *= dv;
          }

          // K^x carries an overall minus sign
          d_blkx[2 * ((bb - e0) * jlen + (jj - j0))]     = -sum_re;
          d_blkx[2 * ((bb - e0) * jlen + (jj - j0)) + 1] = -sum_im;
        } // end of jj
      } // end of bb
      nvtxRangePop();

      // Pull the result block back and scatter into the exchange matrix
      omp_target_memcpy(blk_x, d_blkx, 2*elen*jlen*sizeof(double), 0, 0, host_dev, dev);
      for (b = lidx + e0; b < lidx + eN; b++) {
        long bb = b - lidx;
        for (j = j0; j < jN; j++) {
          long ibs = listibs[(a - lidx) * n_ho + i];
          long jbs = listibs[(b - lidx) * n_ho + j];
          if (ibs < jbs) continue;
          exchange[ibs * n_xton + jbs].re = blk_x[2 * ((bb - e0) * jlen + (j - j0))];
          exchange[ibs * n_xton + jbs].im = blk_x[2 * ((bb - e0) * jlen + (j - j0)) + 1];
        }
      }

      // Print progress (only during the first tile block to avoid duplicates)
      if (odd_rank == 0 && te == 0 && th == 0){
        if ( (cntr == 0) || (0 == cntr % (ncycles/8 + 1)) || (cntr == (ncycles - 1)) ){
            print_progress_bar(cntr, ncycles);
            fflush(0);
        }
        cntr++;
      }
    } // end of ai

      } // end of th
    } // end of te

    // Write the completed exchange matrix for this rank. Format matches the
    // original per-(a,i) output so load_coulomb_mat / restart can still read it.
    for (ai = start; ai < ai_tot; ai += odd_size){
      a = lista[ai];
      i = listi[ai];
      for (b = lidx; b < lidx + n_el; b++){
        for (j = 0; j < n_ho; j++){
          ibs = listibs[(a - lidx) * n_ho + i];
          jbs = listibs[(b - lidx) * n_ho + j];
          if (ibs < jbs) continue;
          fprintf(pf,"%lu %lu %lu %lu %lu %lu %.16g %.16g\n", a, b, i, j, ibs, jbs, \
              exchange[ibs * n_xton + jbs].re, exchange[ibs * n_xton + jbs].im
          );
        }
      }
    }
    fflush(0);

	  fclose(pf);
  	printf("  Done computing exchange mat\n");
	  fflush(0);

    // Cleanup persistent GPU memory
    omp_target_free(d_elec, dev);
    omp_target_free(d_hole, dev);
    omp_target_free(d_pot,  dev);
    omp_target_free(d_blkx, dev);
    omp_target_free(d_lidx, dev);

    free(blk_x);
    free(lista);
    free(listi);

	} // close mpi rank 2

  /************************************************************/
  /****************    SHARE DATA W/ RANKS     ****************/
  /************************************************************/
  MPI_Barrier(MPI_COMM_WORLD);

	// Reduction for even ranks (to rank 0 in even_comm)
  if (rank_parity == 0) {
      // Use MPI_Reduce to sum data from all even ranks into rank 0
      if (mpir == 0){
          MPI_Reduce(MPI_IN_PLACE, direct, 2*sqr(ist->n_xton), MPI_DOUBLE, MPI_SUM, 0, even_comm);
      }else{
          MPI_Reduce(direct, direct, 2*sqr(ist->n_xton), MPI_DOUBLE, MPI_SUM, 0, even_comm);
      }
  }
  if (mpir == 0) printf("Successfully reduced direct mat from even ranks | %s\n", get_time()); fflush(0);

  // Reduction for odd ranks (to rank 1 in odd_comm)
  if (rank_parity == 1) {
      // Use MPI_Reduce to sum data from all odd ranks into rank 1
      if (mpir == 1){
          MPI_Reduce(MPI_IN_PLACE, exchange, 2*sqr(ist->n_xton), MPI_DOUBLE, MPI_SUM, 0, odd_comm);
      }else{
          MPI_Reduce(exchange, exchange, 2*sqr(ist->n_xton), MPI_DOUBLE, MPI_SUM, 0, odd_comm);
      }
  }
  if (mpir == 1) printf("Successfully reduced exchange mat from odd ranks | %s\n", get_time()); fflush(0);

  MPI_Barrier(MPI_COMM_WORLD);

  // If multiple ranks were used to compute the kernel
  // Send exchange data from rank 1 to rank 0
  if (parallel->mpi_size > 1){
    if (mpir == 1) {
        MPI_Send(exchange, 2*sqr(ist->n_xton), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      }
      if (mpir == 0) {
        MPI_Recv(exchange, 2*sqr(ist->n_xton), MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
  } 
  // else{
  //     printf("ERROR: mpi_size < 1 (how are you using MPI code? u gon seg bruh)\n");
  //     exit(EXIT_FAILURE);
  // }

  if (mpir == 0) printf("Successfully sent exchange mat to mpi_rank 0 | %s\n", get_time()); fflush(0);

  MPI_Barrier(MPI_COMM_WORLD);


  // Clean up communicators
  if (even_comm != MPI_COMM_NULL && even_comm != MPI_COMM_WORLD) {
      MPI_Comm_free(&even_comm);
  }
  if (odd_comm != MPI_COMM_NULL && odd_comm != MPI_COMM_WORLD) {
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

long load_coulomb_mat(zomplex* mat, char* fileName, index_st* ist){

  FILE *pf;

  long start;
  long cntr = 0;
  int ieof = 0;
  long a, b, i, j, ibs, jbs;
  long max_st_num = ist->n_holes*ist->n_holes*ist->n_elecs*ist->n_elecs;
  double tmp_re, tmp_im;
  
  // Open the direct/exchange.dat file

  pf = fopen(fileName, "r");

  if (pf == NULL){
      printf("ERROR: could not open file %s\n", fileName);
      fprintf(stderr, "ERROR: could not open file %s\n", fileName);
      exit(EXIT_FAILURE);
  }

  
  // Scan all the lines and load the values into mat
  // Note, order of a,b,i,j doesn't matter because ibs and jbs
  // are the actual indices of the matrices
  while ((ieof != EOF) && (cntr < max_st_num)){
      
      // Scan the file and grab matrix elements
      ieof = fscanf(pf, "%lu %lu %lu %lu %lu %lu %lg %lg", &a, &b, &i, &j, &ibs, &jbs, &tmp_re, &tmp_im);
      
      // Load the matrix elements
      mat[ibs * ist->n_xton + jbs].re = tmp_re;
      mat[ibs * ist->n_xton + jbs].im = tmp_im;

      cntr++;
  }

  // printf("Final value of a = %llu\n", a);
  // printf("Max value of a = %llu\n", max_a);

  start = a;

  return start;
}

/***************************************************************************************/

void build_h0_mat(double *h0mat, double *eval, index_st* ist){
  long a, i, j, ibs;

  FILE *ppsi;
  ibs = 0UL;

  for (a = ist->lumo_idx; a < ist->lumo_idx + ist->n_elecs; a++) {
    for (i = 0; i < ist->n_holes; i++, ibs++) {
        h0mat[ibs*ist->n_xton + ibs] = eval[a] - eval[i];
    }
  }
  
  ppsi = fopen("h0.dat", "w");
  for (i = 0; i < ist->n_xton; i++, fprintf(ppsi,"\n")){
    for (j = 0; j < ist->n_xton; j++){
      fprintf(ppsi,"%.*g ", PR_LEN, h0mat[i*ist->n_xton+j]);
    }
  }
  fclose(ppsi);

  return;
}
/***************************************************************************************/

void build_BSE_mat(zomplex *bsmat, zomplex *direct, zomplex *exchange, index_st* ist){

  FILE *ppsi;
  long ibs, jbs;
  long i, j;

  // Construct the BSE matrix from the exchange and direct kernels
  for (ibs = 0; ibs < ist->n_xton; ibs++){
    for (jbs = 0; jbs <= ibs; jbs++){
      // Symmetrize the matrices
      direct[jbs*ist->n_xton+ibs].re = direct[ibs*ist->n_xton+jbs].re;
      direct[jbs*ist->n_xton+ibs].im = - direct[ibs*ist->n_xton+jbs].im;
      
      exchange[jbs*ist->n_xton+ibs].re = exchange[ibs*ist->n_xton+jbs].re;
      exchange[jbs*ist->n_xton+ibs].im = - exchange[ibs*ist->n_xton+jbs].im;
      
      // Collect values for bsmat
      bsmat[ibs*ist->n_xton+jbs].re = direct[ibs*ist->n_xton+jbs].re + exchange[ibs*ist->n_xton+jbs].re;
      bsmat[jbs*ist->n_xton+ibs].re = direct[jbs*ist->n_xton+ibs].re + exchange[jbs*ist->n_xton+ibs].re;
      
      // Enforce Hermitivity
      if (ibs==jbs){
          bsmat[ibs*ist->n_xton+jbs].im = 0.0;
      } else{
          bsmat[ibs*ist->n_xton+jbs].im = direct[ibs*ist->n_xton+jbs].im + exchange[ibs*ist->n_xton+jbs].im;
          bsmat[jbs*ist->n_xton+ibs].im = direct[jbs*ist->n_xton+ibs].im + exchange[jbs*ist->n_xton+ibs].im;
      }
    }
  }

  
  ppsi = fopen("bsRE.dat", "w");
  for (i = 0; i < ist->n_xton; i++, fprintf(ppsi,"\n")){
    for (j = 0; j < ist->n_xton; j++){
      fprintf(ppsi,"%.*g ", PR_LEN, bsmat[i*ist->n_xton+j].re);
    }
  }
  fclose(ppsi);

  ppsi = fopen("bsIM.dat", "w");
  for (i = 0; i < ist->n_xton; i++, fprintf(ppsi,"\n")){
    for (j = 0; j < ist->n_xton; j++){
      fprintf(ppsi,"%.*g ", PR_LEN, bsmat[i*ist->n_xton+j].im);
    }
  }
  fclose(ppsi);

  return;
}

/*****************************************************************************/
