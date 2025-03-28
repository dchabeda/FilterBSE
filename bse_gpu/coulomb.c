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


    // First, offload the hole states of psi_qp & listibs onto GPU
    nvtxRangePushA("Offloading hole states to device");
    int dev = omp_get_default_device();
    printf("Offloading to GPU device %d\n", dev);
    #pragma omp target enter data map(to: psi_qp[0:nspngr*n_ho], listibs[0:n_ho*n_ho])
    nvtxRangePop();
    // Allocate pot_htree once before ab loop
    nvtxRangePushA("Allocating pot_htree");
    #pragma omp target enter data map(alloc: pot_htree[0:ngrid])
    nvtxRangePop();

    pf = fopen(fileName , "w");
    printf("Starting at ab = %lu on even rank %d\n", start, even_rank); fflush(0);

    // Profile expensive loop
    nvtxRangePushA("Computing BSE direct matrix elements");
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

      nvtxRangePushA("Updating pot Hartree on GPU");
      // Update pot_htree on GPU instead of re-allocating
      #pragma omp target update to(pot_htree[0:ngrid])
      nvtxRangePop();
      // loop over hole states i, j
      nvtxRangePushA("i,j loop of direct");
      
      #pragma omp target teams distribute collapse(2) thread_limit(256) \
      map(to: a, b, i, j, lidx, stlen, n_ho, cplx_idx, ngrid, n_xton, dv, psi_qp, pot_htree, listibs) \
      map(tofrom: direct[0:n_xton * n_xton])
      for (i = 0; i < n_ho; i++) {
        // int jstart = (a == b) ? i: 0;
        // printf("a = %ld b = %ld i = %d jstart = %d\n", a, b, i, jstart);
				for (j = 0; j < n_ho; j++) {
					//get the matrix indicies for {ai,bj}
          long i_st = i * stlen;
          long j_st = j * stlen;
          long ibs = listibs[(a - lidx) * n_ho + i];
          long jbs = listibs[(b - lidx) * n_ho + j];
          
          // Compute only the upper triangle to utilize symmetry
          if (ibs < jbs) continue;

          long   jgr;
          double tmp_re, tmp_im;
          double sum_re, sum_im;
          sum_re = sum_im = 0.0;

          // K^d_{ai,bj}=\int h_d(r) \sum_\sigma psi_{i}(r,\sigma) psi_{j}^{*}(r,\sigma) d^3r
          // Parallelize `jgr` using `#pragma omp parallel for`
          #pragma omp parallel for private(tmp_re, tmp_im) reduction(+:sum_re, sum_im)
          for (jgr = 0; jgr < ngrid; jgr++){
            long jgur = cplx_idx * jgr;
            long jgui = jgur + 1;
            long jgdr = jgur + cplx_idx * ngrid;
            long jgdi = jgdr + 1;

            // Grab pot_htree value at this grid point
            double pot_h_re = pot_htree[jgur];
            double pot_h_im = pot_htree[jgui];

            // Set local values for up spin
            double psi_iur = psi_qp[i_st + jgur];        
            double psi_iui = psi_qp[i_st + jgui];
            double psi_jur = psi_qp[j_st + jgur];        
            double psi_jui = psi_qp[j_st + jgui];
            
            // Set local values for dn spin
            double psi_idr = psi_qp[i_st + jgdr];        
            double psi_idi = psi_qp[i_st + jgdi];
            double psi_jdr = psi_qp[j_st + jgdr];        
            double psi_jdi = psi_qp[j_st + jgdi];
            
            // Perform integrals for up spin
            tmp_re = (psi_jur * psi_iur + psi_jui * psi_iui);
            tmp_im = (psi_jur * psi_iui - psi_jui * psi_iur);
            
            sum_re += (pot_h_re * tmp_re - pot_h_im * tmp_im);
            sum_im += (pot_h_re * tmp_im + pot_h_im * tmp_re);
            
            // Perform integrals for dn spin
            tmp_re = (psi_jdr * psi_idr + psi_jdi * psi_idi);
            tmp_im = (psi_jdr * psi_idi - psi_jdi * psi_idr);
            
            sum_re += (pot_h_re * tmp_re - pot_h_im * tmp_im);
            sum_im += (pot_h_re * tmp_im + pot_h_im * tmp_re);                              
          }

					sum_re *= dv;
					sum_im *= dv;
					
					direct[ibs * n_xton + jbs].re = sum_re;
          direct[ibs * n_xton + jbs].im = sum_im;
				} // end of j
			} // end of i
      nvtxRangePop();

      for (loop_idx = ab; loop_idx < ab + 1; loop_idx++){
        a = lista[loop_idx];
        b = listb[loop_idx];
        for (i = 0; i < n_ho; i++){
          for (j = 0; j < n_ho; j++){
            ibs = listibs[(a - lidx) * n_ho + i];
            jbs = listibs[(b - lidx) * n_ho + j];
            if (ibs < jbs){
              continue;
            }

            fprintf(pf,"%lu %lu %lu %lu %lu %lu %.16g %.16g\n", a, b, i, j, ibs, jbs, \
                direct[ibs * ist->n_xton + jbs].re, direct[ibs * ist->n_xton + jbs].im
            );
          }
        }
      }
      fflush(0);
      

      // Print progress
      if (even_rank == 0){
        if ( (cntr == 0) || (0 == cntr % (ncycles/8+1)) || (cntr == (ncycles - 1)) ){
          print_progress_bar(cntr, ncycles);
          fflush(0);
        }
        cntr++;
      }
		} // end of ab
    
		fclose(pf);
    printf("  Done computing direct mat\n"); 
    fflush(0);
    
    // Free 
    free(lista);
    free(listb);

    // Cleanup GPU memory
    #pragma omp target exit data map(delete: psi_qp[0:nspngr*n_ho], listibs[0:n_ho*n_ho], pot_htree[0:ngrid])  // Free after loop
    nvtxRangePop(); // End marker

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

    pf = fopen(fileName , "w");
    printf("Starting at ai = %lu on odd rank %d\n", start, odd_rank); fflush(0);
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
            
      // loop over electron-hole pairs b, j
      #pragma omp parallel for private(b, j)
      for (b = lidx; b < lidx + n_el; b++) {
        for (j = 0; j < n_ho; j++) {
          long b_st = b * stlen;
          long j_st = j * stlen;
          long ibs = listibs[(a-lidx) * n_ho + i];
          long jbs = listibs[(b-lidx) * n_ho + j];
          
          // Compute only the upper triangle to utilize symmetry
          if (ibs < jbs) continue;
          
          long   jgr;
          double tmp_re, tmp_im;
          double sum_re, sum_im;
          sum_re = sum_im = 0.0;

          // Declare local psi_qp values for reduced mem lookup
          double psi_bur, psi_bui;
          double psi_jur, psi_jui;
          double psi_bdr, psi_bdi;
          double psi_jdr, psi_jdi;
          double pot_h_re, pot_h_im;

          //integrate the effective potential to get K^x_{ai,bj}=\int h_x(r) \sum_\sigma psi_{b}(r,\sigma) psi_{j}^{*}(r,\sigma) d^3r
          for (jgr = 0; jgr < ngrid; jgr++){
            jgur = cplx_idx * jgr;
            jgui = jgur + 1;
            jgdr = jgur + cplx_idx * ngrid;
            jgdi = jgdr + 1;

            // Grab pot_htree value at this grid point
            pot_h_re = pot_htree[jgur];
            pot_h_im = pot_htree[jgui];

            // Set local values for up spin
            psi_bur = psi_qp[b_st + jgur];         psi_bui = psi_qp[b_st + jgui];
            psi_jur = psi_qp[j_st + jgur];         psi_jui = psi_qp[j_st + jgui];
            
            // Set local values for dn spin
            psi_bdr = psi_qp[b_st + jgdr];        psi_bdi = psi_qp[b_st + jgdi];
            psi_jdr = psi_qp[j_st + jgdr];        psi_jdi = psi_qp[j_st + jgdi];
            
            // Perform integrals for up spin
            tmp_re = (psi_jur * psi_bur + psi_jui * psi_bui);
            tmp_im = (psi_jur * psi_bui - psi_jui * psi_bur);
            
            sum_re += (pot_h_re * tmp_re -  pot_h_im * tmp_im);
            sum_im += (pot_h_re * tmp_im +  pot_h_im * tmp_re);
            
            // Handle dn spin
            tmp_re = (psi_jdr * psi_bdr + psi_jdi * psi_bdi);
            tmp_im = (psi_jdr * psi_bdi - psi_jdi * psi_bdr);
            
            sum_re += (pot_h_re * tmp_re -  pot_h_im * tmp_im);
            sum_im += (pot_h_re * tmp_im +  pot_h_im * tmp_re);                              
          }
                
          sum_re *= par->dv;
          sum_im *= par->dv;

          exchange[ibs * n_xton + jbs].re = - 1.0 * sum_re;
          exchange[ibs * n_xton + jbs].im = - 1.0 * sum_im;
        } // end of b
      } // end of j
      
      // Print progress
      if (odd_rank == 0){
        if ( (cntr == 0) || (0 == cntr % (ncycles/8 + 1)) || (cntr == (ncycles - 1)) ){
            print_progress_bar(cntr, ncycles);
            fflush(0);
        }
        cntr++;
      }

      for (loop_idx = ai; loop_idx < ai + 1; loop_idx++){
        // printf("This is ai: %lu\n", ai);
        a = lista[loop_idx];
        i = listi[loop_idx];
        for (b = lidx; b < lidx + n_el; b++){
          for (j = 0; j < n_ho; j++){
            // printf("a %lu b %lu i %lu j %lu\n", a, b, i, j);
            ibs = listibs[(a - lidx) * n_ho + i];
            jbs = listibs[(b - lidx) * n_ho + j];
            // printf("This is ibs = %lu and jbs = %lu\n", ibs, jbs);
            if (ibs < jbs){
              continue;
            }
            // printf("printing exchange %lu %lu %lu %lu\n", loop_idx, b, i, j);
            fprintf(pf,"%lu %lu %lu %lu %lu %lu %.16g %.16g\n", a, b, i, j, ibs, jbs, \
                exchange[ibs * ist->n_xton + jbs].re, exchange[ibs * ist->n_xton + jbs].im
            );
          }
        }
      }
      fflush(0);
      
    } // end of ai

	  fclose(pf);
  	printf("  Done computing exchange mat\n"); 
	  fflush(0);

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
