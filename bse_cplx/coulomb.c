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
	double complex* restrict psi_qp, 
	double complex* restrict pot_bare,
	double complex* restrict pot_screened,
	double complex* restrict direct,
	double complex* restrict exchange,
	index_st*       ist,
	par_st*         par,
	flag_st*        flag,
	parallel_st*    parallel
	){
	
	/************************************************************/
	/*******************  DECLARE VARIABLES   *******************/
	/************************************************************/

  FILE*             pf;

  const int         mpir       = parallel->mpi_rank;

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
  long              jg;
  long              jsg;

  const long        nspngr     = ist->nspinngrid;
  const long        ngrid      = ist->ngrid;
  const long        lidx       = ist->lumo_idx;
  const long        n_el       = ist->n_elecs;
  const long        n_ho       = ist->n_holes;
  const long        n_xton     = ist->n_xton;
  
  long*             listibs;

  const double      dv         = par->dv;

  char*             fileName;
  fileName          = (char*) malloc(30*sizeof(char)+1);
  fileName[30]      = '\0';

  double complex* restrict   rho;
  double complex* restrict   pot_htree;

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
  fftw_plan_loc     planfw;
  fftw_plan_loc     planbw;
  fftw_complex*     fftwpsi;
  long fft_flags    = FFTW_MEASURE;

  // Create FFT structs and plans for Fourier transform
  fftw_init_threads();
  fftw_plan_with_nthreads(ist->nthreads);
  
  fftwpsi = fftw_malloc(sizeof (fftw_complex) * ngrid);
  
  planfw  = fftw_plan_dft_3d(ist->nz, ist->ny ,ist->nx, 
    fftwpsi, fftwpsi, FFTW_FORWARD,fft_flags
  );
  
  planbw  = fftw_plan_dft_3d(ist->nz, ist->ny, ist->nx, 
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

    pf = fopen(fileName , "w");
    printf("Starting at ab = %lu on even rank %d\n", start, even_rank); fflush(0);

    // Profile expensive loop
    nvtxRangePushA("Computing BSE direct matrix elements");
    for (ab = start; ab < ab_tot; ab += even_size) {
      a = lista[ab];
      b = listb[ab];

      // Grab indices of electron-electron states a, b
      a_st = a * nspngr;
      b_st = b * nspngr;

      // Compute hartree potential for a, b density
      // 1) Compute joint density and store in rho
      nvtxRangePushA("Computing ab joint density");

      #pragma omp simd safelen(8) aligned(rho, psi_qp: BYTE_BOUNDARY)
      for (jg = 0; jg < ngrid; jg++){
        rho[jg] = conjmul(psi_qp[a_st + jg], psi_qp[b_st + jg]);
      }
      for (jg = 0; jg < ngrid; jg++){
        jsg = jg + ngrid;
        rho[jg] += conjmul(psi_qp[a_st + jsg], psi_qp[b_st + jsg]);
      }
      
      nvtxRangePop();

      // Compute the hartree potential and store in pot_htree 
      // h_d(r) = \int W(r,r') \rho_{ab}(r') d^3r' via fourier transform
      nvtxRangePushA("Computing hartree pot");
      hartree(rho, pot_screened, pot_htree, ist, planfw, planbw, fftwpsi);            
      nvtxRangePop();
      
      // Copy the "up" component of pot_htree to the "dn" for seamless integration
      memcpy(&pot_htree[ngrid], &pot_htree[0], ngrid * sizeof(pot_htree[0]));
      
      // loop over hole states i, j
      nvtxRangePushA("i,j loop of direct");
      
      #pragma omp parallel for private (j, jsg)
      for (i = 0; i < n_ho; i++) {
        for (j = 0; j < n_ho; j++) {
					//get the matrix indicies for {ai,bj}
          long i_st = i * nspngr;
          long j_st = j * nspngr;
          long ibs = listibs[(a - lidx) * n_ho + i];
          long jbs = listibs[(b - lidx) * n_ho + j];
          
          // Compute only the upper triangle to utilize symmetry
          if (ibs < jbs) continue;

          long   jg;
          double complex sum;
          sum = 0.0 + 0.0 * I;
 
          // K^d_{ai,bj}=\int h_d(r) \sum_\sigma psi_{i}(r,\sigma) psi_{j}^{*}(r,\sigma) d^3r
          
          #pragma omp simd safelen(8) aligned(psi_qp, pot_htree: BYTE_BOUNDARY) reduction(+: sum)
          for (jg = 0; jg < nspngr; jg++){
            sum += pot_htree[jg] * conjmul(psi_qp[j_st + jg], psi_qp[i_st + jg]);                          
          }
          
					sum *= dv;
					
					direct[ibs * n_xton + jbs] = sum;
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
                direct[ibs * ist->n_xton + jbs]
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
    
    if (parallel->mpi_size == 1){
      for (jg = 0; jg < ngrid; jg++) rho[jg] = 0.0 + I*0.0;
    }
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

      a_st = a * nspngr;
      i_st = i * nspngr;

      // 1) Compute joint density and store in rho
      // Up spin
      for (jg = 0; jg < ngrid; jg++){
        rho[jg] = conjmul(psi_qp[a_st + jg], psi_qp[i_st + jg]);
      }
      // Down spin
      for (jg = 0; jg < ngrid; jg++){
        jsg = jg + ngrid;
        // Handle up spin
        rho[jg] += conjmul(psi_qp[a_st + jsg], psi_qp[i_st + jsg]);
      }
      // Compute the hartree potential and store in pot_htree 
      // h_d(r) = \int W(r,r') \rho_{ab}(r') d^3r' via fourier transform
      hartree(rho, pot_bare, pot_htree, ist, planfw, planbw, fftwpsi);

      // Copy "up" component of pot_htree to "dn" for seamless integration
      memcpy(&pot_htree[ngrid], &pot_htree[0], ngrid * sizeof(pot_htree[0]));

      // loop over electron-hole pairs b, j
      #pragma omp parallel for private(b, j)
      for (b = lidx; b < lidx + n_el; b++) {
        for (j = 0; j < n_ho; j++) {
          long b_st = b * nspngr;
          long j_st = j * nspngr;
          long ibs = listibs[(a-lidx) * n_ho + i];
          long jbs = listibs[(b-lidx) * n_ho + j];
          
          // Compute only the upper triangle to utilize symmetry
          if (ibs < jbs) continue;
          
          long   jg;
          long   jsg;
          double complex sum;
          sum = 0.0 + 0.0 * I;

          //integrate the effective potential to get K^x_{ai,bj}=\int h_x(r) \sum_\sigma psi_{b}(r,\sigma) psi_{j}^{*}(r,\sigma) d^3r
          
          for (jg = 0; jg < nspngr; jg++){
            sum += pot_htree[jg] * conjmul(psi_qp[j_st + jg], psi_qp[b_st + jg]);                            
          }
                
          sum *= dv;

          exchange[ibs * n_xton + jbs] = - sum;
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
                exchange[ibs * ist->n_xton + jbs]
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

long load_coulomb_mat(double complex* mat, char* fileName, index_st* ist){

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
      mat[ibs * ist->n_xton + jbs] = CMPLX(tmp_re, tmp_im);

      cntr++;
  }
 
  // printf("Final value of a = %llu\n", a);
  // printf("Max value of a = %llu\n", max_a);

  start = a;

  return start;
}


/*****************************************************************************/
