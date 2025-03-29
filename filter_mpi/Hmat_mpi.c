#include "Hmat.h"

/****************************************************************************/

void diag_H_mpi(
  double*        psitot, 
  double*        pot_local, 
  zomplex*       LS, 
  nlc_st*        nlc, 
  long*          nl, 
  double*        ksqr, 
  double*        eval, 
  index_st*      ist, 
  par_st*        par, 
  flag_st*       flag, 
  parallel_st*   parallel
  ){
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
  *  [planfw] FFTW3 plan for executing 3D backwards DFT              *
  *  [fftwpsi] location to store outcome of Fourier transform        *
  * outputs: void                                                    *
  ********************************************************************/

  FILE *pg;
  
  long ims, jms, jgrid, jgrid_real, jgrid_imag;
  
  const long long mn_states_tot = (int)ist->mn_states_tot, lwk = 3*ist->mn_states_tot;
  long long info;
  
  double *rwork, sumre, sumim; zomplex *tpsi;
  double *H, *work;

  MKL_Complex16 *H_z, *work_z;

  zomplex *psi, *phi;

  const int   mpir   =  parallel->mpi_rank;
  const int   mpi_sz =  parallel->mpi_size;

  const unsigned long stlen = ist->nspinngrid * ist->complex_idx;

  unsigned long start, end;

  // FFT
  long fft_flags = FFTW_MEASURE;
  fftw_plan_loc planfw, planbw;
  fftw_complex *fftwpsi;
  
  // Create FFT structs and plans for Fourier transform
  fftw_init_threads();
  fftw_plan_with_nthreads(par->ham_threads);
  fftwpsi = fftw_malloc(sizeof (fftw_complex) * ist->ngrid);
  planfw = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_FORWARD,fft_flags);
  planbw = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_BACKWARD,fft_flags);
  

  // Allocate memory specific to scalar or complex arrays
  if (0 == flag->isComplex){
    H = (double *) calloc(ist->mn_states_tot*ist->mn_states_tot, sizeof(double));
    work = (double *) calloc(lwk, sizeof(double));
    tpsi = (zomplex *) calloc(ist->mn_states_tot, sizeof(zomplex));
  } 
  if (1 == flag->isComplex){
    H_z = (MKL_Complex16*)calloc(ist->mn_states_tot*ist->mn_states_tot,sizeof(MKL_Complex16));
    work_z = (MKL_Complex16*)calloc(lwk,sizeof(MKL_Complex16));
    rwork = (double*)calloc(3*ist->mn_states_tot,sizeof(double));
    tpsi = (zomplex*)calloc(ist->mn_states_tot,sizeof(zomplex));
  }

  if ((psi = (zomplex*) calloc(ist->nspinngrid, sizeof(psi[0]))) == NULL){
    fprintf(stderr, "\nOUT OF MEMORY: psi in diag_H\n\n"); exit(EXIT_FAILURE);
  }
  if ((phi = (zomplex*) calloc(ist->nspinngrid, sizeof(phi[0]))) == NULL){
    fprintf(stderr, "\nOUT OF MEMORY: phi in diag_H\n\n"); exit(EXIT_FAILURE);
  }
  
  if (mpir == 0) printf("Constructing Hamiltonian matrix\n"); fflush(stdout);

  omp_set_max_active_levels(1);
  omp_set_num_threads(parallel->nthreads);
  /*** calculate H|psi_i> ***/
  
  printf("Determining the workload for mn_tot = %ld on MPI rank %d\n", ist->mn_states_tot, mpir);
  double stride = sqrt( (double) mpir / (double) mpi_sz);
  printf("MPI rank %d sqrt odd/size = %lg\n", mpir, stride);
  
  start = (unsigned long) ((double) ist->mn_states_tot * stride);
  printf("MPI rank %d start = %lu\n", mpir, start);
  end   = (unsigned long) ( (double) ist->mn_states_tot * sqrt((double) (mpir + 1) /(double) mpi_sz ) );
  
  // handle remainder so that all states are assigned to a rank
  if ( mpir == (mpi_sz - 1) ){
      end = ist->mn_states_tot;
  }
  printf("MPI rank %d end = %lu\n", mpir, end);

  for (ims = start; ims < end; ims++){
    
    if (1 == flag->isComplex){
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
        jgrid_real = ist->complex_idx * jgrid;
        jgrid_imag = jgrid_real + 1;

        psi[jgrid].re = psitot[ims * stlen + jgrid_real];
        psi[jgrid].im = psitot[ims * stlen + jgrid_imag];
      }
    }
    else if (0 == flag->isComplex){
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
        psi[jgrid].re = psitot[ims * stlen + jgrid];
        psi[jgrid].im = 0.0;
      }
    }

    memcpy(&phi[0],&psi[0],ist->nspinngrid*sizeof(phi[0]));
    p_hamiltonian(phi, psi, pot_local, LS, nlc, nl, ksqr, ist, par, flag, planfw, planbw, fftwpsi, par->ham_threads);

    /*** calculate <psi_j|H|psi_i> ***/
    #pragma omp parallel for private(jms, jgrid, jgrid_real, jgrid_imag) shared(H, H_z)
    for (jms = 0; jms <= ims; jms++){
      // Utilize Hermitian property of Hamiltonian matrix elements to reduce computation to only lower triangle
      if (0 == flag->isComplex){
        H[ims*ist->mn_states_tot + jms] = H[jms*ist->mn_states_tot + ims] = dotpreal(phi,psitot,jms,ist->nspinngrid,par->dv);
      } 
      if (1 == flag->isComplex){
        H_z[ims*ist->mn_states_tot + jms] = H_z[jms*ist->mn_states_tot + ims] = dotp(phi,psitot,jms,ist->nspinngrid,par->dv);
        H_z[jms*ist->mn_states_tot + ims].imag *= -1;
      }
    }
    
    // Print out progress
    if (mpir == 0){
      if ( (ims == 0) || (0 == (ims % (ist->mn_states_tot/4 + 1))) || (ims == (ist->mn_states_tot - 1))){
        print_progress_bar(ims, ist->mn_states_tot);
      }
    }
    
  }



  // free dynamically allocated memory for psi and phi
  free(psi); free(phi);

  // print out Hmat for debugging purposes
  if (mpir == 0) {
    pg = fopen("hmat.dat" , "w");
    for (ims = 0; ims < ist->mn_states_tot; ims++){
      for (jms = 0; jms < ist->mn_states_tot; jms++){
        if (0 == flag->isComplex){
          fprintf(pg, "%lg ", H[ims*ist->mn_states_tot + jms]);
        }
        if (1 == flag->isComplex){
          fprintf(pg, "%lg+i%lg ", H_z[ims*ist->mn_states_tot + jms].real, H_z[ims*ist->mn_states_tot + jms].imag);
        }
      }
      fprintf(pg, "\n");
    }
    fclose(pg);
  }
  /*** diagonalize the Hamiltonian H ***/
  if (mpir == 0) printf("Diagonalizing Hamiltonian | %s\n", get_time());
  
  // Use real, symmetric diagonalization routine for real wavefunctions
  if (0 == flag->isComplex){
    dsyev_("V", "U", &mn_states_tot, &H[0], &mn_states_tot, &eval[0], &work[0], &lwk, &info);
    if (info) { 
      fprintf(stderr, "error in dsyev_ H\n"); exit(EXIT_FAILURE);
    }
  }
  
  // Use Hermitian diagonalization routine for complex wavefunctions
  if (1 == flag->isComplex){
    zheev_("V","U", &mn_states_tot, &H_z[0], &mn_states_tot, &eval[0], &work_z[0], &lwk, &(rwork[0]), &info);
    if (info){
      fprintf(stderr, "error in zheev_ H\n"); exit(EXIT_FAILURE);
    }
  }

  /*** copy the new function into psitot ***/
  if (mpir == 0) printf("Diagonalization complete! | %s\n", get_time()); fflush(stdout);
  
  // The eigenvectors have been computed in the basis of orthogonalized
  // filtered functions (Phi_filter). In order to obtain them in the grid basis
  // as Psi_grid, we must perform a change of basis.
  // The matrix H is holding the eigenvectors, so we perform
  // (Psi_grid)_a = SUM_i H_ai * (Phi_filter)_i
  // for each grid point
  
  if (mpir == 0) printf("Writing out eigenvectors in the grid basis\n"); 
  
  for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
    jgrid_real = ist->complex_idx * jgrid;
    jgrid_imag = ist->complex_idx * jgrid + 1;
    
    #pragma omp parallel for private (jms)
    for (jms = 0; jms < ist->mn_states_tot; jms++){
      tpsi[jms].re = psitot[ist->complex_idx*jms*ist->nspinngrid + jgrid_real];
      psitot[ist->complex_idx*jms*ist->nspinngrid + jgrid_real] = 0.0;

      if (1 == flag->isComplex){
        tpsi[jms].im = psitot[ist->complex_idx*jms*ist->nspinngrid + jgrid_imag];
        psitot[ist->complex_idx*jms*ist->nspinngrid + jgrid_imag] = 0.0;
      }
      
    }

    #pragma omp parallel for private (jms,sumre,sumim)
    for (jms = 0; jms < ist->mn_states_tot; jms++){

      sumre = sumim = 0.0;    
      for (ims = 0; ims < ist->mn_states_tot; ims++){
        if (0 == flag->isComplex){
          sumre += H[jms*ist->mn_states_tot+ims] * tpsi[ims].re;
        }
        if (1 == flag->isComplex){
          sumre += (H_z[jms*ist->mn_states_tot+ims].real * tpsi[ims].re + H_z[jms*ist->mn_states_tot+ims].imag * tpsi[ims].im);
          sumim += (H_z[jms*ist->mn_states_tot+ims].real * tpsi[ims].im - H_z[jms*ist->mn_states_tot+ims].imag * tpsi[ims].re);
        }
      }

      //#pragma omp critical
      psitot[ist->complex_idx*jms*ist->nspinngrid + jgrid_real] = sumre;
      if (1 == flag->isComplex){
        psitot[ist->complex_idx*jms*ist->nspinngrid+jgrid_imag] = sumim;
      }
    }

    if (!(jgrid % (ist->nspinngrid / 4) ) && (mpir == 0)) {
      printf("\tFinished grid point no. %ld | %s\n",jgrid, get_time()); fflush(stdout);
    }
  }
  

  free(tpsi); 
  if (0 == flag->isComplex){
    free(work); free(H); 
  }
  if (1 == flag->isComplex){
    free(H_z); free(work_z); free(rwork);
  }

  return;
}
  