#include "fd.h"

/*****************************************************************************/

double energy(zomplex *psi, zomplex *phi, double *pot_local, nlc_st *nlc, long *nl, double *ksqr, index_st *ist,
  par_st *par, flag_st *flag, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi){
  /*******************************************************************
  * This function calculates Exp[E] of a filtered state by evaluating*
  * |phi> = H|psi> ~ E|psi>, then projecting with <psi|.             *
  * E = <psi|H|psi>                                                  *
  * inputs:                                                          *
  *  [psi] ngrid-long arr of double/zomplex holding orig. wavefnc    *
  *  [phi] ngrid-long arr to hold |phi> = H|psi>                     *
  *  [pot_local] ngrid-long arr holding the value of the local pot   *
  *  [nlc] nlc struct holding values for computing SO and NL pots    *
  *  [nl] natom-long arr holding the number of NL gridpts per atom   *
  *  [ksqr] ngrid-long arr holding the values of k^2 for KE calc     *
  *  [ist] ptr to counters, indices, and lengths                     *
  *  [par] ptr to par_st holding VBmin, VBmax... params              *
  *  [flag] ptr to flag_st holding job flags                         *
  *  [planfw] FFTW3 plan for executing 3D forward DFT                *
  *  [planfw] FFTW3 plan for executing 3D backwards DFT              *
  *  [fftwpsi] location to store outcome of Fourier transform        *
  * outputs: [double]  E = <psi|H|psi>                               *
  ********************************************************************/

  long i; 
  double ene = 0.0;

  memcpy(&phi[0], &psi[0], ist->nspinngrid*sizeof(phi[0]));
  hamiltonian(phi, psi, pot_local, nlc, nl, ksqr, ist, par, flag, planfw, planbw, fftwpsi);

  for (i = 0; i < ist->nspinngrid; i++) {
    ene += (psi[i].re * phi[i].re + psi[i].im * phi[i].im);
  }
  ene *= par->dv;

  return (ene);
}

/***************************************************************************/

void energy_all(double *psitot, long n_states, double *pot_local, nlc_st *nlc, long *nl,double *ksqr,
  double *ene_filters, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){
  /*******************************************************************
  * This function calculates Exp[E] of all filtered states           *
  * inputs:                                                          *
  *  [psi] ngrid-long arr of double/zomplex to hold orig. wavefnc    *
  *  [phi] ngrid-long arr to hold |phi> = H|psi>                     *
  *  [psims] ms*ngrid-long arr holding all ms wavefuncs              *
  *  [pot_local] ngrid-long arr holding the value of the local pot   *
  *  [nlc] nlc struct holding values for computing SO and NL pots    *
  *  [nl] natom-long arr holding the number of NL gridpts per atom   *
  *  [ksqr] ngrid-long arr holding the values of k^2 for KE calc     *
  *  [ene_filters] ms-long array to store energies of filter states  *
  *  [ist] ptr to counters, indices, and lengths                     *
  *  [par] ptr to par_st holding VBmin, VBmax... params              *
  *  [flag] ptr to flag_st holding job flags                         *
  *  [planfw] FFTW3 plan for executing 3D forward DFT                *
  *  [planfw] FFTW3 plan for executing 3D backwards DFT              *
  *  [fftwpsi] location to store outcome of Fourier transform        *
  * outputs: void                                                    *
  ********************************************************************/
  long jmn;
  
// #pragma omp parallel for private(jmn)
  for (jmn = 0; jmn < n_states; jmn++) {
    
    // Indexes for arrays
    long jgrid, jgrid_real, jgrid_imag, j_state;
    // Arrays for hamiltonian evaluation
    zomplex *psi, *phi;  
    // FFT
    long fft_flags = 0;
    fftw_plan_loc planfw, planbw;
    fftw_complex *fftwpsi;

    if ((psi = (zomplex*)calloc(ist->nspinngrid,sizeof(zomplex)))==NULL){ 
      fprintf(stderr, "\nOUT OF MEMORY: psi in run_filter_cycle\n\n"); exit(EXIT_FAILURE);
    }
    if ((phi = (zomplex*)calloc(ist->nspinngrid,sizeof(zomplex)))==NULL){ 
      fprintf(stderr, "\nOUT OF MEMORY: phi in run_filter_cycle\n\n"); exit(EXIT_FAILURE);
    }
    // Create FFT structs and plans for Fourier transform
    fftwpsi = fftw_malloc(sizeof (fftw_complex) * ist->ngrid);
    planfw = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_FORWARD,fft_flags);
    planbw = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_BACKWARD,fft_flags);
  
    j_state = jmn * ist->complex_idx * ist->nspinngrid;
    
    // copy the wavefunction for state jmn into psi
    if (1 == flag->isComplex){
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
        jgrid_real = ist->complex_idx * jgrid;
        jgrid_imag = ist->complex_idx * jgrid + 1;

        psi[jgrid].re = psitot[j_state + jgrid_real];
        psi[jgrid].im = psitot[j_state + jgrid_imag];
      }
    }
    else{
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
        psi[jgrid].re = psitot[j_state + jgrid];
        psi[jgrid].im = 0.0;
      }
    }
    
    memcpy(&phi[0], &psi[0], ist->nspinngrid*sizeof(phi[0]));
    hamiltonian(phi, psi, pot_local, nlc, nl, ksqr, ist, par, flag, planfw, planbw, fftwpsi);
    
    // Compute the expectation value <psi|H|psi>
    // The quantity H|psi> is stored in phi
    
    ene_filters[jmn] = 0.0;
    for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
      ene_filters[jmn] += (psi[jgrid].re * phi[jgrid].re + psi[jgrid].im * phi[jgrid].im);
    }
    ene_filters[jmn] *= par->dv;
    
    // Free dynamically allocated memory
    free(psi); free(phi);
    fftw_destroy_plan(planfw);
    fftw_destroy_plan(planbw);
    fftw_free(fftwpsi);
  }

  return;
}

/***************************************************************************/

void get_energy_range(zomplex *psi,zomplex *phi,double *pot_local, grid_st *grid, nlc_st *nlc, long *nl, double *ksqr,\
  index_st *ist, par_st *par, parallel_st *parallel, flag_st *flag, fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi){
  /*******************************************************************
  * This function calculates the range of the Hamiltonian spectrum   *
  * by imaginary time evolution. Min. obtained by propagating        *
  * exp[-H]|psi>, max obtained by propagating exp[H]|psi>            *
  * inputs:                                                          *
  *  [psi] ngrid-long arr of double/zomplex to hold orig. wavefnc    *
  *  [phi] ngrid-long arr to hold |phi> = H|psi>                     *
  *  [psims] ms*ngrid-long arr holding all ms wavefuncs              *
  *  [pot_local] ngrid-long arr holding the value of the local pot   *
  *  [nlc] nlc struct holding values for computing SO and NL pots    *
  *  [nl] natom-long arr holding the number of NL gridpts per atom   *
  *  [ksqr] ngrid-long arr holding the values of k^2 for KE calc     *
  *  [ene_filters] ms-long array to store energies of filter states  *
  *  [ist] ptr to counters, indices, and lengths                     *
  *  [par] ptr to par_st holding VBmin, VBmax... params              *
  *  [flag] ptr to flag_st holding job flags                         *
  *  [planfw] FFTW3 plan for executing 3D forward DFT                *
  *  [planfw] FFTW3 plan for executing 3D backwards DFT              *
  *  [fftwpsi] location to store outcome of Fourier transform        *
  * outputs: void                                                    *
  ********************************************************************/
  
  FILE *pf;
  long i, ispn, jgrid; 
  long rand_seed = -874917403;
  double ene_old; 
  double norma, Emin, Emax, tau = 0.05; //tau = 0.025
  long max_iter = 500;
  
  
  if (0 == flag->approxEnergyRange) {
    if (parallel->mpi_rank == 0) printf("Iteratively determining range of Hamiltonian\n");
    // Find E_min
    pf = fopen("Emin-init.dat" , "w");
    // Initialize random state to begin propagation
    for (ispn = 0; ispn < ist->nspin; ispn++)  {
      init_psi(&phi[ispn*ist->ngrid], &rand_seed, grid, ist, par, flag, parallel);
    }
    
    Emin = (ene_old = 0.0) + 10.0; // 
    for (i = 0; (fabs((Emin - ene_old) / Emin) > 1.0e-6) && (i < max_iter) ; i++){
      // Apply the Hamiltonian, shift orig by |phi>, and normalize (equivalent to forward imag. time propagation step)
      memcpy(&psi[0], &phi[0], ist->nspinngrid*sizeof(phi[0]));
      
      hamiltonian(phi, psi, pot_local, nlc, nl,ksqr, ist, par, flag, planfw, planbw, fftwpsi);
      
      for (ispn = 0; ispn < ist->nspin; ispn++) {
        for (jgrid = 0; jgrid < ist->ngrid; jgrid++) {
          phi[ispn*ist->ngrid+jgrid].re = psi[ispn*ist->ngrid+jgrid].re - tau*phi[ispn*ist->ngrid+jgrid].re;
          phi[ispn*ist->ngrid+jgrid].im = psi[ispn*ist->ngrid+jgrid].im - tau*phi[ispn*ist->ngrid+jgrid].im;
        }
      }

      norma = normalize(&phi[0], ist->nspinngrid, ist, par, flag, parallel);
      // set Emin for next iteration
      ene_old = Emin;
      Emin = energy(phi, psi, pot_local, nlc, nl, ksqr, ist, par, flag, planfw, planbw, fftwpsi);
      // print progress
      fprintf(pf, "%ld %.16g %.16g %.16g\n", i, ene_old, Emin, fabs((Emin - ene_old) / Emin)); fflush(pf);

      if ((i > 5) && (Emin - ene_old) > 0){
        if (parallel->mpi_rank == 0) printf("\nWarning: positive step in energy minimization. Check Emin-init.dat\n");
        fprintf(pf, "Warning: positive step in energy minimization\n");
        tau -= 0.025;
      }
    }
    fclose(pf);

    // Find Emax 
    pf = fopen("Emax-init.dat" , "w");
    // Initialize random state to begin propagation
    for (ispn = 0; ispn < ist->nspin; ispn++){
      init_psi(&psi[ispn*ist->ngrid], &rand_seed, grid, ist, par, flag, parallel);
    }
    
    Emax = (ene_old = 0.0) + 0.1;
    for (i = 0; (fabs((Emax-ene_old)/Emax)>1.0e-6) & (i < max_iter); i++){
      // Apply the Hamiltonian and normalize (equivalent to propagation step)
      memcpy(&phi[0], &psi[0], ist->nspinngrid*sizeof(phi[0]));
      hamiltonian(psi, phi, pot_local,nlc, nl, ksqr, ist, par, flag, planfw, planbw, fftwpsi);
      
      norma = normalize(&psi[0], ist->nspinngrid, ist, par, flag, parallel);
      // reset the max energy after the last iteration
      ene_old = Emax;
      Emax = energy(psi,phi,pot_local,nlc,nl,ksqr,ist,par,flag,planfw,planbw,fftwpsi);
      // print progress
      fprintf (pf,"%ld %.16g %.16g %.16g\n", i, ene_old, Emax, fabs((Emax-ene_old)/Emax)); fflush(pf);
    }
    fclose(pf);
  } else if (1 == flag->approxEnergyRange){
    if (parallel->mpi_rank == 0) printf("Approximating energy range of Hamiltonian as [Vmin, Vmax + KE_max]\n");
    Emin = par->Vmin + 0.5;
    Emax = par->Vmax + par->KE_max;
    if (1 == flag->NL){
      Emax += 3.0;
    }
  } else {fprintf(stderr, "ERROR: invalid Hamiltonian energy range strategy selected\n"); exit(EXIT_FAILURE);}


  // Expand the delta E range artificially to help algorithm convergence (less efficient, but more robust).
  Emax *= 1.2;
  Emin -= 0.2 * fabs(Emin);
  
  // Set parameters for calc'ing Filter coefficients
  par->Vmin = Emin;
  par->dE = (Emax - Emin);
  par->dE_1 = 4.0 / par->dE;

  if (parallel->mpi_rank == 0) printf("Emin = %lg, Emax = %lg, dE = %lg\n", Emin, Emax, par->dE);
  fflush(stdout);

  return;
}


/****************************************************************************************/

void calc_sigma_E(double *psitot, double *pot_local, nlc_st *nlc, long *nl, double *ksqr,
  double *sigma_E,index_st *ist,par_st *par,flag_st *flag){
  /*******************************************************************
  * This function calculates the quality of the eigenstates by       *
  * evaluating sigma_E^2 = <psi|H^2|psi> - <psi|H|psi>^2             *
  * inputs:                                                          *
  *  [psi] ngrid-long arr of double/zomplex to hold orig. wavefnc    *
  *  [phi] ngrid-long arr to hold |phi> = H|psi>                     *
  *  [psitot] m*n*ngrid-long arr holding all wavefuncs               *
  *  [pot_local] ngrid-long arr holding the value of the local pot   *
  *  [nlc] nlc struct holding values for computing SO and NL pots    *
  *  [nl] natom-long arr holding the number of NL gridpts per atom   *
  *  [ksqr] ngrid-long arr holding the values of k^2 for KE calc     *
  *  [sigma_E] m*n-long array to store energies of filtered states   *
  *  [ist] ptr to counters, indices, and lengths                     *
  *  [par] ptr to par_st holding VBmin, VBmax... params              *
  *  [flag] ptr to flag_st holding job flags                         *
  *  [planfw] FFTW3 plan for executing 3D forward DFT                *
  *  [planfw] FFTW3 plan for executing 3D backwards DFT              *
  *  [fftwpsi] location to store outcome of Fourier transform        *
  * outputs: void                                                    *
  ********************************************************************/

  long ims;
  
  // Loop over all M*N states
  #pragma omp parallel for private(ims)
  for (ims = 0; ims < ist->mn_states_tot; ims++) {
    long jgrid, jgrid_real, jgrid_imag, jstate;
    double eval, eval2;
    int fft_flags = 0;
    fftw_plan_loc planfw, planbw; 
    fftw_complex *fftwpsi;
    // Arrays for hamiltonian evaluation
    zomplex *psi, *phi; 

    jstate = ist->complex_idx*ims*ist->nspinngrid;
    
    if ((psi = (zomplex*)calloc(ist->nspinngrid,sizeof(zomplex)))==NULL){ 
      fprintf(stderr, "\nOUT OF MEMORY: psi in calc_sigma_E\n\n"); exit(EXIT_FAILURE);
    }
    if ((phi = (zomplex*)calloc(ist->nspinngrid,sizeof(zomplex)))==NULL){ 
      fprintf(stderr, "\nOUT OF MEMORY: phi in calc_sigma_E\n\n"); exit(EXIT_FAILURE);
    }
    fftwpsi = fftw_malloc(sizeof (fftw_complex )*ist->ngrid);
    planfw = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_FORWARD,fft_flags);
    planbw = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_BACKWARD,fft_flags);
  
    // select the current state to compute sigma_E for
    if (1 == flag->isComplex){
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
        jgrid_real = ist->complex_idx * jgrid;
        jgrid_imag = ist->complex_idx * jgrid + 1;

        psi[jgrid].re = psitot[jstate+jgrid_real];
        
        psi[jgrid].im = psitot[jstate+jgrid_imag];
      } 
    }
    else {
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
      psi[jgrid].re = psitot[jstate + jgrid];
      
      psi[jgrid].im = 0.0;
      } 
    }


    memcpy(&phi[0],&psi[0],ist->nspinngrid*sizeof(phi[0]));
    // Apply the Hamiltonian to |psi>: |phi> = H|psi>
    hamiltonian(phi, psi, pot_local, nlc, nl, ksqr, ist, par, flag, planfw, planbw, fftwpsi);
    // Calculate the expectation value of H for wavefunc psi: <psi|H|psi> = <psi|phi> = sum_{jgrid} psi[jgrid] * phi[jgrid] * dv
    eval = 0.0;
    if (1 == flag->isComplex){
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
        jgrid_real = ist->complex_idx*jgrid;
        jgrid_imag = ist->complex_idx*jgrid + 1;

        eval += psitot[jstate+jgrid_real] * phi[jgrid].re;
        eval += psitot[jstate+jgrid_imag] * phi[jgrid].im;
      }
    }
    else {
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
        eval += psitot[jstate+jgrid] * phi[jgrid].re;
      }
    }
    eval *= par->dv;
    
    // Apply the Hamiltonian again onto phi: H|phi> = H^2|psi>
    memcpy(&psi[0], &phi[0], ist->nspinngrid*sizeof(psi[0]));
    hamiltonian(phi, psi, pot_local, nlc, nl, ksqr, ist, par,flag, planfw, planbw, fftwpsi);
    // Calculate the expectation value of H^2: <psi|H^2|psi>
    eval2 = 0.0;
    if (1 == flag->isComplex){
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
        jgrid_real = ist->complex_idx*jgrid;
        jgrid_imag = ist->complex_idx*jgrid + 1;

        eval2 += psitot[jstate+jgrid_real] * phi[jgrid].re;
        eval2 += psitot[jstate+jgrid_imag] * phi[jgrid].im;
      }
    }
    else{
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
        eval2 += psitot[jstate+jgrid] * phi[jgrid].re;
      }
    }
    eval2 *= par->dv;
    // var = <psi|H^2|psi> - <psi|H|psi>^2
    eval2 -= sqr(eval);
    // sigma_E is the sqrt of the variance
    sigma_E[ims] = sqrt(fabs(eval2));

    // Free dynamically allocated memory
    free(psi); free(phi);
    fftw_destroy_plan(planfw);
    fftw_destroy_plan(planbw);
    fftw_free(fftwpsi);
  }

  return;
}

/*****************************************************************************/
