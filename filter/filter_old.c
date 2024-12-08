#include "fd.h"

/*****************************************************************************/

void run_filter_cycle(double *psims, double *pot_local, nlc_st *nlc, long *nl, 
  double *ksqr, zomplex *an, double *zn, double *ene_targets, long thread_id, long jns, 
  grid_st *grid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){
  /*******************************************************************
  * This function runs a filter cycle on m ene_targets with one of   *
  * the initial random states as the starting point.                 *
  * inputs:                                                          *
  *  [psims] arr that will hold all ms filtered wavefunctions        *
  *  [pot_local] ngrid-long arr holding the value of the local pot   *
  *  [nlc] nlc struct holding values for computing SO and NL pots    *
  *  [nl] natom-long arr holding the number of NL gridpts per atom   *
  *  [ksqr] ngrid-long arr holding the values of k^2 for KE calc     *
  *  [an] the Newton interpolation coefficients for the filter func  *
  *  [zn] Chebyshev polynomial support points                        *
  *  [ene_targets] target energies where filter funcs are centered   *
  *  [thread_id] filtering of each random state occurs on one thread *
  *  [jns] index of filter cycle                                     *
  *  [grid] struct holding the grid and grid parameters (dx, xmin...)*
  *  [ist] ptr to counters, indices, and lengths                     *
  *  [par] ptr to par_st holding VBmin, VBmax... params              *
  *  [flag] ptr to flag_st holding job flags                         *
  *  [parallel] holds options for parallelization                    *
  * outputs: void                                                    *
  ********************************************************************/

  FILE *pf; char str[100]; 
  long flags = 0, jms, jgrid, jgrid_real, jgrid_imag;
  zomplex *psi, *phi;  
  fftw_plan_loc planfw, planbw; fftw_complex *fftwpsi;
  double *ene_filters;

  fftwpsi = fftw_malloc(sizeof (fftw_complex )*ist->ngrid);
  if ((psi = (zomplex*)calloc(ist->nspinngrid,sizeof(zomplex)))==NULL)nerror("psi");
  if ((phi = (zomplex*)calloc(ist->nspinngrid,sizeof(zomplex)))==NULL)nerror("phi");
  if ((ene_filters = (double*)calloc(ist->m_states_per_filter,sizeof(double)))==NULL)nerror("ene_filters");
  
  planfw = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_FORWARD,flags);
  planbw = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_BACKWARD,flags);
  
  for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
    jgrid_real = ist->complex_idx * jgrid;
    jgrid_imag = ist->complex_idx * jgrid + 1;

    psi[jgrid].re = psims[jgrid_real];

    if (1 == flag->isComplex){
      psi[jgrid].im = psims[jgrid_imag];
    }
  }
  
  /***********************************************************************/
  /*** filter the states and normalize them ***/
  filter(psi,phi,psims,pot_local,nlc,nl,ksqr,an,zn,thread_id,jns,planfw,planbw,fftwpsi,ist,par,flag,parallel); 
  normalize_all(psims,ist, par, flag, parallel);
  
  /*** calculate and print the energy of the filtered states ***/
  energy_all(psi,phi,psims,pot_local,nlc,nl,ksqr,ene_filters,ist,par,flag,planfw,planbw,fftwpsi);

  sprintf (str,"ene-filt-%ld-%ld.dat",thread_id,jns);
  pf = fopen(str , "w");
  for (jms = 0; jms < ist->m_states_per_filter; jms++){
    fprintf (pf,"%ld %.16g %.16g\n", jms, ene_filters[jms], ene_targets[jms]);
  }
  fclose(pf);

  if (flag->printPsiFilt == 1){
    /*** write the filtered states to a binary file ***/
    sprintf (str,"psi-filt-%ld-%ld.dat",thread_id,jns);
    pf = fopen(str , "w");
    fwrite (psims,sizeof(psims[0]),ist->complex_idx*ist->m_states_per_filter*ist->nspinngrid,pf);
    fclose(pf);

    /*** Human-readable output of psi ***/
    // for (jms = 0; jms < ist->m_states_per_filter; jms++){
    //   sprintf(str, "psi-filt-%ld-%ld.dat", jns, jms);
    //   pf = fopen(str, "w");
    //   for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
    //     fprintf(pf, "%ld %lg %lg\n", jgrid, psims[jms*ist->nspinngrid + jgrid].re, psims[jms*ist->nspinngrid + jgrid].re);
    //   }
    //   fclose(pf);
    // }

    /*** print cube files for the filtered states ***/
    double *rho;
    if ((rho = (double *) calloc(ist->ngrid, sizeof(rho[0]))) == NULL){
      fprintf(stderr, "\nOUT OF MEMORY: filter rho\n\n"); exit(EXIT_FAILURE);
    }
    
    for (jms = 0; jms < ist->m_states_per_filter; jms++){
      for (jgrid = 0; jgrid < ist->ngrid; jgrid++){
        jgrid_real = ist->complex_idx * jgrid;
        jgrid_imag = ist->complex_idx * jgrid + 1;

        rho[jgrid] = sqr(psims[ist->complex_idx*jms*ist->ngrid + jgrid_real]);
        if (1 == flag->isComplex){
          rho[jgrid] += sqr(psims[ist->complex_idx*jms*ist->ngrid + jgrid_imag]);
        }
      }
      sprintf(str, "psi-filt-%ld-%ld.cube", jns, jms);
      write_cube_file(rho, grid, str);
    }
  }

  free(psi); free(ene_filters); free(phi);
  fftw_destroy_plan(planfw);
  fftw_destroy_plan(planbw);
  fftw_free(fftwpsi);
  return;
}

/*****************************************************************************/

void filter(zomplex *psin, zomplex *psim1, double *psims, double *pot_local, nlc_st *nlc, long *nl, double *ksqr, \
  zomplex *an, double *zn, long thread_id, long jns, fftw_plan_loc planfw, fftw_plan_loc planbw,\
  fftw_complex *fftwpsi, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){
  /*******************************************************************
  * This function performs the filter algorithm on m states,         *
  * iteratively applying the Hamiltonian to generate the Newton      *
  * interpolation polynomial approximating a filter func centered    *
  * at ene_target                                                    *
  * inputs:                                                          *
  *  [psin] arr that will hold all ms filtered wavefunctions         *
  *  [psim1] ngrid-long arr holding the value of the local pot       *
  *  [psims] ngrid-long arr holding the value of the local pot       *
  *  [pot_local] ngrid-long arr holding the value of the local pot   *
  *  [nlc] nlc struct holding values for computing SO and NL pots    *
  *  [nl] natom-long arr holding the number of NL gridpts per atom   *
  *  [ksqr] ngrid-long arr holding the values of k^2 for KE calc     *
  *  [an] the Newton interpolation coefficients for the filter func  *
  *  [zn] Chebyshev polynomial support points                        *
  *  [ene_targets] target energies where filter funcs are centered   *
  *  [thread_id] filtering of each random state occurs on one thread *
  *  [jns] index of filter cycle                                     *
  *  [planfw] FFTW3 plan for executing 3D forward DFT                *
  *  [planfw] FFTW3 plan for executing 3D backwards DFT              *
  *  [fftwpsi] location to store outcome of Fourier transform        *
  *  [ist] ptr to counters, indices, and lengths                     *
  *  [par] ptr to par_st holding VBmin, VBmax... params              *
  *  [flag] ptr to flag_st holding job flags                         *
  *  [parallel] holds options for parallelization                    *
  * outputs: void                                                    *
  ********************************************************************/

  FILE *pf; char str[100];
  long jgrid, jgrid_real, jgrid_imag, jc, jms, ncjms, jmsg;

  for (jms = 0; jms < ist->m_states_per_filter; jms++){
    //printf("\njms = %ld\n", jms); fflush(0);
    ncjms = ist->ncheby*jms; 
    jmsg = ist->complex_idx * jms * ist->nspinngrid;

    for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
      jgrid_real = ist->complex_idx * jgrid;
      jgrid_imag = ist->complex_idx * jgrid + 1;

      psims[jmsg+jgrid_real] = an[ncjms+0].re * psin[jgrid].re - an[ncjms+0].im * psin[jgrid].im;
      if (1 == flag->isComplex){
        psims[jmsg+jgrid_imag] = an[ncjms+0].re * psin[jgrid].im + an[ncjms+0].im * psin[jgrid].re;
      }
    }
  }

  sprintf (str,"prop%ld.dat",thread_id);
  pf = fopen(str , "w");
  
  for (jc = 1; jc < ist->ncheby; jc++){
    memcpy(&psim1[0],&psin[0],ist->nspinngrid*sizeof(psim1[0]));
    long i;
    
    hamiltonian(psin,psim1,pot_local,nlc,nl,ksqr,ist,par,flag,planfw,planbw,fftwpsi);

    for (i = 0; i < ist->nspinngrid; i++){
      /*** par->dE_1 = 4.0 / par->dE and therefore I don't multiply by 4 ***/
      psin[i].re = par->dE_1 * psin[i].re - (2.0 + zm1 + par->Vmin * par->dE_1) * psim1[i].re;
      psin[i].im = par->dE_1 * psin[i].im - (2.0 + zm1 + par->Vmin * par->dE_1) * psim1[i].im;
    }
    for (jms = 0; jms < ist->m_states_per_filter; jms++){
      ncjms = ist->ncheby*jms; 
      jmsg = ist->complex_idx * jms * ist->nspinngrid;

      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
        jgrid_real = ist->complex_idx * jgrid;
        jgrid_imag = ist->complex_idx * jgrid + 1;

	      psims[jmsg+jgrid_real] += (an[ncjms+jc].re * psin[jgrid].re - an[ncjms+jc].im * psin[jgrid].im);
	      if (1 == flag->isComplex){
          psims[jmsg+jgrid_imag] += (an[ncjms+jc].re * psin[jgrid].im + an[ncjms+jc].im * psin[jgrid].re);
        }
      }
    }
    if (1 == flag->printNorm) normalize_all(psims, par->dv,ist->m_states_per_filter,ist->nspinngrid,parallel->nthreads, ist->complex_idx, flag->printNorm);
    
    if (!(jc % 100)) {
      fprintf (pf,"%ld %ld\n",jc,jns); fflush(pf);
      }
  }
  //fclose(pf);
  
  return;
}

/*****************************************************************************/

void scale_eigs_for_cheby(zomplex *psim1, zomplex *psin, double *pot_local, nlc_st *nlc, long *nl, 
  double *ksqr, double zm1, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi, 
  index_st *ist, par_st *par, flag_st *flag){
  /*******************************************************************
  * This function scales the spectrum of the filter Hamiltonian      *
  * to be between [-2,2], the convergence radius of the Cheby series *
  * inputs:                                                          *
  *  [psin] arr that will hold all ms filtered wavefunctions         *
  *  [psim1] ngrid-long arr holding the value of the local pot       *
  *  [pot_local] ngrid-long arr holding the value of the local pot   *
  *  [nlc] nlc struct holding values for computing SO and NL pots    *
  *  [nl] natom-long arr holding the number of NL gridpts per atom   *
  *  [ksqr] ngrid-long arr holding the values of k^2 for KE calc     *
  *  [zm1] Chebyshev polynomial support points                       *
  *  [planfw] FFTW3 plan for executing 3D forward DFT                *
  *  [planfw] FFTW3 plan for executing 3D backwards DFT              *
  *  [fftwpsi] location to store outcome of Fourier transform        *
  *  [ist] ptr to counters, indices, and lengths                     *
  *  [par] ptr to par_st holding VBmin, VBmax... params              *
  *  [flag] ptr to flag_st holding job flags                         *
  * outputs: void                                                    *
  ********************************************************************/

  long i;
  
  hamiltonian(psin,psim1,pot_local,nlc,nl,ksqr,ist,par,flag,planfw,planbw,fftwpsi);

  for (i = 0; i < ist->nspinngrid; i++){
    /*** par->dE_1 = 4.0 / par->dE and therefore I don't multiply by 4 ***/
    psin[i].re = par->dE_1 * psin[i].re - (2.0 + zm1 + par->Vmin * par->dE_1) * psim1[i].re;
    psin[i].im = par->dE_1 * psin[i].im - (2.0 + zm1 + par->Vmin * par->dE_1) * psim1[i].im;
  }
  
  return;
}

/*****************************************************************************/
