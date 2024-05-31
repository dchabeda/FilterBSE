#include "fd.h"

/*****************************************************************************/

void diag_H(zomplex *psi,zomplex *phi,double *psitot,double *pot_local,nlc_st *nlc,long *nl,double *ksqr,double *eval,index_st *ist,par_st *par,flag_st *flag,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi){
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
  long long mn_states_tot = (long long)ist->mn_states_tot, info, lwk = 3*ist->mn_states_tot;
  double *rwork, sumre, sumim; zomplex *tpsi;
  double *H, *work;
  MKL_Complex16 *H_z, *work_z;

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

  
  /*omp_set_dynamic(0);
    omp_set_num_threads(ist->nthreads);*/

  /*** calculate H|psi_i> ***/
  pg = fopen("hmat.dat" , "w");
  fprintf (pg,"calculate the H matrix\n"); fflush(pg);
  for (ims = 0; ims < ist->mn_states_tot; ims++){
    //#pragma omp parallel for private(jgrid)
    for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++) {
      jgrid_real = ist->complex_idx * jgrid;
      jgrid_imag = ist->complex_idx * jgrid + 1;

      psi[jgrid].re = psitot[ist->complex_idx*ims*ist->nspinngrid + jgrid_real];
      
      if (1 == flag->isComplex){
        psi[jgrid].im = psitot[ist->complex_idx*ims*ist->nspinngrid+jgrid_imag];
      } else if (0 == flag->isComplex){
        psi[jgrid].im = 0.0;
      }
    }
    memcpy(&phi[0],&psi[0],ist->nspinngrid*sizeof(phi[0]));
    hamiltonian(phi,psi,pot_local,nlc,nl,ksqr,ist,par,flag,planfw,planbw,fftwpsi);

    /*** calculate <psi_j|H|psi_i> ***/
    //#pragma omp parallel for private(jms)
    for (jms = 0; jms < ist->mn_states_tot; jms++){
      if (0 == flag->isComplex){
        H[ims*ist->mn_states_tot + jms] = dotpreal(phi,psitot,jms,ist->nspinngrid,par->dv);
      } 
      if (1 == flag->isComplex){
        H_z[ims*ist->mn_states_tot + jms] = dotp(phi,psitot,jms,ist->nspinngrid,par->dv);
      }
      //fprintf (pg,"%ld %ld %g %g\n",ims,jms,H[ims*ist->mn_states_tot+jms].real,H[ims*ist->mn_states_tot+jms].imag);
    }
    fprintf (pg,"finshed row %ld\n",ims); fflush(pg);
  }

  /*** diagonalize the Hamiltonian H ***/
  fprintf (pg,"diagonalize the new Hamiltonian H\n"); fflush(pg);
  // Use real, symmetric diagonalization routine for real wavefunctions
  if (0 == flag->isComplex){
    dsyev_("V", "U", &mn_states_tot, &H[0], &mn_states_tot, &eval[0], &work[0], &lwk, &info);
    if (info) { 
      fprintf(stderr, "error in dsyev_ H"); exit(EXIT_FAILURE);
    }
  }
  // Use Hermitian diagonalization routine for complex wavefunctions
  if (1 == flag->isComplex){
    zheev_("V","U", &mn_states_tot, &H_z[0], &mn_states_tot, &eval[0], &work_z[0], &lwk, &(rwork[0]), &info);
    if (info){
      fprintf(stderr, "error in zheev_ H"); exit(EXIT_FAILURE);
    }
  }

  /*** copy the new function into psitot ***/
  fprintf (pg, "copy the new function into psitot\n"); fflush(pg);
  // The eigenvectors have been computed in the basis of orthogonalized
  // filtered functions (Phi_filter). In order to obtain them in the grid basis
  // as Psi_grid, we must perform a change of basis.
  // The matrix H is holding the eigenvectors, so we perform
  // (Psi_grid)_a = SUM_i H_ai * (Phi_filter)_i
  // for each grid point

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
    if (!(jgrid % 1000)) fprintf (pg,"finished grid point %ld\n",jgrid); fflush(pg);
  }
  fprintf (pg, "free memory\n");
  fclose(pg);

  free(tpsi); 
  if (0 == flag->isComplex){
    free(work); free(H); 
  }
  if (1 == flag->isComplex){
    free(H_z); free(work_z); free(rwork);
  }

  return;
}

/****************************************************************************/

double dotpreal(zomplex *psi, double *phi, long m, long ngrid, double dv){
  
  long jgrid;
  double sum;
  
  sum = 0.0;
  for (jgrid = 0; jgrid < ngrid; jgrid++){
    sum += (psi[jgrid].re * phi[m*ngrid+jgrid]);
  }
  sum *= dv;

  return(sum);
}

/****************************************************************************/

MKL_Complex16 dotp(zomplex *psi, double *psitot,long m,long ngrid,double dv){
  
  long jgrid, jgrid_real, jgrid_imag;
  MKL_Complex16 sum;
  
  sum.real = sum.imag = 0.0;
  for (jgrid = 0; jgrid < ngrid; jgrid++){
    jgrid_real = 2*jgrid;
    jgrid_imag = 2*jgrid + 1;

    sum.real += (psi[jgrid].re * psitot[2*m*ngrid + jgrid_real] + psi[jgrid].im * psitot[2*m*ngrid + jgrid_imag]);
    sum.imag += (psi[jgrid].re * psitot[2*m*ngrid + jgrid_imag] - psi[jgrid].im * psitot[2*m*ngrid + jgrid_real]);
  }
  sum.real *= dv;
  sum.imag *= dv;

  return(sum);
}

/*****************************************************************************/
