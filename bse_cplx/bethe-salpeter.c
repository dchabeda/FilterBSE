/*****************************************************************************/

#include "bethe-salpeter.h"

/*****************************************************************************/

void bethe_salpeter(
  double complex*  direct, 
  double complex*  exchange,
  double complex*  bsmat, 
  double complex*  bs_coeff, 
  double*          h0mat, 
  double*          xton_ene,
  grid_st*         grid, 
  index_st*        ist, 
  par_st*          par,
  flag_st*         flag,
  parallel_st*     parallel
  ){

  /************************************************************/
	/*******************  DECLARE VARIABLES   *******************/
	/************************************************************/

  FILE*           pf;

  unsigned long   a;
  unsigned long   b;
  unsigned long   i;
  unsigned long   j;
  unsigned long   k;
  unsigned long   l;
  unsigned long   ibs;
  unsigned long   jbs;
  unsigned long   idx;

  const long      n_el   = ist->n_elecs;
  const long      n_ho   = ist->n_holes;
  const long      lidx   = ist->lumo_idx;
  const long      n_xton = ist->n_xton;

  double complex* mat;
  double complex* h;
  double complex  sum;
  double complex  tmp;

  double          start_t;
  double          end_t;

  long*           listibs;

  omp_set_num_threads(parallel->nthreads);

  ALLOCATE(&listibs, n_xton, "listibs in bethe-salpeter");

  for (ibs = 0, a = lidx; a < lidx + n_el; a++) {
    for (i = 0; i < n_ho; i++, ibs++) {
      listibs[(a - lidx) * n_ho + i] = ibs;
    }
  }

  /************************************************************/
	/********************   BUILD BSE "H"    ********************/
	/************************************************************/

  ALLOCATE(&mat, n_xton * n_xton, "mat in bethe-salpeter");
  ALLOCATE(&h,   n_xton * n_xton, "h in bethe-salpeter");
  
  #pragma omp parallel for private(i)
  for (i = 0; i < n_xton*n_xton; i++) {
    h[i] = bs_coeff[i] = h0mat[i] - bsmat[i];
  }

  // Print out the full BSE mat
  pf = fopen("HBSmatRE.dat", "w");
  for (i = 0; i < n_xton; i++, fprintf(pf,"\n")) {
    for (j = 0; j < n_xton; j++) {
      fprintf (pf,"%.*g \t", PR_LEN, creal(h[i*n_xton+j]));
	}
  }
  fclose(pf);

  pf = fopen("HBSmatIM.dat", "w");
  for (i = 0; i < n_xton; i++, fprintf(pf,"\n")) {
    for (j = 0; j < n_xton; j++) {
      fprintf (pf,"%.*g \t", PR_LEN, cimag(h[i*n_xton+j]));
  }
  }
  fclose(pf);

  /************************************************************/
	/********************    DIAG BSE MAT    ********************/
	/************************************************************/

  // Diagonalize the BSE matrix to obtain the coefficients
  start_t = omp_get_wtime();
  diag((int)n_xton, ist->nthreads, bs_coeff, xton_ene);
  end_t = omp_get_wtime();

  printf("Done diagonalizing BSE matrix | %s\n", format_duration(end_t - start_t));
  fflush(0);

  // Prints the coefficients for the first 100 (or n_xton) lowest energy excitonic states
  long numExcStatesToPrint = 100;
  if (n_xton < numExcStatesToPrint) numExcStatesToPrint = n_xton;

  pf = fopen("BSEcoeffRE.dat", "w");
  for (i = 0; i < n_xton; i++, fprintf(pf,"\n")) {
    for (j = 0; j < n_xton; j++) {  
	  fprintf (pf,"%.*g\t", PR_LEN, creal(bs_coeff[i*n_xton+j]));
    }
  }
  fclose(pf);

  pf = fopen("BSEcoeffIM.dat", "w");
  for (i = 0; i < n_xton; i++, fprintf(pf,"\n")) {
    for (j = 0; j < n_xton; j++) {  
	  fprintf (pf,"%.*g\t", PR_LEN, cimag(bs_coeff[i*n_xton+j]));
    }
  }
  fclose(pf);

  
  /************************************************************/
	/******************    GROUND XTON ENE    *******************/
	/************************************************************/

  sum = 0.0 + 0.0*I;

  for (i = 0; i < n_xton; i++) {
    for (j = 0; j < n_xton; j++) {
      tmp = bs_coeff[i * n_xton] * h[j*n_xton + i];
      sum += conjmul(bs_coeff[j * n_xton], tmp);
    }
  }
  
  printf("\nGround state exciton has energy = %.5f a.u. | %.5f eV (%.5f Imag)\n", creal(sum), creal(sum)*AUTOEV, cimag(sum));
  fflush(0);

  /************************************************************/
	/*****************  SOLVE BSE | ALL XTONS  ******************/
	/************************************************************/

  // Compute $mat = h \cdot u$

  if (flag->timingSpecs) start_t = omp_get_wtime();
  #pragma omp parallel for private(l,j,k,sum)
  for (l = 0; l < n_xton; l++) {
    for (j = 0; j < n_xton; j++) {
      sum = 0.0 + 0.0*I;
      for (k = 0; k < n_xton; k++) {
        sum += h[l*n_xton+k] * bs_coeff[k*n_xton+j];
      }
      mat[l*n_xton+j] = sum;
    }
  }
  if (flag->timingSpecs){
    end_t = omp_get_wtime();
    printf("  Step 1: %s\n", format_duration(end_t - start_t));
  }
  
  //compute $u^\dagger \cdot mat = u^\dagger \cdot h \cdot u$
  
  if (flag->timingSpecs) start_t = omp_get_wtime();
  #pragma omp parallel for private(i,j,l,sum)
  for (i = 0; i < n_xton; i++) {
    for (j = 0; j < n_xton; j++) {
      sum = 0.0 + 0.0*I;
      for (l = 0; l < n_xton; l++) {
        sum += conjmul(bs_coeff[l*n_xton+i], mat[l*n_xton+j]);
      }
      h[i*n_xton+j] = sum;
    }
  }
  if (flag->timingSpecs){
    end_t = omp_get_wtime();
    printf("  Step 2: %s\n", format_duration(end_t - start_t));
  }
  
  if (flag->timingSpecs) start_t = omp_get_wtime();
  #pragma omp parallel for private(l,j,k,sum)
  for (l = 0; l < n_xton; l++) {
    for (j = 0; j < n_xton; j++) {
      sum = 0.0 + 0.0*I;
      for (k = 0; k < n_xton; k++) {
        sum += h0mat[l*n_xton + k] * bs_coeff[k*n_xton + j];
      }
      mat[l*n_xton+j] = sum;
    }
  }
  if (flag->timingSpecs){
    end_t = omp_get_wtime();
    printf("  Step 3: %s\n", format_duration(end_t - start_t));
  }
  
  if (flag->timingSpecs) start_t = omp_get_wtime();
  #pragma omp parallel for private(i,j,l,sum)
  for (i = 0; i < n_xton; i++) {
    for (j = 0; j < n_xton; j++) {
      sum = 0.0 + 0.0*I;
      for (l = 0; l < n_xton; l++) {
        sum += conjmul(bs_coeff[l*n_xton+i], mat[l*n_xton+j]);
      }
      h0mat[i*n_xton+j] = creal(sum);
    }
  }
  if (flag->timingSpecs){
    end_t = omp_get_wtime();
    printf("  Step 4: %s\n", format_duration(end_t - start_t));
  }
  
  if (flag->timingSpecs) start_t = omp_get_wtime();
  #pragma omp parallel for private(l,j,k,sum)
  for (l = 0; l < n_xton; l++) {
    for (j = 0; j < n_xton; j++) {
      sum = 0.0 + 0.0*I;
      for (k = 0; k < n_xton; k++) {
        sum += direct[l*n_xton+k] * bs_coeff[k*n_xton + j];
      }
      mat[l*n_xton+j] = sum;
    }
  }
  if (flag->timingSpecs){
    end_t = omp_get_wtime();
    printf("  Step 5: %s\n", format_duration(end_t - start_t));
  }

  //compute $u^\dagger \cdot mat = u^\dagger \cdot h \cdot u$
  
  if (flag->timingSpecs) start_t = omp_get_wtime();
  #pragma omp parallel for private(i,j,l,sum)
  for (i = 0; i < n_xton; i++) {
    for (j = 0; j < n_xton; j++) {
      sum = 0.0 + 0.0*I;
      for (l = 0; l < n_xton; l++) {
        sum += conjmul(bs_coeff[l*n_xton+i], mat[l*n_xton+j]);
      }
      direct[i*n_xton+j] = sum;
    }
  }
  if (flag->timingSpecs){
    end_t = omp_get_wtime();
    printf("  Step 6: %s\n", format_duration(end_t - start_t));
  }
  
  if (flag->timingSpecs) start_t = omp_get_wtime();
  #pragma omp parallel for private(l,j,k,sum)
  for (l = 0; l < n_xton; l++) {
    for (j = 0; j < n_xton; j++) {
      sum = 0.0 + 0.0*I;
      for (k = 0; k < n_xton; k++) {
        sum += exchange[l*n_xton+k] * bs_coeff[k*n_xton+j];
      }
      mat[l*n_xton+j] = sum;
    }
  }
  if (flag->timingSpecs){
    end_t = omp_get_wtime();
    printf("  Step 7: %s\n", format_duration(end_t - start_t));
  }

  //compute $u^\dagger \cdot mat = u^\dagger \cdot h \cdot u$
  
  if (flag->timingSpecs) start_t = omp_get_wtime();
  #pragma omp parallel for private(i,j,l,sum)
  for (i = 0; i < n_xton; i++) {
    for (j = 0; j < n_xton; j++) {
      sum = 0.0 + 0.0*I;
      for (l = 0; l < n_xton; l++) {
        sum += conjmul(bs_coeff[l*n_xton+i], mat[l*n_xton+j]);
      }
      exchange[i*n_xton+j] = sum;
    }
  }
  if (flag->timingSpecs){
    end_t = omp_get_wtime();
    printf("  Step 8: %s\n", format_duration(end_t - start_t));
  }
  fflush(0);
  /************************************************************/
	/******************   PRINT EXCITON.DAT   *******************/
	/************************************************************/
  
  pf = fopen("exciton.dat" , "w");

  fprintf(pf,"#n \t E_n \t <H> \t <H_dir> \t <H_exc> \t <H_0> \t E_B (eV)\n");

  for (i = 0; i < n_xton; i++) {  
    fprintf(pf,"%ld % .12f % .12f % .12f  % .12f  % .12f  % .12f\n", i, 
      xton_ene[i], creal(h[i*n_xton+i]), creal(direct[i*n_xton+i]), creal(exchange[i*n_xton+i]),
     h0mat[i*n_xton+i], (xton_ene[i]-h0mat[i*n_xton+i])*AUTOEV);
  }
  fclose(pf);

  
  free(h); 
  free(mat);
  free(listibs);
  
  return;
}

/*****************************************************************************/
