#include "spin.h"

void qp_spin_frac(
  double complex*  psi_qp, 
  double*          eig_vals,
  grid_st*         grid,
  index_st*        ist,
  par_st*          par,
  flag_st*         flag,
  parallel_st*     parallel
  ){

  /************************************************************/
	/*******************  DECLARE VARIABLES   *******************/
	/************************************************************/

  FILE*                    pf   = fopen("qp_spins.dat", "w");

  const unsigned long      nspngr = ist->nspinngrid;
  const unsigned long      ngrid = ist->ngrid;

  unsigned long            i;
  unsigned long            i_st;

  double *                 pct_spn;
  const double             dv = par->dv;
  
  ALLOCATE(&pct_spn, ist->nspin * ist->n_qp, "per_spn");
  
  /************************************************************/
	/*******************  COMPUTE SPIN FRAC   *******************/
	/************************************************************/

  printf("\nComputing spin composition of quasiparticle spinors | %s\n", get_time()); 
  fflush(stdout);
  
  double start = omp_get_wtime();
  #pragma omp parallel for private(i_st)
  for (i  = 0; i < ist->n_qp; i++){
    i_st = i * nspngr;

    unsigned long  jg;  // jgrid
    unsigned long  jsg; // j(spin)grid

    double         pct = 0;
    
    for (int s = 0; s < 2; s++){
      pct = 0;
      #pragma omp simd safelen(8) aligned(psi_qp: BYTE_BOUNDARY) reduction(+:pct)
      for (jg = 0; jg < ngrid; jg++){
        jsg = jg + s * ngrid;
        pct += cnorm(psi_qp[i_st + jsg]);
      }
      pct_spn[s * ist->n_qp + i] = pct * dv;
    }
  }

  double end = omp_get_wtime();

  printf("\n Time for computing all spins: %lf\n", (end - start));

  
  for (i = 0; i < ist->n_qp; i++){
    fprintf(pf, "\nStats on state%d: (E=%lg)\n", i, eig_vals[i]);
    for (int s = 0; s < 2; s++){
      fprintf(pf, " Spin %d fraction: %f\n", s, pct_spn[s * ist->n_qp + i]);
    }
  }
  
  fclose(pf);

  free(pct_spn);

  return;
}