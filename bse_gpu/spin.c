#include "spin.h"

void qp_spin_frac(
  double *restrict psi_qp, 
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

  const long               stlen = ist->complex_idx * ist->nspinngrid;
  const long               ngrid = ist->ngrid;
  const long               cplx_idx = ist->complex_idx;

  unsigned long            i;
  unsigned long            i_st;

  double *restrict         pct_spn;
  
  ALLOCATE(&pct_spn, ist->nspin * ist->n_qp, "per_up");
  
  /************************************************************/
	/*******************  COMPUTE SPIN FRAC   *******************/
	/************************************************************/

  printf("\nComputing spin composition of quasiparticle spinors | %s\n", get_time()); 
  fflush(stdout);
  
  double start = omp_get_wtime();
  #pragma omp parallel for private(i_st)
  for (i  = 0; i < ist->n_qp; i++){
    i_st = i * stlen;
    
    unsigned long            jg;    // jgrid
    unsigned long            jgr;   // jgrid real
    unsigned long            jgi;   // jgrid imag

    double                   psi_i_r;
    double                   psi_i_i;
    double                   pct = 0;
    
    #pragma omp simd safelen(8) aligned(psi_qp: 32) reduction(+:pct)
    for (int s = 0; s < 2; s++){
      pct = 0;
      for (jg = 0; jg < ngrid; jg++){ 
        jgr = cplx_idx * (jg + s * ngrid);
        jgi = jgr + 1;

        psi_i_r = psi_qp[i_st + jgr];
        psi_i_i = psi_qp[i_st + jgi];
        
        pct += sqr(psi_i_r) + sqr(psi_i_i);
      }
      pct_spn[s * ist->n_qp + i] = pct * par->dv;
      // sprintf(fileName, "rho%d-%d.cube", i, s);
      // write_cube_file(rho, grid, fileName);
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