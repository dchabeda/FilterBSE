#include "mod_bse.h"

/***************************************************************************************/

void mod_bse(
  double complex*  psi_qp,
  double complex*  direct, 
  double complex*  exchange,
  double complex** bsmat, 
  double complex** bs_coeff, 
  double**         h0mat, 
  double**         xton_ene,
  double*          eig_vals,
	xyz_st*          s_mom, 
  xyz_st*          l_mom, 
  double complex*  l2_mom, 
  double complex*  ldots, 
  grid_st*         grid, 
  index_st*        ist, 
  par_st*          par,
  flag_st*         flag,
  parallel_st*     parallel
  ){

  ALLOCATE(bsmat, ist->n_xton * ist->n_xton, "bsmat");
  ALLOCATE(h0mat, ist->n_xton * ist->n_xton, "h0mat");
  ALLOCATE(xton_ene, ist->n_xton, "xton_ene");
  ALLOCATE(bs_coeff, ist->n_xton * ist->n_xton, "bs_coeff");
  
  build_BSE_mat(*bsmat, direct, exchange, ist);
  build_h0_mat(*h0mat, eig_vals, ist);
  
  bethe_salpeter(
    direct, exchange, *bsmat, *bs_coeff, *h0mat, *xton_ene, 
    s_mom, l_mom, l2_mom, ldots, grid, ist, par, flag, parallel
  );
  
  return;
}

/***************************************************************************************/

void build_h0_mat(
  double *h0mat, 
  double *eval, 
  index_st* ist
  ){

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

void build_BSE_mat(
  double complex* bsmat, 
  double complex* direct, 
  double complex* exchange, 
  index_st*       ist){

  FILE *ppsi;
  long ibs, jbs;
  long i, j;
  long ut; // upper triangle
  long lt; // lower triangle

  // Construct the BSE matrix from the exchange and direct kernels
  for (ibs = 0; ibs < ist->n_xton; ibs++){
    for (jbs = 0; jbs <= ibs; jbs++){
      ut = ibs*ist->n_xton + jbs;
      lt = jbs*ist->n_xton + ibs;
      // Symmetrize the matrices
      direct[lt] = conj(direct[ut]);
      
      exchange[lt] = conj(exchange[ut]);
      
      // Collect values for bsmat
      bsmat[ut] = direct[ut] + exchange[ut];
      bsmat[lt] = direct[lt] + exchange[lt];
      
      // Enforce Hermitivity by setting imag part of diag elems to 0.0
      if (ibs==jbs){
          bsmat[ut] = creal(bsmat[ut]) + 0.0*I; // 
      }
    }
  }

  
  ppsi = fopen("bsRE.dat", "w");
  for (i = 0; i < ist->n_xton; i++, fprintf(ppsi,"\n")){
    for (j = 0; j < ist->n_xton; j++){
      fprintf(ppsi,"%.*g ", PR_LEN, creal(bsmat[i*ist->n_xton+j]));
    }
  }
  fclose(ppsi);

  ppsi = fopen("bsIM.dat", "w");
  for (i = 0; i < ist->n_xton; i++, fprintf(ppsi,"\n")){
    for (j = 0; j < ist->n_xton; j++){
      fprintf(ppsi,"%.*g ", PR_LEN, cimag(bsmat[i*ist->n_xton+j]));
    }
  }
  fclose(ppsi);

  return;
}

/*****************************************************************************/
