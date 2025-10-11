#include "mod_init.h"

void mod_init(
  double complex**  psitot,
  double complex**  psi_qp,
  double**          eig_vals,
  double**          sigma_E,
  xyz_st**          R,
  grid_st*          grid,
  double**          gridx,
  double**          gridy,
  double**          gridz,
  index_st*         ist,
  par_st*           par,
  flag_st*          flag,
  parallel_st*      parallel
  ){

  /************************************************************/
	/*******************  DECLARE VARIABLES   *******************/
	/************************************************************/

  long             jg;
  
  const int        mpir = parallel->mpi_rank;

  /************************************************************/
	/********************     READ INPUT     ********************/
	/************************************************************/

  if (mpir == 0){
    write_separation(stdout, "T");
    printf("\n1.\tINITIALIZING JOB | %s\n", get_time());
    write_separation(stdout, "B"); fflush(stdout);
  } 
  
      
  if (0 == flag->initUnsafe){
    if (mpir == 0) printf("\nReading filter output from output.dat:\n"); fflush(stdout);
    read_filter_output("output.dat", psitot, eig_vals, sigma_E, R, grid, gridx, gridy, gridz, ist, par, flag);

  } 
  else if (1 == flag->initUnsafe){
    if (mpir == 0) printf("\nReading UNSAFE input from unsafe_input.par:\n"); fflush(stdout);
    read_unsafe_input(
      psitot, eig_vals, sigma_E, R, grid, gridx, gridy, gridz, 
      ist, par, flag, parallel
    );
  }

  
  /*** Read initial setup from input.par ***/
  if (mpir == 0) printf("\nReading BSE job specifications from input.par:\n"); fflush(stdout);
  
  read_input(flag, grid, ist, par, parallel);

  /************************************************************/
	/*******************   GET QP BASIS IDXS   ******************/
	/************************************************************/

  if (mpir == 0) printf("\nSetting quasiparticle basis indices:\n"); fflush(stdout);
  
  // Allocate memory for the lists of the indices of eigenstates
  // The maximum possible number of hole states is mn_states_tot from filter.
  // We allocate an entire block of that size for both elecs and holes 
  // because we will reallocate after the get_qp_basis_indices function.
  
  ALLOCATE(&(ist->eval_hole_idxs), ist->mn_states_tot, "eval_hole_idxs");
  ALLOCATE(&(ist->eval_elec_idxs), ist->mn_states_tot, "eval_elec_idxs");
  
  get_qp_basis_indices(*eig_vals, *sigma_E, &ist->eval_hole_idxs, &ist->eval_elec_idxs, ist, par, flag, parallel);
  
  
  /************************************************************/
	/********************    BUILD QP BASIS    ******************/
	/************************************************************/

  if (mpir == 0) printf("\nReading quasiparticle basis wavefunctions:\n"); fflush(stdout);
  
  ALLOCATE(psi_qp, ist->n_qp * ist->nspinngrid, "psi_qp");
  
  get_qp_basis(*psi_qp, *psitot, *eig_vals, *sigma_E, ist, par, flag);
  
  ist->n_xton = ist->n_elecs * ist->n_holes;
  
  /************************************************************/
	/********************    SET LOOP IDXS     ******************/
	/************************************************************/
  // ALLOCATE(&ist->jsg_conv, ist->nspinngrid, "jsg");
  
  // #pragma omp simd safelen(4)
  // for (int s = 0; s < 2; s++){
  //   for (jg = 0; jg < ist->ngrid; jg++){ 
  //     long jsg = jg + s * ist->ngrid;
  //     ist->jsg_conv[jsg] = jg;
  //   }
  // }

  free(*psitot);                psitot = NULL;
  free(*sigma_E);               sigma_E = NULL;
  free(ist->eval_elec_idxs);    ist->eval_elec_idxs = NULL;
  free(ist->eval_hole_idxs);    ist->eval_hole_idxs = NULL;
  
  return;
}