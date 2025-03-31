/*****************************************************************************/

#include "basis.h"

/*****************************************************************************/

void get_qp_basis_indices(
  double*         eig_vals, 
  double*         sigma_E, 
  long**          eval_hole_idxs, 
  long**          eval_elec_idxs, 
  index_st*       ist, 
  par_st*         par, 
  flag_st*        flag, 
  parallel_st*    parallel
  ){
  
  /*******************************************************************
  * This function determines the indices of the electron and hole    *
  * states that will be used in the quasiparticle basis. The eigen-  *
  * values and sigma_E values are read from the formatted filter     *
  * output file, output.dat. Based on the user's request for max     *
  * number of states in the VB and CB, the indices will be set to    *
  * read the appropriate wavefunctions from psitot.                  *
  * inputs:                                                          *
  *  [eig_vals] mn_states_tot-long array of eig vals from eval.dat   *
  *  [sigma_E] mn_states_tot-long array of sigma vals from eval.dat  *
  *  [eval_hole_idxs] ptr to arr for the indices of usable h+ states *
  *  [eval_elec_idxs] ptr to arr for the indices of usable e- states *
  *  [ist] ptr to counters, indices, and lengths                     *
  *  [par] ptr to par_st holding VBmin, VBmax... params              *
  *  [flag] ptr to flag_st holding job flags                         *
  * outputs: void                                                    *
  ********************************************************************/
  
  FILE *pf;

  const int mpir = parallel->mpi_rank;
  long i, cntr, eval_homo_idx, eval_lumo_idx, old_n_holes;
  double deltaE;
  
  // Allocate memory for eval_hole_idxs etc if it is NULL
  if (NULL == *eval_hole_idxs){
    printf("\tAllocating new memory for eval_hole_idxs\n");
    *eval_hole_idxs = malloc(500 * sizeof(long));
  }
  if (NULL == *eval_elec_idxs){
    printf("\tAllocating new memory for eval_hole_idxs\n");
    *eval_elec_idxs = malloc(500 * sizeof(long));
  }

  // Set the quasiparticle basis indices to zero
  ist->homo_idx = ist->lumo_idx = ist->n_holes = ist->n_elecs = 0;

  cntr = 0; // We will reorder eig_vals and sigma_E to only contain the converged eigenstates
  // Get the indices of the HOMO, LUMO, and also the number of eigstates in VB/CB 
  for (i = 0; i < ist->mn_states_tot; i++){
    // Find HOMO and total number of VB states
    // If the state has variance < cutoff and energy below fermiE, it is a hole eigenstate
    if (sigma_E[i] < par->sigma_E_cut && eig_vals[i] < par->fermi_E){
      ist->n_holes++; // increment the counter for number of hole states
      ist->homo_idx = cntr; // homo_idx in RAM is the value of cntr for which condition is met
      eval_homo_idx = i; // the value of homo_idx from the eval.dat file
      
      eig_vals[cntr] = eig_vals[i]; // reorder the eig_vals
      sigma_E[cntr] = sigma_E[i]; // reorder the sigma_E
      (*eval_hole_idxs)[cntr] = i; // add this index to the array holding all indices from eval.dat
      cntr++;
    }
    // Find LUMO and total number of CB states
    // If the state has variance < cutoff and energy above fermiE, it is an elec eigenstate
    if (sigma_E[i] < par->sigma_E_cut && eig_vals[i] > par->fermi_E){
      if (0 == ist->lumo_idx){
        eval_lumo_idx = i;
        ist->lumo_idx = ist->homo_idx + 1; // Get first value of i for which condition is met
      }
      ist->n_elecs++; 
      
      eig_vals[cntr] = eig_vals[i];
      sigma_E[cntr] = sigma_E[i];
      (*eval_elec_idxs)[(long) (cntr - ist->n_holes)] = i;
      cntr++;
    }
  }
  old_n_holes = ist->n_holes;

  if (mpir == 0){
    printf("\n\tTotal # of filtered hole eigenstates = %ld\n", ist->n_holes);
    printf("\tTotal # of filtered electron eigenstates = %ld\n", ist->n_elecs);
    printf("\tThe filter eval.dat index of the HOMO state = %ld  LUMO state = %ld\n", eval_homo_idx, eval_lumo_idx);  
  }
  
  // Check how the quasiparticle basis compares to the desired energy range
  // First, in the VB
  deltaE = eig_vals[ist->homo_idx] - eig_vals[ist->homo_idx - ist->n_holes + 1];
  if (mpir == 0){
    if (deltaE < par->delta_E_hole){
      printf("\n\tUnconstrained energy span of holes, %.3lg a.u. (%.2lg eV) < desired span = %.3lg a.u. (%.2lg eV)\n", deltaE, deltaE * 27.211, par->delta_E_hole, par->delta_E_hole * 27.211 );
      printf("\t**Increase size of VB basis states to reach desired result**\n");
    } else {
      printf("\n\tUnconstrained energy span of holes %.3lg a.u. (%.2lg eV) > desired span = %.3lg a.u. (%.2lg eV)\n", deltaE, deltaE * 27.211, par->delta_E_hole, par->delta_E_hole * 27.211);
      printf("\t**Fewer VB basis states would reach the desired result**\n");
    }
  }
  
  // Then in the conduction band
  deltaE = eig_vals[ist->n_elecs + ist->lumo_idx - 1] - eig_vals[ist->lumo_idx];
  if (mpir == 0){
    if (deltaE < par->delta_E_hole){
      printf("\tUnconstrained energy span of elecs, %.3lg a.u. (%.2lg eV) < desired span = %.3lg a.u. (%.2lg eV)\n", deltaE, deltaE * 27.211, par->delta_E_elec, par->delta_E_elec * 27.211);
      printf("\t**Increase size of CB basis states to reach desired result**\n");
    } else {
      printf("\tUnconstrained energy span of elecs, %.3lg a.u. (%.2lg eV) > desired span = %.3lg a.u. (%.2lg eV)\n", deltaE, deltaE * 27.211, par->delta_E_elec, par->delta_E_elec * 27.211);
      printf("\t**Fewer CB basis states would reach the desired result**\n");
    }
  }
  

  // If the max number of hole states was set, then constrain the hole eigenbasis
  // We will also modify the eig_vals and sigma_E array again so they only contain the
  // energies of the selected qp basis states
  
  if ((-1 != ist->max_hole_states) && (ist->n_holes > ist->max_hole_states)){
    
    ist->n_holes = ist->max_hole_states;
    if (mpir == 0) printf("\n\tConstraining hole basis states to maxHoleStates\n\t  new ist->n_holes = %ld\n", ist->n_holes);
    
    // Reorder eig_vals and sigma_E to only ist->max_hole_states 
    // number of hole eigenstates closest to the band edge
    cntr = 0;
    for (i = 0; i < ist->n_holes; i++){
      // Replace the ith element of eig_vals with the value (max_hole_states - i) below the band edge
      eig_vals[i] = eig_vals[ist->homo_idx - ist->n_holes + 1 + i];
      sigma_E[i] = sigma_E[ist->homo_idx - ist->n_holes + 1 + i];
      (*eval_hole_idxs)[i] = (*eval_hole_idxs)[(long) ( ist->homo_idx - ist->n_holes + 1 + i)];
      cntr++;
    }
    // check that the counter has the same value as ist->n_holes
    if ((cntr - ist->n_holes) != 0){
      printf("ERROR: something went wrong reordering eig_vals in get_qp_basis_indices\n");
      exit(EXIT_FAILURE);
    }
    ist->homo_idx = ist->homo_idx - (old_n_holes - ist->n_holes);
    // Determine the new energy span
    deltaE = eig_vals[ist->homo_idx] - eig_vals[ist->homo_idx - ist->n_holes + 1];
    if (mpir == 0){
      if (deltaE < par->delta_E_hole){
        printf("\tConstrained energy span of holes, %.3lg a.u. (%.2lg eV) < desired span = %lg a.u. (%.2lg eV)\n", deltaE, deltaE * 27.211, par->delta_E_hole, par->delta_E_hole * 27.211);
        printf("\t**Increase size of VB basis states to reach desired result**\n");
      } else {
        printf("\n\tConstrained energy span of holes %.3lg a.u. (%.2lg eV) > desired span = %lg a.u. (%.2lg eV)\n", deltaE, deltaE * 27.211, par->delta_E_hole, par->delta_E_hole * 27.211);
        printf("\t**Fewer VB basis states would reach the desired result**\n");
      }
    }
    
    
    ist->lumo_idx = ist->homo_idx + 1; 
    // set the lumo_idx here to handle the case where the next if statement does not evaluate as true
    // in that case, the lumo_idx will be left as the one from the old_n_holes
  }

  // If the max number of elec states was set, then constrain the elec eigenbasis
  if ((-1 != ist->max_elec_states) && (ist->n_elecs > ist->max_elec_states) ){
    
    ist->n_elecs = ist->max_elec_states;
    
    if (mpir == 0) printf("\tConstraining elec basis states to maxElecStates\n\t  new ist->n_elecs = %ld\n", ist->n_elecs);
    
    // Reorder eig_vals and sigma_E to only contain eigenstates
    cntr = ist->homo_idx + 1;
    for (i = old_n_holes; i < old_n_holes + ist->n_elecs; i++){
      eig_vals[cntr] = eig_vals[i];
      sigma_E[cntr] = sigma_E[i];
      // (*eval_elec_idxs)[(long) cntr - ist->n_holes] = (*eval_elec_idxs)[i - ist->n_holes];
      // we have to subtract by n_holes because we start the iterator from n_holes + 1
      cntr++;
    }
 
    // Determine the new energy span
    deltaE = eig_vals[ist->n_elecs + ist->lumo_idx - 1] - eig_vals[ist->lumo_idx];
    if (mpir == 0){
      if (deltaE < par->delta_E_elec){
        printf("\tConstrained energy span of elecs, %.3lg a.u. (%.2lg eV) < desired span = %.3lg a.u. (%.2lg eV)\n", deltaE, deltaE * 27.211, par->delta_E_elec, par->delta_E_elec * 27.211);
        printf("\t**Increase size of CB basis states to reach desired result**\n");
      } 
      else {
        printf("\tConstrained energy span of elecs %.3lg a.u. (%.2lg eV) > desired span = %.3lg a.u. (%.2lg eV)\n", deltaE, deltaE * 27.211, par->delta_E_elec, par->delta_E_elec * 27.211);
        printf("\t**Fewer CB basis states would reach the desired result**\n");
      }
    }

  } 
  // If the user does not specify max_elec or max_hole_states, then use all the quasiparticle states for BSE
  if ((-1 == ist->max_elec_states) && (-1 == ist->max_elec_states)) {
    printf("\n\tAll filtered quasiparticle states will be used for exciton basis\n");
  }

  ist->n_qp = ist->n_holes + ist->n_elecs;

  // Print QP basis info
  if (mpir == 0){
    printf("\n\tSelected # of filtered h+ qp basis states = %ld\n", ist->n_holes);
    printf("\tSelected # of filtered e- qp basis states = %ld\n", ist->n_elecs);
    printf("\tTotal number of quasiparticle states, n_qp = %ld\n", ist->n_qp);
    printf("\tThe BSEeval.par index of the HOMO state = %ld  LUMO state = %ld\n", ist->homo_idx, ist->lumo_idx);
    printf("\tThe HOMO energy = % .6g a.u. % .3f eV\n", eig_vals[ist->homo_idx], eig_vals[ist->homo_idx]*AUTOEV);
    printf("\tThe LUMO energy = % .6g a.u. % .3f eV\n", eig_vals[ist->lumo_idx], eig_vals[ist->lumo_idx]*AUTOEV);
    printf("\tFundamental gap = %.6g a.u. %.3f eV\n", eig_vals[ist->lumo_idx]-eig_vals[ist->homo_idx], (eig_vals[ist->lumo_idx]-eig_vals[ist->homo_idx])*AUTOEV);
  }
  
  // Print BSEeval.par
  if (mpir == 0){
    pf = fopen("BSEeval.par", "w");
    for (i = 0; i < ist->n_qp; i++){
      fprintf(pf, "% .12f %lg\n", eig_vals[i], sigma_E[i]);
    }
    fclose(pf);
  }
  

  return;
}

/*****************************************************************************//*****************************************************************************/

void get_qp_basis(
  double complex* psi_qp, 
  double complex* psitot, 
  double*         eig_vals, 
  double*         sigma_E, 
  index_st*       ist, 
  par_st*         par, 
  flag_st*        flag
  ){
  
  /*******************************************************************
  * This function copies the wavefunctions for the elec and hole     *
  * states for the quasiparticle basis into psi_hole and psi_elec.   *
  * inputs:                                                          *
  *  [eig_vals] mn_states_tot-long array of eig vals from eval.dat   *
  *  [sigma_E] mn_states_tot-long array of sigma vals from eval.dat  *
  *  [eval_hole_idxs] ptr to arr for the indices of usable h+ states *
  *  [eval_elec_idxs] ptr to arr for the indices of usable e- states *
  *  [ist] ptr to counters, indices, and lengths                     *
  *  [par] ptr to par_st holding VBmin, VBmax... params              *
  *  [flag] ptr to flag_st holding job flags                         *
  * outputs: void                                                    *
  ********************************************************************/
 
  long i, cntr, state_idx;
  const unsigned long nspngr = ist->nspinngrid;
  
  // Copy the hole states
  cntr = 0;
  for (i = 0; i < ist->n_holes; i++){
    state_idx = ist->eval_hole_idxs[i];
    memcpy(&psi_qp[i * nspngr], &psitot[state_idx * nspngr], nspngr * sizeof(psitot[0])) ;
    cntr++;
  }

  // Copy the electron states
  for (i = 0; i < ist->n_elecs; i++){
    state_idx = ist->eval_elec_idxs[i];
    memcpy(&psi_qp[cntr * nspngr], &psitot[state_idx * nspngr], nspngr * sizeof(psitot[0]));
    cntr++;
  } 

  // printf("   Finished reading in quasiparticle basis functions\n");

  return;

}

