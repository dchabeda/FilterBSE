#include "mod_portho.h"

void mod_portho(
  double**       psi_rank,
  double**       psitot,
  double*        pot_local,
  double*        eig_vals,
  double*        sigma_E,
  xyz_st*        R,
  zomplex*       LS,
  nlc_st*        nlc, 
  long*          nl, 
  double*        ksqr,
  grid_st*       grid,
  index_st*      ist,
  par_st*        par,
  flag_st*       flag,
  parallel_st*   parallel
  ){

  // Perform distributed memory orthogonalization on filtered states
  // Useful when the total size of the filtered states exceeds device RAM

  if (parallel->mpi_rank == 0){
    write_separation(stdout, "T");
    printf("****  DISTRIBUTED ORTHO  ***  DISTRIBUTED ORTHO  ***  DISTRIBUTED ORTHO   ****");
    write_separation(stdout, "B"); fflush(stdout);
  }
  
  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  const int mpir = parallel->mpi_rank;
  
  long tot_sz = ist->complex_idx*ist->nspinngrid*ist->mn_states_tot;
  printf("tot_sz = %ld\n", tot_sz);

  double init_clock;
  double init_wall;
  
  unsigned long        stlen          = ist->complex_idx * ist->nspinngrid;
  unsigned long        ns_p_rnk       = ist->n_states_per_rank;
  unsigned long        psi_rnk_sz     = ns_p_rnk * stlen;
  unsigned long        psi_rank_mem   = psi_rnk_sz * sizeof(double);
  unsigned long        max_node_mem   = 450UL * 1024 * 1024 * 1024;
  unsigned long        max_st_p_node  = max_node_mem / (stlen * sizeof(double));
  unsigned long        tot_mem;
  unsigned long        start;
  unsigned long        end;
  
  printf("stlen = %lu\n", stlen);
  printf("ns_p_rank = %lu\n", ns_p_rnk);
  printf("psi_rnk_sz = %lu\n", psi_rnk_sz);
  printf("psi_rank_mem = %lu\n", psi_rank_mem);
  printf("max_node_mem = %lu\n", max_node_mem);
  printf("max_st_p_node = %lu\n", max_st_p_node); fflush(0);

  // 450 GiB max per node is a heuristic to not OOM on Perlmutter

  /************************************************************/
  /********************     ALLOC MEM      ********************/
  /************************************************************/

  ALLOCATE(psi_rank, psi_rnk_sz, "psi_rank in portho");
  
  /************************************************************/
  /*******************   READ PSI BLOCKS    *******************/
  /************************************************************/

  start = mpir * ns_p_rnk;
  end   = (mpir + 1) * ns_p_rnk;

  printf("start = %lu\n", start);
  printf("end = %lu\n", end);

  read_psi_from_disk(*psi_rank, start, end, stlen, "psi-filt.dat");

  /************************************************************/
  /******************  BLOCK ORTHOGONALIZE  *******************/
  /************************************************************/

  if (mpir == 0){
    write_separation(stdout, "T");
    printf("\n5. DISTRIBUTED ORTHO | %s\n", get_time()); 
    write_separation(stdout, "B"); fflush(stdout);
  }

  init_clock = (double)clock(); 
  init_wall = (double)time(NULL);

  ns_p_rnk = ortho_cplx((MKL_Complex16*)(*psi_rank), grid->dv, ist, par, flag, parallel);

  normalize_all(*psi_rank, ist->mn_states_tot, ist, par, flag, parallel);
        
  if (mpir == 0){
    printf(
      "\ndone calculating distributed ortho, CPU time (sec) %g, wall run time (sec) %g\n",
      ((double)clock()-init_clock)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-init_wall
    ); 
    fflush(stdout);
  } 

  /************************************************************/
  /******************   BLOCK DIAGONALIZE   *******************/
  /************************************************************/

  if (mpir == 0){
    write_separation(stdout, "T");
    printf("\n6. DISTRIBUTED DIAGONALIZATION | %s\n", get_time()); 
    write_separation(stdout, "B"); fflush(stdout);
  }

  init_clock = (double)clock(); 
  init_wall = (double)time(NULL);

  diag_H(*psi_rank, pot_local, LS, nlc, nl, ksqr, eig_vals, ist, par, flag, parallel);
  normalize_all(*psi_rank, ns_p_rnk, ist, par, flag, parallel);
  

  /************************************************************/
  /*******************     JOIN BLOCKS     ********************/
  /************************************************************/

  // Select min(max_s_p_node, ns_p_rnk / mpi_size) states to join

  if (max_st_p_node < (ns_p_rnk / parallel->mpi_size)){
    printf("max_s_p_node = %lu < ns_p_rnk = %lu\n", max_st_p_node, ns_p_rnk);
    ns_p_rnk = max_st_p_node;
  }
  
  tot_mem = ns_p_rnk * stlen * sizeof(double);

  printf("JOINING BLOCKS:\n");
  printf("%lu states from each rank\n", ns_p_rnk);
  printf("Total memory on rank %d: %lu\n", parallel->mpi_root, tot_mem);

  if (tot_mem > max_node_mem){
    ns_p_rnk -= 10;
    printf("WARNING: initial join would exceed device RAM\n");
    printf("Decreasing ns_p_rnk to %lu\n", ns_p_rnk);
  }

  ist->mn_states_tot = ns_p_rnk;
  tot_sz = parallel->mpi_size * ns_p_rnk * stlen;
  psi_rnk_sz = ns_p_rnk * stlen;

  if (mpir == 0){
    ALLOCATE(psitot, tot_sz, "psitot in portho");
  }

  
  MPI_Gather(*psi_rank, psi_rnk_sz, MPI_DOUBLE, *psitot, psi_rnk_sz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  free(*psi_rank);
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  if (mpir == 0){

    /************************************************************/
    /******************   RE-ORTHOGONALIZE   ********************/
    /************************************************************/

    write_separation(stdout, "T");
    printf("\n6. RE-ORTHOGONALIZING STATES | %s\n", get_time()); 
    write_separation(stdout, "B"); fflush(stdout);
    
    init_clock = (double)clock(); 
    init_wall = (double)time(NULL);

    ist->mn_states_tot = ortho_cplx((MKL_Complex16*)(*psitot), grid->dv, ist, par, flag, parallel);

    normalize_all(*psitot, ist->mn_states_tot, ist, par, flag, parallel);

    printf(
      "\ndone calculating single node ortho, CPU time (sec) %g, wall run time (sec) %g\n",
      ((double)clock()-init_clock)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-init_wall
    ); 
    fflush(stdout);
  
    /************************************************************/
    /******************    RE-DIAGONALIZE    ********************/
    /************************************************************/

    
    write_separation(stdout, "T");
    printf("\n6. RE-DIAGONALIZING HAMILTONIAN | %s\n", get_time()); 
    write_separation(stdout, "B"); fflush(stdout);
    
  
    init_clock = (double)clock(); 
    init_wall = (double)time(NULL);
  
    diag_H(*psitot, pot_local, LS, nlc, nl, ksqr, eig_vals, ist, par, flag, parallel);
    normalize_all(*psitot, ist->mn_states_tot, ist, par, flag, parallel);
    
    printf(
      "\ndone calculating single node diag, CPU time (sec) %g, wall run time (sec) %g\n",
      ((double)clock()-init_clock)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-init_wall
    ); 
    fflush(stdout);

    /************************************************************/
    /********************    CALC SIGMA E    ********************/
    /************************************************************/

    mod_sigma(*psitot, pot_local, eig_vals, sigma_E, grid, LS, nlc, nl,
      ksqr, ist, par, flag, parallel
    );

    /************************************************************/
    /**********************    MOD OUTPUT    ********************/
    /************************************************************/

    mod_output(*psitot, R, eig_vals, sigma_E, grid, ist, par, flag, parallel);
    
    free(*psitot);
  }

  exit(0);
  
}

/************************************************************************/

void read_psi_from_disk(
  double*             psibuf, 
  unsigned long       start, 
  unsigned long       end, 
  unsigned long       stlen, 
  char*               filename
){

  // Read the file "filename" from start to end
  FILE *pf;

  unsigned long j;
  unsigned long offset;
  unsigned long cntr = 0;

  pf = fopen(filename, "r");

  if (pf == NULL){
    fprintf(stderr, "ERROR could not open file %s\n", filename);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  for (j = start; j < end; j++){
    // printf("Reading state %lu\n", j); fflush(0);
    offset = j * stlen * sizeof(psibuf[0]);

    fseek(pf, offset, SEEK_SET);
    
    fread(&psibuf[cntr * stlen], stlen * sizeof(double), 1, pf);
    cntr++;
  }
  
  printf("Done reading states from %s\n", filename); fflush(0);

  fclose(pf);

  return;
}


