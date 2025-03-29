#include "mod_optional_output.h"

void mod_optional_output(
  double*       psitot,
  xyz_st*       R,
  double*       eig_vals,
  double*       sigma_E,
  grid_st*      grid,
  double*       ksqr,
  index_st*     ist,
  par_st*       par,
  flag_st*      flag,
  parallel_st*  parallel){

  if (mpir == 0) {
    write_separation(stdout, "T");
    printf("\nCALCULATING OPTIONAL OUTPUT | %s\n", get_time()); 
    write_separation(stdout, "B"); fflush(stdout);
  }

  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  int           i;
  long          jstate;
  long          jgrid;
  long          jgrid_real;
  long          jgrid_imag;
  long          jspin;
  long          jmn;
  long          jns;

  const int     mpir = parallel->mpi_rank;

  /************************************************************/
  /*******************   PRINT CUBE FILES   *******************/
  /************************************************************/

  
  long i, a, ieof, nval;
  char str[50];
  double evalloc, deloc;
  FILE *pf;

  ist->homo_idx = ist->lumo_idx = 0;
  pf = fopen("eval.dat" , "r");
  for (i = ieof = 0; ieof != EOF; i++){
      ieof = fscanf(pf, "%ld %lg %lg", &a, &evalloc, &deloc);
      if (deloc < par->sigma_E_cut && evalloc < par->fermi_E) ist->homo_idx = i;
      if (i > ist->mn_states_tot){
      if (mpir == 0) printf("No hole states converged to within %lg a.u.\n", par->sigma_E_cut);
      break;
  }
  }
  fclose(pf);

  // nval = i - 1;
  pf = fopen("eval.dat" , "r");
  for (i = 0; i <= ist->homo_idx; i++) {
      fscanf(pf, "%ld %lg %lg", &a, &evalloc, &deloc);
      if (i > ist->mn_states_tot){
      if (mpir == 0) printf("No electron states converged to within %lg a.u.\n", par->sigma_E_cut);
      break;
      }
  }
  for (i = ist->homo_idx+1, ieof = 0; ieof != EOF; i++) {
      fscanf(pf, "%ld %lg %lg", &a, &evalloc, &deloc);
      if (deloc < par->sigma_E_cut) {
      ist->lumo_idx = i;
      break;
      }
      if (i > ist->mn_states_tot){
      if (mpir == 0) printf("No electron states converged to within %lg a.u.\n", par->sigma_E_cut);
      break;
      }
  }
  fclose(pf);

  if (mpir == 0) printf("index of homo, homo_idx = %ld; index of lumo, lumo_idx = %ld\n", ist->homo_idx, ist->lumo_idx); fflush(0);
  // Set the total number of electron and hole states in order to calculate the potential overlap integrals
  ist->total_homo = ist->homo_idx + 1; ist->total_lumo = ist->mn_states_tot - ist->total_homo;
  if (mpir == 0) printf("total_homo = %ld total_lumo = %ld\n", ist->total_homo, ist->total_lumo); fflush(0);


  if (flag->printCubes == 1){
      /*** Write homo and lumo cube files ***/
      
      if (mpir == 0) write_separation(stdout, top);
      if (mpir == 0) printf("\nWRITING CUBE FILES\n"); 
      if (mpir == 0) write_separation(stdout, bottom); fflush(stdout);

      if ((ist->homo_idx == 0) || (ist->lumo_idx == 0)){
      if (mpir == 0) printf("\nDid not converge enough electron or hole states to visualize cube files.\n");
      } else{
      if ((rho = (double *) calloc(ist->ngrid, sizeof(double))) == NULL){
      if (mpir == 0) fprintf(stderr, "\nOUT OF MEMORY: rho\n\n"); exit(EXIT_FAILURE);
      }

      inital_clock_t = (double)clock(); initial_wall_t = (double)time(NULL);

      for (i = 0; (i < ist->total_homo) && (i < ist->ncubes); i++){
      //Spin Up Wavefunction
      sprintf(str,"homo-%ld-Up.cube",i);
      for (jgrid = 0; jgrid < ist->ngrid; jgrid++){
          jgrid_real = ist->complex_idx * jgrid;
          jgrid_imag = ist->complex_idx * jgrid + 1;
          
          rho[jgrid] = sqr(psitot[ist->complex_idx*(ist->homo_idx-i)*ist->nspinngrid + jgrid_real]);
          if (1 == flag->isComplex) rho[jgrid] += sqr(psitot[ist->complex_idx*(ist->homo_idx-i)*ist->nspinngrid + jgrid_imag]);
      }
      write_cube_file(rho, &grid, str);
      //Spin Down Wavefunction
      if (1 == flag->useSpinors){    
          sprintf(str,"homo-%ld-Dn.cube", i);
          for (jgrid = 0; jgrid < ist->ngrid; jgrid++){
          jgrid_real = ist->complex_idx * jgrid;
          jgrid_imag = ist->complex_idx * jgrid + 1;
          
          rho[jgrid] = sqr(psitot[ist->complex_idx*((ist->homo_idx-i)*ist->nspinngrid+ist->ngrid)+jgrid_real]) 
              + sqr(psitot[ist->complex_idx*((ist->homo_idx-i)*ist->nspinngrid+ist->ngrid)+jgrid_imag]);    
          }
          write_cube_file(rho, &grid, str);
      } 
      }

      for (i = 0;  (i < ist->total_lumo) && (i < ist->ncubes); i++){
      sprintf(str,"lumo+%ld-Up.cube",i);
      for (jgrid = 0; jgrid < ist->ngrid; jgrid++){
          jgrid_real = ist->complex_idx * jgrid;
          jgrid_imag = ist->complex_idx * jgrid + 1;
          
          rho[jgrid] = sqr(psitot[ist->complex_idx*(ist->lumo_idx+i)*ist->nspinngrid + jgrid_real]);
          if (1 == flag->isComplex) rho[jgrid] += sqr(psitot[ist->complex_idx*(ist->lumo_idx+i)*ist->nspinngrid + jgrid_imag]);
      }
      write_cube_file(rho, &grid, str);

      if (1 == flag->useSpinors){
          sprintf(str,"lumo+%ld-Dn.cube",i);
          for (jgrid = 0; jgrid < ist->ngrid; jgrid++){
          jgrid_real = ist->complex_idx * jgrid;
          jgrid_imag = ist->complex_idx * jgrid + 1;
          
          rho[jgrid] = sqr(psitot[ist->complex_idx*((ist->lumo_idx+i)*ist->nspinngrid+ist->ngrid)+jgrid_real]) 
              + sqr(psitot[ist->complex_idx*((ist->lumo_idx+i)*ist->nspinngrid+ist->ngrid)+jgrid_imag]);
          }
          write_cube_file(rho, &grid, str);
      }
      }
      free(rho);

      if (mpir == 0) printf("\ndone calculating cubes, CPU time (sec) %g, wall run time (sec) %g\n",
      ((double)clock()-inital_clock_t)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-initial_wall_t);
      }
  }
  if (flag->calcPotOverlap == 1){

      if (mpir == 0) write_separation(stdout, top);
      current_time = time(NULL);
      c_time_string = ctime(&current_time);
      if (mpir == 0) printf("\nCALCULATING POTENTIAL MATRIX ELEMENTS | %s\n", get_time()); fflush(0);
      if (mpir == 0) write_separation(stdout, bottom); fflush(stdout);

      calc_pot_overlap(&psitot[0], pot_local, nlc, nl, eig_vals, &par, &ist, &flag);
  } 
  if (flag->calcSpinAngStat == 1) {
      if (mpir == 0) write_separation(stdout, top);
      current_time = time(NULL);
      c_time_string = ctime(&current_time);
      if (mpir == 0) printf("\nCALCULATING SPIN & ANG. MOM. STATISTICS | %s\n", get_time()); 
      if (mpir == 0) write_separation(stdout, bottom); fflush(stdout);

      calc_angular_exp(psitot, &grid,0, ist->mn_states_tot, &ist, &par, &flag, &parallel, planfw, planbw, fftwpsi);
  
  } 
  if ( (flag->printCubes != 1) && (flag->calcPotOverlap != 1) && (flag->calcSpinAngStat != 1) ) {
      if (mpir == 0) printf("\nNo optional output requested.\n");
  }
  
  return;
}