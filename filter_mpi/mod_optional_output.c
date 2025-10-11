#include "mod_optional_output.h"

void mod_optional_output(
    double *psitot,
    grid_st *grid,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel)
{

  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  const int mpir = parallel->mpi_rank;

  if (mpir == 0)
  {
    write_separation(stdout, "T");
    printf("\nCALCULATING OPTIONAL OUTPUT | %s\n", get_time());
    write_separation(stdout, "B");
    fflush(stdout);
  }

  /************************************************************/
  /*******************   PRINT CUBE FILES   *******************/
  /************************************************************/
  write_separation(stdout, "T");
  printf("\nWRITING CUBE FILES\n");
  write_separation(stdout, "B");
  fflush(stdout);

  if (1 == flag->printCubes)
  {
    print_eigstate_densities(psitot, grid, ist, par, flag, parallel);
  }

  // Set the total number of electron and hole states in order to calculate the potential overlap integrals
  ist->total_homo = ist->homo_idx + 1;
  ist->total_lumo = ist->mn_states_tot - ist->total_homo;

  if (mpir == 0)
  {
    printf("total_homo = %ld total_lumo = %ld\n", ist->total_homo, ist->total_lumo);
    fflush(0);
  }

  if (flag->calcSpinAngStat == 1)
  {
    if (mpir == 0)
      write_separation(stdout, "T");
    if (mpir == 0)
      printf("\nCALCULATING SPIN & ANG. MOM. STATISTICS | %s\n", get_time());
    if (mpir == 0)
      write_separation(stdout, "B");
    fflush(stdout);

    calc_angular_exp(psitot, grid, ist->homo_idx - 4, ist->lumo_idx + 12, ist, par, flag, parallel);
  }
  if ((flag->printCubes != 1) && (flag->calcSpinAngStat != 1))
  {
    if (mpir == 0)
      printf("\nNo optional output requested.\n");
  }

  return;
}

/************************************************************/

void print_eigstate_densities(double *psitot, grid_st *grid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel)
{

  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  long i;
  long a;
  long ieof;
  long jgrid;
  long jgrid_real;
  long jgrid_imag;

  char str[50];
  double evalloc;
  double deloc;
  double *rho;

  FILE *pf;

  ist->homo_idx = ist->lumo_idx = 0;

  const int mpir = parallel->mpi_rank;

  /************************************************************/
  /*******************  FIND HOMO/LUMO IDX   ******************/
  /************************************************************/

  // Find homo_idx
  pf = fopen("eval.dat", "r");
  for (i = 0; i < ist->mn_states_tot; i++)
  {
    ieof = fscanf(pf, "%ld %lg %lg", &a, &evalloc, &deloc);
    if (ieof == EOF)
    {
      break;
    }
    if ((evalloc < par->fermi_E) && (deloc < par->sigma_E_cut))
    {
      ist->homo_idx = i;
    }

    if (i > ist->mn_states_tot)
    {
      printf("No hole states converged to within %lg a.u.\n", par->sigma_E_cut);
      break;
    }
  }
  fclose(pf);

  // Find lumo_idx
  pf = fopen("eval.dat", "r");
  for (i = 0; i <= ist->homo_idx; i++)
  {
    fscanf(pf, "%ld %lg %lg", &a, &evalloc, &deloc);
  }
  for (i = ist->homo_idx + 1, ieof = 0; ieof != EOF; i++)
  {
    fscanf(pf, "%ld %lg %lg", &a, &evalloc, &deloc);

    if (deloc < par->sigma_E_cut)
    {
      ist->lumo_idx = i;
      break;
    }

    if (i > ist->mn_states_tot)
    {
      printf("No electron states converged to within %lg a.u.\n", par->sigma_E_cut);
      break;
    }
  }
  fclose(pf);

  printf("homo_idx = %ld; lumo_idx = %ld\n", ist->homo_idx, ist->lumo_idx);
  fflush(0);

  /*** Write homo and lumo cube files ***/

  if ((ist->homo_idx == 0) || (ist->lumo_idx == 0))
  {
    printf("\nDid not converge enough electron or hole states to visualize cube files.\n");
    return;
  }
  else
  {
    ALLOCATE(&rho, ist->ngrid, "rho in print_eigstate_densities");
  }

  for (i = 0; (i < ist->total_homo) && (i < ist->ncubes); i++)
  {
    // Spin Up Wavefunction
    sprintf(str, "homo-%ld-Up.cube", i);
    for (jgrid = 0; jgrid < ist->ngrid; jgrid++)
    {
      jgrid_real = ist->complex_idx * jgrid;
      jgrid_imag = ist->complex_idx * jgrid + 1;

      rho[jgrid] = sqr(psitot[ist->complex_idx * (ist->homo_idx - i) * ist->nspinngrid + jgrid_real]);
      if (1 == flag->isComplex)
        rho[jgrid] += sqr(psitot[ist->complex_idx * (ist->homo_idx - i) * ist->nspinngrid + jgrid_imag]);
    }
    write_cube_file(rho, grid, str);
    // Spin Down Wavefunction
    if (1 == flag->useSpinors)
    {
      sprintf(str, "homo-%ld-Dn.cube", i);
      for (jgrid = 0; jgrid < ist->ngrid; jgrid++)
      {
        jgrid_real = ist->complex_idx * jgrid;
        jgrid_imag = ist->complex_idx * jgrid + 1;

        rho[jgrid] = sqr(psitot[ist->complex_idx * ((ist->homo_idx - i) * ist->nspinngrid + ist->ngrid) + jgrid_real]) + sqr(psitot[ist->complex_idx * ((ist->homo_idx - i) * ist->nspinngrid + ist->ngrid) + jgrid_imag]);
      }
      write_cube_file(rho, grid, str);
    }
  }

  for (i = 0; (i < ist->total_lumo) && (i < ist->ncubes); i++)
  {
    sprintf(str, "lumo+%ld-Up.cube", i);
    for (jgrid = 0; jgrid < ist->ngrid; jgrid++)
    {
      jgrid_real = ist->complex_idx * jgrid;
      jgrid_imag = ist->complex_idx * jgrid + 1;

      rho[jgrid] = sqr(psitot[ist->complex_idx * (ist->lumo_idx + i) * ist->nspinngrid + jgrid_real]);
      if (1 == flag->isComplex)
        rho[jgrid] += sqr(psitot[ist->complex_idx * (ist->lumo_idx + i) * ist->nspinngrid + jgrid_imag]);
    }
    write_cube_file(rho, grid, str);

    if (1 == flag->useSpinors)
    {
      sprintf(str, "lumo+%ld-Dn.cube", i);
      for (jgrid = 0; jgrid < ist->ngrid; jgrid++)
      {
        jgrid_real = ist->complex_idx * jgrid;
        jgrid_imag = ist->complex_idx * jgrid + 1;

        rho[jgrid] = sqr(psitot[ist->complex_idx * ((ist->lumo_idx + i) * ist->nspinngrid + ist->ngrid) + jgrid_real]) + sqr(psitot[ist->complex_idx * ((ist->lumo_idx + i) * ist->nspinngrid + ist->ngrid) + jgrid_imag]);
      }
      write_cube_file(rho, grid, str);
    }
  }
  free(rho);

  if (mpir == 0)
  {
    printf("\ndone calculating cubes, %s\n", get_time());
  }

  return;
}