#include "mod_optional_output.h"

// Cap on the supercell expansion along any one axis when drawing Bloch cubes.
// The Bloch phase period along an axis is 1/frac unit cells (frac = fractional
// k coordinate); near Gamma this diverges, so we cap it here. A capped axis does
// not tile perfectly but still shows the modulation.
#define MAX_BLOCH_CELLS 16

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

  if (mpir == 0)
  {
    printf("total_homo = %ld total_lumo = %ld\n", ist->total_homo, ist->total_lumo);
    fflush(0);
  }

  if (1 == flag->printCubes)
  {
    print_eigstate_densities(psitot, grid, ist, par, flag, parallel);
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

  ist->total_homo = ist->homo_idx;
  ist->total_lumo = ist->mn_states_tot - ist->total_homo;

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

/************************************************************/
// Smallest supercell factor N in [1, MAX_BLOCH_CELLS] for which the Bloch phase
// N*frac is (near) integer, so e^{ik.r} tiles the supercell. frac ~ 0 (no
// dispersion along this axis) returns 1.

int bloch_supercell_factor(double frac)
{
  frac = fabs(frac);
  if (frac < 1.0e-9)
    return 1;
  for (int N = 1; N <= MAX_BLOCH_CELLS; N++)
  {
    double v = N * frac;
    if (fabs(v - round(v)) < 1.0e-6)
      return N;
  }
  return MAX_BLOCH_CELLS;
}

/************************************************************/
// Build Re(e^{ik.r} u_nk) for one eigenstate on an N1xN2xN3 supercell into rho,
// then write it to k<ik>-<tag><n>.cube. u_nk is the cell-periodic part stored in
// psitot_k (complex, stride stlen); replicating it across cells and multiplying
// by the Bloch phase reconstructs the full wavefunction over the supercell.

static void emit_bloch_cube(
    double *rho, double *psitot_k, long is, long stlen, vector k,
    grid_st *grid, index_st *ist, flag_st *flag,
    int N1, int N2, int N3, int ik, const char *tag, int n)
{
  const double Lx = grid->nx * grid->dx, Ly = grid->ny * grid->dy, Lz = grid->nz * grid->dz;
  const long NXs = (long)N1 * grid->nx, NYs = (long)N2 * grid->ny, NZs = (long)N3 * grid->nz;
  double *base = &psitot_k[is * stlen];

  for (long IZ = 0; IZ < NZs; IZ++)
  {
    long jz = IZ % grid->nz;
    double Z = grid->z[jz] + (IZ / grid->nz) * Lz;
    for (long IY = 0; IY < NYs; IY++)
    {
      long jy = IY % grid->ny;
      double Y = grid->y[jy] + (IY / grid->ny) * Ly;
      for (long IX = 0; IX < NXs; IX++)
      {
        long jx = IX % grid->nx;
        double X = grid->x[jx] + (IX / grid->nx) * Lx;
        long jg = jx + grid->nx * (jy + grid->ny * jz);
        double ure = base[ist->complex_idx * jg];
        double uim = (1 == flag->isComplex) ? base[ist->complex_idx * jg + 1] : 0.0;
        double ph = k.x * X + k.y * Y + k.z * Z;
        rho[IX + NXs * (IY + NYs * IZ)] = ure * cos(ph) - uim * sin(ph);
      }
    }
  }

  char name[64];
  sprintf(name, "k%d-%s%d.cube", ik, tag, n);
  write_cube_super(rho, grid, N1, N2, N3, name);
}

/************************************************************/
// Write .cube files of the frontier orbitals at one k-point, visualizing the
// full Bloch function psi_nk(r) = e^{ik.r} u_nk(r) (its real part) on a supercell
// expanded so the waveform fits. Only states with sigma_E < par->sigma_E_cut are
// drawn; the manifolds are built outward from par->fermi_E (homo down, lumo up),
// ist->ncubes per side. psitot_k holds this k's diagonalized eigenvectors
// (complex, stride complex_idx*nspinngrid); eig_vals/sigma_E are length n_states.

void print_bloch_cubes_k(
    double *psitot_k, double *eig_vals, double *sigma_E, long n_states,
    int ik, vector k, grid_st *grid, index_st *ist, par_st *par, flag_st *flag)
{
  const long stlen = (long)ist->complex_idx * ist->nspinngrid;
  const double cut = par->sigma_E_cut;

  // Frontier indices among converged states (eig_vals ascending from zheev).
  long homo = -1, lumo = -1;
  for (long i = 0; i < n_states; i++)
    if ((eig_vals[i] < par->fermi_E) && (sigma_E[i] < cut))
      homo = i;
  for (long i = 0; i < n_states; i++)
    if ((eig_vals[i] >= par->fermi_E) && (sigma_E[i] < cut))
    {
      lumo = i;
      break;
    }
  if (homo < 0 && lumo < 0)
  {
    printf("ik=%d: no converged frontier states (sigma < %g) to draw cubes\n", ik, cut);
    return;
  }

  // Per-axis supercell expansion from the fractional k coordinate frac = k.L/2pi.
  double Lx = grid->nx * grid->dx, Ly = grid->ny * grid->dy, Lz = grid->nz * grid->dz;
  int N1 = bloch_supercell_factor(k.x * Lx / TWOPI);
  int N2 = bloch_supercell_factor(k.y * Ly / TWOPI);
  int N3 = bloch_supercell_factor(k.z * Lz / TWOPI);
  long NXs = (long)N1 * grid->nx, NYs = (long)N2 * grid->ny, NZs = (long)N3 * grid->nz;

  printf("ik=%d: |k|=%.5f, Bloch supercell %dx%dx%d (homo=%ld lumo=%ld)\n",
         ik, k.mag, N1, N2, N3, homo, lumo);
  fflush(stdout);

  double *rho;
  ALLOCATE(&rho, NXs * NYs * NZs, "rho in print_bloch_cubes_k");

  // Hole manifold: converged states from homo downward.
  int nc = 0;
  for (long i = homo; (i >= 0) && (nc < ist->ncubes); i--)
    if (sigma_E[i] < cut)
      emit_bloch_cube(rho, psitot_k, i, stlen, k, grid, ist, flag, N1, N2, N3, ik, "homo-", nc++);

  // Electron manifold: converged states from lumo upward.
  nc = 0;
  for (long i = lumo; (i >= 0) && (i < n_states) && (nc < ist->ncubes); i++)
    if (sigma_E[i] < cut)
      emit_bloch_cube(rho, psitot_k, i, stlen, k, grid, ist, flag, N1, N2, N3, ik, "lumo+", nc++);

  free(rho);
  return;
}