/****************************************************************************/
/* This file does the printing for the program */

#include "fd.h"

/****************************************************************************/

void write_cube_file(double *rho, grid_st *grid, char *fileName)
{
  /*****************************************************************
   * This function prints out cube files to visualize wavefunctions *
   * inputs: [double *rho] is a pointer to an ngrid long array      *
   *         [grid_st *grid] is a pointer to the grid struct        *
   *         [char *fileName] is a pointer to the output file name  *
   * outputs: void                                                  *
   ******************************************************************/

  FILE *pf, *pConfFile;
  long jgrid, iX, iY, iZ, iYZ, natoms, atomType;
  double x, y, z;
  char line[80], atomSymbol[10];

  pConfFile = fopen("conf.par", "r");
  fscanf(pConfFile, "%ld", &natoms);
  pf = fopen(fileName, "w");
  fprintf(pf, "CUBE FILE\n");
  fprintf(pf, "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n");
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", natoms, grid->xmin, grid->ymin, grid->zmin);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", grid->nx, grid->dx, 0.0, 0.0);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", grid->ny, 0.0, grid->dy, 0.0);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", grid->nz, 0.0, 0.0, grid->dz);
  fgets(line, 80, pConfFile);
  while (fgets(line, 80, pConfFile) != NULL)
  {
    sscanf(line, "%2s %lf %lf %lf %ld %*lf", (char *)&atomSymbol, &x, &y, &z, &atomType);

    // TODO: make sane atom handling....
    if (!strcmp(atomSymbol, "Cd"))
    {
      atomType = 48;
    }
    else if (!strcmp(atomSymbol, "S"))
    {
      atomType = 16;
    }
    else if (!strcmp(atomSymbol, "Se"))
    {
      atomType = 34;
    }
    else if (!strcmp(atomSymbol, "Zn"))
    {
      atomType = 30;
    }
    else if (!strcmp(atomSymbol, "Te"))
    {
      atomType = 52;
    }
    else if (!strcmp(atomSymbol, "C"))
    {
      atomType = 6;
    }
    else if (!strcmp(atomSymbol, "Si"))
    {
      atomType = 14;
    }
    else if (!strcmp(atomSymbol, "Cs"))
    {
      atomType = 55;
    }
    else if (!strcmp(atomSymbol, "Pb"))
    {
      atomType = 82;
    }
    else if (!strcmp(atomSymbol, "I"))
    {
      atomType = 53;
    }
    else
    {
      atomType = 1;
    }
    fprintf(pf, "%5li%12.6f%12.6f%12.6f%12.6f\n", atomType, 0.0, x, y, z);
  }
  for (iX = 0; iX < grid->nx; iX++)
  {
    for (iY = 0; iY < grid->ny; iY++)
    {
      for (iZ = 0; iZ < grid->nz; iZ++)
      {
        iYZ = grid->nx * (grid->ny * iZ + iY);
        jgrid = iYZ + iX;
        fprintf(pf, "%12.5g ", rho[jgrid]);
        if (iZ % 6 == 5)
        {
          fprintf(pf, "\n");
        }
      }
      fprintf(pf, "\n");
    }
  }
  fclose(pConfFile);
  fclose(pf);

  return;
}

/****************************************************************************/
// prints out current time to stdout

void write_current_time(FILE *pf)
{
  time_t startTime;

  startTime = time(NULL);
  fprintf(pf, "%s", ctime(&startTime));

  return;
}

/****************************************************************************/

void write_separation(FILE *pf, char *top_bttm)
{
  /*****************************************************************
   * This function prints asterisk separation lines in stdout       *
   * inputs:                                                        *
   *   [FILE *pf] pointer to output file stream                     *
   *   [char *top_bttm] pointer to char for top/bottom formatting   *
   * outputs: void                                                  *
   ******************************************************************/

  const char top_key[] = "T";
  const char bttm_key[] = "B";

  if (strcmp(top_bttm, top_key) == 0)
  {
    fprintf(pf, "\n\n******************************************************************************\n");
  }
  else if (strcmp(top_bttm, bttm_key) == 0)
  {
    fprintf(pf, "\n******************************************************************************\n");
  }
  else
  {
    fprintf(stderr, "Invalid string supplied to write_separation. Exiting!\n");
    exit(EXIT_FAILURE);
  }

  return;
}

/****************************************************************************/

void write_state_dat(double complex *buf, long n_elems, char *fileName)
{
  FILE *pf;
  long i;

  pf = fopen(fileName, "w");
  for (i = 0; i < n_elems; i++)
  {
    fprintf(pf, "%lg %lg\n", creal(buf[i]), cimag(buf[i]));
  }
  fclose(pf);

  return;
}

/****************************************************************************/

void write_BSE_coeff(double complex *bs_coeff, long n_xton, char *fileName)
{
  /*****************************************************************
   * Writes the full set of BSE exciton coefficients (eigenvectors) *
   * to a single complex file. Layout matches the downstream reader *
   * (filter/carrier_dens.c): row i is the two-particle basis       *
   * component, column j is the exciton, and bs_coeff[i*n_xton + j]  *
   * is component i of exciton j. Each entry is written as a        *
   * "re, im" pair so the reader can parse it with                  *
   * fscanf("%lg, %lg").                                            *
   * inputs: [double complex *bs_coeff] n_xton x n_xton coeff matrix*
   *         [long n_xton] number of two-particle basis states      *
   *         [char *fileName] output file name                      *
   * outputs: void                                                  *
   ******************************************************************/

  FILE *pf;
  long i, j;

  pf = fopen(fileName, "w");
  if (pf == NULL)
  {
    fprintf(stderr, "ERROR: could not open %s for writing "
                    "(is the working directory writable?)\n", fileName);
    return;
  }

  for (i = 0; i < n_xton; i++, fprintf(pf, "\n"))
  {
    for (j = 0; j < n_xton; j++)
    {
      fprintf(pf, "%.*g, %.*g\t", PR_LEN, creal(bs_coeff[i * n_xton + j]),
              PR_LEN, cimag(bs_coeff[i * n_xton + j]));
    }
  }
  fclose(pf);

  return;
}

/****************************************************************************/
// Phase sign for cube-file densities (identical convention to
// filter_mpi/makecube). Returns the sign (+1, -1, or 0) of whichever of the
// real/imaginary components has the larger magnitude. For a real-valued field
// (im == 0) this recovers sign(re). The exciton carrier densities below are
// mathematically real and non-negative, so this always returns +1 for them; it
// is applied to stay bit-faithful to the makecube rendering method.

double density_phase_sign(double re, double im)
{
  double dom = (fabs(re) >= fabs(im)) ? re : im;
  if (dom > 0.0)
    return 1.0;
  if (dom < 0.0)
    return -1.0;
  return 0.0;
}

/****************************************************************************/
// Helper: accumulate one fixed-point carrier density for exciton n onto rho.
//
//   mode = 0 (ELECTRON carrier, integrate out the hole):
//     rho(r_e) = sum_i sum_s | sum_a C^n_{ia} phi_a^s(r_e) |^2
//   mode = 1 (HOLE carrier, integrate out the electron):
//     rho(r_h) = sum_a sum_s | sum_i C^n_{ia} phi_i^s(r_h) |^2
//
// s runs over the two spinor blocks (up = [0,ngrid), dn = [ngrid,2*ngrid)),
// mirroring coulomb.c which always sums both blocks. C^n_{ia} lives at
// bs_coeff[ibs*n_xton + n] with ibs = (a - lumo_idx)*n_holes + i; hole state i
// is psi_qp state i and electron state a is psi_qp state a (stride nspinngrid).
// The result is real and >= 0 by construction (it factors as a sum of |.|^2).

static void accumulate_carrier_density(
    double *rho, double complex *bs_coeff, double complex *psi_qp,
    long n, int mode, index_st *ist)
{
  const long ngrid = ist->ngrid;
  const long nspngr = ist->nspinngrid;
  const long n_xton = ist->n_xton;
  const long lidx = ist->lumo_idx;
  const long n_el = ist->n_elecs;
  const long n_ho = ist->n_holes;

#pragma omp parallel for schedule(static)
  for (long g = 0; g < ngrid; g++)
  {
    double acc = 0.0;
    if (mode == 0) /* electron carrier: sum over holes, contract electrons */
    {
      for (long i = 0; i < n_ho; i++)
      {
        for (long s = 0; s < 2; s++)
        {
          const long off = s * ngrid + g;
          double complex chi = 0.0 + 0.0 * I;
          for (long a = lidx; a < lidx + n_el; a++)
          {
            const long ibs = (a - lidx) * n_ho + i;
            chi += bs_coeff[ibs * n_xton + n] * psi_qp[a * nspngr + off];
          }
          acc += creal(chi) * creal(chi) + cimag(chi) * cimag(chi);
        }
      }
    }
    else /* hole carrier: sum over electrons, contract holes */
    {
      for (long a = lidx; a < lidx + n_el; a++)
      {
        for (long s = 0; s < 2; s++)
        {
          const long off = s * ngrid + g;
          double complex chi = 0.0 + 0.0 * I;
          for (long i = 0; i < n_ho; i++)
          {
            const long ibs = (a - lidx) * n_ho + i;
            chi += bs_coeff[ibs * n_xton + n] * psi_qp[i * nspngr + off];
          }
          acc += creal(chi) * creal(chi) + cimag(chi) * cimag(chi);
        }
      }
    }
    /* Signed density, rendered by the makecube method. acc is real and >= 0,
     * so this equals acc; the transform keeps parity with filter_mpi/makecube
     * (and would sign any residual, rather than squaring it as |rho|^2). */
    rho[g] = sqrt(sqr(acc)) * density_phase_sign(acc, 0.0);
  }
}

/****************************************************************************/
// Writes fixed-point exciton carrier densities as Gaussian cube files, for the
// lowest ist->n_FP_density excitons. For each exciton n it writes both
// rhoElec<n>.cube (electron carrier, hole integrated out) and rhoHole<n>.cube
// (hole carrier, electron integrated out). The signed density is rendered with
// the same method as filter_mpi/makecube (see accumulate_carrier_density).
// Rank-0 post-processing: bs_coeff (full) is on rank 0 and psi_qp is the
// node-shared copy, both readable here.

void write_exciton_densities(
    double complex *bs_coeff, double complex *psi_qp,
    grid_st *grid, index_st *ist, par_st *par,
    flag_st *flag, parallel_st *parallel)
{
  long n_FP = ist->n_FP_density;
  if (n_FP <= 0)
    return;

  /* write_cube_file requires conf.par in the run directory (same config file the
   * rest of bse_cplx reads); check up front so a missing file is a clear message
   * rather than a crash inside the cube writer. */
  if (access("conf.par", F_OK) != 0)
  {
    fprintf(stderr, "WARNING: nExcitonDensities > 0 but conf.par is missing; "
                    "skipping fixed-point density cubes.\n");
    return;
  }

  if (n_FP > ist->n_xton)
  {
    printf("\tnExcitonDensities (%ld) > n_xton (%ld); clamping to %ld\n",
           n_FP, ist->n_xton, ist->n_xton);
    n_FP = ist->n_xton;
  }

  double *rho;
  ALLOCATE(&rho, ist->ngrid, "rho in write_exciton_densities");

  char fileName[64];
  for (long n = 0; n < n_FP; n++)
  {
    accumulate_carrier_density(rho, bs_coeff, psi_qp, n, 0, ist);
    snprintf(fileName, sizeof(fileName), "rhoElec%ld.cube", n);
    write_cube_file(rho, grid, fileName);

    accumulate_carrier_density(rho, bs_coeff, psi_qp, n, 1, ist);
    snprintf(fileName, sizeof(fileName), "rhoHole%ld.cube", n);
    write_cube_file(rho, grid, fileName);

    printf("\tWrote fixed-point densities for exciton %ld (rhoElec%ld.cube, rhoHole%ld.cube) | %s\n",
           n, n, n, get_time());
    fflush(stdout);
  }

  free(rho);

  return;
}

/****************************************************************************/
