#include "fd.h"

/****************************************************************************/
// Map a conf.dat atom label to the atomic number written in a .cube header.
// (Extracted from write_cube_file so the supercell writer can share it.)

long cube_atom_number(char *atomSymbol)
{
  if (!strcmp(atomSymbol, "C1"))
    return 2;
  else if (!strcmp(atomSymbol, "C2"))
    return 3;
  else if (!strcmp(atomSymbol, "C3"))
    return 4;
  else if (!strcmp(atomSymbol, "P1"))
    return 84;
  else if (!strcmp(atomSymbol, "P2"))
    return 85;
  else if (!strcmp(atomSymbol, "P3"))
    return 86;
  else if (!strcmp(atomSymbol, "PC5"))
    return 87;
  else if (!strcmp(atomSymbol, "PC6"))
    return 88;
  else if (!strcmp(atomSymbol, "PA1"))
    return 3;
  else if (!strcmp(atomSymbol, "PA2"))
    return 4;
  else if (!strcmp(atomSymbol, "I0"))
    return 53;
  else if (!strcmp(atomSymbol, "I1"))
    return 53;
  else if (!strcmp(atomSymbol, "Br0"))
    return 35;
  else if (!strcmp(atomSymbol, "Br1"))
    return 35;
  else if ((!strcmp(atomSymbol, "PA3")) || (!strcmp(atomSymbol, "PR1")) ||
           (!strcmp(atomSymbol, "PR2")) || (!strcmp(atomSymbol, "PR3")))
    return 1;
  else
    return assign_atom_number(atomSymbol);
}

/****************************************************************************/
//

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

  if (access("conf.dat", F_OK) == -1)
  {
    printf("ERROR: no conf.dat file exists in directory\n");
    fprintf(stderr, "ERROR: no conf.dat file exists in directory\n");
    exit(EXIT_FAILURE);
  }
  else
  {
    pConfFile = fopen("conf.dat", "r");
  }
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
    sscanf(line, "%3s %lf %lf %lf", (char *)&atomSymbol, &x, &y, &z);

    atomType = cube_atom_number(atomSymbol);
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
        fprintf(pf, "%.5g ", rho[jgrid]);
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
// Write a .cube for an N1 x N2 x N3 supercell of the base grid. rho must hold
// (N1*nx) x (N2*ny) x (N3*nz) values in cube order (x fastest, z slowest):
//   SG = IX + (N1*nx) * (IY + (N2*ny) * IZ).
// Atoms from conf.dat are replicated across all cells, each copy shifted by the
// integer cell offset times the box length (nx*dx, ny*dy, nz*dz). Used to draw
// the full Bloch waveform e^{ik.r} u_nk, whose period exceeds one unit cell.

void write_cube_super(double *rho, grid_st *grid, int N1, int N2, int N3, char *fileName)
{
  FILE *pf, *pConfFile;
  long natoms, ia, m1, m2, m3, iX, iY, iZ;
  long NXs = (long)N1 * grid->nx, NYs = (long)N2 * grid->ny, NZs = (long)N3 * grid->nz;
  double Lx = grid->nx * grid->dx, Ly = grid->ny * grid->dy, Lz = grid->nz * grid->dz;
  char line[80], atomSymbol[10];

  if (access("conf.dat", F_OK) == -1)
  {
    fprintf(stderr, "ERROR: no conf.dat file exists in directory\n");
    exit(EXIT_FAILURE);
  }
  pConfFile = fopen("conf.dat", "r");
  fscanf(pConfFile, "%ld", &natoms);
  fgets(line, 80, pConfFile); // consume rest of the count line

  // Read the base-cell atoms (symbol -> Z, position) once.
  long *atom_Z = (long *)malloc(natoms * sizeof(long));
  double *ax = (double *)malloc(natoms * sizeof(double));
  double *ay = (double *)malloc(natoms * sizeof(double));
  double *az = (double *)malloc(natoms * sizeof(double));
  for (ia = 0; ia < natoms && fgets(line, 80, pConfFile) != NULL; ia++)
  {
    double x, y, z;
    sscanf(line, "%3s %lf %lf %lf", (char *)&atomSymbol, &x, &y, &z);
    atom_Z[ia] = cube_atom_number(atomSymbol);
    ax[ia] = x; ay[ia] = y; az[ia] = z;
  }
  fclose(pConfFile);

  pf = fopen(fileName, "w");
  fprintf(pf, "CUBE FILE (Bloch supercell %dx%dx%d)\n", N1, N2, N3);
  fprintf(pf, "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n");
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", natoms * N1 * N2 * N3,
          grid->xmin, grid->ymin, grid->zmin);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", NXs, grid->dx, 0.0, 0.0);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", NYs, 0.0, grid->dy, 0.0);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", NZs, 0.0, 0.0, grid->dz);

  // Replicate atoms across all cells (cell-major order is irrelevant to viewers).
  for (m1 = 0; m1 < N1; m1++)
    for (m2 = 0; m2 < N2; m2++)
      for (m3 = 0; m3 < N3; m3++)
        for (ia = 0; ia < natoms; ia++)
          fprintf(pf, "%5li%12.6f%12.6f%12.6f%12.6f\n", atom_Z[ia], 0.0,
                  ax[ia] + m1 * Lx, ay[ia] + m2 * Ly, az[ia] + m3 * Lz);

  // Volumetric data, same ordering convention as write_cube_file.
  for (iX = 0; iX < NXs; iX++)
    for (iY = 0; iY < NYs; iY++)
    {
      for (iZ = 0; iZ < NZs; iZ++)
      {
        fprintf(pf, "%.5g ", rho[iX + NXs * (iY + NYs * iZ)]);
        if (iZ % 6 == 5)
          fprintf(pf, "\n");
      }
      fprintf(pf, "\n");
    }

  fclose(pf);
  free(atom_Z); free(ax); free(ay); free(az);
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

  char *top_key;
  top_key = malloc(2 * sizeof(top_key[0]));
  char *bttm_key;
  bttm_key = malloc(2 * sizeof(bttm_key[0]));

  strcpy(top_key, "T");
  strcpy(bttm_key, "B");

  if (0 == strcmp(top_bttm, (const char *)top_key))
  {
    fprintf(pf, "\n\n******************************************************************************\n");
  }
  else if (0 == strcmp(top_bttm, (const char *)bttm_key))
  {
    fprintf(pf, "\n******************************************************************************\n");
  }
  else
  {
    fprintf(stderr, "Invalid string supplied to writeSeparation. Exiting!\n");
    exit(EXIT_FAILURE);
  }

  free(top_key);
  free(bttm_key);

  return;
}

/****************************************************************************/

void write_state_dat(zomplex *psi, long n_elems, char *fileName)
{
  FILE *pf;
  long i;

  pf = fopen(fileName, "w");
  for (i = 0; i < n_elems; i++)
  {
    fprintf(pf, "%lg %lg\n", psi[i].re, psi[i].im);
  }
  fclose(pf);
}

/****************************************************************************/

void write_eval_dat(double *eig_vals, double *sigma_E, long n_elems, char *fileName)
{

  FILE *pf;
  long i;

  pf = fopen(fileName, "w");

  if (pf == NULL)
  {
    fprintf(stderr, "ERROR: disk full, could not open %s\n", fileName);
    exit(EXIT_FAILURE);
  }

  for (i = 0; i < n_elems; i++)
  {
    fprintf(pf, "%ld %.16lg %lg\n", i, eig_vals[i], sigma_E[i]);
  }

  fclose(pf);

  return;
}

/****************************************************************************/

void write_psi_dat(double *psitot, long n_elems, char *fileName)
{

  FILE *pf;

  pf = fopen(fileName, "w");

  fwrite(psitot, sizeof(psitot[0]), n_elems, pf);

  return;
}
/****************************************************************************/

void write_vector_dat(vector *vec, int n_elems, char *fileName)
{

  FILE *pf;
  int i;

  pf = fopen(fileName, "w");

  fprintf(pf, "units of 2pi/a\n");
  for (i = 0; i < n_elems; i++)
  {
    fprintf(pf, "%lf %lf %lf\n", vec[i].x, vec[i].y, vec[i].z);
  }
  fclose(pf);

  return;
}
