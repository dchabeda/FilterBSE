/*****************************************************************************/
// Main file for cube printing utility.
#include "fd.h"
#include "aux.h"
/*****************************************************************************/
int countlines(char *filename);

/*****************************************************************************/
int main(int argc, char *argv[])
{
  FILE *ppsi;
  // custom structs
  grid_st grid;
  par_st par;
  index_st ist;
  parallel_st parallel;
  flag_st flag;
  xyz_st *R;
  atom_info *atom;
  // double arrays
  double *rho;
  double *psi;
  // long int arrays and counters
  long i, jms;
  int j, start, end;
  int jgrid_real, jgrid_imag;
  ist.atom_types = malloc(N_MAX_ATOM_TYPES * sizeof(ist.atom_types[0]));
  time_t currentTime = time(NULL);

  parallel.mpi_rank = 0;
  parallel.mpi_size = 1;

  // command line input parsing
  if (argc != 4)
  {
    printf("\nUsage: ./makecube.x start end filename\n");
    exit(EXIT_FAILURE);
  }
  char psifile[256];

  start = atoi(argv[1]);
  end = atoi(argv[2]);
  // snprintf is bounds-checked (won't overflow psifile if argv[3] is long)
  snprintf(psifile, sizeof(psifile), "%s", argv[3]);

  if (start > end)
  {
    printf("Invaid start (%d), end(%d): start > end\n", start, end);
    exit(EXIT_FAILURE);
  }
  if (start < 0)
  {
    printf("Invaid start (%d): start < 0\n", start);
    exit(EXIT_FAILURE);
  }

  printf("This calculation began at: %s", ctime(&currentTime));
  fflush(stdout);

  /*** read initial setup from input.par ***/
  printf("\nReading job specifications from input.par:\n");
  fflush(0);
  read_input(&flag, &grid, &ist, &par, &parallel);

  /*** allocating memory ***/
  printf("Allocating memory\n");
  fflush(0);
  // the positions of the atoms in the x, y, and z directions
  if ((R = (xyz_st *)calloc(ist.natoms, sizeof(xyz_st))) == NULL)
    nerror("R");
  // the atom specific information
  if ((atom = (atom_info *)calloc(ist.natoms, sizeof(atom_info))) == NULL)
    nerror("atom");

  /*** read the nanocrystal configuration ***/
  printf("\nReading atomic configuration from conf.par:\n");
  fflush(0);
  read_conf(R, atom, &ist, &par, &flag, &parallel);

  /*** initialize the grid ***/
  printf("\nInitializing the grid parameters:\n");
  fflush(0);
  init_grid_params(&grid, R, &ist, &par, &flag, &parallel);

  if ((rho = (double *)calloc(ist.ngrid, sizeof(double))) == NULL)
    nerror("rho");

  // allocate memory for psi
  if ((psi = (double *)calloc(ist.complex_idx * ist.nspinngrid, sizeof(double))) == NULL)
    nerror("psi");

  // open the wavefunction file (name supplied on the command line)
  ppsi = fopen(psifile, "r");
  if (ppsi == NULL)
  {
    printf("Error: could not open wavefunction file '%s'\n", psifile);
    exit(EXIT_FAILURE);
  }

  // count the number of states stored in psifile from its byte size, rather than
  // assuming a fixed "eval.dat" companion (per-k output uses eval-<k>.dat/psi-<k>.dat)
  const long state_size = ist.complex_idx * ist.nspinngrid; // doubles per state
  fseek(ppsi, 0, SEEK_END);
  long file_bytes = ftell(ppsi);
  rewind(ppsi);
  jms = file_bytes / (long)(state_size * sizeof(double));
  printf("%ld total states in %s\n", jms, psifile);
  fflush(0);

  // make sure the requested range is actually present in the file
  if (end >= jms)
  {
    printf("Error: requested end state %d but %s holds only %ld states (valid 0..%ld)\n",
           end, psifile, jms, jms - 1);
    exit(EXIT_FAILURE);
  }

  char filename[20];
  for (j = start; j <= end; j++)
  {
    printf("Reading state %d from %s\n", j, psifile);

    if (fseek(ppsi, j * ist.complex_idx * ist.nspinngrid * sizeof(double), SEEK_SET) != 0)
    {
      printf("Error reading from %s!\n", psifile);
      exit(EXIT_FAILURE);
    }
    fread(&psi[0], sizeof(double), ist.complex_idx * ist.nspinngrid, ppsi);

    if (1 == flag.isComplex)
    {
      for (i = 0; i < ist.ngrid; i++)
      {
        jgrid_real = i * ist.complex_idx;
        jgrid_imag = jgrid_real + 1;
        double re = psi[jgrid_real], im = psi[jgrid_imag];
        // sqrt is not physical, just to make values larger for visualization.
        // Signed by the dominant component so the complex state keeps a phase.
        rho[i] = sqrt(sqr(re) + sqr(im)) * density_phase_sign(re, im);
      }
      sprintf(filename, "rhoUp%i.cube", j);
      write_cube_file(rho, &grid, filename);
    }
    else
    {
      for (i = 0; i < ist.ngrid; i++)
      {
        // sqrt is not physical, just to make values larger for visualization.
        // Real wavefunction: sign(re) recovers its +/- phase (= signed |psi|).
        rho[i] = sqrt(sqr(psi[i])) * density_phase_sign(psi[i], 0.0);
      }
      sprintf(filename, "rhoUp%i.cube", j);
      write_cube_file(rho, &grid, filename);
    }

    if (flag.useSpinors == 1)
    {
      if (1 == flag.isComplex)
      {
        for (i = 0; i < ist.ngrid; i++)
        {
          jgrid_real = i * ist.complex_idx;
          jgrid_imag = jgrid_real + 1;
          double re = psi[ist.ngrid * ist.complex_idx + jgrid_real];
          double im = psi[ist.ngrid * ist.complex_idx + jgrid_imag];
          // sqrt is not physical, just to make values larger for visualization.
          // Signed by the dominant component so the complex state keeps a phase.
          rho[i] = sqrt(sqr(re) + sqr(im)) * density_phase_sign(re, im);
        }
        sprintf(filename, "rhoDn%i.cube", j);
        write_cube_file(rho, &grid, filename);
      }
      else
      {
        for (i = 0; i < ist.ngrid; i++)
        {
          // sqrt is not physical, just to make values larger for visualization.
          // Real wavefunction: sign(re) recovers its +/- phase (= signed |psi|).
          double re = psi[ist.ngrid * ist.complex_idx + i];
          rho[i] = sqrt(sqr(re)) * density_phase_sign(re, 0.0);
        }
        sprintf(filename, "rhoDn%i.cube", j);
        write_cube_file(rho, &grid, filename);
      }
    }
  }
  fclose(ppsi);

  printf("Done with makecube.x\n");

  return 0;
}

/*****************************************************************************/
int countlines(char *filename)
{
  FILE *fp = fopen(filename, "r");
  int lines = 0;
  int ch;
  if (fp == NULL)
  {
    printf("Error: could not open '%s' to count lines\n", filename);
    return -1;
  }
  while (1)
  {
    ch = fgetc(fp);
    if (feof(fp))
    {
      break;
    }
    if (ch == '\n')
    {
      lines++;
    }
  }
  fclose(fp);
  return lines;
}
