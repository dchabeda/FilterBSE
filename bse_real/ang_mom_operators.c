#include "angular.h"

/************************************************************/
// Calculate the three vector components of the action of the L operator on the
// spatial part of the grid (no spin part)
void l_operator(
    double complex *Lxpsi,
    double complex *Lypsi,
    double complex *Lzpsi,
    double complex *psi_qp,
    double *g_vecs,
    grid_st *grid,
    index_st *ist,
    par_st *par,
    fftw_plan_loc planfw,
    fftw_plan_loc planbw,
    fftw_complex *fftwpsi)
{

  double x, y, z;
  long jx, jy, jz, jyz, jxyz;
  double complex px, py, pz;

  // Take derivative along x direction
  p_operator("X", g_vecs, psi_qp, Lxpsi, grid, ist, par, planfw, planbw, fftwpsi);
  // Take derivative along y direction
  p_operator("Y", g_vecs, psi_qp, Lypsi, grid, ist, par, planfw, planbw, fftwpsi);
  // Take derivative along x direction
  p_operator("Z", g_vecs, psi_qp, Lzpsi, grid, ist, par, planfw, planbw, fftwpsi);

  /*** Now do the cross product part at each grid point***/
  omp_set_dynamic(0);
  omp_set_num_threads(ist->nthreads);
#pragma omp parallel for private(jz, jy, jyz, jx, jxyz, px, py, pz, x, y, z)
  for (jz = 0; jz < grid->nz; jz++)
  {
    z = grid->z[jz];
    for (jy = 0; jy < grid->ny; jy++)
    {
      jyz = grid->nx * (grid->ny * jz + jy);
      y = grid->y[jy];
      for (jx = 0; jx < grid->nx; jx++)
      {
        x = grid->x[jx];
        jxyz = jyz + jx;
        // copy over tmp varibales
        px = Lxpsi[jxyz];
        py = Lypsi[jxyz];
        pz = Lzpsi[jxyz];

        Lxpsi[jxyz] = -(y * pz - z * py);

        Lypsi[jxyz] = -(z * px - x * pz);

        Lzpsi[jxyz] = -(x * py - y * px);
      }
    }
  }

  return;
}

/*************************************************************************************************/

void init_g_vecs(double *kindex, double *kx, double *ky, double *kz, grid_st *grid, index_st *ist, par_st *par)
{

  long jx, jy, jz, jyz, jgrid;

  // The k-vectors in each direction are initialized so that
  // they are compatible with FFTW3 definition
  init_kvec(kx, grid->nx, grid->dkx, grid->nx_1 * grid->ny_1 * grid->nz_1);
  init_kvec(ky, grid->ny, grid->dky, grid->nx_1 * grid->ny_1 * grid->nz_1);
  init_kvec(kz, grid->nz, grid->dkz, grid->nx_1 * grid->ny_1 * grid->nz_1);

  // for (jy = 0; jy < grid->ny; jy++)
  // {
  //   printf("ky[%ld] = %lg\n", jy, ky[jy]);
  // }

  for (jz = 0; jz < grid->nz; jz++)
  {
    for (jy = 0; jy < grid->ny; jy++)
    {
      jyz = grid->nx * (grid->ny * jz + jy);
      for (jx = 0; jx < grid->nx; jx++)
      {
        jgrid = jyz + jx;
        kindex[3 * jgrid] = kx[jx];
        kindex[3 * jgrid + 1] = ky[jy];
        kindex[3 * jgrid + 2] = kz[jz];
      }
      // printf("kindex[%ld] = %lg kindex[%ld] = %lg kindex[%ld] = %lg\n", 3 * jgrid, kx[jx], 3 * jgrid + 1, ky[jy], 3 * jgrid + 2, kz[jz]);
    }
  }

  return;
}

/*************************************************************************************************/
// k-vector initialization (works for both even and odd n)
// negative frequencies count backwards from end
// overall transform fw and bw aquires factor of nx*ny*nz so normalize here
void init_kvec(double *k, int n, double dk, double norm)
{
  int npos = (n % 2 == 0) ? n / 2 : (n - 1) / 2; // number of positive freqs
  k[0] = 0.0;

  for (int j = 1; j <= npos; j++)
  {
    double val = j * dk * norm;
    k[j] = val;      // positive frequencies
    k[n - j] = -val; // negative frequencies
  }

  return;
}

/*************************************************************************************************/

void p_operator(char *direc, double *kindex, double complex *psi, double complex *Lpsi, grid_st *grid, index_st *ist, par_st *par, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi)
{
  /*** First use the fft to get the action of the p operator on each axis
   * p = (-i d/dx) imag part from the fft definition***/
  // To compute d/dx, we compute FT^-1[-i*k*FT(psi)]
  // multiplying this by -i, we get p = -i * FT^-1[ -i*k * FT(psi)]
  // p = FT^-1[ -k * FT(psi)]. We already initialized the k-vectors to be -k.
  // so here we just multiply by k.

  long jgrid;
  long dir_offset;

  if (strcmp(direc, "X") == 0)
  {
    // printf("The deriv. is along X direction.\n");
    // printf("Setting dir_offset to 0.\n");
    dir_offset = 0;
  }
  else if (strcmp(direc, "Y") == 0)
  {
    // printf("The deriv. is along Y direction.\n");
    // printf("Setting dir_offset to 1.\n");
    dir_offset = 1;
  }
  else if (strcmp(direc, "Z") == 0)
  {
    // printf("The deriv. is along Z direction.\n");
    // printf("Setting dir_offset to 2.\n");
    dir_offset = 2;
  }
  else
  {
    printf("No direc recognized for derivative.\n");
    exit(EXIT_FAILURE);
  }

  // Copy psi to fftwpsi
  memcpy(&fftwpsi[0], &psi[0], grid->ngrid * sizeof(fftwpsi[0]));

  // FT from r-space to k-space
  fftw_execute(planfw);

  omp_set_dynamic(0);
  omp_set_num_threads(ist->nthreads);
#pragma omp parallel for private(jgrid)
  for (jgrid = 0; jgrid < grid->ngrid; jgrid++)
  {
    // multiply by -k_x/y/z to get p along x/y/z-axis
    fftwpsi[jgrid] *= kindex[3 * jgrid + dir_offset];
  }

  // Inverse FT back to r-space
  fftw_execute(planbw);
  // Copy fftwpsi to psi to store Lx|psi_qp> into |Lxpsi>
  memcpy(&Lpsi[0], &fftwpsi[0], ist->ngrid * sizeof(Lpsi[0]));

  return;
}
