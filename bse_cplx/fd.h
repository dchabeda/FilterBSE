/*****************************************************************************/
#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <complex.h>
#include <malloc.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>
#include "unistd.h"
#include <fftw3.h>
// #include <lapack.h>
#include <mpi.h>
#include <omp.h>
#include <nvToolsExt.h>
#include <mkl.h>

/*****************************************************************************/
typedef struct flag
{
  int SO, NL, LR, useSpinors, isComplex;
  int setSeed;
  int printFPDensity;
  int calcDarkStates, calcSpinAngStat;
  int timingSpecs, saveCheckpoints, restartFromChk, saveOutput;
  int restartCoulomb, coulombDone, calcCoulombOnly, noCalcExciton;
  int initUnsafe;
} flag_st;

typedef struct grid
{
  double *x, *y, *z; // stores the value of the grid points
  double dx, dy, dz, dr, dv, dkx, dky, dkz;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  long nx, ny, nz;
  double nx_1, ny_1, nz_1;
  // Redundancies
  long ngrid;
} grid_st;

struct pzdata
{
  double z;
  double pz;
};

typedef struct st0
{
  double re, im;
} zomplex;

typedef struct st1
{
  double dx, dy, dz, dr, dkx, dky, dkz, dv, epsX, epsY, epsZ;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double KE_max;
  double delta_E_elec, delta_E_hole, sigma_E_cut, fermi_E;
  long psi_file_pos;
  int checkpoint_id;
} par_st;

typedef struct st4
{
  long max_elec_states, max_hole_states, mn_states_tot;
  long nthreads;
  long n_qp, n_xton;
  long nx, ny, nz, ngrid, nspinngrid;
  long n_atom_types, *atom_types;
  long natoms, homo_idx, lumo_idx, n_holes, n_elecs;
  long *eval_hole_idxs, *eval_elec_idxs;
  long *jsg_conv;
  double ngrid_1;
  int n_FP_density;
  int printFPDensity; // 0 = False (default) or 1 = True
  int calcDarkStates; // 0 = False (default) or 1 = True
  int nspin, complex_idx;
} index_st;

typedef struct st5
{
  double complex x, y, z;
} xyz_st;

typedef fftw_plan fftw_plan_loc;

typedef struct parallel
{
  long nthreads;
  int mpi_rank, mpi_size, mpi_root;
} parallel_st;

/*****************************************************************************/

#define rlong(x) (rint(x))

#define cnorm(z) (creal(z) * creal(z) + cimag(z) * cimag(z))
#define conjmul(x, y) ((conj(x)) * (y))

#define sqr(x) ((x) * (x))
#define cube(x) ((x) * (x) * (x))
#define forth(x) ((x) * (x) * (x) * (x))

#define ISWAP(a, b)    \
  {                    \
    double temp = (a); \
    (a) = (b);         \
    (b) = -temp;       \
  }
#define ISWAP2(a, b)   \
  {                    \
    double temp = (a); \
    (a) = 2.0 * (b);   \
    (b) = -2.0 * temp; \
  }

#define AUTOEV 27.2114
#define PIE 3.14159265358979323846
#define TWOPI 6.28318530717958647692
#define FOURPI (4.0 * 3.14159265358979323846)
#define SQRTPI (sqrt(3.14159265358979323846))
#define SVDEPS 1.0e-10
#define EPSR 1.0e-10
#define EPSCHI 1.0e-8
#define DENE 1.0e-10
#define N_MAX_ATOM_TYPES 20
// #define PR_LEN 16
#define PR_LEN 6
#define BYTE_BOUNDARY 32
// Enum to define supported variable types
typedef enum
{
  INT_TYPE,
  LONG_TYPE,
  DOUBLE_TYPE
} VarType;

// memory allocation macro for code readability
#define ALLOCATE(dblptr, length, message) allocate_aligned_memory((void **)(dblptr), (length), sizeof(typeof(**(dblptr))), message)

/*****************************************************************************/

// save.c
void print_input_state(FILE *pf, flag_st *flag, grid_st *grid, par_st *par, index_st *ist, parallel_st *parallel);
void read_filter_output(char *file_name, double **eig_vals, double **sigma_E, xyz_st **R, grid_st *grid, double **gridx, double **gridy, double **gridz, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);

// norm.c
double norma(zomplex *, double, long);
double normalize_zomplex(zomplex *psi, double dr, long ngrid);
void normalize_all(double *, double, long, long);
void norm_vector(double *vector, double dV, long length);

void scalar_product(zomplex *, zomplex *, zomplex *, double, long, long);

double findmaxabsre(zomplex *dwmat, long n);
double findmaxabsim(zomplex *dwmat, long n);

void print_pz_one(double *psi, double *vz, par_st par, index_st ist, char *str);
void print_pz(double *psi, double *sige, double *vz, par_st par, index_st ist);
int z_project(double *vector, double *vz, par_st par, index_st ist, char *fname);
void print_cube(double *pgrid, index_st ist, par_st par, char *fName);
void print_fixed_qp_density(double *psi, double *Cbs, double *vz, index_st ist, par_st par);

// optical.c
void calc_optical_exc(double complex *bs_coeff, double *xton_ene, double *eig_vals, xyz_st *mu, xyz_st *m, index_st *ist, par_st *par);
/*****************************************************************************/

#include "aux.h"
#include "write.h"
