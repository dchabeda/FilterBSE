/*****************************************************************************/
#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <malloc.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>
#include <fftw3.h>
#include <lapack.h>
#include <mpi.h>
#include <nvToolsExt.h>
#include <omp.h>
#include "unistd.h"

/*****************************************************************************/
typedef struct flag {
  int SO, NL, LR, useSpinors, isComplex;
  int setSeed;
  int printFPDensity;
  int calcDarkStates, calcSpinAngStat;
  int timingSpecs, saveCheckpoints, restartFromChk, saveOutput;
  int restartCoulomb, coulombDone, calcCoulombOnly;
  int initUnsafe;
} flag_st;

typedef struct grid {
  double *x, *y, *z; // stores the value of the grid points
  double dx, dy, dz, dr, dv, dkx, dky, dkz;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  long nx, ny, nz;
  double nx_1, ny_1, nz_1;
  // Redundancies
  long ngrid;
} grid_st;

struct pzdata {
    double z;
    double pz;
};

typedef struct st0 {
  double re, im;
} zomplex;

typedef struct st1 {
  double dx, dy, dz, dr, dkx, dky, dkz, dv, epsX, epsY, epsZ;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double KE_max;
  double delta_E_elec, delta_E_hole, sigma_E_cut, fermi_E;
  int checkpoint_id;
} par_st;

typedef struct st4 {
  long max_elec_states, max_hole_states, mn_states_tot;
  long nthreads;
  long n_qp, n_xton;
  long nx, ny, nz, ngrid, nspinngrid;
  long n_atom_types, *atom_types;
  long natoms, homo_idx, lumo_idx, n_holes, n_elecs;
  long *eval_hole_idxs, *eval_elec_idxs;
  long *jgur;
  long *jgui;
  long *jgdr;
  long *jgdi;
  double ngrid_1;
  int n_FP_density;
  int printFPDensity; // 0 = False (default) or 1 = True
  int calcDarkStates; // 0 = False (default) or 1 = True
  int nspin, complex_idx;
} index_st;


typedef struct st5 {
  double x, y, z;
  double x_re, x_im, y_re, y_im, z_re, z_im;
} xyz_st;

typedef fftw_plan fftw_plan_loc;

typedef struct parallel{
  long nthreads;
  int mpi_rank, mpi_size, mpi_root;
} parallel_st;

/*****************************************************************************/

#define rlong(x)      (rint(x)) 

#define sqr(x)       ((x) * (x))
#define cube(x)      ((x) * (x) * (x))
#define forth(x)     ((x) * (x) * (x) * (x))

#define ISWAP(a,b)   {double temp = (a); (a) = (b); (b) = -temp;}
#define ISWAP2(a,b)  {double temp = (a); (a) = 2.0 * (b); (b) = -2.0 * temp;}

#define AUTOEV    27.2114
#define PIE       3.14159265358979323846
#define TWOPI     6.28318530717958647692
#define FOURPI    (4.0*3.14159265358979323846)
#define SQRTPI    (sqrt(3.14159265358979323846))
#define SVDEPS    1.0e-10
#define EPSR      1.0e-10
#define EPSCHI    1.0e-8
#define DENE      1.0e-10
#define N_MAX_ATOM_TYPES 20
#define PR_LEN 16
#define BYTE_BOUNDARY 32
// Enum to define supported variable types
typedef enum { INT_TYPE, LONG_TYPE, DOUBLE_TYPE } VarType;

// memory allocation macro for code readability
#define ALLOCATE(dblptr, length, message) allocate_aligned_memory((void**)(dblptr), (length), sizeof(typeof(**(dblptr))), message)


/*****************************************************************************/

/*#define DEPS  0.02
#define DENERGY 0.01*/


// save.c
void print_input_state(FILE *pf, flag_st *flag, grid_st *grid, par_st *par, index_st *ist, parallel_st *parallel);
void read_filter_output(char *file_name, double **psitot, double **eig_vals, double **sigma_E, xyz_st **R, grid_st *grid, double **gridx, double **gridy, double **gridz, index_st *ist, par_st *par, flag_st *flag);

// write.c
void write_cube_file(double *rho, grid_st *grid, char *fileName);
void write_current_time(FILE *pf);
void write_separation(FILE *pf, char *top_bttm);
void write_state_dat(zomplex *psi, long n_elems, char* fileName);

// norm.c
double norm(zomplex *, double,long);
double normalize_zomplex(zomplex *psi, double dr, long ngrid);
void normalize_all(double *,double,long,long);
void norm_vector(double *vector, double dV, long length);

void scalar_product(zomplex *,zomplex *,zomplex *,double,long,long);

// hartree.c
void hartree(zomplex *rho,zomplex *potq,zomplex *poth,index_st *ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);

void calc_eh_kernel_real(
  double        *psi_qp, 
  zomplex       *pot_bare,
  zomplex       *pot_screened,
  zomplex       *pot_htree,
  zomplex       *bsmat,
  zomplex       *direct,
  zomplex       *exchange,
  double        *h0mat,
  double        *eval,
  index_st      *ist,
  par_st        *par,
  flag_st       *flag,
  fftw_plan_loc *planfw,
  fftw_plan_loc *planbw,
  fftw_complex  *fftwpsi);

void diag(const int n, int nthreads, double _Complex *mat, double *eval);
void bethe_salpeter(zomplex *bsmat, zomplex *direct, zomplex *exchage, zomplex *bs_coeff, double *h0mat, double *xton_ene, zomplex *psi, 
          xyz_st *s_mom, xyz_st *l_mom, zomplex *l2_mom, zomplex *LdotS, grid_st *grid, index_st *ist, par_st *par);

double findmaxabsre(zomplex *dwmat,long n);
double findmaxabsim(zomplex *dwmat,long n);

void print_pz_one(double *psi,double *vz,par_st par,index_st ist,char *str);
void print_pz(double *psi,double *sige,double *vz,par_st par,index_st ist);
int z_project(double *vector, double *vz, par_st par, index_st ist, char *fname);
void print_cube(double *pgrid,index_st ist,par_st par,char *fName);
void print_fixed_qp_density(double *psi, double *Cbs, double *vz, index_st ist, par_st par);


// optical.c
void calc_optical_exc(zomplex *bs_coeff, double *eval, xyz_st *mu, xyz_st *m, index_st *ist, par_st *par);
/*****************************************************************************/

long load_coulomb_mat(zomplex* mat, char* fileName, index_st* ist);
void build_h0_mat(double *h0mat, double *eval, index_st* ist);
void build_BSE_mat(zomplex *bsmat, zomplex *direct, zomplex *exchange, index_st* ist);
