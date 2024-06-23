/*****************************************************************************/

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
#include <omp.h>
#include "mkl.h"
#include "unistd.h"

/*****************************************************************************/
typedef struct flag {
  int SO, NL, LR, useSpinors, isComplex;
  int setSeed;
  int printFPDensity;
  int calcDarkStates, calcSpinAngStat;
  int timingSpecs, saveCheckpoints, restartFromCheckpoint, saveOutput;
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
  double dx, dy, dz, dr, dkx, dky, dkz, dv, epsX, epsY, epsZ, boxl, minr;
  double xmin, xmax, ymin, ymax, zmin, zmax, kxmin, kymin, kzmin;
  double KE_max, gamma, gamma2, Elmin, Elmax, Elumo, Ehomo, Vmin, Vmax;
  double delta_E_elec, delta_E_hole, sigma_E_cut, fermi_E;
  int checkpoint_id;
} par_st;

typedef struct st4 {
  long max_elec_states, max_hole_states, mn_states_tot;
  long n1, n2, n12, natom, nthreads;
  long ms, ms2, niter, nc, npot, npsi;
  long nx, ny, nz, ngrid, nspinngrid, n_atom_types;
  long *atom_types;
  long natoms, nhomo, nlumo, total_homo, total_lumo;
  double nx_1, ny_1, nz_1, ngrid_1;
  int n_FP_density;
  int printFPDensity; // 0 = False (default) or 1 = True
  int calcDarkStates; // 0 = False (default) or 1 = True
  int nspin, complex_idx;
} index_st;


typedef struct st5 {
  double x, y, z;
} xyz_st;

typedef fftw_plan fftw_plan_loc;

typedef struct parallel{
  long nthreads;
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
#define SVDEPS    1.0e-10
#define EPSR      1.0e-10
#define EPSCHI    1.0e-8
#define DENE      1.0e-10
#define N_MAX_ATOM_TYPES 20
/*****************************************************************************/

/*#define DEPS  0.02
#define DENERGY 0.01*/

// init.c
void init(double *potl,double *vx,double *vy,double *vz,double *ksqr,double *rx,double *ry,double *rz,par_st *par,index_st *ist);
void init_size(long, char *argv[],par_st *,index_st *);
void init_pot(zomplex *potq, zomplex *potqx, grid_st *grid, par_st *par,index_st *ist, fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void init_psi(zomplex *psi,double *vx,double *vy,double *vz,index_st ist,par_st par,long *idum);
double screenedcoulomb(double dr, double gamma);
double longerpolate(double r,double dr,double *vr,double *pot,long npot,long n,long j);
long assign_atom_number(char atyp[2]);
void assign_atom_type(char *atype,long j);

// read.c
void read_input(flag_st *flag, grid_st *grid, index_st *ist, par_st *par, parallel_st *parallel);
void read_conf(FILE *pf, xyz_st *R, long n);

// norm.c
double norm(zomplex *, double,long);
double normalize_zomplex(zomplex *psi, double dr, long ngrid);
void normalize(double *vector, double dV, long ngrid);
void normalize_all(double *,double,long,long);
void norm_vector(double *vector, double dV, long length);
double norm_rho(zomplex *rho,double dr,long ngrid);

// write.c
void write_psi(double *,double *,double *,double *,double *,index_st,par_st);
void write_pot(double *,double *,double *,double *,index_st);

void scalar_product(zomplex *,zomplex *,zomplex *,double,long,long);

void nerror(char *);

double rlong(double); 
/*
long int random(void);
void srandom(unsigned long);
*/
double ran_nrc(long *idum);
void Randomize();
double ran();

double get_dot_ligand_size(double *,double *,double *,long);
double get_dot_ligand_size_z(double *rz,long n);

void hartree(zomplex *rho,zomplex *potq,zomplex *poth,index_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);

void single_coulomb_openmp(zomplex       *psi, 
                           zomplex       *potq,
                           zomplex       *potqx,
                           zomplex       *poth,
                           double        *eval,
                           index_st        ist,
                           par_st        par,
                           fftw_plan_loc *planfw,
                           fftw_plan_loc *planbw,
                           fftw_complex  *fftwpsi,
                           zomplex       *bsmat,
                           zomplex       *direct,
                           zomplex       *exchange,
                           double       *h0mat);



void prlong_pz(double *psi,double *sige,double *vz,par_st par,index_st ist);

void diag(const int n, int nthreads, zomplex *mat, double *eval);
void bethe_salpeter(zomplex *bsmat, zomplex *direct, zomplex *exchage, double *h0mat, zomplex *psi, double *vz, zomplex *mux, zomplex *muy, zomplex * muz,
          double *mx, double *my, double *mz,zomplex *sx, zomplex *sy, zomplex *sz,zomplex *lx, zomplex *ly, zomplex *lz, zomplex* lsqr, zomplex* ls, index_st ist, par_st par);
void psi_rnd(zomplex *psi,long ngrid,double dv,long *idum);

double findmaxabsre(zomplex *dwmat,long n);
double findmaxabsim(zomplex *dwmat,long n);

// dipole.c
void dipole(zomplex *psi,zomplex *mux,zomplex *muy,zomplex *muz,double *eval,grid_st *grid, index_st *ist,par_st *par);
void mag_dipole(double *vx, double *vy, double *vz, double *psi, double *mx, double *my, double *mz, 
  double *eval, fftw_plan_loc *planfw,fftw_plan_loc *planbw,fftw_complex *fftwpsi, index_st ist, par_st par);
void rotational_strength(double *rs, double *mux, double *muy, double *muz, double *mx, 
  double *my, double *mz, double *eval, index_st ist);

void prlong_cube(double *pgrid,index_st ist,par_st par);


void print_pz_one(double *psi,double *vz,par_st par,index_st ist,char *str);
void print_pz(double *psi,double *sige,double *vz,par_st par,index_st ist);
int z_project(double *vector, double *vz, par_st par, index_st ist, char *fname);
void print_cube(double *pgrid,index_st ist,par_st par,char *fName);
void print_fixed_qp_density(double *psi, double *Cbs, double *vz, index_st ist, par_st par);

// Functions that write input or output - write.c
void write_cube_file(double *rho, grid_st *grid, char *fileName);
void write_current_time(FILE *pf);
void write_separation(FILE *pf, char *top_bttm);

//angular.c
void spins(zomplex *sx, zomplex *sy, zomplex *sz,zomplex *psi,index_st ist,par_st par);
void angular(zomplex* lx, zomplex* ly, zomplex* lz, zomplex* lsqr, zomplex* ls, zomplex *psi, grid_st *grid, index_st *ist, par_st *par,
  fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi);
void lOpp(zomplex* Lxpsi, zomplex* Lypsi, zomplex* Lzpsi, zomplex* psi, 
  double* vx, double* vy, double* vz,
  fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi,index_st ist, par_st par);

// save.c
void print_input_state(FILE *pf, flag_st *flag, grid_st *grid, par_st *par, index_st *ist, parallel_st *parallel);
void read_filter_output(char *file_name, double **psitot, double **eig_vals, double **sigma_E, xyz_st **R, grid_st *grid, double **gridx, double **gridy, double **gridz, index_st *ist, par_st *par, flag_st *flag);
/*****************************************************************************/
