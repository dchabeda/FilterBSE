/*****************************************************************************/
// Required libraries

#include <stdio.h>
#include <string.h>
#include "unistd.h"
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <time.h>
#include <assert.h>
#include <fftw3.h>
#include <omp.h>
#include <mkl.h>
#include "vector.h"

/*****************************************************************************/
// Application specific structures

typedef struct flag {
  int centerConf, setTargets, setSeed, interpolatePot, useStrain;
  int SO, NL, useSpinors, isComplex;
  int printPsiFilt, printOrtho, printNorm, printCubes;
  int calcPotOverlap, getAllStates, timeHamiltonian, calcSpinAngStat;
  int retryFilter, alreadyTried, saveCheckpoints, restartFromCheckpoint;
} flag_st;

typedef struct index {
  long m_states_per_filter, n_filter_cycles, mn_states_tot;
  long homo_idx, lumo_idx, total_homo, total_lumo;
  long ngrid, nspinngrid, ncheby;
  long natoms, n_atom_types, n_max_atom_types;
  long *atom_types;
  long max_pot_file_len, n_NL_gridpts, n_NL_atoms, nproj;
  int nspin, ncubes, ngeoms;
  int complex_idx;
  int crystal_structure_int, outmost_material_int;
  // Redundancies
  long nx, ny, nz;
  long nthreads;
} index_st;

typedef struct par {
  double Vmin, Vmax, VBmin, VBmax, CBmin, CBmax;
  double scale_surface_Cs;
  long n_targets_VB, n_targets_CB;
  long rand_seed;
  double KE_max, fermi_E, dt, dE, dE_1;
  double R_NLcut2, sigma_E_cut;
  int t_rev_factor;
  int checkpoint_id;
  char crystal_structure[15], outmost_material[15];
  // Redundnacies
  double dv;
} par_st;

typedef struct atom_info {
  long idx;
  char atyp[3];
  int Zval;
  double SO_par, geom_par, LR_par, NL_par[2];
} atom_info;

typedef struct st5 {
  double x, y, z;
} xyz_st;

typedef struct st10 {
  double SO[2], NL1[2], NL2[2];
} coeff_st;

typedef struct grid {
  double *x, *y, *z; // stores the value of the grid points
  double dx, dy, dz, dr, dv, dkx, dky, dkz;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  long nx, ny, nz;
  double nx_1, ny_1, nz_1;
  // Redundancies
  long ngrid;
} grid_st;

typedef struct pot_st {
  double *r, *r_LR;
  double *pseudo, *pseudo_LR;
  double *dr;
  long *file_lens;
  double *a4_params, *a5_params;
} pot_st;

typedef struct zomplex {
  double re, im;
} zomplex;

typedef struct st11 {
  long jxyz;
  zomplex y1[3];
  double proj[5];
  double NL_proj[5];
  int NL_proj_sign[5];
  double r, r2_1, r2, Vr;
} nlc_st;

typedef fftw_plan fftw_plan_loc;

typedef struct parallel{
  long nthreads;
} parallel_st;



/*****************************************************************************/
// Macro definitions

// double macros
#define sqr(x)       ((x) * (x))
#define cube(x)      ((x) * (x) * (x))
#define forth(x)     ((x) * (x) * (x) * (x))
#define ISWAP(a,b)   {double temp = (a); (a) = (b); (b) = -temp;}
#define ISWAP2(a,b)  {double temp = (a); (a) = 2.0 * (b); (b) = -2.0 * temp;}

// complex number macros
#define cplus(a,b,c)  {zomplex tmp; (tmp).re=(a).re+(b).re; (tmp).im=(a).im+(b).im; c=tmp;}
#define cminus(a,b,c) {zomplex tmp; (tmp).re=(a).re-(b).re; (tmp).im=(a).im-(b).im; c=tmp;}
#define cmul(a,b,c)   {zomplex tmp; (tmp).re=(a).re*(b).re-(a).im*(b).im; (tmp).im=(a).re*(b).im+(a).im*(b).re; c=tmp;}
#define cmuls(a,b,c)  {zomplex tmp; (tmp).re=(a).re*(b).re+(a).im*(b).im; (tmp).im=(a).re*(b).im-(a).im*(b).re; c=tmp;}
#define cdev(a,b,c)   {zomplex tmp; double mechane; mechane = (1.0 / ((b).re*(b).re+(b).im*(b).im)); (tmp).re=((a).re*(b).re+(a).im*(b).im) * mechane; (tmp).im = ((a).im * (b).re - (a).re*(b).im)*mechane; c=tmp;}
#define cexp(a,c)     {double texp = exp((a).re); zomplex tmp; (tmp).re=texp * cos((a).im); (tmp).im = texp * sin((a).im); c= tmp;}
#define cexpminx(a,c) {double texp = exp(-(a).re); (c).re=texp * cos((a).im); (c).im = -texp * sin((a).im);}

// constant macros
#define AUTOEV    27.211385
#define PIE       3.14159265358979323846
#define TWOPI     6.28318530717958647692
#define SVDEPS    1.0e-10
#define EPSR      1.0e-10
#define EPSE      1.0e-10
#define EPSR02	  1.0e-2
#define EPSR0     1.0e-4
#define EPSCHI    1.0e-8
#define DENE      1.0e-10
#define EPSDX     1.0e-20
#define ANGTOBOHR 1.889726125
#define PROJ_LEN	1024
#define N_MAX_ATOM_TYPES 20
// physical parameters
#define PbIBondMax  ( (3.3) * (ANGTOBOHR) ) 
#define orthoBondAngle  160.6344 
#define minCubicBondAngle  175.0

/*****************************************************************************/
// Function declarations

//init.c
void init_grid_params(grid_st *grid, xyz_st *R, index_st *ist, par_st *par);
void build_grid_ksqr(double *ksqr, xyz_st *R, grid_st *grid, index_st *ist, par_st *par, flag_st *flag);
void build_local_pot(double *pot_local, pot_st *pot, xyz_st *R, double *ksqr, atom_info *atom, grid_st *grid,
    index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
void set_ene_targets(double *ene_targets, index_st *ist, par_st *par, flag_st *flag);
void init_SO_projectors(double *SO_projectors, grid_st *grid, xyz_st *R, atom_info *atm, index_st *ist, par_st *par);
void init_NL_projectors(nlc_st *nlc, long *nl, double *SO_projectors, grid_st *grid, xyz_st *R, atom_info *atm, index_st *ist, par_st *par, flag_st *flag);
void init_psi(zomplex *psi, long *rand_seed, int isComplex, grid_st *grid, parallel_st *parallel);
double calc_dot_dimension(xyz_st *R, long n, char *dir);
double ret_ideal_bond_len(long natyp_1, long natyp_2, int crystal_structure_int);

//read.c
void read_input(flag_st *flag, grid_st *grid, index_st *ist, par_st *par, parallel_st *parallel);
void read_conf(xyz_st *R, atom_info *atm, index_st *ist, par_st *par, flag_st *flag);
void read_pot(pot_st *pot, xyz_st *R, atom_info *atom, index_st *ist, par_st *par, flag_st *flag);
void interpolate_pot(xyz_st *R, atom_info *atom, index_st *ist, par_st *par);
void calc_geom_par(xyz_st *R,atom_info *atm, index_st *ist );
double calc_bond_angle(long index1,long index2,long index3, xyz_st *R);
long assign_atom_number(char atyp[3]);
void assign_atom_type(char *atype,long j);
long get_number_of_atom_types(atom_info *atm,index_st *ist, long *list);
int assign_crystal_structure(char *crystal_structure);
int assign_outmost_material(char *outmost_material);
double get_ideal_bond_len(long natyp_1, long natyp_2, int crystalStructureInt);

// strain.c
void read_nearest_neighbors(vector *atom_neighbors, double *tetrahedron_vol_ref, long natoms, int crystal_structure, int outmost_material);
void calc_strain_scale(double *strain_scale, vector *atom_neighbors, double *tetrahedron_vol_ref, atom_info *atom, double *a4_params, double *a5_params, long natoms);
double calc_regular_tetrahedron_volume(double bond_length1, double bond_length2, double bond_length3, double bond_length4);

//interpolate.c
double interpolate(double r,double dr,double *vr,double *vr_LR,double *pot,double *pot_LR,long potFileLen,long n,long j, int scale_LR, double scale_LR_par, double strain_factor);

//rand.c
double ran_nrc(long *rand_seed);
void Randomize();

//hamiltonian.c
void kinetic(zomplex *psi, double *ksqr, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi, index_st *ist);
void potential(zomplex *phi, zomplex *psi, double *pot_local, nlc_st *nlc, long *nl, index_st *ist, par_st *par, flag_st *flag);
void spin_orbit_proj_pot(zomplex *phi, zomplex *psi, nlc_st *nlc, long* nl, index_st *ist, par_st *par);
void nonlocal_proj_pot(zomplex *phi, zomplex *psi, nlc_st *nlc, long* nl, index_st *ist, par_st *par);
void hamiltonian(zomplex *phi, zomplex *psi, double *pot_local, nlc_st *nlc, long *nl, double *ksqr,
  index_st *ist, par_st *par, flag_st *flag, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi);
void time_reverse_all(double *psitot, double *dest, index_st *ist, parallel_st *parallel);

//filter.c
void run_filter_cycle(double *psims, double *pot_local, nlc_st *nlc, long *nl, double *ksqr, zomplex *an, 
  double *zn, double *ene_targets, long thread_id, long jns, grid_st *grid, index_st *ist, par_st *par,
  flag_st *flag, parallel_st *parallel);
void filter(zomplex *psin, zomplex *psim1, double *psims, double *pot_local, nlc_st *nlc, long *nl, 
  double *ksqr, zomplex *an, double *zn, long thread_id, long jn, fftw_plan_loc planfw, fftw_plan_loc planbw,
  fftw_complex *fftwpsis, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
void scale_eigs_for_cheby(zomplex *psim1, zomplex *psin, double *pot_local, nlc_st *nlc, long *nl, double *ksqr, double zm1,
  fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi, index_st *ist, par_st *par, flag_st *flag);

//norm.c
double calc_norm(zomplex *, double,long,long);
double normalize(zomplex *, double,long,long);
void normalize_all(double *psi, double dv, long ms, long ngrid, long nthreads, int complex_idx, int printNorm);

//energy.c
double energy(zomplex *psi, zomplex *phi, double *pot_local, nlc_st *nlc, long *nl, double *ksqr, index_st *ist,
  par_st *par, flag_st *flag, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi);
void energy_all(zomplex *psi, zomplex *phi, double *psims, double *pot_local, nlc_st *nlc, long *nl, double *ksqr,
  double *ene, index_st *ist, par_st *par, flag_st *flag, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi);
void get_energy_range(zomplex *psi, zomplex *phi, double *pot_local, grid_st *grid, nlc_st *nlc, long *nl, double *ksqr,
  index_st *ist, par_st *par, parallel_st *parallel, flag_st *flag, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi);
void calc_sigma_E(zomplex *psi, zomplex *phi, double *psitot, double *pot_local, nlc_st *nlc, long *nl, double *ksqr,
  double *eval2, index_st *ist, par_st *par, flag_st *flag, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi);


//coeff.c
void gen_newton_coeff(zomplex *an, double *samp, double *ene_targets, index_st *ist, par_st *par, parallel_st *parallel);
void chebyshev_reordered(double *,double,double,long);
double samp_points_ashkenazy(zomplex *point,double min,double max,long ncheby);
void check_function(zomplex *an,zomplex *samp,index_st *ist,par_st *par, double ene_targets);


//Hmat.c
void diag_H(zomplex *psi,zomplex *phi, double *psitot,double *pot_local,nlc_st *nlc,long *nl,double *ksqr,double *eval,index_st *ist,par_st *par,flag_st *flag,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
MKL_Complex16 dotp(zomplex *psi, double *phi,long m,long ngrid,double dv);
double dotpreal(zomplex *psi,double *phi,long m,long ngrid,double dv);

//nerror.c
void nerror(char *);

//ortho.c
long ortho(double *psitot, double dv, index_st *ist, par_st *par, flag_st *flag);

//projectors.c
void gen_SO_projectors(double dx, double rcut, long nproj,double*  projectors, double* vr);
void gen_nlc_projectors(double dx, double rcut, long nproj, double*  projectors,int* sgnProj, double* vr, atom_info *atm,long jatom);

// pseudopotential.c
void calc_pot_overlap(double *psitot, double *pot_local, nlc_st *nlc, long *nl, double *eval, par_st *par,index_st *ist, flag_st *flag);

//angular.c
void calc_angular_exp(double *psitot, grid_st *grid, int start, int stop, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi);
void apply_J_op( zomplex* Jxpsi, zomplex* Jypsi, zomplex* Jzpsi,zomplex* psi, 
	grid_st* grid, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi,index_st *ist, par_st *par, parallel_st *parallel);
void apply_L_op(zomplex* Lxpsi, zomplex* Lypsi, zomplex* Lzpsi, zomplex* psi, 
	grid_st* grid, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi,index_st *ist, par_st *par, parallel_st *parallel);
void low_pass_filter(zomplex* psi, grid_st *grid, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi,
	index_st *ist, par_st *par, parallel_st *parallel);

//write.c
void write_cube_file(double *rho, grid_st *grid, char *fileName);
void write_separation(FILE *pf, char *top_bttm);

//save.c
void print_input_state(FILE *pf, flag_st *flag, grid_st *grid, par_st *par, index_st *ist, parallel_st *parallel);
void save_job_state(char *file_name, int checkpoint_id, double *psitot, double *pot_local, double *ksqr, zomplex *an, double *zn, double *ene_targets, long *nl,\
    nlc_st *nlc, grid_st *grid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
void restart_from_save(char *file_name, int checkpoint_id, double *psitot, double *pot_local, double *ksqr, zomplex *an, double *zn, double *ene_targets, long *nl,\
    nlc_st *nlc, grid_st *grid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
/*****************************************************************************/
