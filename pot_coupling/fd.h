/*****************************************************************************/
#pragma once

// Standard library includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <malloc.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>
#include <unistd.h>

// External libraries
#include <fftw3.h>
#include <lapack.h>
#include <mpi.h>
#include <nvToolsExt.h>
#include <omp.h>

/*****************************************************************************/

// Complex number struct
typedef struct st0 {
  double re, im;
} zomplex;

// Basic real-space grid
typedef struct grid {
  double *x, *y, *z;
  double dx, dy, dz, dr, dv, dkx, dky, dkz;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  long nx, ny, nz;
  double nx_1, ny_1, nz_1;
  long ngrid;
} grid_st;

// Simulation flags
typedef struct flag {
  int printFPDensity;
  int calcDarkStates;
  int timingSpecs, restartFromChk;
  int restartCoulomb, coulombDone, calcCoulombOnly;
  int initUnsafe;

  // Configuration flags
  int   centerConf;
  int   setTargets;
  int   setSeed;
  int   interpolatePot;
  int   useStrain;

  // Spin and interaction settings
  int   SO;
  int   NL;
  int   LR;
  int   useSpinors;
  int   isComplex;

  // Input control
  int   inputPsiFilt;

  // Output control
  int   printPsiFilt;
  int   printPsiDiag;
  int   printPsiOrtho;
  int   printNorm;
  int   printCubes;
  int   printGaussCubes;

  // Computation settings
  int   calcPotOverlap;
  int   getAllStates;
  int   timeHamiltonian;
  int   calcSpinAngStat;

  // Restart and checkpointing
  int   retryFilter;
  int   alreadyTried;
  int   restartFromOrtho;
  int   restartFromSigma;
  int   saveCheckpoints;
  int   restartFromCheckpoint;
  int   saveOutput;

  // Special modes
  int   calcFilterOnly;
  int   approxEnergyRange;
  int   readProj;
  int   noTimeRev;

  // Performance and parallelization
  int   useFastHam;
  int   useMPIOMP;
  int   MPIOrtho;
  int   MPIDiag;
  int   parallelSigma;

  // Structural settings
  int   periodic;
  int   useGaussianBasis;
  int   readKPath;
} flag_st;

// Parameters
typedef struct par {
  double delta_E_elec, delta_E_hole;
  double Vmin, Vmax, VBmin, VBmax, CBmin, CBmax;
  double scale_surface_Cs;
  long n_targets_VB, n_targets_CB;
  long rand_seed;
  double KE_max, fermi_E, dt, dE, dE_1;
  double R_NLcut2, sigma_E_cut;
  int t_rev_factor, checkpoint_id;
  char crystal_structure[15], outmost_material[15];
  char fft_wisdom_dir[200], fftw_wisdom[250];
  double psi_zero_cut;
  int ham_threads;
  int n_orbitals, n_orbitals_per_atom, n_gauss_per_orbital;
  double R_gint_cut2, pot_cut_rad2, box_z;
  int nb_min, nb_max;
  double dv;
  double dx, dy, dz, dr, dkx, dky, dkz, epsX, epsY, epsZ;
  double xmin, xmax, ymin, ymax, zmin, zmax;

  double phonon_shift;
} par_st;

// Indexing and state management
typedef struct index {
  long m_states_per_filter, n_filter_cycles, mn_states_tot;
  long n_filters_per_rank, n_states_per_rank;
  long n_states_for_ortho, n_states_for_sigma;
  long long psi_rank_size;
  long init_psi_start, init_psi_end;
  long homo_idx, lumo_idx, total_homo, total_lumo;
  long max_elec_states, max_hole_states, n_qp, n_xton, n_holes, n_elecs;
  long *eval_hole_idxs, *eval_elec_idxs;
  long ngrid, nspinngrid, ncheby, natoms, n_atom_types, n_max_atom_types;
  long *atom_types, max_pot_file_len;
  long n_NL_gridpts, n_NL_atoms, nproj;
  int nspin, ncubes, n_gcubes_s, n_gcubes_e, ngeoms, complex_idx;
  double ngrid_1;
  int n_FP_density, printFPDensity, calcDarkStates;
  int crystal_structure_int, outmost_material_int;
  int n_s_ang_mom, n_l_ang_mom, n_j_ang_mom;
  int n_bands, n_G_vecs, n_G_zeros, n_k_pts, nk1, nk2, nk3;
  long nx, ny, nz, nthreads;
} index_st;

// Atom and potential data
typedef struct atom_info {
  long idx;
  char atyp[3];
  int Zval;
  double SO_par, geom_par, LR_par;
  double NL_par[2];
} atom_info;

typedef struct pot_st {
  double *r, *r_LR, *pseudo, *pseudo_LR, *dr;
  long *file_lens;
  double *a4_params, *a5_params;
} pot_st;

// Projector and grid-point data
typedef struct st11 {
  long jxyz;
  zomplex y1[3];
  double proj[5], NL_proj[5];
  int NL_proj_sign[5];
  double r, r2_1, r2, Vr;
} nlc_st;

typedef struct st5 {
  double x, y, z;
  double x_re, x_im, y_re, y_im, z_re, z_im;
} xyz_st;

// pzdata struct (used in projection analysis)
struct pzdata {
  double z;
  double pz;
};

// Parallel environment info
typedef fftw_plan fftw_plan_loc;
typedef struct parallel {
  long nthreads;
  int mpi_rank, mpi_size, mpi_root;
  long   nranks;
  int    nestedOMP;
  int    n_outer_threads;
  int    n_inner_threads;
  long  *jns;
  long  *jms;
} parallel_st;

/************************************************************/
/************************************************************/

typedef struct st10 {
  double  SO[2];
  double  NL1[2];
  double  NL2[2];
} coeff_st;


// Memory and math macros
#define rlong(x)      (rint(x)) 
#define sqr(x)       ((x) * (x))
#define cube(x)      ((x) * (x) * (x))
#define forth(x)     ((x) * (x) * (x) * (x))
#define ISWAP(a,b)   {double temp = (a); (a) = (b); (b) = -temp;}
#define ISWAP2(a,b)  {double temp = (a); (a) = 2.0 * (b); (b) = -2.0 * temp;}

// Physical constants and thresholds
#define AUTOEV    27.2114
#define PIE       3.14159265358979323846
#define TWOPI     6.28318530717958647692
#define FOURPI    (4.0 * PIE)
#define SQRTPI    (sqrt(PIE))
#define SVDEPS    1.0e-10
#define EPSR      1.0e-10
#define EPSCHI    1.0e-8
#define DENE      1.0e-10
#define ANGTOBOHR 1.889726125
#define PROJ_LEN	1024
#define N_MAX_G_VECS 10000
#define EPSDX     1.0e-20
#define EPSR02	  1.0e-2
#define EPSR0     1.0e-4
#define Qa_TO_A   0.0982270230776
// Limits and enums
#define N_MAX_ATOM_TYPES 20
#define PR_LEN 16
#define BYTE_BOUNDARY 32
typedef enum { INT_TYPE, LONG_TYPE, DOUBLE_TYPE } VarType;


// physical parameters
#define PbIBondMax  ( (3.3) * (ANGTOBOHR) ) 
#define orthoBondAngle  160.6344 
#define minCubicBondAngle  175.0
#define HBAR 1.0
#define MASS_E 1.0


// Memory allocation macro
#define ALLOCATE(dblptr, length, message) allocate_aligned_memory((void**)(dblptr), (length), sizeof(typeof(**(dblptr))), message)

/*****************************************************************************/
// Function declarations

// save.c
void print_input_state(FILE *pf, flag_st *flag, grid_st *grid, par_st *par, index_st *ist, parallel_st *parallel);
void read_filter_output(char *file_name, double **psitot, double **eig_vals, double **sigma_E, xyz_st **R, grid_st *grid, double **gridx, double **gridy, double **gridz, index_st *ist, par_st *par, flag_st *flag);

// write.c
void write_cube_file(double *rho, grid_st *grid, char *fileName);
void write_current_time(FILE *pf);
void write_separation(FILE *pf, char *top_bttm);
void write_state_dat(zomplex *psi, long n_elems, char* fileName);

// norm.c
double norm(zomplex *, double, long);
double normalize_zomplex(zomplex *psi, double dr, long ngrid);
void normalize_all(double *, double, long, long);
void norm_vector(double *vector, double dV, long length);
void scalar_product(zomplex *, zomplex *, zomplex *, double, long, long);

// hartree.c
void hartree(zomplex *rho, zomplex *potq, zomplex *poth, index_st *ist, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi);

void calc_eh_kernel_real(
  double *psi_qp, 
  zomplex *pot_bare,
  zomplex *pot_screened,
  zomplex *pot_htree,
  zomplex *bsmat,
  zomplex *direct,
  zomplex *exchange,
  double *h0mat,
  double *eval,
  index_st *ist,
  par_st *par,
  flag_st *flag,
  fftw_plan_loc *planfw,
  fftw_plan_loc *planbw,
  fftw_complex *fftwpsi);

void diag(const int n, int nthreads, double _Complex *mat, double *eval);
void bethe_salpeter(zomplex *bsmat, zomplex *direct, zomplex *exchage, zomplex *bs_coeff, double *h0mat, double *xton_ene, zomplex *psi, 
                    xyz_st *s_mom, xyz_st *l_mom, zomplex *l2_mom, zomplex *LdotS, grid_st *grid, index_st *ist, par_st *par);

double findmaxabsre(zomplex *dwmat, long n);
double findmaxabsim(zomplex *dwmat, long n);

void print_pz_one(double *psi, double *vz, par_st par, index_st ist, char *str);
void print_pz(double *psi, double *sige, double *vz, par_st par, index_st ist);
int z_project(double *vector, double *vz, par_st par, index_st ist, char *fname);
void print_cube(double *pgrid, index_st ist, par_st par, char *fName);
void print_fixed_qp_density(double *psi, double *Cbs, double *vz, index_st ist, par_st par);

// optical.c
void calc_optical_exc(zomplex *bs_coeff, double *eval, xyz_st *mu, xyz_st *m, index_st *ist, par_st *par);

// coulomb.c
long load_coulomb_mat(zomplex* mat, char* fileName, index_st* ist);
void build_h0_mat(double *h0mat, double *eval, index_st* ist);
void build_BSE_mat(zomplex *bsmat, zomplex *direct, zomplex *exchange, index_st* ist);

/*****************************************************************************/
/*****************************************************************************/
// Function declarations

void allocate_memory(void **ptr, size_t length, size_t type_size, char* message);
//init.c
void init_grid_params(grid_st *grid, xyz_st *R, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
void build_grid_ksqr(double *ksqr, xyz_st *R, grid_st *grid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
void set_ene_targets(double *ene_targets, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
void init_psi(zomplex *psi, long *rand_seed, grid_st *grid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
double calc_dot_dimension(xyz_st *R, long n, char *dir);
double ret_ideal_bond_len(long natyp_1, long natyp_2, int crystal_structure_int);

// gauss.c
void gauss_driver(
  gauss_st *gauss, double *pot_local, double *eig_vals, xyz_st *R, atom_info *atom, grid_st *grid,
  index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
void init_gauss_params(gauss_st *gauss, xyz_st *R, atom_info *atom, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
void build_gauss_hamiltonian(double *H_mat, double *T_mat, double *V_mat, index_st *ist, par_st *par);
void proj_gauss_on_grid(MO_st *MO, gauss_st *gauss, grid_st *grid, int start, int end, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
void calc_X_canonical(double *S_mat, double *X, double *U, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
void transform_H(double *H_mat, double *X, double *eig_vals, MO_st *MO, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);

//read.c
void read_input(flag_st *flag, grid_st *grid, index_st *ist, par_st *par, parallel_st *parallel);
void read_conf(xyz_st *R, atom_info *atm, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
void read_pot(pot_st *pot, xyz_st *R, atom_info *atom, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
void read_pot_file(FILE *pf, pot_st *pot, long j, long n, char *req);
void interpolate_pot(xyz_st *R, atom_info *atom, index_st *ist, par_st *par, parallel_st *parallel);
void read_periodic_input(lattice_st *lattice, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
void calc_geom_par(xyz_st *R,atom_info *atm, index_st *ist, parallel_st *parallel);
double calc_bond_angle(long index1,long index2,long index3, xyz_st *R, parallel_st *parallel);
long assign_atom_number(char atyp[4]);
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
double interpolate(double r,double dr,double *vr,double *vr_LR,double *pot,double *pot_LR,long pot_file_len,long n,long j, int scale_LR, double scale_LR_par, double strain_factor, int is_LR);

//rand.c
double ran_nrc(long *rand_seed);
void Randomize(void);


//filter.c
void scale_eigs_for_cheby(zomplex *phi, zomplex *psi, double *pot_local, vector *G_vecs, vector k, grid_st *grid,nlc_st *nlc, long *nl, double zm1,
  fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi, index_st *ist, par_st *par, flag_st *flag);


// integral.c
void overlap_gauss(double *S_mat, gauss_st *gauss, atom_info *atom, index_st *ist, par_st *par, flag_st *flag);
void kinetic_gauss(double *T_mat, gauss_st *gauss, atom_info *atom, index_st *ist, par_st *par, flag_st *flag);
void potential_gauss(double *V_mat, double *pot_local, gauss_st *gauss, grid_st *grid, atom_info *atom, index_st *ist, par_st *par, flag_st *flag);

// norm.c
double calc_norm(zomplex *, double,long,long);
double normalize(zomplex *, long ngrid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
void normalize_all(double *psitot, long n_states, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);


//coeff.c
void chebyshev_reordered(double *,double,double,long);
double samp_points_ashkenazy(zomplex *point,double min,double max,long ncheby);
void check_function(zomplex *an,zomplex *samp,index_st *ist,par_st *par, double ene_targets);


//Hmat.c
MKL_Complex16 dotp(zomplex *psi, double *phi,long m,long ngrid,double dv);
double dotpreal(zomplex *psi,double *phi,long m,long ngrid,double dv);

//nerror.c
void nerror(char *);

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

//save.c
void print_input_state(FILE *pf, flag_st *flag, grid_st *grid, par_st *par, index_st *ist, parallel_st *parallel);
void save_job_state(char *file_name, int checkpoint_id, double *psitot, double *pot_local, double *ksqr, zomplex *an, double *zn, double *ene_targets, long *nl,\
    nlc_st *nlc, grid_st *grid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
void restart_from_save(char *file_name, int checkpoint_id, double *psitot, double *pot_local, double *ksqr, zomplex *an, double *zn, double *ene_targets, long *nl,\
    nlc_st *nlc, grid_st *grid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
void save_output(char *file_name, double *psitot, double *eig_vals, double *sigma_E, xyz_st *R, grid_st *grid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);

// write.c
void write_cube_file(double *rho, grid_st *grid, char *fileName);
void write_separation(FILE *pf, char *top_bttm);
void write_state_dat(zomplex *psi, long n_elems, char* fileName);
void write_vector_dat(vector *vec, int n_elems, char* fileName);
void write_eval_dat(double* eig_vals, double* sigma_E, long n_elems, char* fileName);
void write_psi_dat(double* psitot, long n_elems, char* fileName);
