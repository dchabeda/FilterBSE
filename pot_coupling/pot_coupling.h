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
  int SO, NL, LR, useSpinors, isComplex;
  int printPsiFilt, printOrtho, printNorm, printCubes;
  int calcPotOverlap, getAllStates, timeHamiltonian, calcSpinAngStat;
  int retryFilter, alreadyTried, saveCheckpoints, restartFromCheckpoint, saveOutput;
  int readGrid, useGaussianBasis;

  int printFPDensity;
  int calcDarkStates;
  int timingSpecs;
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

  long max_elec_states, max_hole_states;
  long n_qp, n_xton;
  long n_holes, n_elecs;
  long *eval_hole_idxs, *eval_elec_idxs;
  double ngrid_1;
  int n_FP_density;
  int printFPDensity; // 0 = False (default) or 1 = True
  int calcDarkStates; // 0 = False (default) or 1 = True

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

  double dx, dy, dz, dr, dkx, dky, dkz, epsX, epsY, epsZ;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double delta_E_elec, delta_E_hole;
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

#define N_MAX_ATOM_TYPES 20
#define PROJ_LEN	1024

/*****************************************************************************/
// Function declerations

// save
void print_input_state(FILE *pf, flag_st *flag, grid_st *grid, par_st *par, index_st *ist, parallel_st *parallel);
void read_filter_output(char *file_name, double **psitot, double **eig_vals, double **sigma_E, xyz_st **R, grid_st *grid, double **gridx, double **gridy, double **gridz, index_st *ist, par_st *par, flag_st *flag);

// read
void read_input(flag_st *flag, grid_st *grid, index_st *ist, par_st *par, parallel_st *parallel);
void read_conf(char *file_name, xyz_st *R, atom_info *atm, index_st *ist, par_st *par, flag_st *flag);
void read_pot(pot_st *pot, xyz_st *R, atom_info *atom, index_st *ist, par_st *par, flag_st *flag);
void read_pot_file(FILE *pf, pot_st *pot, long j, long n, char *req);
void interpolate_pot(xyz_st *R, atom_info *atom, index_st *ist, par_st *par);
void calc_geom_par(xyz_st *R,atom_info *atm, index_st *ist );
double calc_bond_angle(long index1,long index2,long index3, xyz_st *R);
long assign_atom_number(char atyp[4]);
void assign_atom_type(char *atype,long j);
long get_number_of_atom_types(atom_info *atm,index_st *ist, long *list);
int assign_crystal_structure(char *crystal_structure);
int assign_outmost_material(char *outmost_material);
double get_ideal_bond_len(long natyp_1, long natyp_2, int crystalStructureInt);

//strain 
void read_nearest_neighbors(vector *atom_neighbors, double *tetrahedron_vol_ref, long natoms, int crystal_structure, int outmost_material);
void calc_strain_scale(double *strain_scale, vector *atom_neighbors, double *tetrahedron_vol_ref, atom_info *atom, double *a4_params, double *a5_params, long natoms);
double calc_regular_tetrahedron_volume(double bond_length1, double bond_length2, double bond_length3, double bond_length4);

// interpolate
double interpolate(double r,double dr,double *vr,double *vr_LR,double *pot,double *pot_LR,long pot_file_len,long n,long j, int scale_LR, double scale_LR_par, double strain_factor, int is_LR);

// potential
void potential(zomplex *phi, zomplex *psi, double *pot_local, nlc_st *nlc, long *nl, index_st *ist, par_st *par, flag_st *flag);
void spin_orbit_proj_pot(zomplex *phi, zomplex *psi, nlc_st *nlc, long* nl, index_st *ist, par_st *par);
void nonlocal_proj_pot(zomplex *phi, zomplex *psi, nlc_st *nlc, long* nl, index_st *ist, par_st *par);

//
double dotpreal(zomplex *psi,double *phi,long m,long ngrid,double dv);

// basis
void get_qp_basis_indices(double *eig_vals, double *sigma_E, long **eval_hole_idxs, long **eval_elec_idxs, index_st *ist, par_st *par, flag_st *flag);
void get_qp_basis(double *psi, double *psitot, double *psi_hole, double *psi_elec, double *eig_vals, double *sigma_E, index_st *ist, par_st *par, flag_st *flag);

// init
void build_local_pot(double *pot_local, pot_st *pot, xyz_st *R, atom_info *atom, xyz_st *grid, grid_st *grid_par,
    index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
void init_SO_projectors(double *SO_projectors, xyz_st *R, atom_info *atom, grid_st *grid, index_st *ist, par_st *par);
void init_NL_projectors(nlc_st *nlc, long *nl, double *SO_projectors, xyz_st *R, atom_info *atom, xyz_st *grid, grid_st *grid_par, index_st *ist, par_st *par, flag_st *flag);

// write
void write_cube_file(double *rho, grid_st *grid, char *fileName);
void write_separation(FILE *pf, char *top_bttm);

// nerror
void nerror(char *);
