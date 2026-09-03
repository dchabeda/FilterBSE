#include "mod_init.h"

void mod_init(
    double complex **psi_qp,
    double **eig_vals,
    double **sigma_E,
    xyz_st **R,
    grid_st *grid,
    double **gridx,
    double **gridy,
    double **gridz,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel)
{

  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  long jg;

  const int mpir = parallel->mpi_rank;

  /************************************************************/
  /********************     READ INPUT     ********************/
  /************************************************************/

  if (mpir == 0)
  {
    write_separation(stdout, "T");
    printf("\n1.\tINITIALIZING JOB | %s\n", get_time());
    write_separation(stdout, "B");
    fflush(stdout);
  }

  if (0 == flag->initUnsafe)
  {
    if (mpir == 0)
      printf("\nReading filter output from output.dat:\n");
    fflush(stdout);
    read_filter_output("output.dat", eig_vals, sigma_E, R, grid, gridx, gridy, gridz, ist, par, flag, parallel);
  }
  else if (1 == flag->initUnsafe)
  {
    if (mpir == 0)
      printf("\nReading UNSAFE input from unsafe_input.par:\n");
    fflush(stdout);
    read_unsafe_input(
        psi_qp, eig_vals, sigma_E, R, grid, gridx, gridy, gridz,
        ist, par, flag, parallel);
  }

  /*** Read initial setup from input.par ***/
  if (mpir == 0)
    printf("\nReading BSE job specifications from input.par:\n");
  fflush(stdout);

  read_input(flag, grid, ist, par, parallel);

  /************************************************************/
  /*******************   GET QP BASIS IDXS   ******************/
  /************************************************************/

  if (mpir == 0)
    printf("\nSetting quasiparticle basis indices:\n");
  fflush(stdout);

  // Allocate memory for the lists of the indices of eigenstates
  // The maximum possible number of hole states is mn_states_tot from filter.
  // We allocate an entire block of that size for both elecs and holes
  // because we will reallocate after the get_qp_basis_indices function.

  ALLOCATE(&(ist->eval_hole_idxs), ist->mn_states_tot, "eval_hole_idxs");
  ALLOCATE(&(ist->eval_elec_idxs), ist->mn_states_tot, "eval_elec_idxs");

  get_qp_basis_indices(*eig_vals, *sigma_E, &ist->eval_hole_idxs, &ist->eval_elec_idxs, ist, par, flag, parallel);

  /************************************************************/
  /********************    BUILD QP BASIS    ******************/
  /************************************************************/

  if (mpir == 0)
    printf("\nReading quasiparticle basis wavefunctions:\n");
  fflush(stdout);

  // Allocate psi_qp in NODE-SHARED memory (MPI-3 shared window): one physical
  // copy per node, shared read-only by all ranks on the node. Each state is
  // ~0.5 GB, so the old per-rank copy (N per node) OOMs at large basis sizes;
  // one copy per node removes that duplication. Only node-rank 0 requests the
  // bytes and reads the wavefunctions from disk; every rank queries the same
  // base pointer. The compute pipeline and GPU tiling are unchanged -- they read
  // psi_qp exactly as before (it is written once here, then read-only).
  {
    const MPI_Aint nelem = (MPI_Aint)ist->n_qp * (MPI_Aint)ist->nspinngrid;
    const MPI_Aint bytes =
        (parallel->node_rank == 0) ? nelem * (MPI_Aint)sizeof(double complex) : 0;
    MPI_Aint seg_bytes;
    int disp_unit;
    double complex *base = NULL;
    MPI_Win_allocate_shared(bytes, (int)sizeof(double complex), MPI_INFO_NULL,
                            parallel->node_comm, &base, &parallel->psi_win);
    // Every rank resolves node-rank 0's segment -> the single shared buffer.
    MPI_Win_shared_query(parallel->psi_win, 0, &seg_bytes, &disp_unit, &base);
    *psi_qp = base;
    MPI_Win_lock_all(MPI_MODE_NOCHECK, parallel->psi_win);
  }

  // Only node-rank 0 fills the shared buffer from disk; all node ranks then
  // synchronize so the freshly written states are visible before any read.
  if (parallel->node_rank == 0)
    get_qp_basis(*psi_qp, *eig_vals, *sigma_E, ist, par, flag);
  MPI_Win_sync(parallel->psi_win);
  MPI_Barrier(parallel->node_comm);
  MPI_Win_sync(parallel->psi_win);

  ist->n_xton = ist->n_elecs * ist->n_holes;

  /************************************************************/
  /********************    SET LOOP IDXS     ******************/
  /************************************************************/
  // ALLOCATE(&ist->jsg_conv, ist->nspinngrid, "jsg");

  // #pragma omp simd safelen(4)
  // for (int s = 0; s < 2; s++){
  //   for (jg = 0; jg < ist->ngrid; jg++){
  //     long jsg = jg + s * ist->ngrid;
  //     ist->jsg_conv[jsg] = jg;
  //   }
  // }

  free(*sigma_E);
  sigma_E = NULL;
  free(ist->eval_elec_idxs);
  ist->eval_elec_idxs = NULL;
  free(ist->eval_hole_idxs);
  ist->eval_hole_idxs = NULL;

  return;
}