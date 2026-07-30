#include "init_periodic.h"

/*****************************************************************************/

void init_periodic(
    lattice_st *lattice,
    vector **G_vecs,
    vector **k_vecs,
    grid_st *grid,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel)
{

  /************************************************************/
  /*******************  DECLARE VARIABLES   *******************/
  /************************************************************/

  const int mpir = parallel->mpi_rank;

  /************************************************************/
  /*******************      READ INPUT      *******************/
  /************************************************************/

  if (mpir == 0)
    printf("\nReading lattice params from periodic_input.par:\n");

  read_periodic_input(lattice, ist, par, flag, parallel);

  /************************************************************/
  /****************   CALC RECIP LAT & G_VECS  ****************/
  /************************************************************/

  if (mpir == 0)
    printf("\nGenerating reciprocal lattice vectors:\n");

  gen_recip_lat_vecs(lattice, ist, par, flag, parallel);

  if (mpir == 0)
    printf("\nGenerating G vectors:\n");

  ALLOCATE(G_vecs, ist->ngrid, "G_vecs");

  gen_G_vecs(G_vecs, grid, ist, par, flag, parallel);

  if (mpir == 0)
    printf("  %d G vectors in the plane-wave basis\n", ist->n_G_vecs);

  write_vector_dat(*G_vecs, ist->n_G_vecs, "G_vecs.dat");

  /************************************************************/
  /*******************      GEN K GRID      *******************/
  /************************************************************/

  if (1 == flag->readKPath)
  {
    if (mpir == 0)
      printf("\nReading k-path from file kpath.par:\n");

    read_k_path(k_vecs, lattice, ist, par, flag, parallel);
  }
  else
  {
    if (mpir == 0)
      printf("\nGenerating %d x %d x %d k grid:\n", ist->nk1, ist->nk2, ist->nk3);

    ALLOCATE(k_vecs, ist->n_k_pts, "k_vecs");

    gen_k_vecs(k_vecs, lattice, ist, par, flag, parallel);
  }

  if (mpir == 0)
    printf("  n_kpts = %d n-G_vecs %d\n", ist->n_k_pts, ist->n_G_vecs);

  return;
}

/*****************************************************************************/

void setup_k_communicators(
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel)
{
  /*******************************************************************
   * Split MPI_COMM_WORLD into k-point groups for the periodic path.  *
   *                                                                  *
   * Each group owns a contiguous block of global k-points and        *
   * gathers its filtered states to its master (k_rank == 0) for      *
   * orthogonalization/diagonalization. Filtering stays MPI-parallel  *
   * over n_filter_cycles within each group.                          *
   *                                                                  *
   * Unified policy (NF = n_filter_cycles, nk = n_k_pts):             *
   *   group_size      = ranks per k-group                            *
   *   n_k_groups      = mpi_size / group_size                        *
   *   n_filters/rank  = NF / group_size                              *
   *   color           = world_rank / group_size  (group id)          *
   *   key             = world_rank % group_size  (rank in group)     *
   *                                                                  *
   *   mpi_size >= nk : group_size = mpi_size/nk, n_k_groups = nk,     *
   *                    each group owns 1 k (NF split over the group). *
   *   mpi_size <  nk : group_size = mpi_size, n_k_groups = 1, the     *
   *                    single group owns all nk k and loops over them *
   *                    (NF split over all ranks).                     *
   ********************************************************************/

  const int mpis = parallel->mpi_size;
  const int mpir = parallel->mpi_rank;
  const int nk = ist->n_k_pts;
  const long NF = ist->n_filter_cycles;

  int group_size;

  if (nk <= 0)
  {
    if (mpir == 0)
      fprintf(stderr, "ERROR: setup_k_communicators called with n_k_pts = %d\n", nk);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  // v1: checkpoint/restart is not supported in the k-grouped periodic path. The
  // checkpoint files carry no k dimension and assume a single global rank-0 block.
  // restartFromOrtho is the exception: it is a k-aware restart that reads per-k
  // psi-filt-[k].dat files (see run_periodic_restart_from_ortho in mod_diag.c).
  // read_input forces restartFromCheckpoint = 1 when restartFromOrtho is set so the
  // driver's switch jumps to its ortho stage, so allow that specific combination.
  if ((0 != flag->restartFromCheckpoint) && (1 != flag->restartFromOrtho))
  {
    if (mpir == 0)
      fprintf(stderr, "ERROR: checkpoint restart (restartFromCheckpoint=%d) is not "
                      "supported with periodic k-grouping\n",
              flag->restartFromCheckpoint);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if (1 == flag->saveCheckpoints && mpir == 0)
    printf("\tNOTE: saveCheckpoints is ignored on the periodic path (no k-aware checkpoints yet)\n");

  // Determine ranks-per-k-group from the two regimes
  if (mpis >= nk)
  {
    if (mpis % nk != 0)
    {
      if (mpir == 0)
        fprintf(stderr,
                "ERROR: mpi_size (%d) must be an integer multiple of n_k_pts (%d) "
                "when mpi_size >= n_k_pts\n",
                mpis, nk);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    group_size = mpis / nk;
  }
  else
  {
    group_size = mpis;
  }

  // NF must distribute evenly over the ranks of a k-group
  if (NF % group_size != 0)
  {
    if (mpir == 0)
      fprintf(stderr,
              "ERROR: ranks-per-k-group (%d) must divide n_filter_cycles (%ld)\n",
              group_size, NF);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  parallel->n_k_groups = mpis / group_size;
  parallel->k_size = group_size; // confirmed below from k_comm

  // Split the world communicator into k-groups
  const int color = mpir / group_size;
  const int key = mpir % group_size;
  parallel->k_color = color;

  MPI_Comm_split(MPI_COMM_WORLD, color, key, &parallel->k_comm);
  MPI_Comm_rank(parallel->k_comm, &parallel->k_rank);
  MPI_Comm_size(parallel->k_comm, &parallel->k_size);

  // Block-distribute the nk global k-points over the n_k_groups groups
  const int n_groups = parallel->n_k_groups;
  const int base = nk / n_groups;
  const int rem = nk % n_groups;
  parallel->n_k_local = base + (color < rem ? 1 : 0);
  parallel->k_global_start = color * base + (color < rem ? color : rem);

  if ((nk % n_groups != 0) && (mpir == 0))
    printf("\tWARNING: n_k_pts (%d) not divisible by n_k_groups (%d); "
           "k-points distributed unevenly (%d or %d per group)\n",
           nk, n_groups, base + 1, base);

  ALLOCATE(&(parallel->k_local_to_global), (parallel->n_k_local > 0 ? parallel->n_k_local : 1),
           "k_local_to_global");
  for (int ikl = 0; ikl < parallel->n_k_local; ikl++)
    parallel->k_local_to_global[ikl] = parallel->k_global_start + ikl;

  // Recompute the per-rank state counts for the grouped layout. read_input set these
  // assuming a single global pool of NF cycles over mpi_size ranks; under k-grouping
  // the NF cycles are split over the group_size ranks of each k-group instead.
  const long m = ist->m_states_per_filter;
  ist->n_filters_per_rank = NF / group_size;
  ist->n_states_per_rank = ist->n_filters_per_rank * m;
  ist->mn_states_per_k = NF * m; // states per k gathered across the group, before ortho
  ist->psi_rank_size = ist->n_states_per_rank * ist->nspinngrid * ist->complex_idx *
                       parallel->n_k_local;

  if (mpir == 0)
  {
    printf("\nk-point MPI grouping:\n");
    printf("\tmpi_size = %d, n_k_pts = %d, n_filter_cycles = %ld\n", mpis, nk, NF);
    printf("\tn_k_groups = %d, ranks per group (k_size) = %d\n", n_groups, group_size);
    printf("\tn_filters_per_rank = %ld, n_states_per_rank = %ld, mn_states_per_k = %ld\n",
           ist->n_filters_per_rank, ist->n_states_per_rank, ist->mn_states_per_k);
    fflush(stdout);
  }

  return;
}

/*****************************************************************************/

void gen_recip_lat_vecs(
    lattice_st *lattice,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel)
{

  // Calculate the reciprocal lattice vectors, b_i,
  // from the lattice vectors, a_i

  vector a1, a2, a3;
  vector a2xa3, a3xa1, a1xa2;
  double V_lat;

  a1 = lattice->a1;
  a2 = lattice->a2;
  a3 = lattice->a3;

  // b_1 = (2pi) * (a_2 x a_3) / [a_1 . (a_2 x a_3)]
  a2xa3 = retCrossProduct(a2, a3);
  V_lat = retDotProduct(a1, a2xa3);

  lattice->V_lat = V_lat;
  lattice->b1 = retScaledVector(a2xa3, TWOPI / V_lat);

  // b_2 = (2pi) * (a_3 x a_1) / [a_1 . (a_2 x a_3)]
  a3xa1 = retCrossProduct(a3, a1);
  lattice->b2 = retScaledVector(a3xa1, TWOPI / V_lat);

  // b_3 = (2pi) * (a_1 x a_2) / [a_1 . (a_2 x a_3)]
  a1xa2 = retCrossProduct(a1, a2);
  lattice->b3 = retScaledVector(a1xa2, TWOPI / V_lat);

  if (parallel->mpi_rank == 0)
  {
    printf("\tDirect lattice volume, V_lat = %f\n", V_lat);
    printf("\tb1 = %.4f %.4f %.4f\n", lattice->b1.x, lattice->b1.y, lattice->b1.z);
    printf("\tb2 = %.4f %.4f %.4f\n", lattice->b2.x, lattice->b2.y, lattice->b2.z);
    printf("\tb3 = %.4f %.4f %.4f\n", lattice->b3.x, lattice->b3.y, lattice->b3.z);
  }

  return;
}

/*****************************************************************************/

void gen_G_vecs(
    vector **G_vecs,
    grid_st *grid,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel)
{

  FILE *pf;
  double *gx, *gy, *gz;
  double G_max;
  int jx, jy, jz, jyz, jxyz;
  double KE_pref = sqr(HBAR) / (2 * MASS_E);
  // double gx_min, gy_min, gz_min;
  double scale_mag;

  if ((gx = (double *)calloc(grid->nx, sizeof(double))) == NULL)
    nerror("ksqrx");
  if ((gy = (double *)calloc(grid->ny, sizeof(double))) == NULL)
    nerror("ksqry");
  if ((gz = (double *)calloc(grid->nz, sizeof(double))) == NULL)
    nerror("ksqrz");

  // The kinetic energy is 0.5*|(k+G)|^2
  // Here, we generate and store the system-dependent vectors G
  init_kvec(gx, grid->nx, grid->dkx, 1);
  init_kvec(gy, grid->ny, grid->dky, 1);
  init_kvec(gz, grid->nz, grid->dkz, 1);
  // for (gx[0] = 0.0, jx = 1; jx <= grid->nx / 2; jx++)
  // {
  //   gx[jx] = (double)(jx)*grid->dkx * grid->ngrid_1;
  //   gx[grid->nx - jx] = -1.00 * gx[jx];
  // }

  // for (gy[0] = 0.0, jy = 1; jy <= grid->nx / 2; jy++)
  // {
  //   gy[jy] = (double)(jy)*grid->dky * grid->ngrid_1;
  //   gy[grid->ny - jy] = -1.00 * gy[jy];
  // }

  // for (gz[0] = 0.0, jz = 1; jz <= grid->nz / 2; jz++)
  // {
  //   gz[jz] = (double)(jz)*grid->dkz * grid->ngrid_1;
  //   gz[grid->nz - jz] = -1.00 * gz[jz];
  // }

  // printf("These are the grid and g vectors\n");
  // printf("\nnx gridx  g_x\n");
  // for (jx = 0; jx < grid->nx; jx++)
  // {
  //   printf("%d %lf %lf\n", jx, grid->x[jx], gx[jx]);
  // }
  // printf("\nny gridy  g_y\n");
  // for (jy = 0; jy < grid->ny; jy++)
  // {
  //   printf("%d %lf %lf\n", jy, grid->y[jy], gy[jy]);
  // }
  // printf("\nnz gridz  g_z\n");
  // for (jz = 0; jz < grid->nz; jz++)
  // {
  //   printf("%d %lf %lf\n", jz, grid->z[jz], gz[jz]);
  // }

  pf = fopen("G_vecs.dat", "w");
  // G_max = par->KE_max / KE_pref; // * grid->ngrid_1

  for (jz = 0; jz < grid->nz; jz++)
  {
    for (jy = 0; jy < grid->ny; jy++)
    {
      jyz = grid->nx * (grid->ny * jz + jy);
      for (jx = 0; jx < grid->nx; jx++)
      {
        jxyz = jyz + jx;
        (*G_vecs)[jxyz].x = gx[jx];
        (*G_vecs)[jxyz].y = gy[jy];
        (*G_vecs)[jxyz].z = gz[jz];
        (*G_vecs)[jxyz].mag = sqrt(sqr(gx[jx]) + sqr(gy[jy]) + sqr(gz[jz]));

        // if ((*G_vecs)[jxyz].mag > G_max)
        // {
        //   scale_mag = G_max / (*G_vecs)[jxyz].mag;
        //   (*G_vecs)[jxyz] = retScaledVector((*G_vecs)[jxyz], scale_mag);
        // }
        fprintf(pf, "%d %.4lg %.4lg %.4lg %.4lg\n", jxyz, (*G_vecs)[jxyz].x, (*G_vecs)[jxyz].y, (*G_vecs)[jxyz].z, (*G_vecs)[jxyz].mag);
      }
    }
  }

  ist->n_G_vecs = ist->ngrid;
  // printf("There are %d G vectors\n", ist->n_G_vecs);
  // printf("These are the G_vecs:\n");
  // for (i = 0; i < ist->n_G_vecs; i++){
  //     printf("  %lf %lf %lf mag = %lf\n", G_vecs[i].x, G_vecs[i].y, G_vecs[i].z, G_vecs[i].mag);
  // }

  return;
}

/*****************************************************************************/

void gen_k_vecs(
    vector **k_vecs,
    lattice_st *lattice,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel)
{
  // Generate an nk1 x nk2 x nk3 kgrid for the calculation
  int i;
  int n1, n2, n3;
  double k1_scale = TWOPI / lattice->a; // / ist->ngrid;
  double k2_scale = TWOPI / lattice->b; // / ist->ngrid;
  double k3_scale = TWOPI / lattice->c; // / ist->ngrid;
  // Initialize to 0.0 so that a singleton dimension (nk == 1, n = 0) yields k = 0
  // rather than 0 * (uninitialized) = NaN.
  double dk1 = 0.0, dk2 = 0.0, dk3 = 0.0;
  vector k;

  // Sample only the irreducible HALF of the Brillouin zone along each axis.
  // Time-reversal symmetry gives E(k) = E(-k), so the path [0, 0.5] (Gamma ->
  // zone boundary, in fractional reciprocal-lattice units) already determines
  // the full zone -- no need to sample [0, 1). Each axis is a linspace from 0.0
  // to 0.5 INCLUSIVE over nk points, so the spacing dk = 0.5/(nk-1) adapts to
  // the requested number of points (np.linspace, not np.arange). A singleton
  // axis (nk == 1) stays at Gamma (dk = 0, k = 0).
  if (ist->nk1 > 1)
  {
    dk1 = 0.5 / (double)(ist->nk1 - 1);
  }
  if (ist->nk2 > 1)
  {
    dk2 = 0.5 / (double)(ist->nk2 - 1);
  }
  if (ist->nk3 > 1)
  {
    dk3 = 0.5 / (double)(ist->nk3 - 1);
  }

  // Generate the nk1 x nk2 x nk3 kgrid (each axis spans fractional 0.0 -> 0.5)
  i = 0;
  k.x = k.y = k.z = k.mag = 0.0;
  for (n1 = 0; n1 < ist->nk1; n1++)
  {
    for (n2 = 0; n2 < ist->nk2; n2++)
    {
      for (n3 = 0; n3 < ist->nk3; n3++)
      {
        k.x = k1_scale * n1 * dk1;
        k.y = k2_scale * n2 * dk2;
        k.z = k3_scale * n3 * dk3;
        k.mag = retVectorMagnitude(k);
        (*k_vecs)[i] = k;
        i++;
      }
    }
  }

  if (i != ist->n_k_pts)
  {
    printf("ERROR: invalid k point generation!\n");
    fprintf(stderr, "ERROR: invalid k point generation!\n");
    exit(EXIT_FAILURE);
  }

  if (parallel->mpi_rank == 0)
  {
    printf("There are %d k vectors\n", ist->n_k_pts);
    printf("These are the k_vecs:\n");
    for (i = 0; i < ist->n_k_pts; i++)
    {
      printf("  %lf %lf %lf mag = %lf\n", (*k_vecs)[i].x, (*k_vecs)[i].y, (*k_vecs)[i].z, (*k_vecs)[i].mag);
    }
  }
  return;
}

/*****************************************************************************/

void read_k_path(
    vector **k_vecs,
    lattice_st *lattice,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel)
{
  // Read a k path from kpath.par in order to generate a bandstructure plot
  FILE *pf;
  int n_hsym, n_step, n_diff_k, n_k_pts = 0;
  int i, ik, d;
  int *n_steps_k;
  double kx, ky, kz, dk;
  vector *hsym_k;
  vector k_diff, u, lat_scale;

  // Open k-path.par and read in the high-symmetry k-point list
  if (access("kpath.par", F_OK) != -1)
  {
    pf = fopen("kpath.par", "r");
    fscanf(pf, "%d", &n_hsym);

    // Allocate memory for the high symmetry k-path
    if ((hsym_k = (vector *)calloc(n_hsym, sizeof(hsym_k[0]))) == NULL)
    {
      printf("ERROR: allocating memory for hsym_k\n");
      fprintf(stderr, "ERROR: allocating memory for hsym_k\n");
      exit(EXIT_FAILURE);
    }
    // Allocate memory for the number of steps between k-points
    if ((n_steps_k = (int *)calloc(n_hsym, sizeof(int))) == NULL)
    {
      printf("ERROR: allocating memory for n_steps_k\n");
      fprintf(stderr, "ERROR: allocating memory for n_steps_k\n");
      exit(EXIT_FAILURE);
    }

    // Read in the k_path vectors and number of steps betweenn kpoints
    i = 0;
    for (i = 0; i < n_hsym; i++)
    {
      fscanf(pf, "%lf %lf %lf %d", &kx, &ky, &kz, &n_step);
      hsym_k[i].x = kx;
      hsym_k[i].y = ky;
      hsym_k[i].z = kz;
      n_steps_k[i] = n_step;
      n_k_pts = n_k_pts + n_step;
      // printf("%d %f %f %f %d\n", i, hsym_k[i].x, hsym_k[i].y, hsym_k[i].z, n_steps_k[i]);
    }
  }
  else
  {
    printf("ERROR: no file kpath.par in directory!\n");
    fprintf(stderr, "ERROR: no file kpath.par in directory!\n");
    exit(EXIT_FAILURE);
  }
  // Set the number of k-points for the calculation
  ist->n_k_pts = n_k_pts;
  // Scale all the k vectors by their lattice parameters
  lat_scale.x = TWOPI / lattice->a;
  lat_scale.y = TWOPI / lattice->b;
  lat_scale.z = TWOPI / lattice->c;
  for (i = 0; i < n_hsym; i++)
  {
    hsym_k[i] = retElementWiseVectorMultiplication(hsym_k[i], lat_scale);
  }

  // printf("These are the scaled hsym points:\n");
  // for (i = 0; i < n_hsym; i++){
  //     printf("  %lf %lf %lf mag = %lf\n", hsym_k[i].x, hsym_k[i].y, hsym_k[i].z, hsym_k[i].mag);
  // }

  // Allocate memory for the k_vecs
  if ((*(k_vecs) = (vector *)calloc(ist->n_k_pts, sizeof(vector))) == NULL)
  {
    printf("ERROR allocating memory for k_vecs\n");
    fprintf(stderr, "ERROR allocating memory for k_vecs\n");
    exit(EXIT_FAILURE);
  }

  // We construct the k_vectors by following the path between high-sym k_points
  // That path is defined by the difference vectors between each of the k_points
  // u = k2 - k1

  n_diff_k = n_hsym - 1;
  i = 0;
  for (ik = 0; ik < n_diff_k; ik++)
  {
    k_diff = retSubtractedVectors(hsym_k[ik + 1], hsym_k[ik]);
    n_step = n_steps_k[ik];
    dk = k_diff.mag / n_step;
    k_diff = retScaledVector(k_diff, dk);
    for (d = 0; d < n_step; d++)
    {
      u = retScaledVector(k_diff, d);
      u = retAddedVectors(hsym_k[ik], u);
      (*k_vecs)[i] = u;
      i++;
    }
  }

  // Add on the final k_point that would be left behind from the loop
  (*k_vecs)[ist->n_k_pts - 1] = hsym_k[n_hsym - 1];

  // printf("There are %d k vectors\n", ist->n_k_pts);
  // printf("These are the k_vecs:\n");
  // for (i = 0; i < ist->n_k_pts; i++){
  //     printf("  %lf %lf %lf mag = %lf\n", (*k_vecs)[i].x, (*k_vecs)[i].y, (*k_vecs)[i].z, (*k_vecs)[i].mag);
  // }

  return;
}
