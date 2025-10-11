#include "init.h"

/*****************************************************************************/

void init_grid_params(grid_st *grid, xyz_st *R, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel)
{
  /*****************************************************************
   * This function initializes the parameters for building the grid *
   * The input geometry size is computed to ensure the number of    *
   * requested grid points is sufficient. Final grid parameters are *
   * updated and stored in *grid and redundantly in *ist-> The k^2  *
   * grid is also initialized.                                      *
   * inputs:                                                        *
   *   [grid_st *grid] pointer to grid struct, holds the grid       *
   *   [xyz_st *R] pointer to double arrays containing atom coords  *
   *   [index_st *ist] pointer to counters, indices, and lengths    *
   * outputs: void                                                  *
   ******************************************************************/

  const int mpir = parallel->mpi_rank;

  long ntmp;
  double xd, yd, zd;
  double boxl_x, boxl_y, boxl_z;

  char *X;
  char *Y;
  char *Z;
  X = malloc(2 * sizeof(X[0]));
  Y = malloc(2 * sizeof(Y[0]));
  Z = malloc(2 * sizeof(Z[0]));
  strcpy(X, "X");
  strcpy(Y, "Y");
  strcpy(Z, "Z");

  xd = rint(0.5 * calc_dot_dimension(R, ist->natoms, X) + 5.0);
  yd = rint(0.5 * calc_dot_dimension(R, ist->natoms, Y) + 5.0);
  zd = rint(0.5 * calc_dot_dimension(R, ist->natoms, Z) + 5.0);

  if (mpir == 0)
    printf("\tMin. required box dimension for each direction (Bohr):\n");
  if (mpir == 0)
    printf("\t-----------------------------------------------------\n");
  if (mpir == 0)
    printf("\txd = %.2lf yd = %.2lf zd = %.2lf\n", 2.0 * xd, 2.0 * yd, 2.0 * zd);

  /***initial parameters for the pot reduce mass, etc. in the x direction ***/
  grid->xmin = -xd;
  grid->xmax = xd - grid->dx; // There is an off-by-one in the grid due to the loop not including the final element. Create odd grid to have same grid points on each side.
  ntmp = (long)((grid->xmax - grid->xmin) / grid->dx);
  if (ntmp > grid->nx)
  {
    if (mpir == 0)
      printf("\tinput nx insufficient; updating parameter.\n");
    grid->nx = ntmp;

    // Do not let auto-grid generator construct odd grid
    // If nx is odd, add 1 more grid point to make it even
    if (grid->nx % 2)
    {
      grid->nx += 1;
    }
  }

  // Handle even and odd grid points so that there is always a (0.0, 0.0, 0.0) point
  if (grid->nx % 2 == 0)
  {
    boxl_x = (double)(grid->nx) * grid->dx;
    grid->xmin = -boxl_x / 2.0;
    grid->xmax = ((double)(grid->nx - 1) * grid->dx) / 2.0;
    grid->dkx = TWOPI / boxl_x; // reciprocal lattice vector length
  }
  else
  {
    // odd number of grid points
    boxl_x = (double)(grid->nx - 1) * grid->dx;
    grid->xmin = -boxl_x / 2.0;
    grid->xmax = boxl_x / 2.0;
    grid->dkx = TWOPI / boxl_x; // reciprocal lattice vector length
  }
  if (mpir == 0)
  {
    printf("\tThe x_min = %lg and x_max %lg\n", grid->xmin, grid->xmax);
  }
  /***initial parameters for the pot reduce mass, etc. in the y direction ***/
  grid->ymin = -yd;
  grid->ymax = yd - grid->dy;
  ntmp = (long)((grid->ymax - grid->ymin) / grid->dy);
  if (ntmp > grid->ny)
  {
    if (mpir == 0)
      printf("\tinput ny insufficient; updating parameter.\n");
    grid->ny = ntmp;

    // If nx is odd, add 1 more grid point to make it even
    if (grid->ny % 2)
    {
      grid->ny += 1;
    }
  }

  // Handle even and odd gridd points so that there is always a (0.0, 0.0, 0.0) point
  if (grid->ny % 2 == 0)
  {
    boxl_y = (double)(grid->ny) * grid->dy;
    grid->ymin = -boxl_y / 2.0;
    grid->ymax = ((double)(grid->ny - 1) * grid->dy) / 2.0;
    grid->dky = TWOPI / boxl_y; // reciprocal lattice vector length
  }
  else
  {
    // odd number of grid points
    boxl_y = (double)(grid->ny - 1) * grid->dy;
    grid->ymin = -boxl_y / 2.0;
    grid->ymax = boxl_y / 2.0;
    grid->dky = TWOPI / boxl_y; // reciprocal lattice vector length
  }

  if (mpir == 0)
  {
    printf("\tThe y_min = %lg and y_max %lg\n", grid->ymin, grid->ymax);
  }

  /***initial parameters for the pot reduce mass, etc. in the z direction ***/
  // Handle possible periodic z direction
  if (0 == flag->periodic)
  {
    grid->zmin = -zd;
    grid->zmax = zd - grid->dz;
    ntmp = (long)((grid->zmax - grid->zmin) / grid->dz);
    if ((ntmp > grid->nz) && (0 == flag->periodic))
    {
      if (mpir == 0)
        printf("\tinput nz insufficient; updating parameter.\n");
      grid->nz = ntmp;

      // If nx is odd, add 1 more grid point to make it even
      if (grid->nz % 2)
      {
        grid->nz += 1;
      }
    }

    // Handle even and odd gridd points so that there is always a (0.0, 0.0, 0.0) point
    if (grid->nz % 2 == 0)
    {
      boxl_z = (double)(grid->nz) * grid->dz;
      grid->zmin = -boxl_z / 2.0;
      grid->zmax = ((double)(grid->nz - 1) * grid->dz) / 2.0;
      grid->dkz = TWOPI / boxl_z; // reciprocal lattice vector length
    }
    else
    {
      // odd number of grid points
      boxl_z = (double)(grid->nz - 1) * grid->dz;
      grid->zmin = -boxl_z / 2.0;
      grid->zmax = boxl_z / 2.0;
      grid->dkz = TWOPI / boxl_z; // reciprocal lattice vector length
    }
    if (mpir == 0)
      printf("\tThe z_min = %lg and z_max %lg\n", grid->zmin, grid->zmax);
  }
  else if ((1 == flag->periodic) && (0.0 != par->box_z))
  {
    // Generate the box in the periodic z direction
    // Determine if the number of grid points chosen is sufficient for the max grid spacing
    if (mpir == 0)
      printf("\tUsing periodic box of length %lg Bohr in z dim\n", par->box_z);
    double dz;
    dz = par->box_z / grid->nz;
    if (dz < grid->dz)
    {
      grid->dz = dz;
    }
    else
    {
      // Increase the number of gridpts to have a grid spacing below grid->dz
      grid->nz = (int)(par->box_z / grid->dz) + 1;
      grid->dz = par->box_z / grid->nz;
      if (mpir == 0)
        printf("\tModify periodic box grid: nz = %ld dz = %lg\n", grid->nz, grid->dz);
    }

    // Handle even and odd gridd points so that there is always a (0.0, 0.0, 0.0) point
    if (grid->nz % 2 == 0)
    {
      boxl_z = (double)(grid->nz) * grid->dz;
      grid->zmin = -boxl_z / 2.0;
      grid->zmax = ((double)(grid->nz - 1) * grid->dz) / 2.0;
      grid->dkz = TWOPI / boxl_z; // reciprocal lattice vector length
    }
    else
    {
      // odd number of grid points
      boxl_z = (double)(grid->nz - 1) * grid->dz;
      grid->zmin = -boxl_z / 2.0;
      grid->zmax = boxl_z / 2.0;
      grid->dkz = TWOPI / boxl_z; // reciprocal lattice vector length
    }

    if (mpir == 0)
    {
      printf("\t-L/2 = %lg L/2 = %lg\n", grid->zmin, grid->zmax);
    }
  }
  else
  {
    if (mpir == 0)
    {
      printf("ERROR: periodic flag set but box_z = 0.0!\n");
      fprintf(stderr, "ERROR: periodic flag set but box_z = 0!\n");
    }
    exit(EXIT_FAILURE);
  }

  if ((grid->nx % 2) || (grid->ny % 2) || (grid->nz % 2))
  {
    if (mpir == 0)
    {
      printf("\n\nWARNING: ODD # OF GRID POINTS -> ENSURE (0.0, 0.0, 0.0) POINT EXISTS\n");
      printf("WARNING: OR GRID IS INCOMPATIBLE WITH HARTREE POT IN BSE\n\n");
    }
  }

  grid->nx_1 = 1.0 / (double)(grid->nx);
  grid->ny_1 = 1.0 / (double)(grid->ny);
  grid->nz_1 = 1.0 / (double)(grid->nz);
  grid->ngrid_1 = grid->nx_1 * grid->ny_1 * grid->nz_1;
  grid->dv = par->dv = grid->dx * grid->dy * grid->dz;
  grid->dr = sqrt(sqr(grid->dx) + sqr(grid->dy) + sqr(grid->dz));

  ist->ngrid = grid->ngrid = grid->nx * grid->ny * grid->nz;
  ist->nspinngrid = ist->nspin * ist->ngrid;
  ist->nx = grid->nx;
  ist->ny = grid->ny;
  ist->nz = grid->nz;
  ist->psi_rank_size = ist->n_states_per_rank * ist->nspinngrid * ist->complex_idx;

  if (mpir == 0)
  {
    printf("\n\tFinal grid parameters:\n");
    printf("\t----------------------\n");
    printf("\txmin = %lf ymin = %lf zmin = %lf\n", grid->xmin, grid->ymin, grid->zmin);
    printf("\tdx = %g dy = %g dz = %g dv = %g dr = %g\n", grid->dx, grid->dy, grid->dz, grid->dv, grid->dr);
    printf("\tnx = %ld  ny = %ld  nz = %ld\n", ist->nx, ist->ny, ist->nz);
    printf("\tngrid = %ld, nspin = %d, nspinngrid = %ld\n", ist->ngrid, ist->nspin, ist->nspinngrid);
  }

  free(X);
  free(Y);
  free(Z);
  fflush(stdout);

  return;
}

/*****************************************************************************/

void build_grid_ksqr(double *ksqr, xyz_st *R, grid_st *grid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel)
{
  /*******************************************************************
   * This function builds the real- and k-space grids.                *
   * inputs:                                                          *
   *  [ksqr] array of k^2 values to store k^2 grid                    *
   *  [grid] struct holding the grid and grid parameters (dx, xmin...)*
   *  [ist] ptr to counters, indices, and lengths                     *
   *  [par] ptr to par_st holding VBmin, VBmax... params              *
   *  [flag_st *par] ptr to flag_st holding job flags                 *
   * outputs: void                                                    *
   ********************************************************************/

  FILE *pf;
  long jx, jy, jz, jyz, jxyz, jtmp, jatom;
  double dx, dy, dz, *ksqrx, *ksqry, *ksqrz;
  double KE_max;

  // ****** ****** ****** ****** ****** ******
  // Building the grid
  // ****** ****** ****** ****** ****** ******
  if (parallel->mpi_rank == 0)
    printf("\tBuilding the grid...\n");

  for (jx = 0, dx = grid->xmin; jx < grid->nx; jx++, dx += grid->dx)
    grid->x[jx] = dx;
  for (jy = 0, dy = grid->ymin; jy < grid->ny; jy++, dy += grid->dy)
    grid->y[jy] = dy;
  for (jz = 0, dz = grid->zmin; jz < grid->nz; jz++, dz += grid->dz)
    grid->z[jz] = dz;

  FILE *pgrid;
  pgrid = fopen("grid.dat", "w");
  for (jx = 0; jx < ist->nx; jx++)
  {
    for (jy = 0; jy < ist->ny; jy++)
    {
      for (jz = 0; jz < ist->nz; jz++)
      {
        fprintf(pgrid, "%.6lg %.6lg %.6lg\n", grid->x[jx], grid->y[jy], grid->z[jz]);
      }
    }
  }
  fclose(pgrid);
  // ****** ****** ****** ****** ****** ******
  // Initializing the k^2 vectors
  // ****** ****** ****** ****** ****** ******
  if ((ksqrx = (double *)calloc(grid->nx, sizeof(double))) == NULL)
    nerror("ksqrx");
  if ((ksqry = (double *)calloc(grid->ny, sizeof(double))) == NULL)
    nerror("ksqry");
  if ((ksqrz = (double *)calloc(grid->nz, sizeof(double))) == NULL)
    nerror("ksqrz");

  // hold extra factor of 0.5 for the kinetic energy
  //  The kinetic energy is 0.5k^2
  init_kvec(ksqrx, grid->nx, grid->dkx, 1);
  init_kvec(ksqry, grid->ny, grid->dky, 1);
  init_kvec(ksqrz, grid->nz, grid->dkz, 1);
  for (jx = 1; jx < grid->nx; jx++)
  {
    ksqrx[jx] *= 0.5 * ksqrx[jx] * grid->ngrid_1;
  }
  for (jy = 1; jy < grid->ny; jy++)
  {
    ksqry[jy] *= 0.5 * ksqry[jy] * grid->ngrid_1;
  }
  for (jz = 1; jz < grid->nz; jz++)
  {
    ksqrz[jz] *= 0.5 * ksqrz[jz] * grid->ngrid_1;
  }

  pf = fopen("ksqr.dat", "w");
  KE_max = par->KE_max * (grid->ny_1 * grid->nx_1 * grid->nz_1);
  for (jz = 0; jz < grid->nz; jz++)
  {
    for (jy = 0; jy < grid->ny; jy++)
    {
      jyz = grid->nx * (grid->ny * jz + jy);
      for (jx = 0; jx < grid->nx; jx++)
      {
        jxyz = jyz + jx;
        ksqr[jxyz] = ksqrx[jx] + ksqry[jy] + ksqrz[jz];
        if (ksqr[jxyz] > KE_max)
          ksqr[jxyz] = KE_max; // KE cutoff at KE_max
        fprintf(pf, "%ld %lg\n", jxyz, ksqr[jxyz]);
      }
    }
  }
  fclose(pf);

  if ((1 == flag->NL) || (1 == flag->SO))
  {
    // Count the max number of grid points within Rnlcut of an atom***/
    if (parallel->mpi_rank == 0)
      printf("\tCount max no. grid points in Rnlcut of an atom\n");
    ist->n_NL_gridpts = 0;
    for (jatom = 0; jatom < ist->n_NL_atoms; jatom++)
    {
      for (jtmp = 0, jz = 0; jz < grid->nz; jz++)
      {
        dz = grid->z[jz] - R[jatom].z;

        // Implement periodic boundary conditions
        if (1 == flag->periodic)
        {
          // if dz > l/2, compute distance to periodic replica
          if (dz > par->box_z / 2)
            dz -= par->box_z;
          if (dz < -par->box_z / 2)
            dz += par->box_z;
        }

        for (jy = 0; jy < grid->ny; jy++)
        {
          dy = grid->y[jy] - R[jatom].y;
          for (jx = 0; jx < grid->nx; jx++)
          {
            dx = grid->x[jx] - R[jatom].x;

            if (dx * dx + dy * dy + dz * dz < par->R_NLcut2)
            {
              jtmp++;
            }
          }
        }
      }
      if (jtmp > ist->n_NL_gridpts)
        ist->n_NL_gridpts = jtmp;
    }
    if (parallel->mpi_rank == 0)
      printf("\tn_NL_gridpts = %ld\n", ist->n_NL_gridpts);
    fflush(stdout);
  }

  free(ksqrx);
  free(ksqry);
  free(ksqrz);
}

/*****************************************************************************/

void set_ene_targets(double *ene_targets, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel)
{
  /*****************************************************************
   * This function sets the energy targets for the filter cycles    *
   * These energy targets are the centers of the Gaussian filter    *
   * functions that determine the range of the Hamiltonian spectrum *
   * that the random initial state will converge toward.            *
   * inputs:                                                        *
   *   [double *ene_targets] ptr to double array;stores ene targets *
   *   [index_st *ist] ptr to counters, indices, and lengths        *
   *   [par_st *par] ptr to par_st holding VBmin, VBmax... params   *
   *   [flag_st *par] ptr to flag_st holding job flags              *
   * outputs: void                                                  *
   ******************************************************************/

  FILE *pf;
  double del;
  long jx;

  /*** setting the energy targets ***/
  // If the energy targets were set for VB and CB separately,
  // check that nVB + nCB is equal to the total number of energy targets (states per filter)
  if ((par->n_targets_VB + par->n_targets_CB) != ist->m_states_per_filter)
  {
    fprintf(stderr, "ERROR: n_targets_VB + n_targets_CB not equal to total m_states_per_filter!\n");
    exit(EXIT_FAILURE);
  }

  if (flag->setTargets != 1)
  {
    // If the ene targets were not set in input, place half in VB/CB
    par->n_targets_CB = (long)(ist->m_states_per_filter / 2);
    par->n_targets_VB = ist->m_states_per_filter - par->n_targets_CB;
  }
  // If either band has less than 1 ene target, crash
  if (par->n_targets_CB < 1 || par->n_targets_VB < 1)
  {
    fprintf(stderr, "ERROR: less than one energy target set in VB or CB\n");
    fflush(0);
    exit(EXIT_FAILURE);
  }

  // Compute the spacing between VB ene targets
  del = (par->VBmax - par->VBmin) / (double)par->n_targets_VB;
  if (parallel->mpi_rank == 0)
    printf("\tSpacing between states in VB: %lg a.u.\n", del);
  // Set the VB ene targets so that the VBmax energy is included
  for (jx = 0; jx < par->n_targets_VB; jx++)
  {
    ene_targets[par->n_targets_VB - jx - 1] = par->VBmax - (double)(jx)*del;
  }

  // Compute the spacing between CB ene targets
  del = (par->CBmax - par->CBmin) / (double)par->n_targets_CB;
  if (parallel->mpi_rank == 0)
    printf("\tSpacing between states in CB: %lg a.u.\n", del);
  // Set the CB ene targets so the CBmin energy is included
  for (jx = par->n_targets_VB; jx < ist->m_states_per_filter; jx++)
  {
    ene_targets[jx] = par->CBmin + (double)(jx - par->n_targets_VB) * del;
  }

  if (parallel->mpi_rank == 0)
    printf("\n\tEnergy target list:\n");
  if (parallel->mpi_rank == 0)
    printf("\t-------------------\n");
  // Print out the energy targets
  pf = fopen("ene_targets.dat", "w");
  for (jx = 0; jx < ist->m_states_per_filter; jx++)
  {
    fprintf(pf, "%g\n", ene_targets[jx]);
    if (parallel->mpi_rank == 0)
      printf("\t%g\n", ene_targets[jx]);
  }
  fclose(pf);

  return;
}

/*****************************************************************************/
void build_local_pot(double *pot_local, pot_st *pot, xyz_st *R, atom_info *atom,
                     grid_st *grid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel)
{
  /*******************************************************************
   * This function calculates the local potential at all gridpoints.  *
   * inputs:                                                          *
   *  [pot_local] ngrid long array, holds the local pot               *
   *  [pot_for_atom] pot_file_lens-long array of atom pseudopot. vals *
   *  [r_pot_file] pot_file_lens-long array of atom pseudopot. radii  *
   *  [pot_file_lens] array of pseudopot. file lens for each atom     *
   *  [dr] array of pseudopot. file radii spacings for each atom      *
   *  [R] array of atom cartesian coordinates                         *
   *  [ksqr] array of k^2 values to store k^2 grid                    *
   *  [atom] struct holding per-atom paramaters (Zval, SO_par, etc)   *
   *  [grid] struct holding the grid and grid parameters (dx, xmin...)*
   *  [ist] ptr to counters, indices, and lengths                     *
   *  [par] ptr to par_st holding VBmin, VBmax... params              *
   *  [flag] ptr to flag_st holding job flags                         *
   *  [parallel] holds options for parallelization                    *
   * outputs: void                                                    *
   ********************************************************************/

  // Declare local variables from structs (for readability/documentation)
  const int mpir = parallel->mpi_rank;
  const int nx = grid->nx;
  const int ny = grid->ny;
  const int nz = grid->nz;
  const int nat = ist->natoms;
  const int mpl = ist->max_pot_file_len;
  const int ipot = flag->interpolatePot;
  const int uS = flag->useStrain;
  const int pbc = flag->periodic;
  const int csi = ist->crystal_structure_int;
  const int omi = ist->outmost_material_int;

  const double ssCs = par->scale_surface_Cs;
  const double pcr2 = par->pot_cut_rad2;
  const double boxz = par->box_z;

  long jx, jy, jz, jyz, jxyz, jat;
  double dx, dy, dz, dist, dist2, sum;
  double strainF = 1.0;        // Default: No effect on potential
  int scaleLR = (ssCs != 1.0); // Default: if scale surface != 1.0, then scale LR

  // Notify if surface Cs atom charges are scaled
  if (scaleLR && mpir == 0)
  {
    printf("Surface Cs atom charges being scaled by %lg\n", ssCs);
  }

  // Allocate memory for strain parameters if needed
  vector *neigh_list = NULL;
  double *s_scale = NULL;
  double *vol_ref = NULL;

  if (uS)
  {
    neigh_list = (vector *)calloc(4 * nat, sizeof(vector));
    s_scale = (double *)calloc(nat, sizeof(double));
    vol_ref = (double *)calloc(nat, sizeof(double));
    pot->a4_params = (double *)calloc(ist->ngeoms * ist->n_atom_types, sizeof(double));
    pot->a5_params = (double *)calloc(ist->ngeoms * ist->n_atom_types, sizeof(double));

    if (!neigh_list || !s_scale || !vol_ref || !pot->a4_params || !pot->a5_params)
    {
      fprintf(stderr, "OUT OF MEMORY allocating strain\n");
      exit(EXIT_FAILURE);
    }
  }

  // Read atomic pseudopotentials
  if (mpir == 0)
    printf("\tReading atomic pseudopotentials...\n");

  read_pot(pot, R, atom, ist, par, flag, parallel);

  // Calculate strain_scale if needed
  if (uS)
  {
    read_nearest_neighbors(neigh_list, vol_ref, nat, csi, omi);
    calc_strain_scale(s_scale, neigh_list, vol_ref, atom, pot->a4_params, pot->a5_params, nat);
  }

  // Construct pseudopotential on grid
  if (mpir == 0)
    printf("\tConstructing total pseudopotential on the grid...\n");

  omp_set_dynamic(0);
  omp_set_num_threads(parallel->nthreads);

#pragma omp parallel for private(dx, dy, dz, dist, dist2, jy, jx, jyz, jxyz, sum, jat, strainF)
  for (jz = 0; jz < nz; jz++)
  {
    for (jy = 0; jy < ny; jy++)
    {
      jyz = nx * (ny * jz + jy);
      for (jx = 0; jx < nx; jx++)
      {
        jxyz = jyz + jx;
        sum = 0.0;

        for (jat = 0; jat < nat; jat++)
        {
          dx = grid->x[jx] - R[jat].x;
          dy = grid->y[jy] - R[jat].y;
          dz = grid->z[jz] - R[jat].z;

          // Apply periodic boundary conditions
          if (pbc)
          {
            if (dz > boxz / 2)
              dz -= boxz;
            if (dz < -boxz / 2)
              dz += boxz;
          }

          dist2 = dx * dx + dy * dy + dz * dz;

          // Skip grid points beyond potential cutoff radius
          // Does not skip any grid points if long range flag is on
          // if ((flag->LR == 0) && (dist2 > pcr2))
          // {
          //   continue;
          // }
          // Apply strain factor if needed
          strainF = uS ? s_scale[jat] : 1.0;
          dist = sqrt(dist2);

          // Interpolation handling
          int atyp_idx = atom[jat].idx;
          double gpar = atom[jat].geom_par;

          if (ipot)
          {
            if ((jxyz == 0) && (jat == 0) && (mpir == 0))
            {
              printf("\tComputing interpolated cubic/ortho potential\n");
            }

            // geom1 (cubic)
            sum += (1.0 - gpar) * interpolate(
                                      dist, pot->dr[2 * atyp_idx], pot->r, pot->r_LR, pot->pseudo, pot->pseudo_LR,
                                      mpl, pot->file_lens[2 * atyp_idx], 2 * atyp_idx,
                                      scaleLR, atom[jat].LR_par, strainF, flag->LR);
            // geom2 (orthorhombic)
            sum += gpar * interpolate(
                              dist, pot->dr[2 * atyp_idx + 1], pot->r, pot->r_LR, pot->pseudo, pot->pseudo_LR,
                              mpl, pot->file_lens[2 * atyp_idx + 1], 2 * atyp_idx + 1,
                              scaleLR, atom[jat].LR_par, strainF, flag->LR);
          }
          else
          {
            if (jxyz == 0 && jat == 0 && parallel->mpi_rank == 0)
              printf("\tComputing potential without interpolating over cubic/ortho parameters\n\n");

            sum += interpolate(
                dist, pot->dr[atyp_idx], pot->r, pot->r_LR, pot->pseudo, pot->pseudo_LR,
                mpl, pot->file_lens[atyp_idx], atyp_idx,
                scaleLR, atom[jat].LR_par, strainF, flag->LR);
          }
        }
        pot_local[jxyz] = pot_local[ist->ngrid + jxyz] = sum;
      }
    }
  }

  // Compute potential minimum
  par->Vmin = 1.0e10;
  par->Vmax = -1.0e10;
  for (jxyz = 0; jxyz < ist->ngrid; jxyz++)
  {
    if (par->Vmax < pot_local[jxyz])
      par->Vmax = pot_local[jxyz];
    if (par->Vmin > pot_local[jxyz])
      par->Vmin = pot_local[jxyz];
  }

  if (parallel->mpi_rank == 0)
    printf("\tVmin = %g Vmax = %g dV = %g \n", par->Vmin, par->Vmax, par->Vmax - par->Vmin);

  if (1 == flag->useStrain)
  {
    free(neigh_list);
    free(vol_ref);
    free(s_scale);
    free(pot->a4_params);
    free(pot->a5_params);
  }

  return;
}

/*****************************************************************************/
void init_SO_projectors(double *SO_projectors, grid_st *grid, xyz_st *R, atom_info *atom, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel)
{
  /*******************************************************************
   * This function initializes the spin-orbit pot projectors.         *
   * It also calculates the number of grid points within R_NL_cut2    *
   * inputs:                                                          *
   *  [SO_projectors] PROJ_LEN*nproj-long array holding projector vals*
   *  [grid] struct holding the grid and grid parameters              *
   *  [R] array of atom cartesian coordinates                         *
   *  [atom] struct holding per-atom paramaters (Zval, SO_par, etc)   *
   *  [ist] ptr to counters, indices, and lengths                     *
   *  [par] ptr to par_st holding VBmin, VBmax... params              *
   * outputs: void                                                    *
   ********************************************************************/

  FILE *pf;
  int rpoint;
  double dr, *vr;
  double proj;
  long N = PROJ_LEN;
  long jatom;
  int iproj;
  char fileName[50];
  // Compute a grid on which to defining the spherical harmonics and ang.mom. radial functions
  // name the grid vr, allocate memory
  if ((vr = (double *)calloc(N, sizeof(double))) == NULL)
  {
    nerror("mem_vr");
  }

  dr = sqrt(par->R_NLcut2) / ((double)N);
  for (rpoint = 0; rpoint < N; rpoint++)
  {
    vr[rpoint] = (double)rpoint * dr;
  }

  if (flag->readProj == 0)
  {
    // gen projectors on the fly
    gen_SO_projectors(grid->dx, sqrt(par->R_NLcut2), ist->nproj, SO_projectors, vr);
  }
  else if (flag->readProj == 1)
  {
    for (iproj = 0; iproj < ist->nproj; iproj++)
    {
      sprintf(fileName, "projectorSO_%d.dat", iproj);
      pf = fopen(fileName, "r");
      for (rpoint = 0; rpoint < N; rpoint++)
      {
        fscanf(pf, "%lg %lg", &(vr[0]), &proj);
        SO_projectors[iproj * N + rpoint] = proj;
      }
      fclose(pf);
    }
  }
  else
  {
    printf("ERROR: readProj value %d not recognized\n", flag->readProj);
    fprintf(stderr, "ERROR: readProj value %d not recognized\n", flag->readProj);
    exit(EXIT_FAILURE);
  }

  // Print the spin-orbit parameters for each atom for use in the coupling code
  for (jatom = 0; jatom < ist->n_NL_atoms; jatom++)
  {
    if (atom[jatom].SO_par != 0)
    {
      sprintf(&fileName[0], "SO_proj_const_%ld.dat", jatom);
      pf = fopen(fileName, "w");
      fprintf(pf, "%.10f", atom[jatom].SO_par);
      fclose(pf);
    }
  }

  free(vr);
  return;
}

/*****************************************************************************/
void init_NL_projectors(nlc_st *nlc, long *nl, double *SO_projectors, grid_st *grid, xyz_st *R, atom_info *atom, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel)
{
  /*******************************************************************
   * This function initializes the non local pot projectors.          *
   * inputs:                                                          *
   *  [nlc] struct holding nlc projector val @ each gridpt per-atom   *
   *  [nl] natom long array holding number of NL gridpts around atom  *
   *  [SO_projectors] PROJ_LEN*nproj-long array holding projector vals*
   *  [grid] struct holding the grid and grid parameters              *
   *  [R] array of atom cartesian coordinates                         *
   *  [atom] struct holding per-atom paramaters (Zval, SO_par, etc)   *
   *  [ist] ptr to counters, indices, and lengths                     *
   *  [par] ptr to par_st holding VBmin, VBmax... params              *
   * outputs: void                                                    *
   ********************************************************************/

  FILE *pf;
  long jatom;
  double dr, *vr, dr_proj;
  long rpoint, N = PROJ_LEN;

  // TODO NOTE: the current method for computing the NL potential relies on the SO_projectors being defined.
  // If there is no spin-orbit coupling used in the calculation, then the SO_projectors should be set to
  // the identity.

  // Useful constants
  double tmp1 = 0.5 * sqrt(3.0 / PIE);
  double tmp2 = 0.5 * sqrt(3.0 / TWOPI);

  // Compute a grid on which to defining the spherical harmonics and ang.mom. radial functions
  // name the grid vr, allocate memoRy
  if ((vr = (double *)calloc(N, sizeof(double))) == NULL)
  {
    nerror("nlc_vr");
  }

  dr = sqrt(par->R_NLcut2) / ((double)N);
  for (rpoint = 0; rpoint < N; rpoint++)
  {
    vr[rpoint] = (double)rpoint * dr;
  }
  dr_proj = vr[1];

  // Non-local piece
  // Find all the grid points within par->Rnlcut of each atom and calculate
  // r, r2, y1[1+m], proj(r), etc. at the grid points and store the results in nlc
  /*** for the nonlocal potential ***/

  // 2. Calculate r, r2, y1[1+m], proj(r), etc. at the grid points
  omp_set_dynamic(0);
  omp_set_num_threads(parallel->nthreads);
#pragma omp parallel for private(jatom)
  for (jatom = 0; jatom < ist->n_NL_atoms; jatom++)
  {
    long jx, jy, jz, jyz, jxyz;
    int iproj, rpoint, *sgnProj;
    double dx, dy, dz, dxeps, dyeps, dzeps, dr_1, dr2;
    double *nlcprojectors;

    if ((nlcprojectors = (double *)calloc(N * ist->nproj, sizeof(double))) == NULL)
    {
      nerror("nlc_projector");
    }
    if ((sgnProj = (int *)calloc(ist->nproj, sizeof(int))) == NULL)
    {
      nerror("nlc_sgnProj");
    }

    if (0 == flag->readProj)
    {
      // gen projectors on the fly

      // generate the nonlocal part for each atom
      gen_nlc_projectors(grid->dx, sqrt(par->R_NLcut2), ist->nproj, nlcprojectors, sgnProj, vr, atom, jatom);
    }
    else if (1 == flag->readProj)
    {
      char projNL_file[50], sgnNL_file[50];
      FILE *pNL, *pSgn;
      double eig, scratch;
      int tmp;

      // Read projector info from files
      // NL projectors
      for (iproj = 0; iproj < ist->nproj; iproj++)
      {
        sprintf(projNL_file, "projectorNL_%ld_%d.dat", jatom, iproj);
        pNL = fopen(projNL_file, "r");
        if (pNL != NULL)
        {
          for (rpoint = 0; rpoint < PROJ_LEN; rpoint++)
          {
            fscanf(pNL, "%lf %lf\n", &scratch, &nlcprojectors[PROJ_LEN * iproj + rpoint]);
          }
          fclose(pNL);
        }
        else
        {
          if (55 != atom[jatom].Zval)
          {
            // If any atom not a Cs atom is missing its file, this is BAD
            printf("Cannot open %s file %s\n", atom[jatom].atyp, projNL_file);
            fprintf(stderr, "Cannot open %s file %s\n", atom[jatom].atyp, projNL_file);
            exit(EXIT_FAILURE);
          }
        }
      }

      // NL_projector signs
      sprintf(sgnNL_file, "NL_Proj_Eigs%ld-sorted.dat", jatom);
      pSgn = fopen(sgnNL_file, "r");
      if (pSgn != NULL)
      {
        for (iproj = 0; iproj < ist->nproj; iproj++)
        {
          fscanf(pSgn, "%d %lf %lf\n", &tmp, &scratch, &eig);
          if (eig < 0.0)
          {
            sgnProj[iproj] = -1;
          }
          else
          {
            sgnProj[iproj] = 1;
          }
        }
        fclose(pSgn);
      }
      else
      {
        if (55 != atom[jatom].Zval)
        {
          // If any atom not a Cs atom is missing its file, this is BAD
          printf("Cannot open %s file %s\n", atom[jatom].atyp, sgnNL_file);
          fprintf(stderr, "Cannot open %s file %s\n", atom[jatom].atyp, sgnNL_file);
          exit(EXIT_FAILURE);
        }
      }
    }

    nl[jatom] = 0;
    for (jz = 0; jz < grid->nz; jz++)
    {
      dz = grid->z[jz] - R[jatom].z;

      // Implement periodic boundary conditions
      if (1 == flag->periodic)
      {
        // if dz > l/2, compute distance to periodic replica
        if (dz > par->box_z / 2)
          dz -= par->box_z;
        if (dz < -par->box_z / 2)
          dz += par->box_z;
      }

      dzeps = dz + EPSDX;
      for (jy = 0; jy < grid->ny; jy++)
      {
        jyz = grid->nx * (grid->ny * jz + jy);
        dy = grid->y[jy] - R[jatom].y;
        dyeps = dy + EPSDX;
        for (jx = 0; jx < grid->nx; jx++)
        {
          jxyz = jyz + jx;
          dx = grid->x[jx] - R[jatom].x;
          dxeps = dx + EPSDX;
          dr2 = dx * dx + dy * dy + dz * dz;
          if (dr2 < par->R_NLcut2)
          {
            nlc[jatom * ist->n_NL_gridpts + nl[jatom]].jxyz = jxyz;

            nlc[jatom * ist->n_NL_gridpts + nl[jatom]].r = sqrt(dr2);

            dr_1 = 1.0 / sqrt(dx * dx + dy * dy + dzeps * dzeps);
            nlc[jatom * ist->n_NL_gridpts + nl[jatom]].y1[1].re = tmp1 * dzeps * dr_1;
            nlc[jatom * ist->n_NL_gridpts + nl[jatom]].y1[1].im = 0.0;

            dr_1 = 1.0 / sqrt(dxeps * dxeps + dy * dy + dz * dz);
            nlc[jatom * ist->n_NL_gridpts + nl[jatom]].y1[2].re = -tmp2 * dxeps * dr_1;
            nlc[jatom * ist->n_NL_gridpts + nl[jatom]].y1[0].re = tmp2 * dxeps * dr_1;

            dr_1 = 1.0 / sqrt(dx * dx + dyeps * dyeps + dz * dz);
            nlc[jatom * ist->n_NL_gridpts + nl[jatom]].y1[2].im = -tmp2 * dyeps * dr_1;
            nlc[jatom * ist->n_NL_gridpts + nl[jatom]].y1[0].im = -tmp2 * dyeps * dr_1;

            // write projectors to nlc struct and scale projectors by the SO scaling for this atom
            for (iproj = 0; iproj < ist->nproj; iproj++)
            {
              nlc[jatom * ist->n_NL_gridpts + nl[jatom]].proj[iproj] =
                  interpolate(sqrt(dr2), dr_proj, vr, NULL, &SO_projectors[N * iproj], NULL, 0, N, 0, 0, 1.0, 1.0, 0);
              // scale projectors by the SO scaling for this atom
              nlc[jatom * ist->n_NL_gridpts + nl[jatom]].proj[iproj] *= sqrt(atom[jatom].SO_par);
              // fprintf(pf, "%li %i %f %f %f\n", jatom, iproj,nlc[jatom*ist->n_NL_gridpts + nl[jatom]].r,nlc[jatom*ist->n_NL_gridpts + nl[jatom]].proj[iproj], sqrt(atom[jatom].SO_par) );

              nlc[jatom * ist->n_NL_gridpts + nl[jatom]].NL_proj[iproj] =
                  interpolate(sqrt(dr2), dr_proj, vr, NULL, &nlcprojectors[N * iproj], NULL, 0, N, 0, 0, 1.0, 1.0, 0);
              nlc[jatom * ist->n_NL_gridpts + nl[jatom]].NL_proj_sign[iproj] = sgnProj[iproj];
            }

            if (dr2 > EPSDX)
            {
              nlc[jatom * ist->n_NL_gridpts + nl[jatom]].r2_1 = sqr(dr_1);
              nlc[jatom * ist->n_NL_gridpts + nl[jatom]].r2 = dr2;
            }
            else
            {
              nlc[jatom * ist->n_NL_gridpts + nl[jatom]].r2_1 = 0.0;
              nlc[jatom * ist->n_NL_gridpts + nl[jatom]].r2 = 1.0 / EPSDX;
            }
            nl[jatom]++;
          }
        }
      }
    }
    free(nlcprojectors);
    free(sgnProj);
  }

  pf = fopen("list_NL_grid.dat", "w");
  for (jatom = 0; jatom < ist->n_NL_atoms; jatom++)
  {
    fprintf(pf, "%ld %ld\n", jatom, nl[jatom]);
  }
  fclose(pf);

  free(vr);

  return;
}

/***************************************************************************/
void init_filter_states(
    double *psi_rank,
    zomplex *psi,
    grid_st *grid,
    long *rand_seed,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel)
{

  FILE *pseed;
  FILE *ppsi;
  char str[20];

  int jspin;
  long jns;
  long jms;
  long jstate;
  long jgrid;
  long jgrid_real;
  long jgrid_imag;
  long ns_block;
  long stlen = ist->complex_idx * ist->nspinngrid;

  // If the user wants to read in starting wavefunctions for the filter
  // Then read them in and do not randomly generate states
  // There must be at least ist->mn_states_tot states available in the file
  // psi-init.dat

  if (1 == flag->inputPsiFilt)
  {
    int n_every;
    n_every = (1 == flag->useSpinors) ? 2 : 1;

    ppsi = fopen("psi-init.dat", "r");

    read_psi_dat(ppsi, psi_rank, ist->init_psi_start, ist->init_psi_end, stlen, n_every);

    // For debugging, print out cube files of the psi's that were read in:
    printf("Writing cube output\n");
    fflush(0);
    double *rho = (double *)calloc(ist->ngrid, sizeof(double));
    printf("Allocated rho\n");
    fflush(0);
    for (jns = 0; jns < ist->mn_states_tot; jns++)
    {
      printf("Reading state %ld\n", jns);
      fflush(0);
      jstate = jns * stlen;
      for (jgrid = 0; jgrid < ist->ngrid; jgrid++)
      {
        jgrid_real = ist->complex_idx * jgrid;
        jgrid_imag = jgrid_real + 1;

        rho[jgrid] = sqr(psi_rank[jstate + jgrid_real]) + sqr(psi_rank[jstate + jgrid_imag]);
      }
      sprintf(str, "psi-init-%ld.cube", jns);
      write_cube_file(rho, grid, str);
    }
    printf("Done writing cubes\n");
    fflush(0);

    free(rho);

    fclose(ppsi);
    return;
  }

  // Print out the random seeds for debugging purposes
  sprintf(str, "seed-%d.dat", parallel->mpi_rank);
  pseed = fopen(str, "w");
  fprintf(pseed, "idx  seed\n");

  // Initialize n_filters_per_rank random states on each mpi rank
  for (jns = 0; jns < ist->n_filters_per_rank; jns++)
  {

    // Initialize a random wavefunction for filtering
    for (jspin = 0; jspin < ist->nspin; jspin++)
    {
      fprintf(pseed, "%ld %ld\n", ist->nspin * jns + jspin, *rand_seed);
      init_psi(&psi[jspin * ist->ngrid], rand_seed, grid, ist, par, flag, parallel);
    }

    ns_block = jns * (ist->m_states_per_filter * stlen);

    // Add the initial random psi to the start of the ns_block for psi_rank
    if (1 == flag->isComplex)
    {
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++)
      {
        // handle indexing of real and imaginary components if complex
        // ist->complex_idx = (flag->isComplex + 1) = 2 when complex valued functions are in use, 1 if real
        jgrid_real = ist->complex_idx * jgrid;
        jgrid_imag = jgrid_real + 1;

        psi_rank[ns_block + jgrid_real] = psi[jgrid].re;
        psi_rank[ns_block + jgrid_imag] = psi[jgrid].im;
        // the imaginary components will be stored one double away (16 bytes, or IMAG_IDX) in memory from the real component
      }
    }
    else
    {
      for (jgrid = 0; jgrid < ist->nspinngrid; jgrid++)
      {
        // if the wavefunction is real-valued, then only the real component is stored
        psi_rank[ns_block + jgrid] = psi[jgrid].re;
      }
    }
  }

  fclose(pseed);

  return;
}

/*****************************************************************************/
void init_psi(zomplex *psi, long *rand_seed, grid_st *grid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel)
{
  /*******************************************************************
   * This function initializes a random, normalized wavefnc on grid   *
   * inputs:                                                          *
   *  [psi] ngrid-long array of double/zomplex for storing wavefnc    *
   *  [rand_seed] ptr to value of random seed; value reset on return  *
   *  [grid] struct containing grid dimension values                  *
   *  [parallel] struct holding parallelization options (for norm)    *
   * outputs: void                                                    *
   ********************************************************************/

  long jx, jy, jz, jzy, jxyz;
  long randint = (*rand_seed);
  // Loop over entire grid to set new values at all grid points
  for (jz = 0; jz < grid->nz; jz++)
  {
    for (jy = 0; jy < grid->ny; jy++)
    {
      jzy = grid->nx * (grid->ny * jz + jy);
      for (jx = 0; jx < grid->nx; jx++)
      {
        jxyz = jzy + jx;
        // Initialize wavefunction value at this grid point to a
        // random number between [-1.0,1.0]
        // ran_nrc generates random between [0.0,1.0] and resets the seed
        psi[jxyz].re = (-1.0 + 2.0 * ran_nrc(&randint));
        // If using complex-valued functions, then initialize a random value for imag component
        if (1 == flag->isComplex)
        {
          psi[jxyz].im = (-1.0 + 2.0 * ran_nrc(&randint));
        }
        else if (0 == flag->isComplex)
        {

          // otherwise set imaginary component to 0.0
          psi[jxyz].im = 0.0;
        }
      }
    }
  }

  // normalize this wavefunction and set the value of rand_seed to the new
  // seed so the next wavefunction is different.
  normalize(psi, ist->ngrid, ist, par, flag, parallel);

  (*rand_seed) = randint;

  return;
}

/***************************************************************************/

double calc_dot_dimension(xyz_st *R, long n_atoms, char *dir)
{
  /*******************************************************************
   * This function calculates the size of the NC along one dimension  *
   * inputs:                                                          *
   *  [R] array of coordinates along one dimension                    *
   *  [n_atoms] number of atoms in geometry                           *
   * outputs: [double] distance between atoms in a.u.                 *
   ********************************************************************/

  long i, j;
  double dr2, dis2;
  char *X;
  X = malloc(2 * sizeof(X[0]));
  char *Y;
  Y = malloc(2 * sizeof(Y[0]));
  char *Z;
  Z = malloc(2 * sizeof(Z[0]));

  strcpy(X, "X");
  strcpy(Y, "Y");
  strcpy(Z, "Z");

  if (0 == strcmp(dir, (const char *)X))
  {
    dr2 = 0.0;
    for (i = 0; i < n_atoms - 1; i++)
    {
      for (j = i + 1; j < n_atoms; j++)
      {
        dis2 = sqr(R[i].x - R[j].x);
        if (dis2 > dr2)
          dr2 = dis2;
      }
    }
  }

  if (0 == strcmp(dir, "Y"))
  {
    dr2 = 0.0;
    for (i = 0; i < n_atoms - 1; i++)
    {
      for (j = i + 1; j < n_atoms; j++)
      {
        dis2 = sqr(R[i].y - R[j].y);
        if (dis2 > dr2)
          dr2 = dis2;
      }
    }
  }

  if (0 == strcmp(dir, "Z"))
  {
    dr2 = 0.0;
    for (i = 0; i < n_atoms - 1; i++)
    {
      for (j = i + 1; j < n_atoms; j++)
      {
        dis2 = sqr(R[i].z - R[j].z);
        if (dis2 > dr2)
          dr2 = dis2;
      }
    }
  }

  free(X);
  free(Y);
  free(Z);
  return sqrt(dr2);
}

/*****************************************************************************/

// k-vector initialization (works for both even and odd n)
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
}

/*****************************************************************************/
