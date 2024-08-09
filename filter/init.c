#include "fd.h"
#include "vector.h"

/*****************************************************************************/

void init_grid_params(grid_st *grid, xyz_st *R, index_st *ist, par_st *par, flag_st *flag){
  /*****************************************************************
  * This function initializes the parameters for building the grid *
  * The input geometry size is computed to ensure the number of    *
  * requested grid points is sufficient. Final grid parameters are *
  * updated and stored in *grid and redundantly in *ist. The k^2   *
  * grid is also initialized.                                      *
  * inputs:                                                        *
  *   [grid_st *grid] pointer to grid struct, holds the grid       *
  *   [xyz_st *R] pointer to double arrays containing atom coords  *
  *   [index_st *ist] pointer to counters, indices, and lengths    *
  * outputs: void                                                  *
  ******************************************************************/

  long ntmp; 
  double xd, yd, zd;
  char *X; X = malloc(2*sizeof(X[0])); strcpy(X, "X");
  char *Y; Y = malloc(2*sizeof(Y[0])); strcpy(Y, "Y");
  char *Z; Z = malloc(2*sizeof(Z[0])); strcpy(Z, "Z");

  xd = rint(0.5 * calc_dot_dimension(R, ist->natoms, X) + 5.0);
  yd = rint(0.5 * calc_dot_dimension(R, ist->natoms, Y) + 5.0);
  zd = rint(0.5 * calc_dot_dimension(R, ist->natoms, Z) + 5.0);
  
  printf("\tMin. required box dimension for each direction (Bohr):\n");
  printf("\t-----------------------------------------------------\n");
  printf("\txd = %g yd = %g zd = %g\n", 2*xd, 2*yd, 2*zd);

  /***initial parameters for the pot reduce mass, etc. in the x direction ***/
  grid->xmin = -xd;
  grid->xmax = xd;
  printf("\tThe x_max = %lg and x_min %lg\n", grid->xmax, grid->xmin);
  ntmp  = (long)((grid->xmax - grid->xmin) / grid->dx);
  if (ntmp > grid->nx){
    printf("\tinput nx insufficient; updating parameter.\n");
    grid->nx = ntmp;
  }
  grid->xmin = -((double)(grid->nx) * grid->dx) / 2.0;
  grid->xmax = ((double)(grid->nx) * grid->dx) / 2.0;
  
  grid->dkx = TWOPI / ((double)grid->nx * grid->dx); // reciprocal lattice vector length
  
  /***initial parameters for the pot reduce mass, etc. in the y direction ***/
  grid->ymin = -yd;
  grid->ymax = yd;
  printf("\tThe y_max = %lg and y_min %lg\n", grid->ymax, grid->ymin);
  ntmp  = (long)((grid->ymax - grid->ymin) / grid->dy);
  if (ntmp > grid->ny){
    printf("\tinput ny insufficient; updating parameter.\n");
    grid->ny = ntmp;
  }
  grid->ymin = -((double)(grid->ny) * grid->dy) / 2.0;
  grid->ymax = ((double)(grid->ny) * grid->dy) / 2.0;

  grid->dky = TWOPI / ((double)grid->ny * grid->dy);

  /***initial parameters for the pot reduce mass, etc. in the z direction ***/
  grid->zmin = -zd;
  grid->zmax = zd;
  ntmp  = (long)((grid->zmax - grid->zmin) / grid->dz);
  printf("\tThe z_max = %lg and z_min %lg\n", grid->zmax, grid->zmin);
  if (ntmp > grid->nz){
    printf("\tinput nz insufficient; updating parameter.\n");
    grid->nz = ntmp;
  }
  grid->zmin = -((double)(grid->nz) * grid->dz) / 2.0;
  grid->zmax = ((double)(grid->nz) * grid->dz) / 2.0;

  grid->dkz = TWOPI / ((double)grid->nz * grid->dz);
  
  grid->nx_1 = 1.0 / (double)(grid->nx);
  grid->ny_1 = 1.0 / (double)(grid->ny);
  grid->nz_1 = 1.0 / (double)(grid->nz);
  grid->dv = par->dv = grid->dx * grid->dy * grid->dz;
  grid->dr = sqrt(sqr(grid->dx) + sqr(grid->dy) + sqr(grid->dz));

  ist->ngrid = grid->ngrid = grid->nx * grid->ny * grid->nz;
  ist->nspinngrid = ist->nspin * ist->ngrid;
  ist->nx = grid->nx; ist->ny = grid->ny; ist->nz = grid->nz;


  printf("\n\tFinal grid parameters:\n");
  printf("\t----------------------\n");
  printf("\txd = %.1f yd = %.1f zd = %.1f\n", 2*grid->xmax, 2*grid->ymax, 2*grid->zmax);
  printf("\tdx = %g dy = %g dz = %g dv = %g dr = %g\n", grid->dx, grid->dy, grid->dz, grid->dv, grid->dr);
  printf("\tnx = %ld  ny = %ld  nz = %ld\n", ist->nx, ist->ny, ist->nz);
  printf("\tngrid = %ld, nspin = %d, nspinngrid = %ld\n", ist->ngrid, ist->nspin, ist->nspinngrid);
  

  fflush(stdout);

  return;
}

/*****************************************************************************/

void build_grid_ksqr(double *ksqr, xyz_st *R, grid_st *grid, index_st *ist, par_st *par, flag_st *flag){
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

  // ****** ****** ****** ****** ****** ****** 
  // Building the grid 
  // ****** ****** ****** ****** ****** ******
  printf("\tBuilding the grid...\n");
  
  for (jx = 0, dx = grid->xmin; jx < grid->nx; jx++, dx += grid->dx) grid->x[jx] = dx;
  for (jy = 0, dy = grid->ymin; jy < grid->ny; jy++, dy += grid->dy) grid->y[jy] = dy;
  for (jz = 0, dz = grid->zmin; jz < grid->nz; jz++, dz += grid->dz) grid->z[jz] = dz;
  
  // ****** ****** ****** ****** ****** ****** 
  // Initializing the k^2 vectors
  // ****** ****** ****** ****** ****** ******
  if ((ksqrx  = (double*)calloc(grid->nx,sizeof(double)))==NULL)nerror("ksqrx");
  if ((ksqry  = (double*)calloc(grid->ny,sizeof(double)))==NULL)nerror("ksqry");
  if ((ksqrz  = (double*)calloc(grid->nz,sizeof(double)))==NULL)nerror("ksqrz");

  //hold extra factor of 0.5 for the kinetic energy
  // The kinetic energy is 0.5k^2, which is symmetric, so only loop over half grid points
  for (ksqrx[0] = 0.0, jx = 1; jx <= grid->nx / 2; jx++)
    ksqrx[jx] = (ksqrx[grid->nx-jx] = 0.5 * sqr((double)(jx) * grid->dkx) *
		grid->nx_1 * grid->ny_1 * grid->nz_1);

  for (ksqry[0] = 0.0, jy = 1; jy <= grid->ny / 2; jy++)
    ksqry[jy] = (ksqry[grid->ny-jy] = 0.5 * sqr((double)(jy) * grid->dky) *
		grid->ny_1 * grid->nx_1 * grid->nz_1);

  for (ksqrz[0] = 0.0, jz = 1; jz <= grid->nz / 2; jz++)
    ksqrz[jz] = (ksqrz[grid->nz-jz] = 0.5 * sqr((double)(jz) * grid->dkz) *
		grid->nz_1 * grid->nx_1 * grid->ny_1);

  pf = fopen("ksqr.dat", "w");
  par->KE_max *= (grid->ny_1 * grid->nx_1 * grid->nz_1);
  for (jz = 0; jz < grid->nz; jz++){
    for (jy = 0; jy < grid->ny; jy++){
      jyz = grid->nx * (grid->ny * jz + jy);
      for (jx = 0; jx < grid->nx; jx++){
        jxyz = jyz + jx;
        ksqr[jxyz] = ksqrx[jx] + ksqry[jy] + ksqrz[jz];
        if (ksqr[jxyz] > par->KE_max) ksqr[jxyz] = par->KE_max; // KE cutoff at KE_max
        fprintf(pf, "%ld %lg\n", jxyz, ksqr[jxyz]);
      }
    }
  }
  fclose(pf);

  if ((1 == flag->NL) || (1 == flag->SO)){
    // Count the max number of grid points within Rnlcut of an atom***/
    printf("\tCount max no. grid points in Rnlcut of an atom\n");
    ist->n_NL_gridpts = 0;
    for (jatom = 0; jatom < ist->n_NL_atoms; jatom++) {
      for (jtmp =0, jz = 0; jz < grid->nz; jz++) {
        for (jy = 0; jy < grid->ny; jy++) {
          for (jx = 0; jx < grid->nx; jx++) {
            dx = grid->x[jx] - R[jatom].x;
            dy = grid->y[jy] - R[jatom].y;
            dz = grid->z[jz] - R[jatom].z;
            if (dx*dx + dy*dy + dz*dz < par->R_NLcut2) {
              jtmp++; 
            }
          }
        }
      }
      if (jtmp > ist->n_NL_gridpts) ist->n_NL_gridpts = jtmp;
    }
    printf("\tn_NL_gridpts = %ld\n", ist->n_NL_gridpts);
    fflush(stdout); 
  }

  free(ksqrx); free(ksqry);  free(ksqrz);
}

/*****************************************************************************/

void set_ene_targets(double *ene_targets, index_st *ist, par_st *par, flag_st *flag){
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
  
  if (flag->setTargets != 1){
    // If the energy targets were set for VB and CB separately,
    // check that nVB + nCB is equal to the total number of energy targets (states per filter)
    if (par->n_targets_VB + par->n_targets_CB != ist->m_states_per_filter){
      fprintf(stderr, "ERROR: n_targets_VB + n_targets_CB not equal to total m_states_per_filter!\n");
      exit(EXIT_FAILURE);
    }
    // If the ene targets were not set in input, place half in VB/CB
    par->n_targets_CB = (long)(ist->m_states_per_filter / 2);
    par->n_targets_VB = ist->m_states_per_filter - par->n_targets_CB;
  }
  // If either band has less than 1 ene target, crash
  if(par->n_targets_CB < 1 || par->n_targets_VB < 1){
    fprintf(stderr, "ERROR: less than one energy target set in VB or CB\n"); fflush(0);
    exit(EXIT_FAILURE);
  }

  // Compute the spacing between VB ene targets
  del = (par->VBmax - par->VBmin)/(double)par->n_targets_VB;
  printf("\tSpacing between states in VB: %lg a.u.\n", del);
  // Set the VB ene targets so that the VBmax energy is included
  for (jx = 0; jx < par->n_targets_VB; jx++) {
    ene_targets[par->n_targets_VB - jx - 1] = par->VBmax - (double)(jx) * del;
  }
  
  // Compute the spacing between CB ene targets
  del = (par->CBmax - par->CBmin)/(double)par->n_targets_CB;
  printf("\tSpacing between states in CB: %lg a.u.\n", del);
  // Set the CB ene targets so the CBmin energy is included
  for (jx = par->n_targets_VB; jx < ist->m_states_per_filter; jx++) {
    ene_targets[jx] = par->CBmin + (double)(jx-par->n_targets_VB) * del;
  }
  
  printf("\n\tEnergy target list:\n");
  printf("\t-------------------\n");
  // Print out the energy targets
  pf = fopen("ene_targets.dat","w");
  for (jx = 0; jx < ist->m_states_per_filter; jx++) {
    fprintf(pf, "%g\n", ene_targets[jx]);
    printf("\t%g\n", ene_targets[jx]);
  }
  fclose(pf);

  return;
}


/*****************************************************************************/
void build_local_pot(double *pot_local, pot_st *pot, xyz_st *R, atom_info *atom,
    grid_st *grid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){
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

  long jx, jy, jz, jyz, jxyz, jatom;
  double del, dx, dy, dz;
  double sum;
  vector *atom_neighbor_list;
  double *strain_scale;
  double *vol_ref;
  double strain_factor = 1.0; // Default: when the strain factor is 1.0, it has no effect on the potential
  int scale_LR = 0;

  // turn on the scale LR flag if surface Cs atoms will be scaled
  if (1.0 != par->scale_surface_Cs){
    printf("Surface Cs atom charges being scaled by %lg\n", par->scale_surface_Cs);
    scale_LR = 1;
  }
  // Allocate memory for strain parameters if strain-dependent potentials are requested
  if (1 == flag->useStrain){
    if ((pot->a4_params = (double *) calloc(ist->ngeoms * ist->n_atom_types, sizeof(pot->a4_params[0]))) == NULL){
    fprintf(stderr, "\nOUT OF MEMORY: a4params\n\n"); exit(EXIT_FAILURE);
    }
    if ((pot->a5_params = (double *) calloc(ist->ngeoms * ist->n_atom_types, sizeof(pot->a5_params[0]))) == NULL){
    fprintf(stderr, "\nOUT OF MEMORY: a4params\n\n"); exit(EXIT_FAILURE);
    }
    if ((atom_neighbor_list = (vector *) calloc(4 * ist->natoms, sizeof(vector))) == NULL){
      fprintf(stderr, "OUT OF MEMORY: atom_neighbor_list\n");
      exit(EXIT_FAILURE);
    }
    if ((strain_scale = (double *) calloc(ist->natoms, sizeof(double))) == NULL){
      fprintf(stderr, "OUT OF MEMORY: strain_scale\n");
      exit(EXIT_FAILURE);
    }
    if ((vol_ref = (double *) calloc(ist->natoms, sizeof(double))) == NULL){
      fprintf(stderr, "OUT OF MEMORY: vol_ref\n");
      exit(EXIT_FAILURE);
    }
  }

  // ****** ****** ****** ****** ****** ****** 
  // Read atomic pseudopotentials
  // ****** ****** ****** ****** ****** ******
  printf("\tReading atomic pseudopotentials...\n"); 
  read_pot(pot, R, atom, ist, par, flag);
  
  // ****** ****** ****** ****** ****** ****** 
  // Calculate strain_scale for atomic pots
  // ****** ****** ****** ****** ****** ******
  if (1 == flag->useStrain){
    read_nearest_neighbors(atom_neighbor_list, vol_ref, ist->natoms, ist->crystal_structure_int, ist->outmost_material_int);
    calc_strain_scale(strain_scale, atom_neighbor_list, vol_ref, atom, pot->a4_params, pot->a5_params, ist->natoms);
  }

  // ****** ****** ****** ****** ****** ****** 
  // Construct pseudopotential on grid
  // ****** ****** ****** ****** ****** ******
  printf("\tConstructing total pseudopotential on the grid...\n");

  omp_set_dynamic(0);
  omp_set_num_threads(parallel->nthreads);
#pragma omp parallel for private(dx,dy,dz,del,jy,jx,jyz,jxyz,sum,jatom)
  for (jz = 0; jz < grid->nz; jz++) {
    for (jy = 0; jy < grid->ny; jy++) {
      jyz = grid->nx * (grid->ny * jz + jy);
      for (jx = 0; jx < grid->nx; jx++) {
      	jxyz = jyz + jx;
      	for (sum = 0.0, jatom = 0; jatom < ist->natoms; jatom++){
      	  dx = grid->x[jx] - R[jatom].x;
      	  dy = grid->y[jy] - R[jatom].y;
      	  dz = grid->z[jz] - R[jatom].z;
      	  del = sqrt(dx * dx + dy * dy + dz * dz);

          // If strain dependent terms are requested, then set strain_factor to be the strain term for this atom
          if (1 == flag->useStrain){
            strain_factor = strain_scale[jatom];
          }
          // If potential interpolation is requested, then the potential needs to be the weighted average of the two geometries
          if (flag->interpolatePot == 1){
            if ((jxyz == 0)&& (jatom == 0)) printf("\tComputing interpolated cubic/ortho potential\n"); 
            //cubic part of the function
            sum += (1.0-atom[jatom].geom_par)*interpolate(del,pot->dr[2*atom[jatom].idx],pot->r,pot->r_LR,pot->pseudo,pot->pseudo_LR,ist->max_pot_file_len,
                pot->file_lens[2*atom[jatom].idx],2*atom[jatom].idx,scale_LR,atom[jatom].LR_par, strain_factor, flag->LR);
            //ortho part of the function
            sum += (atom[jatom].geom_par)*interpolate(del,pot->dr[2*atom[jatom].idx+1],pot->r,pot->r_LR,pot->pseudo,pot->pseudo_LR,ist->max_pot_file_len,
                pot->file_lens[2*atom[jatom].idx+1],2*atom[jatom].idx+1,scale_LR,atom[jatom].LR_par, strain_factor, flag->LR);
          } 
          // Default route without potential interpolation between geometries
          else {
            if ((jxyz == 0) && (jatom == 0)) printf("\tComputing potential without interpolating over cubic/ortho parameters\n\n"); 
      	    sum += interpolate(del,pot->dr[atom[jatom].idx],pot->r,pot->r_LR,pot->pseudo,pot->pseudo_LR,ist->max_pot_file_len,
                pot->file_lens[atom[jatom].idx],atom[jatom].idx,scale_LR,atom[jatom].LR_par, strain_factor, flag->LR);
          }
      	}
      	pot_local[jxyz] = sum;
      }
    }
  }

  
  // Compute potential minimum
  par->Vmin = 1.0e10;
  par->Vmax = -1.0e10;
  for (jxyz= 0; jxyz < ist->ngrid; jxyz++){
    if (par->Vmax < pot_local[jxyz]) par->Vmax = pot_local[jxyz];
    if (par->Vmin > pot_local[jxyz]) par->Vmin = pot_local[jxyz];
  }

  printf("\tVmin = %g Vmax = %g dV = %g \n", par->Vmin, par->Vmax, par->Vmax-par->Vmin);

  if (1 == flag->useStrain){
    free(atom_neighbor_list); free(vol_ref); free(strain_scale);
    free(pot->a4_params); free(pot->a5_params);
  }
  return;
}

/*****************************************************************************/
void init_SO_projectors(double *SO_projectors, xyz_st *R, atom_info *atom, grid_st *grid, index_st *ist, par_st *par){
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
  
  int rpoint;
  double dr, dr_proj, *vr; 
  long N = PROJ_LEN;

  // Compute a grid on which to defining the spherical harmonics and ang.mom. radial functions
  // name the grid vr, allocate memory
  if ((vr = (double*) calloc(N, sizeof(double)))==NULL){nerror("mem_vr");}

  dr = sqrt(par->R_NLcut2) / ((double) N) ;
  for ( rpoint = 0; rpoint < N; rpoint++){
		vr[rpoint] = (double) rpoint * dr ;
	}
  dr_proj = vr[1];

  //gen projectors on the fly
  gen_SO_projectors(grid->dx, sqrt(par->R_NLcut2), ist->nproj, SO_projectors, vr); 
  
  //printf("Projector dr = %f\n",dr_proj); fflush(0);
  
  printf("\tSO projectors generated.\n");

  free(vr);
  return;
}


/*****************************************************************************/
void init_NL_projectors(nlc_st *nlc,long *nl, double *SO_projectors, xyz_st *R,atom_info *atom, grid_st *grid, index_st *ist,par_st *par, flag_st *flag){
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
  long rpoint;
  double dr; 
  double *vr, dr_proj; 
  long N = PROJ_LEN;

  // TODO NOTE: the current method for computing the NL potential relies on the SO_projectors being defined.
  // If there is no spin-orbit coupling used in the calculation, then the SO_projectors should be set to
  // the identity.

  // Useful constants
  double tmp1 = 0.5 * sqrt(3.0 / PIE);
  double tmp2 = 0.5 * sqrt(3.0 / TWOPI);

  // Compute a grid on which to defining the spherical harmonics and ang.mom. radial functions
  // name the grid vr, allocate memoRy
  if ((vr = (double*) calloc(N, sizeof(double)))==NULL){nerror("nlc_vr");}

  dr = sqrt(par->R_NLcut2) / ((double) N) ;
  for ( rpoint = 0; rpoint < N; rpoint++){
		vr[rpoint] = (double) rpoint * dr ;
	}
  dr_proj = vr[1];
  // printf("Projector dr = %f\n",dr_proj); fflush(0);
  
  // Non-local piece
  // Find all the grid points within par->Rnlcut of each atom and calculate
  // r, r2, y1[1+m], proj(r), etc. at the grid points and store the results in nlc
  /*** for the nonlocal potential ***/

  // 2. Calculate r, r2, y1[1+m], proj(r), etc. at the grid points
#pragma omp parallel for private(jatom, vr)
  for (jatom = 0; jatom < ist->n_NL_atoms; jatom++) {
    long jx, jy, jz, jyz, jxyz;
    int iproj, *sgnProj;
    double dx, dy, dz, dxeps, dyeps, dzeps, dr_1, dr2;
    double *nlcprojectors;
    
    //gen projectors on the fly
    if ((nlcprojectors = (double*) calloc(N * ist->nproj, sizeof(double)))==NULL){nerror("nlc_projector");}
    if ((sgnProj = (int*) calloc(ist->nproj, sizeof(int)))==NULL){nerror("nlc_sgnProj");}
  
    //generate the nonlocal part for each atom
    gen_nlc_projectors(grid->dx, sqrt(par->R_NLcut2), ist->nproj, nlcprojectors, sgnProj, vr, atom, jatom);
    // printf("Exited gen_nlc_projectors %ld\n", jatom); fflush(0);

    // Print the spin-orbit parameters for each atom for use in the coupling code
    char fileName[50];
    if (atom[jatom].SO_par != 0) {
        sprintf(&fileName[0], "SO_proj_const_%ld.dat",jatom);
        pf = fopen(fileName, "w");
        fprintf(pf, "%.10f", atom[jatom].SO_par);
        fclose(pf);
    }
    
    nl[jatom] = 0;
    for (jz = 0; jz < grid->nz; jz++) {
      dz = grid->z[jz] - R[jatom].z;
      dzeps = dz + EPSDX;
      for (jy = 0; jy < grid->ny; jy++) {
        jyz = grid->nx * (grid->ny * jz + jy);
        dy = grid->y[jy] - R[jatom].y;
        dyeps = dy + EPSDX;
        for (jx = 0; jx < grid->nx; jx++) {
          jxyz = jyz + jx;
          dx = grid->x[jx] - R[jatom].x;
          dxeps = dx + EPSDX;
          dr2 = dx * dx + dy * dy + dz * dz;
          if (dr2 < par->R_NLcut2) {
            nlc[jatom*ist->n_NL_gridpts + nl[jatom]].jxyz = jxyz;

            nlc[jatom*ist->n_NL_gridpts + nl[jatom]].r  = sqrt(dr2);

            dr_1 = 1.0 / sqrt(dx * dx + dy * dy + dzeps * dzeps);
            nlc[jatom*ist->n_NL_gridpts + nl[jatom]].y1[1].re = tmp1 * dzeps * dr_1;
            nlc[jatom*ist->n_NL_gridpts + nl[jatom]].y1[1].im = 0.0;

            dr_1 = 1.0 / sqrt(dxeps * dxeps + dy * dy + dz * dz);
            nlc[jatom*ist->n_NL_gridpts + nl[jatom]].y1[2].re  = -tmp2 * dxeps * dr_1;
            nlc[jatom*ist->n_NL_gridpts + nl[jatom]].y1[0].re = tmp2 * dxeps * dr_1;

            dr_1 = 1.0 / sqrt(dx * dx + dyeps * dyeps + dz * dz);
            nlc[jatom*ist->n_NL_gridpts + nl[jatom]].y1[2].im  = -tmp2 * dyeps * dr_1;
            nlc[jatom*ist->n_NL_gridpts + nl[jatom]].y1[0].im = -tmp2 * dyeps * dr_1;

            //write projectors to nlc struct and scale projectors by the SO scaling for this atom
            for (iproj = 0; iproj< ist->nproj; iproj++){ 
                nlc[jatom*ist->n_NL_gridpts + nl[jatom]].proj[iproj] = 
                  interpolate(sqrt(dr2),dr_proj,vr,NULL, &SO_projectors[N*iproj],NULL,0, N,0,0,1.0, 1.0, 0);
                //scale projectors by the SO scaling for this atom
                nlc[jatom*ist->n_NL_gridpts + nl[jatom]].proj[iproj] *= sqrt(atom[jatom].SO_par);
                //fprintf(pf, "%li %i %f %f %f\n", jatom, iproj,nlc[jatom*ist->n_NL_gridpts + nl[jatom]].r,nlc[jatom*ist->n_NL_gridpts + nl[jatom]].proj[iproj], sqrt(atom[jatom].SO_par) );

                if (flag->NL == 1){
                  nlc[jatom*ist->n_NL_gridpts + nl[jatom]].NL_proj[iproj] =
                  interpolate(sqrt(dr2),dr_proj,vr,NULL, &nlcprojectors[N*iproj],NULL,0, N,0,0,1.0,1.0,0);
                  nlc[jatom*ist->n_NL_gridpts + nl[jatom]].NL_proj_sign[iproj] = sgnProj[iproj];
                }
            }

            if (dr2 > EPSDX) {
              nlc[jatom*ist->n_NL_gridpts + nl[jatom]].r2_1 = sqr(dr_1);
              nlc[jatom*ist->n_NL_gridpts + nl[jatom]].r2 = dr2;
            }
            else {
              nlc[jatom*ist->n_NL_gridpts + nl[jatom]].r2_1 = 0.0;
              nlc[jatom*ist->n_NL_gridpts + nl[jatom]].r2 = 1.0 / EPSDX;
            }
            nl[jatom]++;
          }
        }
      }
    }

    free(nlcprojectors);
    free(sgnProj);
  }

  pf = fopen("list_NL_grid.dat" , "w");
  for (jatom = 0; jatom < ist->n_NL_atoms; jatom++) {
    fprintf(pf, "%ld %ld\n", jatom, nl[jatom]);
  }
  fclose(pf);
  
  printf("\tNL projectors generated.\n");

  free(vr);
  
  return;
}


/*****************************************************************************/
void init_psi(zomplex *psi, long *rand_seed, int isComplex, grid_st *grid, parallel_st *parallel){
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
  for (jz = 0; jz < grid->nz; jz++){
    for (jy = 0; jy < grid->ny; jy++){
      jzy = grid->nx * (grid->ny * jz + jy);
      for (jx = 0; jx < grid->nx; jx++) {
        jxyz = jzy + jx;
        // Initialize wavefunction value at this grid point to a
        // random number between [-1.0,1.0] 
        // ran_nrc generates random between [0.0,1.0] and resets the seed
        psi[jxyz].re = (-1.0 + 2.0 * ran_nrc(&randint));

        // If using complex-valued functions, then initialize a random value for imag component
        if (1 == isComplex){
          psi[jxyz].im = (-1.0 + 2.0 * ran_nrc(&randint));
        } else if (0 == isComplex){
          // otherwise set imaginary component to 0.0
          psi[jxyz].im = 0.0;
        }
        
      }
    }
  }
  
  // normalize this wavefunction and set the value of rand_seed to the new
  // seed so the next wavefunction is different.
  
  normalize(psi, grid->dv, grid->ngrid, parallel->nthreads);
  (*rand_seed) = randint;
  
  
  return;
}

/***************************************************************************/

/***************************************************************************/

double calc_dot_dimension(xyz_st *R, long n_atoms, char *dir){
  /*******************************************************************
  * This function calculates the size of the NC along one dimension  *
  * inputs:                                                          *
  *  [R] array of coordinates along one dimension                    *
  *  [n_atoms] number of atoms in geometry                           *
  * outputs: [double] distance between atoms in a.u.                 *
  ********************************************************************/

  long i, j;
  double dr2, dis2;
  char *X; X = malloc(2*sizeof(X[0]));
  char *Y; Y = malloc(2*sizeof(Y[0]));
  char *Z; Z = malloc(2*sizeof(Z[0]));

  strcpy(X, "X"); strcpy(Y, "Y"); strcpy(Z, "Z");
  
  if (0 == strcmp(dir, (const char *)X ) ){
    dr2 = 0.0;
    for (i = 0; i < n_atoms-1; i++){
      for (j = i+1; j < n_atoms; j++){
        dis2 = sqr(R[i].x - R[j].x);
        if (dis2 > dr2) dr2 = dis2;
      }
    }
  }

  if (0 == strcmp(dir, "Y") ){
    dr2 = 0.0;
    for (i = 0; i < n_atoms-1; i++){
      for (j = i+1; j < n_atoms; j++){
        dis2 = sqr(R[i].y - R[j].y);
        if (dis2 > dr2) dr2 = dis2;
      }
    }
  }

  if (0 == strcmp(dir, "Z") ){
    dr2 = 0.0;
    for (i = 0; i < n_atoms-1; i++){
      for (j = i+1; j < n_atoms; j++){
        dis2 = sqr(R[i].z - R[j].z);
        if (dis2 > dr2) dr2 = dis2;
      }
    }
  }

  return sqrt(dr2);
}

/*****************************************************************************/
