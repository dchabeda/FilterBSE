#include "init_periodic.h"

/*****************************************************************************/

void init_periodic(
  lattice_st*    lattice,
  vector*        G_vecs,
  vector*        k_vecs,
  grid_st*       grid, 
  index_st*      ist, 
  par_st*        par, 
  flag_st*       flag, 
  parallel_st*   parallel){

  const int mpir = parallel->mpi_rank;
  
  /*** Read input for the periodic lattice ***/
  if (mpir == 0) printf("\nReading lattice params from periodic_input.par:\n");
  read_periodic_input(lattice, ist, par, flag, parallel);

  /*** calculate the reciprocal vectors for the lattice ***/
  if (mpir == 0) printf("\nGenerating reciprocal lattice vectors:\n");
  
  gen_recip_lat_vecs(&lattice, &ist, &par, &flag, &parallel);

  /*** generate the G vectors for plane-wave basis ***/
  if (mpir == 0) printf("\nGenerating G vectors:\n");
  
  ALLOCATE(&G_vecs, ist->ngrid, "G_vecs");

  gen_G_vecs(G_vecs, grid, ist, par, flag, parallel);
  
  printf("  %d G vectors in the plane-wave basis\n", ist.n_G_vecs);
  write_vector_dat(G_vecs, ist.n_G_vecs, "G_vecs.dat");


  /*** generate the k grid for computing energies/bandstructure ***/
  if (1 == flag.readKPath){
    if (mpir == 0) printf("\nReading k-path from file kpath.par:\n");
    
    read_k_path(&k_vecs, &lattice, &ist, &par, &flag, &parallel);
    
    printf("  Successfully generated %d k vectors\n", ist.n_k_pts);
  } 
  else{
    if (mpir == 0) printf("\nGenerating %d x %d x %d k grid:\n", ist.nk1, ist.nk2, ist.nk3);
    
    if((k_vecs = calloc(ist.n_k_pts, sizeof(k_vecs[0]))) == NULL){
      if (mpir == 0) fprintf(stderr, "\nOUT OF MEMORY: k_vecs\n\n"); exit(EXIT_FAILURE);
    }
    
    gen_k_vecs(k_vecs, &lattice, &ist, &par, &flag, &parallel);
    printf("  Successfully generated %d k vectors\n", ist.n_k_pts);
  }
  printf("n_kpts = %d n-G_vecs %d\n", ist.n_k_pts, ist.n_G_vecs);

  return;
}

/*****************************************************************************/

void gen_recip_lat_vecs(
  lattice_st*    lattice, 
  index_st*      ist, 
  par_st*        par, 
  flag_st*       flag, 
  parallel_st*   parallel){
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
  lattice->b1 = retScaledVector(a2xa3, TWOPI/V_lat);

  // b_2 = (2pi) * (a_3 x a_1) / [a_1 . (a_2 x a_3)]
  a3xa1 = retCrossProduct(a3, a1);
  lattice->b2 = retScaledVector(a3xa1, TWOPI/V_lat);

  // b_3 = (2pi) * (a_1 x a_2) / [a_1 . (a_2 x a_3)]
  a1xa2 = retCrossProduct(a1, a2);
  lattice->b3 = retScaledVector(a1xa2, TWOPI/V_lat);

  printf("\tDirect lattice volume, V_lat = %f\n", V_lat);
  printf("\tb1 = %.4f %.4f %.4f\n", lattice->b1.x, lattice->b1.y, lattice->b1.z);
  printf("\tb2 = %.4f %.4f %.4f\n", lattice->b2.x, lattice->b2.y, lattice->b2.z);
  printf("\tb3 = %.4f %.4f %.4f\n", lattice->b3.x, lattice->b3.y, lattice->b3.z);
  
  return;

}

/*****************************************************************************/

void gen_G_vecs(
  vector*        G_vecs, 
  grid_st*       grid, 
  index_st*      ist, 
  par_st*        par, 
  flag_st*       flag, 
  parallel_st*   parallel
  ){
  FILE *pf;
  double *gx, *gy, *gz;
  double G_max;
  int jx, jy, jz, jyz, jxyz;
  double KE_pref = sqr(HBAR) / (2*MASS_E);
  double gx_min, gy_min, gz_min;
  double scale_mag;

  if ((gx  = (double*)calloc(grid->nx,sizeof(double)))==NULL)nerror("ksqrx");
  if ((gy  = (double*)calloc(grid->ny,sizeof(double)))==NULL)nerror("ksqry");
  if ((gz  = (double*)calloc(grid->nz,sizeof(double)))==NULL)nerror("ksqrz");

  // The kinetic energy is 0.5*|(k+G)|^2
  // Here, we generate and store the system-dependent vectors G
  gx_min = - (TWOPI / grid->xmin) * grid->ngrid_1;
  for (gx[0] = 0.0, jx = 1; jx <= grid->nx / 2; jx++){
    gx[grid->nx-jx] = -1.00 * (gx[jx] = gx_min + (double)(jx) * grid->dkx * grid->ngrid_1);
  }

  gy_min = - (TWOPI / grid->ymin) * grid->ngrid_1;
  for (gy[0] = 0.0, jy = 1; jy <= grid->nx / 2; jy++){
    gy[grid->ny-jy] = -1.00 * (gy[jy] = gy_min + (double)(jy) * grid->dky * grid->ngrid_1);
  }

  gz_min = - (TWOPI / grid->zmin) * grid->ngrid_1;
  for (gz[0] = 0.0, jz = 1; jz <= grid->nz / 2; jz++){
    gz[grid->nz-jz] = -1.00 * (gz[jz] = gz_min + (double)(jz) * grid->dkz * grid->ngrid_1);
  }

  printf("These are the grid and g vectors\n");
  printf("\nnx gridx  g_x\n");
  for (jx = 0; jx < grid->nx; jx++){
      printf("%d %lf %lf\n", jx, grid->x[jx], gx[jx]);
  }
  printf("\nny gridy  g_y\n");
  for (jy = 0; jy < grid->ny; jy++){
      printf("%d %lf %lf\n", jy, grid->y[jy], gy[jy]);
  }
  printf("\nnz gridz  g_z\n");
  for (jz = 0; jz < grid->nz; jz++){
      printf("%d %lf %lf\n", jz, grid->z[jz], gz[jz]);
  }

  pf = fopen("ksqr.dat", "w");
  G_max = par->KE_max * (grid->ny_1 * grid->nx_1 * grid->nz_1) / KE_pref;
  
  for (jz = 0; jz < grid->nz; jz++){
    for (jy = 0; jy < grid->ny; jy++){
      jyz = grid->nx * (grid->ny * jz + jy);
      for (jx = 0; jx < grid->nx; jx++){
        jxyz = jyz + jx;
        G_vecs[jxyz].x = gx[jx];
        G_vecs[jxyz].y = gy[jy];
        G_vecs[jxyz].z = gz[jz];
        G_vecs[jxyz].mag = sqr(gx[jx]) + sqr(gy[jy]) + sqr(gz[jz]);
        if (G_vecs[jxyz].mag > G_max){
          scale_mag = G_max / G_vecs[jxyz].mag;
          G_vecs[jxyz] = retScaledVector(G_vecs[jxyz], scale_mag);
        } 
        fprintf(pf, "%d %.4lg %.4lg %.4lg %.4lg\n", jxyz, G_vecs[jxyz].x, G_vecs[jxyz].y, G_vecs[jxyz].z, G_vecs[jxyz].mag);
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
  vector*        k_vecs, 
  lattice_st*    lattice, 
  index_st*      ist, 
  par_st*        par, 
  flag_st*       flag, 
  parallel_st*   parallel){
  // Generate an nk1 x nk2 x nk3 kgrid for the calculation
  int i;
  int n1, n2, n3;
  double k1_scale = TWOPI / lattice->a;
  double k2_scale = TWOPI / lattice->b;
  double k3_scale = TWOPI / lattice->c;
  double dk1, dk2, dk3;
  vector k;
  
  // Generate the k grid spacing
  
  if (ist->nk1 > 1){
    dk1 = 1.0 / (double) ist->nk1;
  }
  if (ist->nk2 > 1){
    dk2 = 1.0 / (double) ist->nk2;
  }
  if (ist->nk3 > 1){
    dk3 = 1.0 / (double) ist->nk3;
  }
  
  // Generate the nk1 x nk2 x nk3 kgrid
  i = 0;
  k.x = k.y = k.z = 0.0;
  for (n1 = 0; n1 < ist->nk1; n1++){
    for (n2 = 0; n2 < ist->nk2; n2++){
      for (n3 = 0; n3 < ist->nk3; n3++){
        k.x = k1_scale * n1 * dk1;
        k.y = k2_scale * n2*dk2;
        k.z = k3_scale * n3 * dk3;
        k_vecs[i] = k;
        i++;
      }
    }
  }

  if (i != ist->n_k_pts){
      printf("ERROR: invalid k point generation!\n");
      fprintf(stderr, "ERROR: invalid k point generation!\n");
      exit(EXIT_FAILURE);
  }

  // printf("There are %d k vectors\n", ist->n_k_pts);
  // printf("These are the k_vecs:\n");
  // for (i = 0; i < ist->n_k_pts; i++){
  //     printf("  %lf %lf %lf mag = %lf\n", k_vecs[i].x, k_vecs[i].y, k_vecs[i].z, k_vecs[i].mag);
  // }

  return;
}
  
/*****************************************************************************/

void read_k_path(
  vector**       k_vecs, 
  lattice_st*    lattice, 
  index_st*      ist, 
  par_st*        par, 
  flag_st*       flag, 
  parallel_st*   parallel){
  // Read a k path from kpath.par in order to generate a bandstructure plot
  FILE *pf;
  int n_hsym, n_step, n_diff_k, n_k_pts = 0;
  int i, ik, d;
  int *n_steps_k;
  double kx, ky, kz, dk;
  vector *hsym_k;
  vector k_diff, u, lat_scale;

  // Open k-path.par and read in the high-symmetry k-point list
  if( access( "kpath.par", F_OK) != -1 ) {
    pf = fopen("kpath.par", "r");
    fscanf(pf, "%d", &n_hsym);

    // Allocate memory for the high symmetry k-path
    if ((hsym_k = (vector*) calloc(n_hsym, sizeof(hsym_k[0]))) == NULL){
      printf("ERROR: allocating memory for hsym_k\n");
      fprintf(stderr, "ERROR: allocating memory for hsym_k\n");
      exit(EXIT_FAILURE);
    }
    // Allocate memory for the number of steps between k-points
    if ((n_steps_k = (int*) calloc(n_hsym, sizeof(int))) == NULL){
      printf("ERROR: allocating memory for n_steps_k\n");
      fprintf(stderr, "ERROR: allocating memory for n_steps_k\n");
      exit(EXIT_FAILURE);
    }
    
    // Read in the k_path vectors and number of steps betweenn kpoints
    i = 0;
    for (i = 0; i < n_hsym; i++){
      fscanf(pf, "%lf %lf %lf %d", &kx, &ky, &kz, &n_step);
      hsym_k[i].x = kx;
      hsym_k[i].y = ky;
      hsym_k[i].z = kz;
      n_steps_k[i] = n_step;
      n_k_pts = n_k_pts + n_step;
      // printf("%d %f %f %f %d\n", i, hsym_k[i].x, hsym_k[i].y, hsym_k[i].z, n_steps_k[i]);
    }
  } 
  else{
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
  for (i = 0; i < n_hsym; i++){
      hsym_k[i] = retElementWiseVectorMultiplication(hsym_k[i], lat_scale);
  }
  
  // printf("These are the scaled hsym points:\n");
  // for (i = 0; i < n_hsym; i++){
  //     printf("  %lf %lf %lf mag = %lf\n", hsym_k[i].x, hsym_k[i].y, hsym_k[i].z, hsym_k[i].mag);
  // }

  // Allocate memory for the k_vecs
  if ((*(k_vecs) = (vector*) calloc(ist->n_k_pts, sizeof(vector))) == NULL){
      printf("ERROR allocating memory for k_vecs\n");
      fprintf(stderr, "ERROR allocating memory for k_vecs\n");
      exit(EXIT_FAILURE);
  }
  
  
  // We construct the k_vectors by following the path between high-sym k_points
  // That path is defined by the difference vectors between each of the k_points
  // u = k2 - k1
  
  n_diff_k = n_hsym - 1;
  i = 0;
  for (ik = 0; ik < n_diff_k; ik++){
      k_diff = retSubtractedVectors(hsym_k[ik+1], hsym_k[ik]);
      n_step = n_steps_k[ik];
      dk = k_diff.mag / n_step;
      k_diff = retScaledVector(k_diff, dk);
      for (d = 0; d < n_step; d++){
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
