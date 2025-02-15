#include "fd.h"

/*****************************************************************************/

void gen_recip_lat_vecs(lattice_st *lattice, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){
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

void gen_G_vecs(vector *G_vecs, grid_st *grid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){
  FILE *pf;
  double *gx, *gy, *gz;
  double G_max;
  int jx, jy, jz, jyz, jxyz;
  int i, cntr, idx;
  double KE_pref = 0.5;
  double gx_min, gy_min, gz_min;
  double gmag2;

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

  // printf("These are the grid and g vectors\n");
  // printf("\nnx gridx  g_x\n"); fflush(0);
  // for (jx = 0; jx < grid->nx; jx++){
  //     printf("%d %lf %lf\n", jx, grid->x[jx], gx[jx]);
  // }
  // printf("\nny gridy  g_y\n");
  // for (jy = 0; jy < grid->ny; jy++){
  //     printf("%d %lf %lf\n", jy, grid->y[jy], gy[jy]);
  // }
  // printf("\nnz gridz  g_z\n");
  // for (jz = 0; jz < grid->nz; jz++){
  //     printf("%d %lf %lf\n", jz, grid->z[jz], gz[jz]);
  // }

  pf = fopen("Gsqr.dat", "w");
  
  G_max = par->KE_max * sqr(grid->ngrid_1) / KE_pref;
  
  printf("The G_max = %lg\n", G_max);

  cntr = 0;
  idx = 0;
  for (jz = 0; jz < grid->nz; jz++){
    for (jy = 0; jy < grid->ny; jy++){
      jyz = grid->nx * (grid->ny * jz + jy);
      for (jx = 0; jx < grid->nx; jx++){
        jxyz = jyz + jx;
        gmag2 = sqr(gx[jx]) + sqr(gy[jy]) + sqr(gz[jz]);

        if (gmag2 < G_max){
          G_vecs[jxyz].x = gx[jx];
          G_vecs[jxyz].y = gy[jy];
          G_vecs[jxyz].z = gz[jz];
          G_vecs[jxyz].mag = sqr(gx[jx]) + sqr(gy[jy]) + sqr(gz[jz]);

          cntr++;

          fprintf(pf, "%d %.4lg %.4lg %.4lg %.4lg\n", 
          jxyz, G_vecs[jxyz].x, G_vecs[jxyz].y, G_vecs[jxyz].z, G_vecs[jxyz].mag);
        } 
        else{
          G_vecs[jxyz].x = 0.0;
          G_vecs[jxyz].y = 0.0;
          G_vecs[jxyz].z = 0.0;
          G_vecs[jxyz].mag = 0.0;

          idx++;
        }
      }
    }
  }
  
  ist->n_G_vecs = cntr;
  ist->n_G_zeros = idx;
  // printf("There are %d G vectors\n", ist->n_G_vecs);
  // printf("These are the G_vecs:\n");
  // for (i = 0; i < ist->n_G_vecs; i++){
  //   idx = grid->g_idx[i];
  //   printf("  %lf %lf %lf mag = %lf\n", G_vecs[idx].x, G_vecs[idx].y, G_vecs[idx].z, G_vecs[idx].mag);
  // }

  return;
}

/*****************************************************************************/

void gen_k_vecs(vector *k_vecs, lattice_st *lattice, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){
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
                k.mag = sqrt(sqr(k.x) + sqr(k.y) + sqr(k.z));
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

    printf("There are %d k vectors\n", ist->n_k_pts);
    printf("These are the k_vecs:\n");
    for (i = 0; i < ist->n_k_pts; i++){
        printf("  %lf %lf %lf mag = %lf\n", k_vecs[i].x, k_vecs[i].y, k_vecs[i].z, k_vecs[i].mag);
    }

    return;
}

/*****************************************************************************/

void read_k_path(vector **k_vecs, lattice_st *lattice, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){
  // Read a k path from kpath.par in order to generate a bandstructure plot
  FILE *pf;
  int n_hsym, n_step, n_diff_k, n_k_pts = 0;
  int i, ik, d;
  int *n_steps_k;
  double kx, ky, kz, dk;
  vector *hsym_k;
  vector k_diff, u, lat_scale;
  printf("Openeed read_k_path:\n");
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
        printf("%d %f %f %f %d\n", i, hsym_k[i].x, hsym_k[i].y, hsym_k[i].z, n_steps_k[i]);
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
  for (i = 0; i < n_hsym; i++){
      printf("  %lf %lf %lf mag = %lf\n", hsym_k[i].x, hsym_k[i].y, hsym_k[i].z, hsym_k[i].mag);
  }

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
  printf("These are the k_vecs:\n");
  for (i = 0; i < ist->n_k_pts; i++){
      printf("  %lf %lf %lf mag = %lf\n", (*k_vecs)[i].x, (*k_vecs)[i].y, (*k_vecs)[i].z, (*k_vecs)[i].mag);
  }

  return;
}


// void periodic_driver(double *pot_local, nlc_st *nlc, long *nl, grid_st *grid, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){

//   FILE *pf;
//   lattice_st lattice;
//   vector *G_vecs;
//   vector *k_vecs;
//   vector k;
//   double *H_mat;
//   double **eig_vals;
//   long long H_size;
//   int i, ik, ig;
//   char fileName[50];

//   /*** read the lattice params from periodic_input.par ***/
//   //
//   //
//   read_periodic_input(&lattice, ist, par, flag, parallel);
//   //
//   //
  
//   /*** calculate the reciprocal vectors for the lattice ***/
//   if (parallel->mpi_rank == 0) printf("\nGenerating reciprocal lattice vectors:\n");
//   gen_recip_lat_vecs(&lattice, ist, par, flag, parallel);
  
//   /*** generate the G vectors for plane-wave basis ***/
//   if (parallel->mpi_rank == 0) printf("\nGenerating G vectors:\n");
  
//   // Allocate memory for the G vectors
//   if ((G_vecs = (vector *) calloc(ist->ngrid, sizeof(vector))) == NULL){
//     if (parallel->mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: G_vecs\n\n"); exit(EXIT_FAILURE);
//   }

//   // Allocate memory for the indices where psi(g) = 0 due to cutoff
//   if (( grid->g_idx = (long *) calloc(ist->ngrid, sizeof(long))) == NULL){
//     if (parallel->mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: g_idx\n\n"); exit(EXIT_FAILURE);
//   }
//   if ((grid->g_zeros = (long *) calloc(ist->ngrid, sizeof(long))) == NULL){
//     if (parallel->mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: g_zeros\n\n"); exit(EXIT_FAILURE);
//   }

//   //
//   //
//   gen_G_vecs(G_vecs, grid, ist, par, flag, parallel);
//   //
//   //
//   printf("  %d G vectors in the plane-wave basis\n", ist->n_G_vecs);
//   write_vector_dat(G_vecs, ist->n_G_vecs, "G_vecs.dat");
  
//   /*** generate the k grid for computing energies/bandstructure ***/
//   if (1 == flag->readKPath){
//     if (parallel->mpi_rank == 0) printf("\nReading k-path from file kpath.dat:\n");
//     //
//     //
//     read_k_path(&k_vecs, &lattice, ist, par, flag, parallel);
//     //
//     //
//     printf("  Successfully generated %d k vectors\n", ist->n_k_pts);
//   } 
//   else{
//     if (parallel->mpi_rank == 0) printf("\nGenerating %d x %d x %d k grid:\n", ist->nk1, ist->nk2, ist->nk3);
//     if((k_vecs = calloc(ist->n_k_pts, sizeof(k_vecs[0]))) == NULL){
//       if (parallel->mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: k_vecs\n\n"); exit(EXIT_FAILURE);
//     }
//     //
//     //
//     gen_k_vecs(k_vecs, &lattice, ist, par, flag, parallel);
//     //
//     //
//     printf("  Successfully generated %d k vectors\n", ist->n_k_pts);
//   }
//   printf("n_kpts = %d n-G_vecs %d\n", ist->n_k_pts, ist->n_G_vecs);
//   //
//   //
//   // write_vector_dat(k_vecs, ist->n_k_pts, "kpoints.dat");
//   //
//   //

//   H_size = ist->n_k_pts * ( sqr(ist->nspin) * sqr(ist->n_G_vecs));
//   // printf("Hsize = %lld\n", H_size);
//   // Allocate memory for the Hamiltonian and eigen values
//   if ((H_mat = (double *) calloc( H_size , sizeof(double))) == NULL){
//     if (parallel->mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: H_mat\n\n"); exit(EXIT_FAILURE);
//   }

//   // Allocate memory for the eigen values
//   if ((eig_vals = (double **) calloc( ist->n_k_pts , sizeof(double *))) == NULL){
//     if (parallel->mpi_rank == 0) fprintf(stderr, "\nOUT OF MEMORY: H_mat\n\n"); exit(EXIT_FAILURE);
//   }

//   for (ik = 0; ik < ist->n_k_pts; ik++){
//     eig_vals[ik] = (double *)calloc(ist->n_G_vecs, sizeof(double));
//   }

  
//   for (ik = 0; ik < ist->n_k_pts; ik++){
//     // Initialize FT
//     fftw_init_threads();
//     fftw_plan_with_nthreads(parallel->nthreads);
//     fftw_plan_loc planfw, planbw; fftw_complex *fftwpsi; 
//     long fft_flags=0;
    
//     fftwpsi = fftw_malloc(sizeof(fftw_complex)*ist->ngrid);
//     /*** initialization for the fast Fourier transform ***/
//     planfw = fftw_plan_dft_3d(grid->nz, grid->ny, grid->nx, fftwpsi, fftwpsi, FFTW_FORWARD, fft_flags);
//     planbw = fftw_plan_dft_3d(grid->nz, grid->ny, grid->nx, fftwpsi, fftwpsi, FFTW_BACKWARD, fft_flags);
    
//     // Grab the k vector for this iteration
//     k = k_vecs[ik];

//     /* Generate the Hamiltonian in G space */
//     //
//     //
//     printf("ik = %d\n", ik);
//     H_size = ik * (sqr(ist->nspin) * sqr(ist->n_G_vecs));
//     g_hamiltonian(&H_mat[H_size], G_vecs, k, pot_local, nlc, nl, grid, ist, par, flag, parallel, planfw, fftwpsi);
//     //
//     //

//     // Diagonalize the H matrix to obtain E(k)
//     //
//     //
//     diag_mat(&H_mat[H_size], eig_vals[ik], ist->n_G_vecs);
//     //
//     //
//   }
  
//   /* Print the eigenvalues and eigenvectors to file*/
//   for (ik = 0; ik < ist->n_k_pts; ik++){
//     sprintf(fileName, "eigs-k-%d.dat\n", ik);
//     pf = fopen(fileName, "w");
//     fprintf(pf, "k = %lg %lg %lg\n", k_vecs[ik].x, k_vecs[ik].y, k_vecs[ik].z);
//     for (ig = 0; ig < ist->n_G_vecs; ig++){
//       fprintf(pf, "%d %lg\n", ig, eig_vals[ik][ig]);
//     }
//     fclose(pf);
//   }

//   pf = fopen("bandstructure.dat", "w");
//   /* Print the bandstructure to file*/
//   for (ig = par->nb_min; ig < par->nb_max; ig++){
//     for (ik = 0; ik < ist->n_k_pts; ik++){
//       fprintf(pf, "%d %.6lg %.6lg %.6lg %.8lg\n", ig, k_vecs[ik].x, k_vecs[ik].y, k_vecs[ik].z, eig_vals[ik][ig]);
//     }
//   }
//   fclose(pf);

//   /* Fourier transform the G-space wavefunctions to real space for visualization*/

//   /* End job */
//   free(H_mat); free(G_vecs); free(k_vecs); free(eig_vals);
//   return;
// }