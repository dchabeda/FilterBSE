#include "gauss.h"

/*******************************************************************************************/

void mod_gauss(
  gauss_st *gauss, double *pot_local, double * eig_vals, xyz_st *R, atom_info *atom, grid_st *grid,
  index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){
  // Driver for computing Filter in a Gaussian basis
  double *S_mat;
  double *T_mat;
  double *V_mat;
  double *H_mat;
  int a;
  double *X, *U, *C;
  MO_st *MO;
  
  // Allocate memory for the Gauss struct
  gauss = (gauss_st*) calloc(par->n_orbitals, sizeof(gauss_st));
  for (a = 0; a < par->n_orbitals; a++){
    gauss[a].coeff = (double*) calloc(par->n_gauss_per_orbital, sizeof(double));
    gauss[a].exp = (double*) calloc(par->n_gauss_per_orbital, sizeof(double));
    gauss[a].type = (int*) calloc(par->n_gauss_per_orbital, sizeof(int));
  }

  // Allocate memory for the S, T, and V arrays
  if ((S_mat = (double*) calloc(par->n_orbitals * par->n_orbitals, sizeof(S_mat[0]))) == NULL){
    fprintf(stderr, "ERROR: allocating memory for S_mat in main.c\n");
    exit(EXIT_FAILURE);
  }
  if ((T_mat = (double*) calloc(par->n_orbitals * par->n_orbitals, sizeof(T_mat[0]))) == NULL){
    fprintf(stderr, "ERROR: allocating memory for S_mat in main.c\n");
    exit(EXIT_FAILURE);
  }
  if ((V_mat = (double*) calloc(par->n_orbitals * par->n_orbitals, sizeof(V_mat[0]))) == NULL){
    fprintf(stderr, "ERROR: allocating memory for S_mat in main.c\n");
    exit(EXIT_FAILURE);
  }
  if ((H_mat = (double*) calloc(par->n_orbitals * par->n_orbitals, sizeof(H_mat[0]))) == NULL){
    fprintf(stderr, "ERROR: allocating memory for S_mat in main.c\n");
    exit(EXIT_FAILURE);
  }
  if ((X = (double*) calloc(par->n_orbitals * par->n_orbitals, sizeof(X[0]))) == NULL){
    fprintf(stderr, "ERROR: allocating memory for X in main.c\n");
    exit(EXIT_FAILURE);
  }
  if ((U = (double*) calloc(par->n_orbitals * par->n_orbitals, sizeof(U[0]))) == NULL){
    fprintf(stderr, "ERROR: allocating memory for U in main.c\n");
    exit(EXIT_FAILURE);
  }
  // Allocate memory for the orbital coefficient matrix, C
  if ((C = (double*) calloc(par->n_orbitals * par->n_orbitals, sizeof(double)))==NULL){
    fprintf(stderr, "ERROR: allocating memory for C in main.c\n");
    exit(EXIT_FAILURE);
  }
  if ((MO = (MO_st*)calloc(par->n_orbitals, sizeof(MO_st)))==NULL)nerror("Could not allocate MOs\n");
    for (a = 0; a < par->n_orbitals; a++){
      if ((MO[a].coeff = (double*)calloc(par->n_orbitals, sizeof(double)))==NULL)nerror("Could not allocate AO coefficients to MOs\n");
  }

  //
  //
  init_gauss_params(gauss, R, atom, ist, par, flag, parallel);
  //
  //
  overlap_gauss(S_mat, gauss, atom, ist, par, flag);
  //
  //
  kinetic_gauss(T_mat, gauss, atom, ist, par, flag);
  //
  //
  potential_gauss(V_mat, pot_local, gauss, grid, atom, ist, par, flag);
  //
  //
  build_gauss_hamiltonian(H_mat, T_mat, V_mat, ist, par);
  //
  //
  calc_X_canonical(S_mat, X, U, ist, par, flag, parallel);
  //
  //
  transform_H(H_mat, X, eig_vals, MO, ist, par, flag, parallel);
  //
  //
  // Print out the eigenvalues
  for (a = 0; a < par->n_orbitals; a++){
    printf("\tE%d = %lg\n", a, eig_vals[a]);
  }
  //
  //
  if (1 == flag->printGaussCubes){
    proj_gauss_on_grid(MO, gauss, grid, ist->n_gcubes_s, ist->n_gcubes_e, ist, par, flag, parallel);
  }
  //
  //
  return;
}

void init_gauss_params( 
  gauss_st     *gauss,   xyz_st *R, 
  atom_info    *atom,    index_st *ist, 
  par_st       *par,     flag_st *flag, 
  parallel_st  *parallel)
{
  // This function initializes the parameters necessary to build an atom-centered Gaussian Basis

  // The gaussian functions are represented by two parameters; coeff, and wid
  // We will read in these two parameters for each Gaussian on atom type alpha
  // and store then in the arrays
  
  //Each atom will have a seperate file containing the Gaussian parameters for that atom type
  // e.g. gaussCd.par, gaussSe.par
  // We will loop over all of the atom types, read the corresponding files, and populate the
  // gauss-> coeff and gauss->wid arrays with the appropriate values

  // Set variables 
  char atyp[3];
  char g_filnm[50];
  int atyp_idx, a, i, n_orbitals_loc;
  double **coeff_loc, **exp_loc;
  int **type_loc;
  double min_exp = 1e10;
  long jatom, jatm_orb, jatm_orb_loc;
  double sigma, g_norm;
  FILE *pf;

  // 
  // 
  // We need to allocate n_atom_types amount of space for the switches
  int *alrdy_read;
  alrdy_read = (int*) calloc(ist->n_atom_types, sizeof(int));
  for (i = 0; i < ist->n_atom_types; i++){
    // Initialize all the switches to zero
    alrdy_read[i] = 0;
  }
  // We also need to allocate memory for local copies of the coeffs and exponents
  n_orbitals_loc = ist->n_atom_types*par->n_orbitals_per_atom;
  coeff_loc = (double **) calloc(n_orbitals_loc, sizeof(double*));
  exp_loc = (double **) calloc(n_orbitals_loc, sizeof(double*));
  type_loc  = (int **) calloc(n_orbitals_loc, sizeof(int*));

  for (a = 0; a < n_orbitals_loc; a++){
    coeff_loc[a] = (double*) calloc(par->n_gauss_per_orbital, sizeof(double));
    exp_loc[a] = (double*) calloc(par->n_gauss_per_orbital, sizeof(double));
    type_loc[a]  = (int*) calloc(par->n_gauss_per_orbital, sizeof(int));
  }
  // 
  // 

  for (jatom = 0; jatom < ist->natoms; jatom++){
    strcpy(&atyp, atom[jatom].atyp);
    atyp_idx = atom[jatom].idx;
    // atype now contains the atomic symbol. Ex. if iatm = 48, atype = Cd.
    // atyp_idx contains the index of this atom type
    
    // Check if we have already encountered this atom type and read its params
    // If 0, then the file has NOT been read, and we should read it
    // If 1, then the file has been read, and we should read in the local copy of params
    if ( alrdy_read[atyp_idx] == 0 ){
      double coeff, exp, ieof;
      int type_val;  // Variable to hold the orbital type read from file
      // Make the filename for this particular atom type
      sprintf (g_filnm, "gauss%c%c%c", atyp[0], atyp[1], atyp[2]);
      strcat(g_filnm, ".par");
      printf("\n\tGaussian param file name = %s\n", g_filnm);

      // open the ligand potential file if it exists
      pf = fopen(g_filnm, "r");
      if (pf != NULL) {
        // Loop over each basis function
        for (a = 0; a < par->n_orbitals_per_atom; a++){
          jatm_orb = jatom * par->n_orbitals_per_atom + a;
          jatm_orb_loc = atyp_idx * par->n_orbitals_per_atom + a; 
          
          // Before reading in Gauss params, associate the coordinates of the atom for this basis function
          gauss[jatm_orb].Rx = R[jatom].x;
          gauss[jatm_orb].Ry = R[jatom].y;
          gauss[jatm_orb].Rz = R[jatom].z;
          
          // Grab n_gauss_per_orbital number of coeffs and widths
          for (i = 0; i < par->n_gauss_per_orbital; i++){
            ieof = fscanf(pf, "%d %lg %lg", &type_val, &coeff, &exp);
            
            // Compute the normalization constant for this gaussian
            sigma = 1 / sqrt(1 / (2 * exp));
            // printf("This is sigma = %lg\n", sigma);
            g_norm = pow(1/TWOPI, 1.5) * pow(sigma, 1.5) ;
            printf("This is g_norm = %lg\n", g_norm);
            
            gauss[jatm_orb].type[i]  = type_val;
            gauss[jatm_orb].coeff[i] = g_norm * coeff;
            gauss[jatm_orb].exp[i] = exp;
            
            // Store local copies of the coeffs and widths
            type_loc[jatm_orb_loc][i] = type_val;
            coeff_loc[jatm_orb_loc][i] = g_norm * coeff;
            exp_loc[jatm_orb_loc][i] = exp;

            // Check if this is the minimum exponent (widest Gaussian function)
            if (exp <= min_exp){
              min_exp = exp;
            }

            if (EOF == ieof){
              printf("Error: reached EOF prematurely reading %s\n", g_filnm);
              exit(EXIT_FAILURE);
            }
          }
        }
        // Flip the switch so that we do not return to read this file
        alrdy_read[atyp_idx] = 1;
        fclose(pf);
      } 
      else {
        // File not found
        fprintf(stderr, "Gaussian param file %s not in directory!\n", g_filnm);
        exit(EXIT_FAILURE);
      }
    }

    else if ( alrdy_read[atyp_idx] == 1 ){
      // If the file has already been read, then we still need to place its params 
      // into the gauss struct
      for (a = 0; a < par->n_orbitals_per_atom; a++){
        jatm_orb = jatom * par->n_orbitals_per_atom + a;
        jatm_orb_loc = atyp_idx * par->n_orbitals_per_atom + a; 
        
        // Before reading in Gauss params, associate the coordinates of the atom for this basis function
        gauss[jatm_orb].Rx = R[jatom].x;
        gauss[jatm_orb].Ry = R[jatom].y;
        gauss[jatm_orb].Rz = R[jatom].z;
        
        // Grab n_gauss_per_orbital number of coeffs and widths
        for (i = 0; i < par->n_gauss_per_orbital; i++){
          // Read in local copies of the polarizations, coeffs and widths
          gauss[jatm_orb].type[i]  = type_loc[jatm_orb_loc][i];
          gauss[jatm_orb].coeff[i] = coeff_loc[jatm_orb_loc][i];
          gauss[jatm_orb].exp[i] = exp_loc[jatm_orb_loc][i];
        }
      }
    }
    else {
      printf("ERROR: unknown alredy_read flag %d!\n", alrdy_read[atyp_idx]);
      fprintf(stderr, "ERROR: unknown alredy_read flag %d!\n", alrdy_read[atyp_idx]);
      exit(EXIT_FAILURE);
    }
  }
  
  if (flag->useGaussianBasis == 1){
    // ln(0.01) lets the Gaussian decay to 0.01 percent of its original value
    par->R_gint_cut2 = 4.605 / min_exp;
    printf("\tMin. Gaussian exp = %.4lg -> R_gint_cut2 = %lg\n", min_exp, par->R_gint_cut2);
  }

  // Debug: Print the values read
  for (jatom = 0; jatom < ist->natoms; jatom++){
    printf("Atom type number %ld\n", jatom);
    for (a = 0; a < par->n_orbitals_per_atom; a++){
      jatm_orb = jatom * par->n_orbitals_per_atom + a;
      printf("  Basis function %d\n", a);
      printf("    Rx = %g Ry = %g Rz = %g\n", gauss[jatm_orb].Rx, gauss[jatm_orb].Ry, gauss[jatm_orb].Rz);
      for (i = 0; i < par->n_gauss_per_orbital; i++){
        printf("    g%dpol = %d coeff = %g exp = %g\n", i, gauss[jatm_orb].type[i], gauss[jatm_orb].coeff[i], gauss[jatm_orb].exp[i]);
      }
    }
  }

  free(type_loc);
  free(coeff_loc);
  free(exp_loc);
}

/*******************************************************************************************/

void build_gauss_hamiltonian(double *H_mat, double *T_mat, double *V_mat, index_st *ist, par_st *par) {

  FILE *pf;
  int a, b, idx;

  pf = fopen("gauss_H_mat.dat", "w");

  for (a = 0; a < par->n_orbitals; a++){
    for (b = 0; b < par->n_orbitals; b++){
      idx = a * par->n_orbitals + b;
      
      H_mat[idx] = T_mat[idx] + V_mat[idx];
      fprintf(pf, "<%d|H|%d> = %.16lg\n", a, b, H_mat[idx]);
    }
  }

  fclose(pf);

//   for (a = 0; a < par->n_orbitals; a++){
//     printf("\n\t");
//     for (b = 0; b < par->n_orbitals; b++){
//       idx = a * par->n_orbitals + b;
//       printf("%lg  ", H_mat[idx]);
//     }
//   }
//   printf("\n");

  return;
}

/*******************************************************************************************/

void proj_gauss_on_grid(MO_st *MO, gauss_st *gauss, grid_st *grid, int start, int end, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){
  // project the Gaussian functions onto real space

  long jx, jy, jz, jyz, jxyz;
  int a, i;
  int n_state, H_idx;
  double i_c, i_e;
  double x, y, z;
  double Rx_a, Ry_a, Rz_a;
  double R2; // distance sqr between orbital and grid point
  double sum;
  double eig_coeff;
  double *rho;
  char fileName[50];

  // Create a flattened array containing the matrix element indices
  // i = 0;
  // for (a = 0; a < par->n_orbitals; a++){
  //   for (b = 0; b < par->n_orbitals; b++){
  //     idx = a * par->n_orbitals + b;
  //     mat_elems[i] = a;
  //     i++;
  //   }
  // }
  
  rho = calloc(ist->ngrid, sizeof(double));

  for (n_state = start; n_state < end; n_state++){
    for (jz = 0; jz < grid->nz; jz++){
      z = grid->z[jz];
      for (jy = 0; jy < grid->ny; jy++){
        y = grid->y[jy];
        jyz = grid->nx * (grid->ny * jz + jy);
        for (jx = 0; jx < grid->nx; jx++){
          x = grid->x[jx];
          jxyz = jyz + jx;

          sum = 0.0;
          // Sum the contributions to this grid point from all orbitals
          for (a = 0; a < par->n_orbitals; a++){
            // Get the vector to the center of the orbital
            Rx_a = gauss[a].Rx;
            Ry_a = gauss[a].Ry;
            Rz_a = gauss[a].Rz;
            R2 = sqr(x - Rx_a) + sqr(y - Ry_a) + sqr(z - Rz_a);
            // If the orbital is far from the grid point, do not perform sum
            if (R2 > par->R_gint_cut2){
              continue;
            }

            // Sum contributions from each Gaussian function making orbital a
            for (i = 0; i < par->n_gauss_per_orbital; i++){
              i_c = gauss[a].coeff[i];
              i_e = gauss[a].exp[i];
              // Compute the value of the Gaussian function at this point in real space
              // Add it to the total value from this orbital
              sum += i_c * exp(- i_e * R2);
            }

            // Get the eig vector coefficient for this orbital in state n
            // <a|H|n>
            eig_coeff = MO[n_state].coeff[a];

            sum *= eig_coeff;
          }
          // update the value of rho at this grid point;
          rho[jxyz] = sum;
        }
      }
    }
    // Write a cube file of the orbitals for each state
    sprintf(fileName, "orbital-%d.cube", n_state);
    write_cube_file(rho, grid, fileName);
  }

  free(rho);

  return;
}

/*******************************************************************************************/

void calc_X_canonical(double *S_mat, double *X, double *U, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){

  // Diagonalize the basis to obtain the unitary matrix U, canonical transformation X
  int i, j;
  long long M, N, K, LDA, LWORK, INFO;
  double *A, *s, *WORK, *s_rt_inv;

  N = LDA = par->n_orbitals;
  
  LWORK = 3 * N;

  if ((s = (double*) malloc(N * sizeof(s[0]))) == NULL) nerror("Could not allocate s\n");
  if ((WORK = (double*) malloc(LWORK * sizeof(WORK[0]))) == NULL) nerror("Could not allocate WORK\n");
  if ((A = (double*) malloc(N * N * sizeof(A[0])))==NULL) nerror("Could not allocate A\n");

  // Copy the S_mat matrix into A
  memcpy(&A[0], &S_mat[0], N * N * sizeof(double));

  // Diagonalize the S_mat matrix
  // 
  // 
  dsyev_("V", "U", &N, &A[0], &LDA, &s[0], &WORK[0], &LWORK, &INFO);
  // 
  // 
  if (INFO) { 
    nerror("error in dsyev_ S_mat\n");
  } 
  // printf("Done diagonalizing S_mat matrix\n"); fflush(0);

  // Print the eigenvalues and Eigenvectors
  // printf("Eigenvalues:\n");
  // for (i = 0; i < par->n_orbitals; i++){
  //     printf("% .6f\n", s[i]);
  // }
  // printf("Eigenvectors:\n");
  // for (i = 0; i < par->n_orbitals; i++){
  //     printf("\n");
  //     for (j = 0; j < par->n_orbitals; j++){
  //         printf(" % .6f ", A[i*par->n_orbitals + j]);
  //     }
  // }
  // printf("\n");

  // Make the inverse square root of the diagonalized S_mat matrix, s^-1/2
  if ((s_rt_inv = (double*)malloc(N * N * sizeof(double))) == NULL)nerror("Could not allocate inverse S_mat matrix\n");

  for (i = 0; i < par->n_orbitals; i++){
      for (j = 0; j < par->n_orbitals; j++){
          if (i == j) s_rt_inv[i * par->n_orbitals + j] = 1 / sqrt(s[i]);
          else s_rt_inv[i * par->n_orbitals + j] = 0.0;
      }
  }

  // Print the inv sqrt diag mat and Eigenvectors
  // printf("\nMatrix s^-1/2:\n");
  // for (i = 0; i < par->n_orbitals; i++){
  //     printf("\n");
  //     for (j = 0; j < par->n_orbitals; j++){
  //         printf(" % .6f ", s_rt_inv[i*par->n_orbitals + j]);
  //     }
  // }
  // printf("\n");

  // Compute X_canonical = U * s^-1/2.
  M = N = K = par->n_orbitals;
  matmul(M, N, K, A, s_rt_inv, X);    

  // printf("\nX_canonical:\n");
  // for (i = 0; i < par->n_orbitals; i++){
  //     printf("\n");
  //     for (j = 0; j < par->n_orbitals; j++){
  //         printf(" % .6f ", X[i*par->n_orbitals + j]);
  //     }
  // }
  // printf("\n");

  free(s); free(WORK); free(A); free(s_rt_inv);

  return;
}

/*******************************************************************************************/

void transform_H(double *H_mat, double *X, double *eig_vals, MO_st *MO, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel){

  // ***** ***** ***** ***** ***** ***** ***** *****
  //
  // ***** ***** ***** ***** ***** ***** ***** *****
  // Transform the H matrix by the canonical transformation 
  // so it is in the basis of orthonormal states.
  // H' = X^T * H * X
  FILE *pf;
  int i, j;
  long long N = par->n_orbitals;
  double *Hp;

  pf = stdout;

  if ((Hp = (double*) calloc(par->n_orbitals * par->n_orbitals, sizeof(Hp[0]))) == NULL){
    fprintf(stderr, "ERROR: allocating memory for S_mat in main.c\n");
    exit(EXIT_FAILURE);
  }

  //
  // 
  trans_mat(par->n_orbitals, X, H_mat, Hp);
  // 
  // 

  // Print H'
  // fprintf(pf, "\nHp = ");
  // for (i = 0; i < par->n_orbitals; i++){
  //     fprintf(pf, "\n");
  //     for (j = 0; j < par->n_orbitals; j++){
  //         fprintf(pf, " % .6f ", Hp[i*par->n_orbitals + j]);
  //     }
  // }
  // fprintf(pf, "\n\n");

  // Diagonalize the H' matrix
  // 
  // 
  diag_mat(Hp, eig_vals, par->n_orbitals);
  // 
  // 
  
  // The matrix Hp now contains the MO coefficients. Store them in the MO struct
  for (i = 0; i < par->n_orbitals; i++){
      for (j = 0; j < par->n_orbitals; j++){
          MO[i].coeff[j] = Hp[j*par->n_orbitals + i]; 
          // This looks backwards, but the i-th MO is the i-th COLUMN of Cp. The j-th AO coefficient is the j-th ROW of Cp
      }
  }
  // Print the MOs
  // printf("The first MO is Psi = %.6f Phi_1 + %.6f Phi_2\n", MO[0].coeff[0], MO[0].coeff[1] );
  // printf("The second MO is Psi = %.6f Phi_1 + %.6f Phi_2\n", MO[1].coeff[0], MO[1].coeff[1] );

  // Print the eigenvalues (orbital energies) and eigenvectors (orbital coefficients)
  // fprintf(pf, "Orbital energies:\n");
  // for (i = 0; i < par->n_orbitals; i++){
  //     fprintf(pf, "% .6f\n", eig_vals[i]);
  // }
  // fprintf(pf, "\n");

  // fprintf(pf, "Canonical orbital coefficients:\n");
  // for (i = 0; i < par->n_orbitals; i++){
  //     fprintf(pf, "\n");
  //     for (j = 0; j < par->n_orbitals; j++){
  //         fprintf(pf, " % .6f ", Hp[i*par->n_orbitals + j]);
  //     }
  // }
  // fprintf(pf, "\n\n");

  // Compute the orbital coefficients in terms of the original basis functions
  // C = X * C'
  double *C = (double *) calloc(N * N , sizeof(C[0]));
  matmul(N, N, N, X, Hp, C); 

  // fprintf(pf, "Original-orbital coefficients:\n");
  
  // for (i = 0; i < N; i++){
  //     fprintf(pf, "\n");
  //     for (j = 0; j < N; j++){
  //         fprintf(pf, " % .6f ", C[i*N + j]);  
  //     }
  // }
  // fprintf(pf, "\n\n");

}
