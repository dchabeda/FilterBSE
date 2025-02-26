#include "integral.h"


/*******************************************************************************************/

void overlap_gauss(double *S_mat, gauss_st *gauss, atom_info *atom, index_st *ist, par_st *par, flag_st *flag) {
  FILE *pf;
  int a, b, i, j;
  double i_c, j_c;
  xyz_st R_a;
  xyz_st R_b;
  double R2;
  double d_ij, alpha, beta, rsqr, S_ij, ab, apb, sum;
  int ix, iy, iz;
  int jx, jy, jz;
  
  pf = fopen("gauss_overlap.dat", "w");
  // Compute overlap matrix elements for basis functions 
  // defined as a sum of Gaussians
  for (a = 0; a < par->n_orbitals; a++) {
    R_a.x = gauss[a].Rx;
    R_a.y = gauss[a].Ry;
    R_a.z = gauss[a].Rz;
    for (b = 0; b < par->n_orbitals; b++) {
      R_b.x = gauss[b].Rx;
      R_b.y = gauss[b].Ry;
      R_b.z = gauss[b].Rz;

      // Compute the distance squared between orbitals a and b
      R2 = sqr(R_b.x - R_a.x) + sqr(R_b.y - R_a.y) + sqr(R_b.z - R_a.z);
      if (R2 > par->R_gint_cut2){
        // If the atoms are too far apart, set <a|S|b> = 0 and skip calculation
        S_mat[a * par->n_orbitals + b] = 0.0;
        continue;
      }

      // Initialize Overlap matrix calculation
      sum = 0.0;
      for (i = 0; i < par->n_gauss_per_orbital; i++) {
        for (j = 0; j < par->n_gauss_per_orbital; j++) {
          i_c = gauss[a].coeff[i];
          j_c = gauss[b].coeff[j]; 
          
          d_ij = i_c * j_c;
          alpha = gauss[a].exp[i];
          beta = gauss[b].exp[j];
          apb = alpha + beta;
          // ab = alpha * beta;
          // apb = alpha + beta;
          // rsqr = (sqr(gauss[b].Rx - gauss[a].Rx) +
          //         sqr(gauss[b].Ry - gauss[a].Ry) +
          //         sqr(gauss[b].Rz - gauss[a].Rz));
          // S_ij = pow(PIE/apb , 1.5) * exp( - (ab/apb) * rsqr);

          get_gauss_polarization(gauss[a].type[i], &ix, &iy, &iz);
          get_gauss_polarization(gauss[b].type[j], &jx, &jy, &jz);
  
          S_ij = S(R_a, R_b, alpha, beta, ix, iy, iz, jx, jy, jz);

          sum += d_ij * S_ij;     
        }
      }
      // Store the result in the overlap array
      S_mat[a * par->n_orbitals + b] = sum; // Ensure overlap is declared before this line
      fprintf(pf, "<%d|%d> = %.16lg\n", a, b, sum);
    }
  }
  
  fclose(pf);

  return;
    
}

/*******************************************************************************************/

void get_gauss_polarization(int polz, int* ix, int *iy, int *iz){
  
  if      (polz == 0){ (*ix) = (*iy) = (*iz) = 0;}
  else if (polz == 1){ (*ix) = 1; (*iy) = (*iz) = 0;}
  else if (polz == 2){ (*iy) = 1; (*ix) = (*iz) = 0;}
  else if (polz == 3){ (*iz) = 1; (*ix) = (*iy) = 0;}
  else    {printf("ERROR bad polarization in Gaussian!\n"); exit(EXIT_FAILURE);}

  return;
}

/*******************************************************************************************/

double S(
  xyz_st R_a, xyz_st R_b, 
  double a, double b, 
  int ix, int iy, int iz,
  int jx, int jy, int jz
  ){

  double Sij_x, Sij_y, Sij_z;
  
  double apb = a + b;
  double S_ij;

  Sij_x = E_hermite_gauss(ix, jx, 0, R_a.x - R_b.x, a, b);
  Sij_y = E_hermite_gauss(iy, jy, 0, R_a.y - R_b.y, a, b);
  Sij_z = E_hermite_gauss(iz, jz, 0, R_a.z - R_b.z, a, b);
  
  S_ij = Sij_x * Sij_y * Sij_z * pow(PIE/apb, 1.5);
  
  return S_ij;
}

/*******************************************************************************************/

void kinetic_gauss(double *T_mat, gauss_st *gauss, atom_info *atom, index_st *ist, par_st *par, flag_st *flag) {
  FILE *pf;
  int a, b, i, j;
  double i_c, j_c;
  xyz_st R_a;
  xyz_st R_b;
  double R2;
  double d_ij, alpha, beta, rsqr, T_ij, ab, apb, sum;
  int ix, iy, iz;
  int jx, jy, jz;

  pf = fopen("gauss_kinetic.dat", "w");
  // Compute KE matrix elements for basis functions defined as a sum of Gaussians
  for (a = 0; a < par->n_orbitals; a++) {
    R_a.x = gauss[a].Rx;
    R_a.y = gauss[a].Ry;
    R_a.z = gauss[a].Rz;
    for (b = 0; b < par->n_orbitals; b++) {
      R_b.x = gauss[b].Rx;
      R_b.y = gauss[b].Ry;
      R_b.z = gauss[b].Rz;

      // Compute the distance squared between orbitals a and b
      R2 = sqr(R_b.x - R_a.x) + sqr(R_b.y - R_a.y) + sqr(R_b.z - R_a.z);
      if (R2 > par->R_gint_cut2){
        // If the atoms are too far apart, set <a|S|b> = 0 and skip calculation
        T_mat[a * par->n_orbitals + b] = 0.0;
        continue;
      }

      // Initialize KE matrix calculation
      sum = 0.0;
      for (i = 0; i < par->n_gauss_per_orbital; i++) {
        for (j = 0; j < par->n_gauss_per_orbital; j++) {
          i_c = gauss[a].coeff[i];
          j_c = gauss[b].coeff[j];
          
          d_ij = i_c * j_c;
          alpha = gauss[a].exp[i];
          beta = gauss[b].exp[j];
          // ab = alpha * beta;
          // apb = alpha + beta;
          // rsqr = (
          //   sqr(gauss[b].Rx - gauss[a].Rx) +
          //   sqr(gauss[b].Ry - gauss[a].Ry) +
          //   sqr(gauss[b].Rz - gauss[a].Rz)
          // );
          // T_ij = (ab/apb) * (3. - 2. * ab/apb * rsqr) * pow(PIE/apb , 1.5) * exp( - (ab/apb) * rsqr);

          get_gauss_polarization(gauss[a].type[i], &ix, &iy, &iz);
          get_gauss_polarization(gauss[b].type[j], &jx, &jy, &jz);

          T_ij = T(R_a, R_b, alpha, beta, ix, iy, iz, jx, jy, jz);

          sum += d_ij * T_ij;     
        }
      }
      // Store the result in the overlap array
      T_mat[a * par->n_orbitals + b] = sum; // Ensure overlap is declared before this line
      fprintf(pf, "<%d|%d> = %.16lg\n", a, b, sum);
    }
  }

  fclose(pf);
  
  return;
    
}

/*******************************************************************************************/

double T(
  xyz_st R_a, xyz_st R_b, 
  double a, double b, 
  int ix, int iy, int iz,
  int jx, int jy, int jz
  ){

  double T0, T1, T2;
  
  double apb = a + b;
  double T_ij;

  T0 = b * (2 * (jx+jy+jz) + 3) * S(R_a, R_b, a, b, ix, iy, iz, jx, jy, jz);
  
  T1 = -2 * sqr(b) * (
    S(R_a, R_b, a, b, ix, iy, iz, jx+2, jy, jz) +
    S(R_a, R_b, a, b, ix, iy, iz, jx, jy+2, jz) +
    S(R_a, R_b, a, b, ix, iy, iz, jx, jy, jz+2)
  );

  T2 = -0.5 * (
    jx*(jx-1) * S(R_a, R_b, a, b, ix, iy, iz, jx-2, jy, jz) +
    jy*(jy-1) * S(R_a, R_b, a, b, ix, iy, iz, jx, jy-2, jz) +
    jz*(jz-1) * S(R_a, R_b, a, b, ix, iy, iz, jx, jy, jz-2)
  );

  T_ij = T0 + T1 + T2;

  return T_ij;
}

/*******************************************************************************************/

void potential_gauss(double *V_mat, double *pot_local, gauss_st *gauss, grid_st *grid, atom_info *atom, index_st *ist, par_st *par, flag_st *flag) {
  FILE *pf;
  int a, b, i, j;
  double i_c, j_c, i_e, j_e;
  double x, y, z;
  double Rx_a, Ry_a, Rz_a, Rx_b, Ry_b, Rz_b;
  double R_cut_orb2 = sqr(2.5) * par->R_gint_cut2;
  double R2; // distance sqr between two orbitals
  xyz_st orb_com; // orbital center of mass
  double d_ij, dist_a_sqr, dist_b_sqr, sum;
  double gauss_i, gauss_j;
  double pot_r;
  long jx, jy, jz, jyz, jxyz;
  double x_max, x_min, y_max, y_min, z_max, z_min;
  double r_thresh, dthresh, r_cut1;
  int n_grid; // Just to see how many grid points we actually loop over.

  pf = fopen("gauss_potential.dat", "w");

  for (jxyz = 0; jxyz < ist->ngrid; jxyz++){
      pot_local[jxyz] = 1.0;
  }

  // Compute KE matrix elements for basis functions defined as a sum of Gaussians
  for (a = 0; a < par->n_orbitals; a++) {
    Rx_a = gauss[a].Rx;
    Ry_a = gauss[a].Ry;
    Rz_a = gauss[a].Rz;
    for (b = 0; b < par->n_orbitals; b++) {
      // printf("\n%d %d\n", a, b);
      // Initialize KE matrix calculation
      Rx_b = gauss[b].Rx;
      Ry_b = gauss[b].Ry;
      Rz_b = gauss[b].Rz;

      // Compute distance squared between orbitals a and b
      if (a == b){
        R2 = 0.0;
      }
      else {
        R2 = sqr(Rx_b - Rx_a) + sqr(Ry_b - Ry_a) + sqr(Rz_b - Rz_a);
      }

      // If orbitals are further apart than R_cut_orb2, do not integrate (no overlap)
      if (R2 > R_cut_orb2){
        // If the atoms are too far apart to compute the overlap, then make <a|v(r)|b> = 0.0
        V_mat[a * par->n_orbitals + b] = 0.0;
        continue;
      }

      // Compute the orbital COM
      orb_com.x = 0.5 * (Rx_a + Rx_b);
      orb_com.y = 0.5 * (Ry_a + Ry_b);
      orb_com.z = 0.5 * (Rz_a + Rz_b);
      // define the cutoff thresholds for the grids
      r_cut1 = 0.5 * R2;
      dthresh = r_cut1 + 3.0;
      r_thresh = r_cut1 + sqrt(par->R_gint_cut2) + dthresh;
      x_max = orb_com.x + r_thresh;
      x_min = orb_com.x - r_thresh;
      y_max = orb_com.y + r_thresh;
      y_min = orb_com.y - r_thresh;
      z_max = orb_com.z + r_thresh;
      z_min = orb_com.z - r_thresh;
      
      sum = 0.0;
      for (i = 0; i < par->n_gauss_per_orbital; i++) {
        for (j = 0; j < par->n_gauss_per_orbital; j++) {
          i_c = gauss[a].coeff[i];
          j_c = gauss[b].coeff[j];
          d_ij = i_c * j_c;

          i_e = gauss[a].exp[i];
          j_e = gauss[b].exp[j];

          n_grid = 0;
          // Compute numerical integral between the Gauss i and j on the real space grid
          for (jz = 0; jz < grid->nz; jz++){
            z = grid->z[jz];
            // Implement cutoff in z direction
            if ((z > z_max) || (z < z_min) ){
                // We won't compute these grid pts; they are outside our orb densities
                continue;
            }
            for (jy = 0; jy < grid->ny; jy++){
              y = grid->y[jy];
              if ((y > y_max) || (y < y_min) ){
                // We won't compute these grid pts; they are outside our orb densities
                continue;
              } 
              jyz = grid->nx * (grid->ny * jz + jy);
              for (jx = 0; jx < grid->nx; jx++){
                x = grid->x[jx];
                if ((x > x_max) || (x < x_min) ){
                    // We won't compute these grid pts; they are outside our orb densities
                    continue;
                }
                n_grid++;
                jxyz = jyz + jx;
                // Computing the contribution to the integrand of grid point jgrid
                dist_a_sqr = ( sqr(x - Rx_a) + sqr(y - Ry_a) + sqr(z - Rz_a) );
                gauss_i = exp(-i_e * (dist_a_sqr) );

                dist_b_sqr = ( sqr(x - Rx_b) + sqr(y - Ry_b) + sqr(z - Rz_b) );
                gauss_j = exp(-j_e * (dist_b_sqr) );

                pot_r = pot_local[jxyz];

                sum += d_ij * gauss_i * pot_r * gauss_j;
              }
            }
          }
        }
      }
      // printf("%d %d grid points looped over = %d\n", a, b, n_grid);
      // Store the result in the overlap array
      sum *= par->dv;
      V_mat[a * par->n_orbitals + b] = sum; // Ensure overlap is declared before this line
      fprintf(pf, "<%d|%d> = %.16lg\n", a, b, sum);
    }
  }
  
  return;   
}

/*******************************************************************************************/

double E_hermite_gauss(int i, int j, int t, double Qx, double a, double b){

  // Recursive definition of Hermite Gaussian coefficients.
  // Returns a double.
  // a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)
  // b: orbital exponent on Gaussian 'b' (e.g. beta in the text)
  // i,j: orbital angular momentum number on Gaussian 'a' and 'b'
  // t: number nodes in Hermite (depends on type of integral,
  // e.g. always zero for overlap integrals)
  // Qx: distance between origins of Gaussian 'a' and 'b'
  
  double    p = a + b;
  double    q = (a * b) / p;

  if ( (t < 0) || (t > (i + j)) ){
    // out of bounds for t
    return 0.0;
  }
  else if ( (i == j) && (j == t) && (t == 0) ){
    // base case
    return exp( - q * sqr(Qx)); // K_AB
  }
  else if (j == 0){
    // decrement index i
    return (1/(2*p))*E_hermite_gauss(i-1,j,t-1,Qx,a,b) -\
    (q*Qx/a)*E_hermite_gauss(i-1,j,t,Qx,a,b) +\
    (t+1)*E_hermite_gauss(i-1,j,t+1,Qx,a,b);
  }
  else {
    // Decrement index j
    return (1/(2*p))*E_hermite_gauss(i,j-1,t-1,Qx,a,b) + \
      (q*Qx/b)*E_hermite_gauss(i,j-1,t,Qx,a,b) +\
      (t+1)*E_hermite_gauss(i,j-1,t+1,Qx,a,b);
  }


}
