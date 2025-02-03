#include "fd.h"
#include <math.h> // For pow() and exp()

/*******************************************************************************************/

void overlap_gauss(double *S_mat, gauss_st *gauss, atom_info *atom, index_st *ist, par_st *par, flag_st *flag) {
  FILE *pf;
  int a, b, i, j;
  double i_c, j_c;
  double d_ij, alpha, beta, rsqr, S_ij, ab, apb, sum;

  pf = fopen("gauss_overlap.dat", "w");
  // Compute overlap matrix elements for basis functions 
  // defined as a sum of Gaussians
  for (a = 0; a < par->n_orbitals; a++) {
    for (b = 0; b < par->n_orbitals; b++) {
      // Initialize Overlap matrix calculation
      sum = 0.0;
      for (i = 0; i < par->n_gauss_per_orbital; i++) {
        for (j = 0; j < par->n_gauss_per_orbital; j++) {
          i_c = gauss[a].coeff[i];
          j_c = gauss[b].coeff[j]; 
          
          d_ij = i_c * j_c;
          alpha = gauss[a].exp[i];
          beta = gauss[b].exp[j];
          ab = alpha * beta;
          apb = alpha + beta;
          rsqr = (sqr(gauss[b].Rx - gauss[a].Rx) +
                  sqr(gauss[b].Ry - gauss[a].Ry) +
                  sqr(gauss[b].Rz - gauss[a].Rz));
          S_ij = pow(PIE/apb , 1.5) * exp( - (ab/apb) * rsqr);
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

void kinetic_gauss(double *T_mat, gauss_st *gauss, atom_info *atom, index_st *ist, par_st *par, flag_st *flag) {
  FILE *pf;
  int a, b, i, j;
  double i_c, j_c;
  double d_ij, alpha, beta, rsqr, S_ij, ab, apb, sum;

  pf = fopen("gauss_kinetic.dat", "w");
  // Compute KE matrix elements for basis functions defined as a sum of Gaussians
  for (a = 0; a < par->n_orbitals; a++) {
    for (b = 0; b < par->n_orbitals; b++) {
      // Initialize KE matrix calculation
      sum = 0.0;
      for (i = 0; i < par->n_gauss_per_orbital; i++) {
        for (j = 0; j < par->n_gauss_per_orbital; j++) {
          i_c = gauss[a].coeff[i];
          j_c = gauss[b].coeff[j];
          
          d_ij = i_c * j_c;
          alpha = gauss[a].exp[i];
          beta = gauss[b].exp[j];
          ab = alpha * beta;
          apb = alpha + beta;
          rsqr = (sqr(gauss[b].Rx - gauss[a].Rx) +
                  sqr(gauss[b].Ry - gauss[a].Ry) +
                  sqr(gauss[b].Rz - gauss[a].Rz));
          S_ij = (ab/apb) * (3. - 2. * ab/apb * rsqr) * pow(PIE/apb , 1.5) * exp( - (ab/apb) * rsqr);
          sum += d_ij * S_ij;     
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

  // for (jxyz = 0; jxyz < ist->ngrid; jxyz++){
  //     pot_local[jxyz] = 1.0;
  // }

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

