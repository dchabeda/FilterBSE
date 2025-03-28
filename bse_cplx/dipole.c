/****************************************************************************/

#include "dipole.h"

/****************************************************************************/

void calc_elec_dipole(
  xyz_st*         elec_dip,
  double complex* psi_qp,
  double*         eig_vals,
  grid_st*        grid,
  index_st*       ist, 
  par_st*         par, 
  flag_st*        flag
  ){
  /*******************************************************************
  * This function computes the electric transition dipole matrix     *
  * matrix elements.                                                 *
  * inputs:                                                          *
  *  [elec_dip] array to hold matrix elems in x, y, z direction  *
  *  [psi_qp] array holding all qp_basis states                      *
  *  [eig_vals] array holding the quasiparticle orbital energies     *
  *  [grid] grid_st instance holding values of all grid points       *
  *  [ist] ptr to counters, indices, and lengths                     *
  *  [par] ptr to par_st holding VBmin, VBmax... params              *
  * outputs: void                                                    *
  ********************************************************************/

  /************************************************************/
	/*******************  DECLARE VARIABLES   *******************/
	/************************************************************/

  FILE*             pf; 

  //                Indices
  unsigned long     i;
  unsigned long     a;
  unsigned long     i_st;
  unsigned long     a_st;
  unsigned long     idx;

  //                Scalars
  double            z;  
  double            y;  
  double            x;  
  double            E_ia;  
  double            mu2;  

  //                Constants
  const long        nspngr    = ist->nspinngrid;
  const long        ngrid     = ist->ngrid;
  const long        lidx      = ist->lumo_idx;
  const long        n_el      = ist->n_elecs;
  const long        n_ho      = ist->n_holes;

  const double      dv        = par->dv;
  
  // Output will be written to these files
  pf = fopen("OS0.dat" , "w"); 
  
  
  fprintf(pf, "i  a   sqrt(mu2)     Ea-Ei 	  f_osc       mu_x.re     mu_x.im     mu_y.re     mu_y.im     mu_z.re     mu_z.im");

  
  /************************************************************/
	/*******************   CALC ELEC DIPOLE   *******************/
	/************************************************************/

  double start = omp_get_wtime();

  nvtxRangePushA("Calc elec dipole");
  #pragma omp parallel collapse(2) for private(a, i_st, a_st, x, y, z, idx)
  for (i = 0; i < n_ho; i++){
    for (a = lidx; a < lidx + n_el; a++) {
      // nvtxRangePushA("loop over a");
      i_st = i * nspngr;
      a_st = a * nspngr;

      long              jx;  
      long              jy;
      long              jyz;
      long              jz;  
      long              jg;
      long              jsg;
      
      double complex mu_x;
      double complex mu_y;
      double complex mu_z;
      double complex tmp; 

      mu_x = mu_y = mu_z = 0.0 + 0.0*I;

      // If using spinors, integrate cplx psi over two spins
      // Loop over the grid
      for (jz = 0; jz < ist->nz; jz++) {
        z = grid->z[jz];
        for (jy = 0; jy < ist->ny; jy++) {
          y = grid->y[jy];
          jyz = ist->nx * (ist->ny * jz + jy);
          #pragma omp simd safelen(8) aligned(psi_qp, grid->x: BYTE_BOUNDARY) reduction(+: mu_x, mu_y, mu_z)
          for (jx = 0; jx < ist->nx; jx++) { 
            x = grid->x[jx];
            jg = jyz + jx;
              for (int s = 0; s < 2; s++){
              jsg = jg + s * ngrid;

              // MATRIX ELEM
              tmp = conjmul(psi_qp[i_st + jsg], psi_qp[a_st + jsg]);
              mu_x += tmp * x;
              mu_y += tmp * y;
              mu_z += tmp * z;
            }
          }
        }
      }
      // Place these values into alloc'd memory
      idx = i * n_el + (a - lidx);

      elec_dip[idx].x = dv * mu_x;
      elec_dip[idx].y = dv * mu_y; 
      elec_dip[idx].z = dv * mu_z;
    }
  }

  double end = omp_get_wtime();
  printf("Speed of computing dipole: %f\n", (end - start)/(n_el * n_ho));
  nvtxRangePop();


  /************************************************************/
	/*******************     PRINT VALUES     *******************/
	/************************************************************/

  for (i = 0; i < n_ho; i++){
    for (a = lidx; a < lidx + n_el; a++){
      idx = i * n_el + (a - lidx);

      E_ia = eig_vals[a] - eig_vals[i];
      mu2 = cnorm(elec_dip[idx].x) + cnorm(elec_dip[idx].y) + cnorm(elec_dip[idx].z);

      fprintf(pf,"\n%ld % ld  %.8f %.12f % .8f % .8f % .8f % .8f % .8f % .8f % .8f", 
        i, a, sqrt(mu2), E_ia, (2.0/3.0) * E_ia * mu2, 
        elec_dip[idx].x,
        elec_dip[idx].y,
        elec_dip[idx].z
      );
    }
  }
  fclose(pf); 
  
  return;
}

/****************************************************************************/

void calc_mag_dipole(
  xyz_st*          mag_dip, 
  double complex*  psi_qp, 
  double*          eig_vals, 
  grid_st*         grid, 
  index_st*        ist, 
  par_st*          par, 
  flag_st*         flag
  ){
  // This function calculates the magnetic dipole matrix elements between the
  // single-particle electron (a) and hole (i) states: <psi_a|m|psi_i>
  // where m = -1/2*L = -1/2 * (r x p) where x is the cross product.

  /************************************************************/
	/*******************  DECLARE VARIABLES   *******************/
	/************************************************************/

  //                File pointers  
  FILE*             pf;

  //                Indices and states  
  long              i;
  long              a;
  long              i_st;  
  long              a_st; 
  long              idx;
  long              jg;

  //                Constants
  const long        nspngr    = ist->nspinngrid;
  const long        ngrid     = ist->ngrid;
  const long        lidx      = ist->lumo_idx;
  const long        n_el      = ist->n_elecs;
  const long        n_ho      = ist->n_holes;

  const double      dv        = par->dv;  

  //                Double values
  double            E_ia;
  double            ms;

  //                k-space grid vectors  
  double*           gx;  
  double*           gy;  
  double*           gz;  
  double*           g_vecs;  

  //                Complex arrays
  double complex*   Lxpsi;  
  double complex*   Lypsi;  
  double complex*   Lzpsi;

  // FFTW  
  fftw_plan_loc     planfw;  
  fftw_plan_loc     planbw;  
  fftw_complex*     fftwpsi;  
  

  // Allocate memory for multithreaded FFT
	fftw_init_threads();
	fftw_plan_with_nthreads(ist->nthreads);

	fftwpsi = fftw_malloc(ngrid * sizeof(fftw_complex));
  planfw = fftw_plan_dft_3d(grid->nz, grid->ny, grid->nx, fftwpsi, fftwpsi, FFTW_FORWARD, 0);
  planbw = fftw_plan_dft_3d(grid->nz, grid->ny, grid->nx, fftwpsi, fftwpsi, FFTW_BACKWARD, 0);

  /************************************************************/
	/*********************     INIT G VECS    *******************/
	/************************************************************/

	// G vectors for derivatives in k space
  ALLOCATE(&gx, grid->nx, "gx in mag_dipole");
  ALLOCATE(&gy, grid->ny, "gy in mag_dipole");
  ALLOCATE(&gz, grid->nz, "gz in mag_dipole");
  ALLOCATE(&g_vecs, 3 * ngrid, "g_vecs in mag_dipole");
	
  // Initializing the G vectors
  init_g_vecs(g_vecs, gx, gy, gz, grid, ist, par);

  // Allocate memory
  ALLOCATE(&Lxpsi, nspngr, "Lxpsi in mag_dipole");
  ALLOCATE(&Lypsi, nspngr, "Lypsi in mag_dipole");
  ALLOCATE(&Lzpsi, nspngr, "Lzpsi in mag_dipole");
  
  // Output will be written to this file
  pf = fopen("M0.dat" , "w"); 
  fprintf(pf, "  i   a    sqrt(ms)       Ea-Ei         m_x.re      m_x.im      m_y.re      m_y.im      m_z.re      m_z.im\n");
  
  /************************************************************/
	/********************     CALC MAG DIP    *******************/
	/************************************************************/
  nvtxRangePushA("Calc mag dipole");
  for (i = 0; i < n_ho; i++){
    nvtxRangePushA("loop over i");
    i_st = i * nspngr;

    nvtxRangePushA("Compute L by FFTs");
    //spin up part
		l_operator(&Lxpsi[0], &Lypsi[0], &Lzpsi[0], &psi_qp[i_st], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
    //spin dn part
    l_operator(&Lxpsi[ngrid], &Lypsi[ngrid], &Lzpsi[ngrid], &psi_qp[i_st + ngrid], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
    nvtxRangePop();

    #pragma omp parallel for private(a, a_st, jg, idx)
    for (a = lidx; a < lidx + n_el; a++) {
      nvtxRangePushA("loop over a");
      a_st = a * nspngr;
      
      double complex m_x; // <a|m|i>_x
      double complex m_y; // <a|m|i>_y
      double complex m_z; // <a|m|i>_z

      m_x = m_y = m_z = 0.0 + 0.0*I;

      #pragma omp simd safelen(8) aligned(psi_qp, Lxpsi, Lypsi, Lzpsi: BYTE_BOUNDARY) reduction(+:m_x, m_y, m_z)
      for (jg = 0; jg < nspngr; jg++) { 
        m_x += conjmul(psi_qp[a_st + jg], Lxpsi[jg]);
        m_y += conjmul(psi_qp[a_st + jg], Lxpsi[jg]);
        m_z += conjmul(psi_qp[a_st + jg], Lxpsi[jg]);
      }
      

      m_x *= 0.5 * dv;
      m_y *= 0.5 * dv;
      m_z *= 0.5 * dv;

      idx = i * n_el + (a - lidx);

      mag_dip[idx].x = m_x;
      mag_dip[idx].y = m_y;
      mag_dip[idx].z = m_z;
      nvtxRangePop();
    }
    nvtxRangePop();
  }
  nvtxRangePop();

  /************************************************************/
	/*******************     PRINT VALUES     *******************/
	/************************************************************/

  for (i = 0; i < n_ho; i++){
    for (a = lidx; a < lidx + n_el; a++){
      idx = i*n_el+(a-lidx);

      ms = cnorm(mag_dip[idx].x) + cnorm(mag_dip[idx].y) + cnorm(mag_dip[idx].z);

      E_ia = eig_vals[a] - eig_vals[i];

      fprintf(pf,"%3ld %3ld  %+.8f %+.12f %+.8f %+.8f %+.8f %+.8f %+.8f %+.8f\n", i, a, sqrt(ms), E_ia,
        mag_dip[idx].x,
        mag_dip[idx].y,
        mag_dip[idx].z
      );
    }
  }
  fclose(pf);

  // Free dynamically allocated memory
  free(Lxpsi);
  free(Lypsi);
  free(Lzpsi);

  fftw_free(fftwpsi);
	fftw_destroy_plan(planfw);
	fftw_destroy_plan(planbw);

  free(gx);
  free(gy);
  free(gz);
  free(g_vecs);

  return;
}

// /****************************************************************************/

// void rotational_strength(double *rs, double *mux, double *muy, double *muz, double *mx, 
//   double *my, double *mz, double *eig_vals, index_st ist) {
//   FILE *pf;
//   int i, a, index;

//   pf = fopen("rs0.dat", "w");
//   for (i = 0; i < ist->n_holes; i++) {
//     for (a = ist->lumo_idx; a < ist->lumo_idx+n_el; a++) {
//       index = i*n_el + (a-ist->lumo_idx);
//       rs[index] = mux[index]*mx[index] + muy[index]*my[index] + muz[index]*mz[index];
//       fprintf(pf, "%d %d %.12f  %.16f\n", i, a, eig_vals[a] - eig_vals[i], rs[index]);
//     }
//   }
//   fclose(pf);

//   return;
// }

/****************************************************************************/
