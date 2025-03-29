/****************************************************************************/

#include "dipole.h"

/****************************************************************************/

void calc_elec_dipole(
  xyz_st*         elec_dip,
  double*         psi_qp,
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
  long              i;
  long              a;
  long              i_st;
  long              a_st;
  long              idx;

  //                Scalars
  double            z;  
  double            y;  
  double            x;  
  double            bohr_freq;  
  double            osc_strength;  

  //                Constants
  const long        nspngr    = ist->nspinngrid;
  const long        ngrid     = ist->ngrid;
  const long        cplx_idx  = ist->complex_idx;
  const long        nspin     = ist->nspin;
  const long        stlen     = nspngr * cplx_idx; // state length
  const long        lidx      = ist->lumo_idx;
  const long        n_el      = ist->n_elecs;
  const long        n_ho      = ist->n_holes;

  const double      dv        = par->dv;
  
  // Output will be written to these files
  pf = fopen("OS0.dat" , "w"); 
  
  if (flag->isComplex){
    fprintf(pf, "i  a   sqrt(mu2)     Ea-Ei 	  f_osc       mu_x.re     mu_x.im     mu_y.re     mu_y.im     mu_z.re     mu_z.im");
  } 
  
  else{
    fprintf(pf, "i  a   sqrt(mu2)     Ea-Ei 	  f_osc       mu_x.re     mu_y.re     mu_z.re");
  }
  
  /************************************************************/
	/*******************   CALC ELEC DIPOLE   *******************/
	/************************************************************/

  nvtxRangePushA("Calc elec dipole");
  #pragma omp parallel for private(a, i_st, a_st, x, y, z, idx, osc_strength, bohr_freq)
  for (i = 0; i < n_ho; i++){
    nvtxRangePushA("loop over i");
    for (a = lidx; a < lidx + n_el; a++) {
      nvtxRangePushA("loop over a");
      long              jx;  
      long              jy;  
      long              jz;  
      long              jyz;
      long              jxyz;

      long              jgur;   // jgrid_up_real  
      long              jgui;   // jgrid up imag
      long              jgdr;   // jgrid_dn_real 
      long              jgdi;   // jgrid_dn_imag 
      
      i_st = i * stlen;
      a_st = a * stlen;
      
      double mu_x_re, mu_x_im;
      double mu_y_re, mu_y_im;
      double mu_z_re, mu_z_im;
      double tmp_re, tmp_im; 

      mu_x_re = mu_x_im = 0.0;
      mu_y_re = mu_y_im = 0.0;
      mu_z_re = mu_z_im = 0.0;

      double psi_iur, psi_iui;
      double psi_idr, psi_idi;
      double psi_aur, psi_aui;
      double psi_adr, psi_adi;

      // If using spinors, integrate cplx psi over two spins
      if (2 == nspin){
        // Loop over the grid
        #pragma omp simd
        for (jz = 0; jz < ist->nz; jz++) {
          z = grid->z[jz];
          for (jy = 0; jy < ist->ny; jy++) {
            y = grid->y[jy];
            jyz = ist->nx * (ist->ny * jz + jy);
            for (jx = 0; jx < ist->nx; jx++) {
              x = grid->x[jx];

              // Index this iteration of the loop
              jxyz = jyz + jx;
              jgur = cplx_idx * jxyz;
              jgui = jgur + 1;
              jgdr = jgur + cplx_idx * ngrid;
              jgdi = jgdr + 1;

              // Load appropriate values of psi
              psi_aur = psi_qp[a_st + jgur];
              psi_iur = psi_qp[i_st + jgur];
              psi_aui = psi_qp[a_st + jgui];
              psi_iui = psi_qp[i_st + jgui];
              psi_adr = psi_qp[a_st + jgdr];
              psi_idr = psi_qp[i_st + jgdr];
              psi_adi = psi_qp[a_st + jgdi];
              psi_idi = psi_qp[i_st + jgdi];

              // REAL PART OF MATRIX ELEM
              tmp_re = psi_aur * psi_iur + psi_aui * psi_iui
                     + psi_adr * psi_idr + psi_adi * psi_idi;

              // IMAG PART OF MATRIX ELEM
              tmp_im = psi_aui * psi_iur + psi_adi * psi_idr
                     - psi_aur * psi_iui - psi_adr * psi_idi;
            
              mu_x_re += tmp_re * x; mu_y_re += tmp_re * y; mu_z_re += tmp_re * z;
              mu_x_im += tmp_im * x; mu_y_im += tmp_im * y; mu_z_im += tmp_im * z;
            }
          }
        }
        // Place these values into alloc'd memory
        idx = i * n_el + (a - lidx);

        elec_dip[idx].x_re = dv * mu_x_re;  elec_dip[idx].x_im = dv * mu_x_im;  
        elec_dip[idx].y_re = dv * mu_y_re;  elec_dip[idx].y_im = dv * mu_y_im;  
        elec_dip[idx].z_re = dv * mu_z_re;  elec_dip[idx].z_im = dv * mu_z_im;  
      }
      else if (1 == nspin){
        // Loop over the grid
        #pragma omp simd
        for (jz = 0; jz < ist->nz; jz++) {
          z = grid->z[jz];
          for (jy = 0; jy < ist->ny; jy++) {
            y = grid->y[jy];
            jyz = ist->nx * (ist->ny * jz + jy);
            for (jx = 0; jx < ist->nx; jx++) {
              x = grid->x[jx];

              // Index this iteration of the loop
              jxyz = jyz + jx;
              jgur = cplx_idx * jxyz; // psi is real, cplx_idx = 1

              // Load appropriate values of psi
              psi_aur = psi_qp[a_st + jgur];
              psi_iur = psi_qp[i_st + jgur];

              // REAL MATRIX ELEM
              tmp_re = psi_aur * psi_iur;

              mu_x_re += tmp_re * x; 
              mu_y_re += tmp_re * y; 
              mu_z_re += tmp_re * z;
            }
          }
        }

        // Place these values into alloc'd memory
        idx = i * n_el + (a - lidx);

        elec_dip[idx].x_re = dv * mu_x_re;  
        elec_dip[idx].y_re = dv * mu_y_re; 
        elec_dip[idx].z_re = dv * mu_z_re;
      }
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
      idx = i * n_el + (a - lidx);

      bohr_freq = eig_vals[a] - eig_vals[i];

      osc_strength = (
        sqr(elec_dip[idx].x_re) + sqr(elec_dip[idx].y_re) + sqr(elec_dip[idx].z_re)
        + sqr(elec_dip[idx].x_im) + sqr(elec_dip[idx].y_im) + sqr(elec_dip[idx].z_im)
      );

      if (1 == flag->useSpinors){
        fprintf(pf,"\n%ld % ld  %.8f %.12f % .8f % .8f % .8f % .8f % .8f % .8f % .8f", 
          i, a, sqrt(osc_strength), bohr_freq, (2.0/3.0)*bohr_freq*osc_strength, 
          elec_dip[idx].x_re, elec_dip[idx].x_im,
          elec_dip[idx].y_re, elec_dip[idx].y_re,
          elec_dip[idx].z_re, elec_dip[idx].z_im
        );
      } 
      
      else{
        fprintf(pf,"\n%ld % ld  %.8f %.12f % .8f % .8f % .8f % .8f", i, a, sqrt(osc_strength), bohr_freq, (2.0/3.0)*bohr_freq*osc_strength, 
          elec_dip[idx].x_re, elec_dip[idx].y_re, elec_dip[idx].z_re
        );
      }

    }
  }
  fclose(pf); 

  return;
}

/****************************************************************************/

void calc_mag_dipole(
  xyz_st*     mag_dip, 
  double*     psi_qp, 
  double*     eig_vals, 
  grid_st*    grid, 
  index_st*   ist, 
  par_st*     par, 
  flag_st*    flag
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

  long              jgr; 
  long              jsgr; 
  long              jgur;
  long              jgui;
  long              jgdr;
  long              jgdi;  

  //                Constants
  const long        nspngr    = ist->nspinngrid;
  const long        ngrid     = ist->ngrid;
  const long        cplx_idx  = ist->complex_idx;
  const long        stlen     = nspngr * cplx_idx; // state length
  const long        lidx      = ist->lumo_idx;
  const long        n_el      = ist->n_elecs;
  const long        n_ho      = ist->n_holes;

  const double      dv        = par->dv;  

  //                Double values
  double            bohr_freq;
  double            ms;

  //                k-space grid vectors  
  double*           gx;  
  double*           gy;  
  double*           gz;  
  double*           g_vecs;  

  //                Complex arrays
  zomplex*          Lxpsi;  
  zomplex*          Lypsi;  
  zomplex*          Lzpsi;

  // FFTW  
  fftw_plan_loc     planfw;  
  fftw_plan_loc     planbw;  
  fftw_complex*     fftwpsi;  
  

  // Allocate memory for multithreaded FFT
	fftw_init_threads();
	fftw_plan_with_nthreads(ist->nthreads);

	fftwpsi = fftw_malloc(ist->ngrid*sizeof(fftw_complex));
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
	Lxpsi = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));
	Lypsi = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));
	Lzpsi = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));
  
  // Output will be written to this file
  pf = fopen("M0.dat" , "w"); 
  fprintf(pf, "  i   a    sqrt(ms)       Ea-Ei         m_x.re      m_x.im      m_y.re      m_y.im      m_z.re      m_z.im\n");
  
  /************************************************************/
	/********************     CALC MAG DIP    *******************/
	/************************************************************/
  nvtxRangePushA("Calc mag dipole");
  for (i = 0; i < n_ho; i++){
    nvtxRangePushA("loop over i");
    i_st = i * stlen;

    nvtxRangePushA("Compute L by FFTs");
    //spin up part
		l_operator(&Lxpsi[0], &Lypsi[0], &Lzpsi[0], &psi_qp[i_st], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
		if (1 == flag->useSpinors){
      //spin dn part
		  l_operator(&Lxpsi[ngrid], &Lypsi[ngrid], &Lzpsi[ngrid], &psi_qp[i_st + cplx_idx * ngrid], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
    }
    nvtxRangePop();
    #pragma omp parallel for private(a, a_st, jgr, jsgr, idx)
    for (a = lidx; a < lidx + n_el; a++) {
      nvtxRangePushA("loop over a");
      a_st = a * stlen;
      
      double m_x_re, m_x_im;
      double m_y_re, m_y_im;
      double m_z_re, m_z_im;

      m_x_re = m_x_im = 0.0;
      m_y_re = m_y_im = 0.0;
      m_z_re = m_z_im = 0.0;

      // Declare variables for vals in loop
      double psi_aur, psi_aui;
      double psi_adr, psi_adi;
      double Lxur, Lxui;
      double Lyur, Lyui;
      double Lzur, Lzui;
      double Lxdr, Lxdi;
      double Lydr, Lydi;
      double Lzdr, Lzdi;

      #pragma omp simd
      for (jgr = 0; jgr < ngrid; jgr++) {
        // Handle indexing
        jsgr = jgr + ngrid;
        jgur = cplx_idx * jgr;
        jgui = jgur + 1;
        jgdr = jgur + cplx_idx * ngrid;
        jgdi = jgdr + 1;

        // Load values of psi_a
        psi_aur = psi_qp[a_st + jgur];
        psi_aui = psi_qp[a_st + jgui];
        psi_adr = psi_qp[a_st + jgdr];
        psi_adi = psi_qp[a_st + jgdi];

        // Load values of L|i>
        Lxur = Lxpsi[jgr].re;    Lxui = Lxpsi[jgr].im; // up
        Lyur = Lypsi[jgr].re;    Lyui = Lypsi[jgr].im;
        Lzur = Lzpsi[jgr].re;    Lzui = Lzpsi[jgr].im;
        Lxdr = Lxpsi[jsgr].re;   Lxdi = Lxpsi[jsgr].im; // down
        Lydr = Lypsi[jsgr].re;   Lydi = Lypsi[jsgr].im;
        Lzdr = Lzpsi[jsgr].re;   Lzdi = Lzpsi[jsgr].im;

        // Handle up spin
        m_x_re += psi_aur * Lxur - psi_aui * Lxui;
        m_x_im += psi_aur * Lxui - psi_aui * Lxur;

        m_y_re += psi_aur * Lyur - psi_aui * Lyui;
        m_y_im += psi_aur * Lyui - psi_aui * Lyur;

        m_z_re += psi_aur * Lzur - psi_aui * Lzui;
        m_z_im += psi_aur * Lzui - psi_aui * Lzur;

        // Handle down spin
        m_x_re += psi_adr * Lxdr - psi_adi * Lxdi;
        m_x_im += psi_adr * Lxdi - psi_adi * Lxdr;

        m_y_re += psi_adr * Lydr - psi_adi * Lydi;
        m_y_im += psi_adr * Lydi - psi_adi * Lydr;

        m_z_re += psi_adr * Lzdr - psi_adi * Lzdi;
        m_z_im += psi_adr * Lzdi - psi_adi * Lzdr;
      }
      nvtxRangePop();

      m_x_re *= dv;   m_y_re *= dv;   m_z_re *= dv;
      m_x_im *= dv;   m_y_im *= dv;   m_z_im *= dv;

      idx = i * n_el + (a - lidx);

      mag_dip[idx].x_re =  0.5 * m_x_re;
      mag_dip[idx].y_re =  0.5 * m_y_re;
      mag_dip[idx].z_re =  0.5 * m_z_re;

      mag_dip[idx].x_im =  0.5 * m_x_im;
      mag_dip[idx].y_im =  0.5 * m_y_im;
      mag_dip[idx].z_im =  0.5 * m_z_im;
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

      ms = ( 
        sqr(mag_dip[idx].x_re) + sqr(mag_dip[idx].x_im) 
        + sqr(mag_dip[idx].y_re) + sqr(mag_dip[idx].y_im)
        + sqr(mag_dip[idx].z_re) + sqr(mag_dip[idx].z_im) 
      );

      bohr_freq = eig_vals[a] - eig_vals[i];

      fprintf(pf,"%3ld %3ld  %+.8f %+.12f %+.8f %+.8f %+.8f %+.8f %+.8f %+.8f\n", i, a, sqrt(ms), bohr_freq,
        mag_dip[idx].x_re,
        mag_dip[idx].x_im,
        mag_dip[idx].y_re,
        mag_dip[idx].y_im,
        mag_dip[idx].z_re,
        mag_dip[idx].z_im
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
