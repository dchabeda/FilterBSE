/****************************************************************************/

#include "angular.h"

/****************************************************************************/

void calc_spin_mtrx(
	xyz_st*      s_mom, 
	double*      psi_qp, 
	grid_st*     grid, 
	index_st*    ist, 
	par_st*      par
	){
	
  /*******************************************************************
  * This function computes the value of the spin projection operator *
  * between the electron and hole quasiparticle states. The spin     *
  * projection operators are defined by the Pauli matrices as:       *
  * Sx = 1/2 (|up><down| + |down><up|)                               *
  * Sy = 1/2 (-i|up><down| + i|down><up|)                            *
  * Sz = 1/2 (|up><up| - |down><down|)                               *
  * and we compute the matrix elements <i|Sx|j> as                   *
  * sum_r psi_i(r, s) * Sx * psi_j(r, s)                             *
  * inputs:                                                          *
  *  [s_mom] array to hold components of the spin projection elems   *
  *  [psi_qp] array holding all qp_basis states                      *
  *  [grid] grid_st instance holding values of all grid points       *
  *  [ist] ptr to counters, indices, and lengths                     *
  *  [par] ptr to par_st holding VBmin, VBmax... params              *
  *  [flag] ptr to flag_st holding job flags                         *
  *  [planfw] FFTW3 plan for executing 3D forward DFT                *
  *  [planfw] FFTW3 plan for executing 3D backwards DFT              *
  *  [fftwpsi] location to store outcome of Fourier transform        *
  * outputs: void                                                    *
  ********************************************************************/
 
  /************************************************************/
	/*******************  DECLARE VARIABLES   *******************/
	/************************************************************/
	
  FILE*            pfx;
  FILE*            pfy;
  FILE*            pfz;

	long             i;
	long             j;
	long             a;
	long             b;
	long             i_st;
	long             j_st;
	long             a_st;
	long             b_st;
	long             jgr;
  long             jgur;
  long             jgui;
  long             jgdr;
  long             jgdi;
	long             idx;
  	
	const long       n_ho   = ist->n_holes;
	const long       n_el   = ist->n_elecs;
	const long       lidx   = ist->lumo_idx;
	const long       ngrid  = ist->ngrid;
	const long       nspngr = ist->nspinngrid;
  const long       cplx_idx = ist->complex_idx;
  const long       stlen = cplx_idx * nspngr;
	
	const double     dv     = grid->dv;
  
  pfx = fopen("sx.dat", "w"); 
  pfy = fopen("sy.dat", "w"); 
  pfz = fopen("sz.dat", "w");
  
  /************************************************************/
	/*******************    CALC HOLE SPINS   *******************/
	/************************************************************/

  nvtxRangePushA("Calc qp hole spin mat elems");
	for(i = 0; i < ist->n_holes; i++){
		for (j = 0; j < ist->n_holes; j++){
      nvtxRangePushA("loop ij");
      i_st = i * stlen;
      j_st = j * stlen;

      double sx_re, sx_im;
      double sy_re, sy_im;
      double sz_re, sz_im;

      sx_re = sx_im = sy_re = 0.0;
      sy_im = sz_re = sz_im = 0.0;

      double psi_iur, psi_iui;
			double psi_idr, psi_idi;
      double psi_jur, psi_jui;
			double psi_jdr, psi_jdi;

      #pragma omp simd
			for (jgr = 0; jgr < ist->ngrid; jgr++) {
        // Handle indexing
        jgur = cplx_idx * jgr;
        jgui = jgur + 1;
        jgdr = jgur + cplx_idx * ngrid;
        jgdi = jgdr + 1;
        
        // Load values for this iteration
        // Set local values for up spin
        psi_iur = psi_qp[i_st + jgur];      psi_iui = psi_qp[i_st + jgui];
        psi_jur = psi_qp[j_st + jgur];      psi_jui = psi_qp[j_st + jgui];
        
        // Set local values for dn spin
        psi_idr = psi_qp[i_st + jgdr];      psi_idi = psi_qp[i_st + jgdi];
        psi_jdr = psi_qp[j_st + jgdr];      psi_jdi = psi_qp[j_st + jgdi];

        // S_x component
        // <j|r,dn> * <r,up|i>
        sx_re += psi_jdr * psi_iur + psi_jdi * psi_iui;
        sx_im += psi_jdr * psi_iui - psi_jdi * psi_iur;
        // <j|r,up> * <r,dn|i>
        sx_re += psi_jur * psi_idr + psi_jui * psi_idi;
        sx_im += psi_jur * psi_idi - psi_jui * psi_idr;
        
        // S_y component
        // i*<j|r,dn> * <r,up|i>
        sy_im += psi_jdr * psi_iur + psi_jdi * psi_iui;
        sy_re -= psi_jdr * psi_iui - psi_jdi * psi_iur;
        // -i*<j|r,up> * <r,dn|i>
        sy_im -= psi_jur * psi_idr + psi_jui * psi_idi;
        sy_re += psi_jur * psi_idi - psi_jui * psi_idr;
        
        // S_z component
        // <j|r,up> * <r,up|i>
        sz_re += psi_jur * psi_iur + psi_jui * psi_iui;
        sz_im += psi_jur * psi_iui - psi_jui * psi_iur;
			  // - <j|r,dn> * <r,dn|i>
			  sz_re -= psi_jdr * psi_idr + psi_jdi * psi_idi;
        sz_im -= psi_jdr * psi_idi - psi_jdi * psi_idr;
		  }
			
			// multiply all by (1/2)*dV and complex conjugate bc holes are time reversed
			s_mom[i*ist->n_holes+j].x_re = -0.5 * sx_re * dv;
			s_mom[i*ist->n_holes+j].x_im =  0.5 * sx_im * dv;
			s_mom[i*ist->n_holes+j].y_re = -0.5 * sy_re * dv;
			s_mom[i*ist->n_holes+j].y_im =  0.5 * sy_im * dv;
			s_mom[i*ist->n_holes+j].z_re = -0.5 * sz_re * dv;
			s_mom[i*ist->n_holes+j].z_im =  0.5 * sz_im * dv;
      nvtxRangePop();
		}
	}
  nvtxRangePop();
  /************************************************************/
	/*******************   PRINT HOLE OUTPUT  *******************/
	/************************************************************/

  for (i = 0; i < n_ho; i++){
    for (j = 0; j < n_ho; j++){
      idx = i * n_ho + j;
      fprintf (pfx,"%ld %ld % .10g % .10g\n", i, j, s_mom[idx].x_re, s_mom[idx].x_im); 
      fprintf (pfy,"%ld %ld % .10g % .10g\n", i, j, s_mom[idx].y_re, s_mom[idx].y_im);
      fprintf (pfz,"%ld %ld % .10g % .10g\n", i, j, s_mom[idx].z_re, s_mom[idx].z_im);
    }
  }
  fflush(pfx);
  fflush(pfy);
  fflush(pfz);

  /************************************************************/
	/*******************    CALC ELEC SPINS   *******************/
	/************************************************************/

  nvtxRangePushA("Calc qp elec spins");
	for (a = lidx; a < lidx + n_el; a++){
		for (b = lidx; b < lidx + n_el; b++){
      nvtxRangePop();
      a_st = a * stlen;
      b_st = b * stlen;

			double sx_re, sx_im;
      double sy_re, sy_im;
      double sz_re, sz_im;

      sx_re = sx_im = sy_re = 0.0;
      sy_im = sz_re = sz_im = 0.0;

      double psi_aur, psi_aui;
			double psi_adr, psi_adi;
      double psi_bur, psi_bui;
			double psi_bdr, psi_bdi;
      
      #pragma omp simd
			for (jgr = 0; jgr < ngrid; jgr++) {
        // Handle indexing
        jgur = cplx_idx * jgr;
        jgui = jgur + 1;
        jgdr = jgur + cplx_idx * ngrid;
        jgdi = jgdr + 1;
        
        // Load values for this iteration
        // Set local values for up spin
        psi_aur = psi_qp[a_st + jgur];      psi_aui = psi_qp[a_st + jgui];
        psi_bur = psi_qp[b_st + jgur];      psi_bui = psi_qp[b_st + jgui];
        
        // Set local values for dn spin
        psi_adr = psi_qp[a_st + jgdr];      psi_adi = psi_qp[a_st + jgdi];
        psi_bdr = psi_qp[b_st + jgdr];      psi_bdi = psi_qp[b_st + jgdi];
        
        //Spin x part
        //<b|r,dn> * <r,up|a>
        sx_re += psi_bdr * psi_aur + psi_bdi * psi_aui;
        sx_im += psi_bdr * psi_aui - psi_bdi * psi_aur;
        //<b|r,up> * <r,dn|a>
        sx_re += psi_bur * psi_adr + psi_bui * psi_adi;
        sx_im += psi_bur * psi_adi - psi_bui * psi_adr;
        
        //Spin y part
        //i*<b|r,dn> * <r,up|a>
        sy_im += psi_bdr * psi_aur + psi_bdi * psi_aui;
        sy_re -= psi_bdr * psi_aui - psi_bdi * psi_aur;
        //-i*<b|r,up> * <r,dn|a>
        sy_im -= psi_bur * psi_adr + psi_bui * psi_adi;
        sy_re += psi_bur * psi_adi - psi_bui * psi_adr;
        
        //Spin z part
        //<b|r,up> * <r,up|a>
        sz_re += psi_bur * psi_aur + psi_bui * psi_aui;
        sz_im += psi_bur * psi_aui - psi_bui * psi_aur;
        //<b|r,dn> * <r,dn|a>
        sz_re -= psi_bdr * psi_adr + psi_bdi * psi_adi;
        sz_im -= psi_bdr * psi_adi - psi_bdi * psi_adr;
			}

      idx = sqr(n_ho) + (a - lidx)*n_el + (b - lidx);
      
			//divide all by 2
			s_mom[idx].x_re = 0.5 * sx_re * dv;
			s_mom[idx].x_im = 0.5 * sx_im * dv;
			s_mom[idx].y_re = 0.5 * sy_re * dv;
			s_mom[idx].y_im = 0.5 * sy_im * dv;
			s_mom[idx].z_re = 0.5 * sz_re * dv;
			s_mom[idx].z_im = 0.5 * sz_im * dv;
      nvtxRangePop();
		}
	}
  nvtxRangePop();
  /************************************************************/
	/*******************   PRINT ELEC OUTPUT  *******************/
	/************************************************************/

  for (a = lidx; a < lidx + n_el; a++){
    for (b = lidx; b < lidx + n_el; b++){
      idx = sqr(n_ho) + (a-lidx)*n_el + (b-lidx);
      fprintf (pfx, "%ld %ld % .10g % .10g\n", a, b, s_mom[idx].x_re, s_mom[idx].x_im);
      fprintf (pfy, "%ld %ld % .10g % .10g\n", a, b, s_mom[idx].y_re, s_mom[idx].y_im);
      fprintf (pfz, "%ld %ld % .10g % .10g\n", a, b, s_mom[idx].z_re, s_mom[idx].z_im);
    }
  }
  fflush(pfx);
  fflush(pfy);
  fflush(pfz);

  // Free memory allocated to FILE pointers

	fclose(pfx);
	fclose(pfy);
	fclose(pfz);

	return;
}

/****************************************************************************/

void calc_ang_mom_mtrx(
	xyz_st*        l_mom, 
	zomplex*       l2_mom,
	zomplex*       LdotS,
	double*        psi_qp,
	grid_st*       grid,
	index_st*      ist,
	par_st*        par){
	
	/************************************************************/
	/*******************  DECLARE VARIABLES   *******************/
	/************************************************************/
	
	long             i;
	long             j;
	long             a;
	long             b;
	long             i_st;
	long             j_st;
	long             a_st;
	long             b_st;
	long             jgr;
  long             jgur;
  long             jgui;
  long             jgdr;
  long             jgdi;
	long             jsgr;
	long             idx;
  	
	const long       n_ho   = ist->n_holes;
	const long       n_el   = ist->n_elecs;
	const long       lidx   = ist->lumo_idx;
	const long       ngrid  = ist->ngrid;
	const long       nspngr = ist->nspinngrid;
  const long       cplx_idx = ist->complex_idx;
  const long       stlen  = cplx_idx * nspngr;
	
	const double     dv     = grid->dv;

	double*          gx;
	double*          gy;
	double*          gz;
	double*          g_vecs;

	fftw_plan_loc    planfw;
	fftw_plan_loc    planbw; 
	fftw_complex*    fftwpsi;

	/************************************************************/
	/*******************    INITIALIZE FFT    *******************/
	/************************************************************/
	
	// Note: parallelize FFT to achieve fast performance; 
	// do not parallelize over states! - Daniel C 3.13.25
	// Best strategy:
	//    - Use hybrid MPI/OpenMP to distribute initial states i/a 
	//    - Compute FFT with threads
	//    - Perform integrals on GPU

	fftw_init_threads();
	fftw_plan_with_nthreads(ist->nthreads);

	// Allocate memory for multithreaded FFT
	fftwpsi    = fftw_malloc(ist->ngrid*sizeof(fftw_complex));
  planfw     = fftw_plan_dft_3d(grid->nz, grid->ny, grid->nx, fftwpsi, fftwpsi, FFTW_FORWARD, 0);
	planbw     = fftw_plan_dft_3d(grid->nz, grid->ny, grid->nx, fftwpsi, fftwpsi, FFTW_BACKWARD, 0);

	/************************************************************/
	/*******************    INIT G vecs   *******************/
	/************************************************************/

	// G vectors for derivatives in k space
  ALLOCATE(&gx, grid->nx, "gx");
  ALLOCATE(&gy, grid->ny, "gy");
  ALLOCATE(&gz, grid->nz, "gz");
	ALLOCATE(&g_vecs, 3 * ngrid, "g_vecs");

  // Initializing the G vectors
  init_g_vecs(g_vecs, gx, gy, gz, grid, ist, par);

	/************************************************************/
	/*******************    ALLOC L OP MEM    *******************/
	/************************************************************/

	zomplex* Lxpsi = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));
	zomplex* Lypsi = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));
	zomplex* Lzpsi = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));

	zomplex* Lx2psi = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));
	zomplex* Ly2psi = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));
	zomplex* Lz2psi = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));
	zomplex* temp1 = (zomplex*) calloc(ist->ngrid, sizeof(zomplex));
	zomplex* temp2 = (zomplex*) calloc(ist->ngrid, sizeof(zomplex));
	
	FILE* pfx = fopen("lx.dat", "w"); 
	FILE* pfy = fopen("ly.dat", "w"); 
	FILE* pfz = fopen("lz.dat", "w"); 
	FILE* pfsqr = fopen("lsqr.dat", "w");
	FILE* pfls = fopen("ls.dat", "w");

	/************************************************************/
	/******************    CALC HOLE <L>     ********************/
	/************************************************************/

  printf("Hole States:\n");
	printf("i       j         Lx.re           Ly.re           Lz.re          L^2.re           L.S.re\n"); 
	fflush(0);
	
  nvtxRangePushA("Calc qp hole L elems");
	// Compute <j|Lx|i>, <j|Ly|i>, <j|Lz|i>, <j|Lx^2|i>, <j|Ly^2|i>, <j|Lz^2|i>
	for (i = 0; i < n_ho; i++){
    nvtxRangePushA("loop over i");
		i_st = i * stlen;
		
    nvtxRangePushA("Compute L by FFTs");
		// These computations require FFT with small (max: 500^3) grids -> performed on CPU
		// 1) Compute Lx|i>, Ly|i>, Lz|i>
		//spin up part
		l_operator(&Lxpsi[0], &Lypsi[0], &Lzpsi[0], &psi_qp[i_st], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
		//spin dn part
		l_operator(&Lxpsi[ngrid], &Lypsi[ngrid], &Lzpsi[ngrid], &psi_qp[i_st + cplx_idx * ngrid], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);

		// Compute Lx^2|i>
		//spin up part
		l_operator(&Lx2psi[0], &temp1[0], &temp2[0], &Lxpsi[0], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
		//spin dn part
		l_operator(&Lx2psi[ngrid], &temp1[0], &temp2[0], &Lxpsi[ngrid], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);

		// Compute Ly^2|i>
		//spin up part
		l_operator(&temp1[0], &Ly2psi[0], &temp2[0], &Lypsi[0], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
		//spin dn part
		l_operator(&temp1[0],&Ly2psi[ngrid],&temp2[0], &Lypsi[ngrid], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);

		// Compute Lz^2|i>
		//spin up part
		l_operator(&temp1[0], &temp2[0], &Lz2psi[0], &Lzpsi[0], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
		//spin dn part
		l_operator(&temp1[0], &temp2[0], &Lz2psi[ngrid], &Lzpsi[ngrid], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
    nvtxRangePop();
		#pragma omp parallel for private(j_st, jgr, jsgr)
		for (j = i; j < n_ho; j++){
      nvtxRangePushA("loop over a");
			j_st = j * stlen;

			// Declare variables on the GPU

			double psi_jur;
			double psi_jui;
			double psi_jdr;
			double psi_jdi;

			double lx_re,    lx_im;
			double ly_re,    ly_im;
			double lz_re,    lz_im;
			double lsqr_re,  lsqr_im;
			double ldots_re, ldots_im;

      double Lxur, Lxui;
      double Lyur, Lyui;
      double Lzur, Lzui;
      double Lxdr, Lxdi;
      double Lydr, Lydi;
      double Lzdr, Lzdi;
			
			lx_re    =  lx_im     = 0.0;
			ly_re    =  ly_im     = 0.0;
			lz_re    =  lz_im     = 0.0;
			lsqr_re  =  lsqr_im   = 0.0;
			ldots_re =  ldots_im  = 0.00;
			
			// Perform integration
			#pragma omp simd
			for(jgr = 0; jgr < ngrid; jgr++){
        // Handle indexing
        jgur = cplx_idx * jgr;
        jgui = jgur + 1;
        jgdr = jgur + cplx_idx * ngrid;
        jgdi = jgdr + 1;
        jsgr = jgr + ngrid;
				
				// Load values of psi_j
				psi_jur = psi_qp[j_st + jgur];
				psi_jui = psi_qp[j_st + jgui];
				psi_jdr = psi_qp[j_st + jgdr];
				psi_jdi = psi_qp[j_st + jgdi];

				// Load values of L|psi>
				Lxur = Lxpsi[jgr].re;    Lxui = Lxpsi[jgr].im; // up
        Lyur = Lypsi[jgr].re;    Lyui = Lypsi[jgr].im;
        Lzur = Lzpsi[jgr].re;    Lzui = Lzpsi[jgr].im;
        Lxdr = Lxpsi[jsgr].re;   Lxdi = Lxpsi[jsgr].im; // down
        Lydr = Lypsi[jsgr].re;   Lydi = Lypsi[jsgr].im;
        Lzdr = Lzpsi[jsgr].re;   Lzdi = Lzpsi[jsgr].im;

        // Load values of L^2|psi>
        // Lx2psi[jgr].re

				// Handle up spin
				// Holes are like time-reversed electrons, so -^*
				// -<j|L_x|i>^*
				lx_re -= psi_jur * Lxur + psi_jui * Lxui;
				lx_im += psi_jur * Lxui - psi_jui * Lxur;
				// -<j|L_y|i>^*
				ly_re -= psi_jur * Lyur + psi_jui * Lyui;
				ly_im += psi_jur * Lyui - psi_jui * Lyur;
				// -<j|L_z|i>^*
				lz_re -= psi_jur * Lzur + psi_jui * Lzui;
				lz_im += psi_jur * Lzui - psi_jui * Lzur;
				
				// <j|L_x^2|i>^*
				lsqr_re += psi_jur * Lx2psi[jgr].re + psi_jui * Lx2psi[jgr].im;
				lsqr_im -= psi_jur * Lx2psi[jgr].im - psi_jui * Lx2psi[jgr].re;
				// <j|L_y^2|i>^*
				lsqr_re += psi_jur * Ly2psi[jgr].re + psi_jui * Ly2psi[jgr].im;
				lsqr_im -= psi_jur * Ly2psi[jgr].im - psi_jui * Ly2psi[jgr].re;
				// <j|L_z^2|i>^*
				lsqr_re += psi_jur * Lz2psi[jgr].re + psi_jui * Lz2psi[jgr].im;
				lsqr_im -= psi_jur * Lz2psi[jgr].im - psi_jui * Lz2psi[jgr].re;

				// Handle down spin
				// Holes are like time-reversed electrons, so -^*
				// -<j|L_x|i>^*
				lx_re -= psi_jdr * Lxdr + psi_jdi * Lxdi;
				lx_im += psi_jdr * Lxdi - psi_jdi * Lxdr;
				// -<j|L_y|i>^*
				ly_re -= psi_jdr * Lydr + psi_jdi * Lydi;
				ly_im += psi_jdr * Lydi - psi_jdi * Lydr;
				// -<j|L_z|i>^*
				lz_re -= psi_jdr * Lzdr + psi_jdi * Lzdi;
				lz_im += psi_jdr * Lzdi - psi_jdi * Lzdr;
				
				// <j|L_x^2|i>^*
				lsqr_re += psi_jdr * Lx2psi[jsgr].re + psi_jdi * Lx2psi[jsgr].im;
				lsqr_im -= psi_jdr * Lx2psi[jsgr].im - psi_jdi * Lx2psi[jsgr].re;
				// <j|L_y^2|i>^*
				lsqr_re += psi_jdr * Ly2psi[jsgr].re + psi_jdi * Ly2psi[jsgr].im;
				lsqr_im -= psi_jdr * Ly2psi[jsgr].im - psi_jdi * Ly2psi[jsgr].re;
				// <j|L_z^2|i>^*
				lsqr_re += psi_jdr * Lz2psi[jsgr].re + psi_jdi * Lz2psi[jsgr].im;
				lsqr_im -= psi_jdr * Lz2psi[jsgr].im - psi_jdi * Lz2psi[jsgr].re;

				// Compute LdotS, the dot product of orbital and spin angular momentum
				//sx.lx
				//<j|r,dn> * <r,up|Lxi>
				ldots_re += psi_jdr * Lxur + psi_jdi * Lxui;
				ldots_im += psi_jdr * Lxui - psi_jdi * Lxur;
				//<j|r,up> * <r,dn|Lxi>
				ldots_re += psi_jur * Lxdr + psi_jui * Lxdi;
				ldots_im += psi_jur * Lxdi - psi_jui * Lxdr;

				//sy.ly
				//i*<j|r,dn> * <r,up|Lyi>
				ldots_im += psi_jdr * Lyur + psi_jdi * Lyui;
				ldots_re -= psi_jdr * Lyui - psi_jdi * Lyur;
				//-i*<j|r,up> * <r,dn|Lyi>
				ldots_im -= psi_jur * Lydr + psi_jui * Lydi;
				ldots_re += psi_jur * Lydi - psi_jui * Lydr;

        //sz.lz
        //<j|r,up> * <r,up|Lzi>
				ldots_re += psi_jur * Lzur + psi_jui * Lzui;
				ldots_im += psi_jur * Lzui - psi_jui * Lzur;
				//<j|r,dn> * <r,dn|Lzi>
				ldots_re -= psi_jdr * Lzdr + psi_jdi * Lzdi;
				ldots_im -= psi_jdr * Lzdi - psi_jdi * Lzur;
			}
			
			// normalize and store in allocated array
			l_mom [i*n_ho + j].x_re = lx_re   * dv;  l_mom [i*n_ho + j].x_im = lx_im   * dv;
			l_mom [i*n_ho + j].y_re = ly_re   * dv;  l_mom [i*n_ho + j].y_im = ly_im   * dv;
			l_mom [i*n_ho + j].z_re = lz_re   * dv;  l_mom [i*n_ho + j].z_im = lz_im   * dv;
			l2_mom[i*n_ho + j].re   = lsqr_re * dv;  l2_mom[i*n_ho + j].im   = lsqr_im * dv;
			
			//normalize and complex conjugate LdotS
			LdotS[i*n_ho+j].re =  0.5 * ldots_re * dv; 
			LdotS[i*n_ho+j].im = -0.5 * ldots_im * dv;
      nvtxRangePop();
		}
    nvtxRangePop();
	}
  nvtxRangePop();

  /************************************************************/
	/*******************   FILL LOWER TRIANG  *******************/
	/************************************************************/

  for (i = 0; i < n_ho; i++){
    for (j = i + 1; j < n_ho; j++){
      // Populate the lower triangle with the complex conj. of upper tri
      
      long lt = j*n_ho + i;          // lower triangle
      long ut = i*n_ho + j;  // upper triangle
      // transpose
      l_mom [lt].x_re =   l_mom [ut].x_re;
      l_mom [lt].y_re =   l_mom [ut].y_re;
      l_mom [lt].z_re =   l_mom [ut].z_re;
      l2_mom[lt].re   =   l2_mom[ut].re;
      LdotS [lt].re   =   LdotS [ut].re;
      // imag part gets conjugated ^*
      l_mom [lt].x_im = - l_mom [ut].x_im;
      l_mom [lt].y_im = - l_mom [ut].y_im;
      l_mom [lt].z_im = - l_mom [ut].z_im;
      l2_mom[lt].im   = - l2_mom[ut].im;
      LdotS [lt].im   = - LdotS [ut].im;
      
    }
  }

	/************************************************************/
	/*******************   PRINT HOLE OUTPUT  *******************/
	/************************************************************/

  
	for (i = 0; i < n_ho; i++){
		for (j = 0; j < n_ho; j++){
			fprintf(pfx,   "%ld\t%ld\t%lf\t%lf\n", i, j, l_mom [i*n_ho+j].x_re, l_mom [i*n_ho+j].x_im);
			fprintf(pfy,   "%ld\t%ld\t%lf\t%lf\n", i, j, l_mom [i*n_ho+j].y_re, l_mom [i*n_ho+j].y_im);
			fprintf(pfz,   "%ld\t%ld\t%lf\t%lf\n", i, j, l_mom [i*n_ho+j].z_re, l_mom [i*n_ho+j].z_im);
			fprintf(pfsqr, "%ld\t%ld\t%lf\t%lf\n", i, j, l2_mom[i*n_ho+j].re,   l2_mom[i*n_ho+j].im);
			fprintf(pfls,  "%ld\t%ld\t%lf\t%lf\n", i, j, LdotS [i*n_ho+j].re,   LdotS [i*n_ho+j].im); 
			
			if(i==j){
				printf("%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\n",
					i, j, l_mom[i*n_ho+j].x_re, l_mom[i*n_ho+j].y_re, l_mom[i*n_ho +j].z_re,
					l2_mom[i*n_ho+j].re, LdotS[i*n_ho+j].re
				);
			}
		}
	}
	fflush(0);

	/************************************************************/
	/******************    CALC ELEC <L>     ********************/
	/************************************************************/

	printf("Electron States:\n");
	printf("a       b         Lx.re           Ly.re           Lz.re          L^2.re           L.S.re\n");
	fflush(0);
  nvtxRangePushA("Calc qp elec L elems");
	for (a = lidx; a < lidx + n_el; a++){
    nvtxRangePushA("loop over a");
		a_st = a * stlen;
		
    nvtxRangePushA("Compute L by FFTs");
		//spin up part
		l_operator(Lxpsi, Lypsi, Lzpsi, &psi_qp[a_st], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
		//spin dn part
		l_operator(&Lxpsi[ngrid], &Lypsi[ngrid], &Lzpsi[ngrid], &psi_qp[a_st + cplx_idx * ngrid], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
		
		//Lxsqr parts
		//spin up part
		l_operator(Lx2psi, &temp1[0], &temp2[0], Lxpsi, g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
		//spin dn part
		l_operator(&Lx2psi[ngrid], &temp1[0], &temp2[0], &Lxpsi[ngrid], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);

		//Lysqr parts
		//spin up part
		l_operator(&temp1[0], Ly2psi, &temp2[0], Lypsi, g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
		//spin dn part
		l_operator(&temp1[0], &Ly2psi[ngrid], &temp2[0], &Lypsi[ngrid], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);

		//Lzsqr parts
		//spin up part
		l_operator(&temp1[0], &temp2[0], &Lz2psi[0], &Lzpsi[0], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
		//spin dn part
		l_operator(&temp1[0], &temp2[0], &Lz2psi[ngrid], &Lzpsi[ngrid], g_vecs, grid, ist, par, planfw, planbw, fftwpsi);
		nvtxRangePop();

		#pragma omp parallel for private(b, b_st, jgr, jsgr, idx)
		for (b = a; b < lidx + n_el; b++){
      nvtxRangePushA("loop over b");
			b_st = b * stlen;

			// Declare variables on the GPU

			double psi_bur;
			double psi_bui;
			double psi_bdr;
			double psi_bdi;

			double lx_re,    lx_im;
			double ly_re,    ly_im;
			double lz_re,    lz_im;
			double lsqr_re,  lsqr_im;
			double ldots_re, ldots_im;

      double Lxur, Lxui;
      double Lyur, Lyui;
      double Lzur, Lzui;
      double Lxdr, Lxdi;
      double Lydr, Lydi;
      double Lzdr, Lzdi;
			
			lx_re    =  lx_im     = 0.0;
			ly_re    =  ly_im     = 0.0;
			lz_re    =  lz_im     = 0.0;
			lsqr_re  =  lsqr_im   = 0.0;
			ldots_re =  ldots_im  = 0.00;
			
			idx = sqr(n_ho) + (a - lidx) * n_el + (b - lidx);
      
      #pragma omp simd
			for(jgr = 0; jgr < ngrid; jgr++){
        // Handle indexing
        jgur = cplx_idx * jgr;
        jgui = jgur + 1;
        jgdr = jgur + cplx_idx * ngrid;
        jgdi = jgdr + 1;
        jsgr = jgr + ngrid;
				
				// Load values of psi_b
				psi_bur = psi_qp[b_st + jgur];
				psi_bui = psi_qp[b_st + jgui];
				psi_bdr = psi_qp[b_st + jgdr];
				psi_bdi = psi_qp[b_st + jgdi];

        // Load values of L|psi>
				Lxur = Lxpsi[jgr].re;    Lxui = Lxpsi[jgr].im; // up
        Lyur = Lypsi[jgr].re;    Lyui = Lypsi[jgr].im;
        Lzur = Lzpsi[jgr].re;    Lzui = Lzpsi[jgr].im;
        Lxdr = Lxpsi[jsgr].re;   Lxdi = Lxpsi[jsgr].im; // down
        Lydr = Lypsi[jsgr].re;   Lydi = Lypsi[jsgr].im;
        Lzdr = Lzpsi[jsgr].re;   Lzdi = Lzpsi[jsgr].im;
				
				// Handle up spin
				// <b|L_x|a>
				lx_re += psi_bur * Lxur + psi_bui * Lxui;
				lx_im += psi_bur * Lxui - psi_bui * Lxur;
				// <b|L_y|a>
				ly_re += psi_bur * Lyur + psi_bui * Lyui;
				ly_im += psi_bur * Lyui - psi_bui * Lyur;
				// <b|L_z|a>
				lz_re += psi_bur * Lzur + psi_bui * Lzui;
				lz_im += psi_bur * Lzui - psi_bui * Lzur;
				
				// <b|L_x^2|a>
				lsqr_re += psi_bur * Lx2psi[jgr].re + psi_bui * Lx2psi[jgr].im;
				lsqr_im += psi_bur * Lx2psi[jgr].im - psi_bui * Lx2psi[jgr].re;
				// <b|L_y^2|a>
				lsqr_re += psi_bur * Ly2psi[jgr].re + psi_bui * Ly2psi[jgr].im;
				lsqr_im += psi_bur * Ly2psi[jgr].im - psi_bui * Ly2psi[jgr].re;
				// <b|L_z^2|a>
				lsqr_re += psi_bur * Lz2psi[jgr].re + psi_bui * Lz2psi[jgr].im;
				lsqr_im += psi_bur * Lz2psi[jgr].im - psi_bui * Lz2psi[jgr].re;

				// Handle down spin
				// <b|L_x|a>
				lx_re += psi_bdr * Lxdr + psi_bdi * Lxdi;
				lx_im += psi_bdr * Lxdi - psi_bdi * Lxdr;
				// <b|L_y|a>
				ly_re += psi_bdr * Lydr + psi_bdi * Lydi;
				ly_im += psi_bdr * Lydi - psi_bdi * Lydr;
				// <b|L_z|a>
				lz_re += psi_bdr * Lzdr + psi_bdi * Lzdi;
				lz_im += psi_bdr * Lzdi - psi_bdi * Lzdr;
				
				// <b|L_x^2|a>
				lsqr_re += psi_bdr * Lx2psi[jsgr].re + psi_bdi * Lx2psi[jsgr].im;
				lsqr_im += psi_bdr * Lx2psi[jsgr].im - psi_bdi * Lx2psi[jsgr].re;
				// <b|L_y^2|a>
				lsqr_re += psi_bdr * Ly2psi[jsgr].re + psi_bdi * Ly2psi[jsgr].im;
				lsqr_im += psi_bdr * Ly2psi[jsgr].im - psi_bdi * Ly2psi[jsgr].re;
				// <b|L_z^2|a>
				lsqr_re += psi_bdr * Lz2psi[jsgr].re + psi_bdi * Lz2psi[jsgr].im;
				lsqr_im += psi_bdr * Lz2psi[jsgr].im - psi_bdi * Lz2psi[jsgr].re;

				// Compute LdotS for the electron states
				//s_mom.xlx
				//<b|r,dn> * <r,up|Lxa>
				ldots_re += psi_bdr * Lxur + psi_bdi * Lxui;
				ldots_im += psi_bdr * Lxui - psi_bdi * Lxur;
				//<b|r,up> * <r,dn|Lxa>
				ldots_re += psi_bur * Lxdr + psi_bui * Lxdi;
				ldots_im += psi_bur * Lxdi - psi_bui * Lxdr;

				//s_mom.yly
				//i*<j|r,dn> * <r,up|Lyi>
				ldots_im += psi_bdr * Lyur + psi_bdi * Lyui;
				ldots_re -= psi_bdr * Lyui - psi_bdi * Lyur;
				//-i*<j|r,up> * <r,dn|Lyi>
				ldots_im -= psi_bur * Lydr + psi_bui * Lydi;
				ldots_re += psi_bur * Lydi - psi_bui * Lydr;

        //s_mom.zlz
        //<j|r,up> * <r,up|Lzi>
				ldots_re += psi_bur * Lzur + psi_bui * Lzui;
				ldots_im += psi_bur * Lzui - psi_bui * Lzur;
				//<j|r,dn> * <r,dn|Lzi>
				ldots_re -= psi_bdr * Lzdr + psi_bdi * Lzdi;
				ldots_im -= psi_bdr * Lzdi - psi_bdi * Lzdr;
			}
			
			// Normalize and store in alloc'd arrays
			l_mom [idx].x_re = lx_re   * dv; l_mom [idx].x_im = lx_im   * dv;
			l_mom [idx].y_re = ly_re   * dv; l_mom [idx].y_im = ly_im   * dv;
			l_mom [idx].z_re = lz_re   * dv; l_mom [idx].z_im = lz_im   * dv;
			l2_mom[idx].re   = lsqr_re * dv; l2_mom[idx].im   = lsqr_im * dv;
			LdotS[idx].re = 0.5 * ldots_re * dv;
			LdotS[idx].im = 0.5 * ldots_im * dv;
      nvtxRangePop();
		}
    nvtxRangePop();
	}
  nvtxRangePop();

  /************************************************************/
	/*******************   FILL LOWER TRIANG  *******************/
	/************************************************************/

  for (a = lidx; a < lidx + n_el; a++){
    for (b = a + 1; b < lidx + n_el; b++){
      // Populate the lower triangle with the complex conj. of upper tri
      
      long lt = sqr(n_ho) + (b-lidx)*n_el + (a-lidx);          // lower triangle
      long ut = sqr(n_ho) + (a-lidx)*n_el + (b-lidx);  // upper triangle
      // transpose
      l_mom [lt].x_re =   l_mom [ut].x_re;
      l_mom [lt].y_re =   l_mom [ut].y_re;
      l_mom [lt].z_re =   l_mom [ut].z_re;
      l2_mom[lt].re   =   l2_mom[ut].re;
      LdotS [lt].re   =   LdotS [ut].re;
      // imag part gets conjugated ^*
      l_mom [lt].x_im = - l_mom [ut].x_im;
      l_mom [lt].y_im = - l_mom [ut].y_im;
      l_mom [lt].z_im = - l_mom [ut].z_im;
      l2_mom[lt].im   = - l2_mom[ut].im;
      LdotS [lt].im   = - LdotS [ut].im;
    }
  }

	/************************************************************/
	/*******************   PRINT HOLE OUTPUT  *******************/
	/************************************************************/

	for (a = lidx; a < lidx + n_el; a++){
		for (b = lidx; b < lidx + n_el; b++){
			idx = sqr(n_ho) + (a - lidx) * n_el + (b - lidx);
			fprintf(pfx,"%ld\t%ld\t%lf\t%lf\n",   a, b, l_mom [idx].x_re, l_mom[idx].x_im);
			fprintf(pfy,"%ld\t%ld\t%lf\t%lf\n",   a, b, l_mom [idx].y_re, l_mom[idx].y_im);
			fprintf(pfz,"%ld\t%ld\t%lf\t%lf\n",   a, b, l_mom [idx].z_re, l_mom[idx].z_im);
			fprintf(pfsqr,"%ld\t%ld\t%lf\t%lf\n", a, b, l2_mom[idx].re,   l2_mom[idx].im);
			fprintf(pfls,"%ld\t%ld\t%lf\t%lf\n",  a, b, LdotS [idx].re,   LdotS[idx].im);
			
			if(a==b){
				printf("%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\n",
					a, b, l_mom[idx].x_re, l_mom[idx].y_re, l_mom[idx].z_re,
					l2_mom[idx].re, LdotS[idx].re
				);
			}
		}
	}
	fflush(0);
	
	/************************************************************/
	/*******************   FREE DYNAMIC MEM   *******************/
	/************************************************************/

	fclose(pfx);
	fclose(pfy);
	fclose(pfz);
	fclose(pfsqr);
	fclose(pfls);
	
	free(gx);
	free(gy);
	free(gz);
	free(g_vecs);

	free(Lxpsi);
	free(Lypsi);
	free(Lzpsi);
	free(Lx2psi);
	free(Ly2psi);
	free(Lz2psi);
	free(temp1);
	free(temp2);

	fftw_free(fftwpsi);
	fftw_destroy_plan(planfw);
	fftw_destroy_plan(planbw);

	return;
}


/************************************************************/
//Calculate the three vector components of the action of the L operator on the 
//spatial part of the grid (no spin part)
void l_operator(
	zomplex* Lxpsi, 
	zomplex* Lypsi, 
	zomplex* Lzpsi, 
	zomplex* psi_qp,
	double*  g_vecs,
	grid_st *grid, 
	index_st *ist, 
	par_st *par, 
	fftw_plan_loc planfw, 
	fftw_plan_loc planbw, 
	fftw_complex* fftwpsi){
	
	double x, y, z;
	long jx, jy, jz, jyz, jxyz;
	zomplex px, py, pz;


	// Take derivative along x direction
	p_operator("X", g_vecs, psi_qp, Lxpsi, grid, ist, par, planfw, planbw, fftwpsi);
	// Take derivative along y direction
	p_operator("Y", g_vecs, psi_qp, Lypsi, grid, ist, par, planfw, planbw, fftwpsi);
	// Take derivative along x direction
	p_operator("Z", g_vecs, psi_qp, Lzpsi, grid, ist, par, planfw, planbw, fftwpsi);

	/*** Now do the cross product part at each grid point***/
	omp_set_dynamic(0);
  	omp_set_num_threads(ist->nthreads);
	#pragma omp parallel for private (jz,jy,jyz,jx,jxyz,px,py,pz,x,y,z)
  	for (jz = 0; jz < grid->nz; jz++) { 
  		z = grid->z[jz];
  		for (jy = 0;jy<grid->ny;jy++){
  			jyz = grid->nx * (grid->ny * jz + jy);
  			y = grid->y[jy];
  			for(jx = 0; jx<grid->nx;jx++){
  				x = grid->x[jx];
				jxyz = jyz + jx;
				//copy over tmp varibales
				px.re = Lxpsi[jxyz].re; px.im = Lxpsi[jxyz].im;
				py.re = Lypsi[jxyz].re; py.im = Lypsi[jxyz].im;
				pz.re = Lzpsi[jxyz].re; pz.im = Lzpsi[jxyz].im;
				
				Lxpsi[jxyz].re = (y*pz.re - z*py.re);
				Lxpsi[jxyz].im = (y*pz.im - z*py.im);

				Lypsi[jxyz].re = (z*px.re - x*pz.re);
				Lypsi[jxyz].im = (z*px.im - x*pz.im);

				Lzpsi[jxyz].re = (x*py.re - y*px.re);
				Lzpsi[jxyz].im = (x*py.im - y*px.im);
			}
		}
	}


	return;
}

/*************************************************************************************************/

void init_g_vecs(double *kindex, double *kx, double *ky, double *kz, grid_st *grid, index_st *ist, par_st *par){

	long jx, jy, jz, jyz, jgrid;
	// The k-vectors in each direction are initialized so that
	// they are stored as - k. e.g. The kx values of positive 
	// x values are made negative; simplifies p_operator function
	for (kx[0] = 0.0, jx = 1; jx <= grid->nx / 2; jx++){
    	kx[grid->nx-jx] = -1.00 * (kx[jx] = (double)(jx) * grid->dkx * 
    		grid->nx_1 * grid->ny_1 * grid->nz_1);
  	}
	// for (jx = 0; jx < grid->nx; jx++){
	// 	printf("x[%ld] = %g , kx = %g\n", jx, grid->x[jx], kx[jx]);
	// }

  	for (ky[0] = 0.0, jy = 1; jy <= grid->ny / 2; jy++){
    	ky[grid->ny-jy] = -1.00 * (ky[jy] = (double)(jy) * grid->dky *
			grid->nx_1 * grid->ny_1 * grid->nz_1);
  	}

  	for (kz[0] = 0.0, jz = 1; jz <= grid->nz / 2; jz++){
    	kz[grid->nz-jz] = -1.00 * (kz[jz] = (double)(jz) * grid->dkz *
			grid->nx_1 * grid->ny_1 * grid->nz_1);
  	}
	
	for (jz = 0; jz < grid->nz; jz++) {
    	for (jy = 0; jy < grid->ny; jy++) {
      	jyz = grid->nx * (grid->ny * jz + jy);
      	for (jx = 0; jx < grid->nx; jx++) {   
        	jgrid = jyz + jx;
        	kindex[3*jgrid]   = kx[jx]; 
        	kindex[3*jgrid+1] = ky[jy];
        	kindex[3*jgrid+2] = kz[jz];
      } 
    }
  }

  return;

}

/*************************************************************************************************/

void p_operator(char* direc, double *kindex, zomplex *psi, zomplex *Lpsi, grid_st *grid, index_st *ist, par_st *par, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi){
	/*** First use the fft to get the action of the p operator on each axis 
	* p = (-i d/dx) imag part from the fft definition***/
	// To compute d/dx, we compute FT^-1[-i*k*FT(psi)]
	// multiplying this by -i, we get p = -i * FT^-1[ -i*k * FT(psi)]
	// p = FT^-1[ -k * FT(psi)]. We already initialized the k-vectors to be -k.
	// so here we just multiply by k.
	
	long jgrid;
	long dir_offset;

	if (strcmp(direc, "X") == 0){
		// printf("The deriv. is along X direction.\n");
		// printf("Setting dir_offset to 0.\n");
		dir_offset = 0;
	} else if (strcmp(direc, "Y") == 0){
		// printf("The deriv. is along Y direction.\n");
		// printf("Setting dir_offset to 1.\n");
		dir_offset = 1;
	} else if (strcmp(direc, "Z") == 0){
		// printf("The deriv. is along Z direction.\n");
		// printf("Setting dir_offset to 2.\n");
		dir_offset = 2;
	} else {
		printf("No direc recognized for derivative.\n");
		exit(EXIT_FAILURE);
	}

	// Copy psi to fftwpsi
  	memcpy(&fftwpsi[0], &psi[0], grid->ngrid*sizeof(fftwpsi[0]));
  	
	// FT from r-space to k-space
  	fftw_execute(planfw);
  	
  	omp_set_dynamic(0);
  	omp_set_num_threads(ist->nthreads);
  	#pragma omp parallel for private (jgrid)
  	for (jgrid = 0; jgrid < grid->ngrid; jgrid++) {
		//multiply by -k_x/y/z to get p along x/y/z-axis
		fftwpsi[jgrid][0] *= kindex[3*jgrid + dir_offset];
		fftwpsi[jgrid][1] *= kindex[3*jgrid + dir_offset];
	}

	// Inverse FT back to r-space
	fftw_execute(planbw);
	// Copy fftwpsi to psi to store Lx|psi_qp> into |Lxpsi>
	memcpy(&Lpsi[0], &fftwpsi[0], ist->ngrid*sizeof(Lpsi[0]));

	return;
	
}