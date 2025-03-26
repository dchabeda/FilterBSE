/****************************************************************************/

#include "fd.h"

/****************************************************************************/

void calc_spin_mtrx(xyz_st *s_mom, double *psi_qp, grid_st *grid, index_st *ist, par_st *par){
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
  *  [s_mom] array to hold components of the spin projection mat elems   *
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
 
  long a, i, b, j;
  long astate_idx, istate_idx, bstate_idx, jstate_idx;
  long jgrid, jgridup_re, jgridup_im, jgriddn_re, jgriddn_im; 
  long mtrx_idx;
  FILE *pfx, *pfy, *pfz;

  pfx = fopen("sx.dat", "w"); pfy = fopen("sy.dat", "w"); pfz = fopen("sz.dat", "w");
  //calculate spin matrix elements between all occupied (hole) orbitals
	for(i = 0; i < ist->n_holes; i++){
        istate_idx = i * ist->nspinngrid * ist->complex_idx;
		for (j = 0; j < ist->n_holes; j++){
      		zomplex sx, sy, sz;
      		sx.re = sx.im = sy.re = sy.im = sz.re = sz.im = 0.0;

			for (jgrid = 0; jgrid < ist->ngrid; jgrid++) {
        jgridup_re = ist->complex_idx * jgrid;
        jgridup_im = jgridup_re + 1;
        jgriddn_re = jgridup_re + ist->ngrid*ist->complex_idx;
        jgriddn_im = jgriddn_re + 1;
        jstate_idx = j * ist->nspinngrid * ist->complex_idx;
        
        //Spin x part
        //<j|r,dn> * <r,up|i>
        sx.re += psi_qp[jstate_idx+jgriddn_re] * psi_qp[istate_idx+jgridup_re]
                + psi_qp[jstate_idx+jgriddn_im] * psi_qp[istate_idx+jgridup_im];
        sx.im += psi_qp[jstate_idx+jgriddn_re] * psi_qp[istate_idx+jgridup_im]
                - psi_qp[jstate_idx+jgriddn_im] * psi_qp[istate_idx+jgridup_re];
        //<j|r,up> * <r,dn|i>
        sx.re += psi_qp[jstate_idx+jgridup_re] * psi_qp[istate_idx+jgriddn_re]
                + psi_qp[jstate_idx+jgridup_im] * psi_qp[istate_idx+jgriddn_im];
        sx.im += psi_qp[jstate_idx+jgridup_re] * psi_qp[istate_idx+jgriddn_im]
                - psi_qp[jstate_idx+jgridup_im] * psi_qp[istate_idx+jgriddn_re];
        
        //Spin y part
        //i*<j|r,dn> * <r,up|i>
        sy.im += psi_qp[jstate_idx+jgriddn_re] * psi_qp[istate_idx+jgridup_re]
                + psi_qp[jstate_idx+jgriddn_im] * psi_qp[istate_idx+jgridup_im];
        sy.re -= psi_qp[jstate_idx+jgriddn_re] * psi_qp[istate_idx+jgridup_im]
                - psi_qp[jstate_idx+jgriddn_im] * psi_qp[istate_idx+jgridup_re];
        //-i*<j|r,up> * <r,dn|i>
        sy.im -= psi_qp[jstate_idx+jgridup_re] * psi_qp[istate_idx+jgriddn_re]
                + psi_qp[jstate_idx+jgridup_im] * psi_qp[istate_idx+jgriddn_im];
        sy.re += psi_qp[jstate_idx+jgridup_re] * psi_qp[istate_idx+jgriddn_im]
                - psi_qp[jstate_idx+jgridup_im] * psi_qp[istate_idx+jgriddn_re];
        
        //Spin z part
        //<j|r,up> * <r,up|i>
        sz.re += psi_qp[jstate_idx+jgridup_re] * psi_qp[istate_idx+jgridup_re]
                + psi_qp[jstate_idx+jgridup_im] * psi_qp[istate_idx+jgridup_im];
        sz.im += psi_qp[jstate_idx+jgridup_re] * psi_qp[istate_idx+jgridup_im]
                - psi_qp[jstate_idx+jgridup_im] * psi_qp[istate_idx+jgridup_re];
        //<j|r,dn> * <r,dn|i>
        sz.re -= psi_qp[jstate_idx+jgriddn_re] * psi_qp[istate_idx+jgriddn_re]
                + psi_qp[jstate_idx+jgriddn_im] * psi_qp[istate_idx+jgriddn_im];
        sz.im -= psi_qp[jstate_idx+jgriddn_re] * psi_qp[istate_idx+jgriddn_im]
                - psi_qp[jstate_idx+jgriddn_im] * psi_qp[istate_idx+jgriddn_re];
		  }
			
			//multiply all by (1/2)*dV and complex conjugate bc holes are time reversed
			s_mom[i*ist->n_holes+j].x_re = -0.5 * sx.re * grid->dv;
			s_mom[i*ist->n_holes+j].x_im =  0.5 * sx.im * grid->dv;
			s_mom[i*ist->n_holes+j].y_re = -0.5 * sy.re * grid->dv;
			s_mom[i*ist->n_holes+j].y_im =  0.5 * sy.im * grid->dv;
			s_mom[i*ist->n_holes+j].z_re = -0.5 * sz.re * grid->dv;
			s_mom[i*ist->n_holes+j].z_im =  0.5 * sz.im * grid->dv;

			fprintf (pfx,"%ld %ld % .10g % .10g\n",i,j,s_mom[i*ist->n_holes+j].x_re, s_mom[i*ist->n_holes+j].x_im); fflush(pfx);
      fprintf (pfy,"%ld %ld % .10g % .10g\n",i,j,s_mom[i*ist->n_holes+j].y_re, s_mom[i*ist->n_holes+j].y_im); fflush(pfy);
      fprintf (pfz,"%ld %ld % .10g % .10g\n",i,j,s_mom[i*ist->n_holes+j].z_re, s_mom[i*ist->n_holes+j].z_im); fflush(pfz);
		}
	}

  //calculate spin matrix elements between all unoccupied (electron) orbitals
	for (a = ist->lumo_idx; a < ist->lumo_idx+ist->n_elecs; a++){
		for (b = ist->lumo_idx; b < ist->lumo_idx+ist->n_elecs; b++){
			zomplex sx, sy, sz;
      sx.re = sx.im = sy.re = sy.im = sz.re = sz.im = 0.0;

			for (jgrid = 0; jgrid < ist->ngrid; jgrid++) {
        jgridup_re = ist->complex_idx * jgrid;
        jgridup_im = jgridup_re + 1;
        jgriddn_re = jgridup_re + ist->ngrid*ist->complex_idx;
        jgriddn_im = jgriddn_re + 1;
        astate_idx = a * ist->nspinngrid * ist->complex_idx;
        bstate_idx = b * ist->nspinngrid * ist->complex_idx;
						
        //Spin x part
        //<b|r,dn> * <r,up|a>
        sx.re += psi_qp[bstate_idx+jgriddn_re] * psi_qp[astate_idx+jgridup_re]
                + psi_qp[bstate_idx+jgriddn_im] * psi_qp[astate_idx+jgridup_im];
        sx.im += psi_qp[bstate_idx+jgriddn_re] * psi_qp[astate_idx+jgridup_im]
                - psi_qp[bstate_idx+jgriddn_im] * psi_qp[astate_idx+jgridup_re];
        //<b|r,up> * <r,dn|a>
        sx.re += psi_qp[bstate_idx+jgridup_re] * psi_qp[astate_idx+jgriddn_re]
                + psi_qp[bstate_idx+jgridup_im] * psi_qp[astate_idx+jgriddn_im];
        sx.im += psi_qp[bstate_idx+jgridup_re] * psi_qp[astate_idx+jgriddn_im]
                - psi_qp[bstate_idx+jgridup_im] * psi_qp[astate_idx+jgriddn_re];
        
        //Spin y part
        //i*<b|r,dn> * <r,up|a>
        sy.im += psi_qp[bstate_idx+jgriddn_re] * psi_qp[astate_idx+jgridup_re]
                + psi_qp[bstate_idx+jgriddn_im] * psi_qp[astate_idx+jgridup_im];
        sy.re -= psi_qp[bstate_idx+jgriddn_re] * psi_qp[astate_idx+jgridup_im]
                - psi_qp[bstate_idx+jgriddn_im] * psi_qp[astate_idx+jgridup_re];
        //-i*<b|r,up> * <r,dn|a>
        sy.im -= psi_qp[bstate_idx+jgridup_re] * psi_qp[astate_idx+jgriddn_re]
                + psi_qp[bstate_idx+jgridup_im] * psi_qp[astate_idx+jgriddn_im];
        sy.re += psi_qp[bstate_idx+jgridup_re] * psi_qp[astate_idx+jgriddn_im]
                - psi_qp[bstate_idx+jgridup_im] * psi_qp[astate_idx+jgriddn_re];
        
        //Spin z part
        //<b|r,up> * <r,up|a>
        sz.re += psi_qp[bstate_idx+jgridup_re] * psi_qp[astate_idx+jgridup_re]
                + psi_qp[bstate_idx+jgridup_im] * psi_qp[astate_idx+jgridup_im];
        sz.im += psi_qp[bstate_idx+jgridup_re] * psi_qp[astate_idx+jgridup_im]
                - psi_qp[bstate_idx+jgridup_im] * psi_qp[astate_idx+jgridup_re];
        //<b|r,dn> * <r,dn|a>
        sz.re -= psi_qp[bstate_idx+jgriddn_re] * psi_qp[astate_idx+jgriddn_re]
                + psi_qp[bstate_idx+jgriddn_im] * psi_qp[astate_idx+jgriddn_im];
        sz.im -= psi_qp[bstate_idx+jgriddn_re] * psi_qp[astate_idx+jgriddn_im]
                - psi_qp[bstate_idx+jgriddn_im] * psi_qp[astate_idx+jgriddn_re];
			}

      mtrx_idx = sqr(ist->n_holes) + (a-ist->lumo_idx)*ist->n_elecs + (b-ist->lumo_idx);
      
			//divide all by 2
			s_mom[mtrx_idx].x_re = 0.5 * sx.re * grid->dv;
			s_mom[mtrx_idx].x_im = 0.5 * sx.im * grid->dv;
			s_mom[mtrx_idx].y_re = 0.5 * sy.re * grid->dv;
			s_mom[mtrx_idx].y_im = 0.5 * sy.im * grid->dv;
			s_mom[mtrx_idx].z_re = 0.5 * sz.re * grid->dv;
			s_mom[mtrx_idx].z_im = 0.5 * sz.im * grid->dv;
      
			fprintf (pfx,"%ld %ld % .10g % .10g\n",a,b,s_mom[mtrx_idx].x_re, s_mom[mtrx_idx].x_im);
      fprintf (pfy,"%ld %ld % .10g % .10g\n",a,b,s_mom[mtrx_idx].y_re, s_mom[mtrx_idx].y_im);
      fprintf (pfz,"%ld %ld % .10g % .10g\n",a,b,s_mom[mtrx_idx].z_re, s_mom[mtrx_idx].z_im);
		}
	}
  // Free memory allocated to FILE pointers
	fclose(pfx);
	fclose(pfy);
	fclose(pfz);

	return;
}

/****************************************************************************/

void calc_ang_mom_mtrx(xyz_st *l_mom, zomplex *l2_mom, zomplex *LdotS, double *psi_qp, grid_st *grid, index_st *ist, par_st *par){
	long istate_idx, jstate_idx, astate_idx, bstate_idx; 
	long jdridup_re, jdridup_im, jgriddn_re, jdriddn_im;
  	long i,j,a,b,jgrid,index, jgridup, jgriddn;
	fftw_plan_loc planfw, planbw; fftw_complex *fftwpsi;

	// Allocate memory for multithreaded FFT
	fftw_init_threads();
	fftw_plan_with_nthreads(ist->nthreads);

	fftwpsi = fftw_malloc(ist->ngrid*sizeof(fftw_complex));
    planfw = fftw_plan_dft_3d(grid->nz, grid->ny, grid->nx, fftwpsi, fftwpsi, FFTW_FORWARD, 0);
    planbw = fftw_plan_dft_3d(grid->nz, grid->ny, grid->nx, fftwpsi, fftwpsi, FFTW_BACKWARD, 0);


	zomplex *psi1 = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));;
	
	zomplex* Lxpsi = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));
	zomplex* Lypsi = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));
	zomplex* Lzpsi = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));

	zomplex* Lxsqrpsi = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));
	zomplex* Lysqrpsi = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));
	zomplex* Lzsqrpsi = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));
	zomplex* temp1 = (zomplex*) calloc(ist->ngrid, sizeof(zomplex));
	zomplex* temp2 = (zomplex*) calloc(ist->ngrid, sizeof(zomplex));
	
	FILE* pfx = fopen("lx.dat", "w"); 
	FILE* pfy = fopen("ly.dat", "w"); 
	FILE* pfz = fopen("lz.dat", "w"); 
	FILE* pfsqr = fopen("lsqr.dat", "w");
	FILE* pfls = fopen("ls.dat", "w");
	
  	printf("Hole States:\n");
	printf("i       j         Lx.re           Ly.re           Lz.re          L^2.re           L.S.re\n");
	for (i = 0; i < ist->n_holes; i++){
		istate_idx = i * ist->nspinngrid * ist->complex_idx;
		memcpy(psi1, &psi_qp[istate_idx], ist->nspinngrid*ist->complex_idx*sizeof(psi_qp[0]));
		
		//spin up part
		l_operator(&Lxpsi[0],&Lypsi[0],&Lzpsi[0], &psi1[0], grid, ist, par, planfw, planbw, fftwpsi);
		//spin dn part
		l_operator(&Lxpsi[ist->ngrid],&Lypsi[ist->ngrid],&Lzpsi[ist->ngrid], &psi1[ist->ngrid], grid, ist, par, planfw, planbw, fftwpsi);

		//Lxsqr parts
		//spin up part
		l_operator(&Lxsqrpsi[0],&temp1[0],&temp2[0], &Lxpsi[0], grid, ist, par, planfw, planbw, fftwpsi);
		//spin dn part
		l_operator(&Lxsqrpsi[ist->ngrid],&temp1[0],&temp2[0], &Lxpsi[ist->ngrid], grid, ist, par, planfw, planbw, fftwpsi);

		//Lysqr parts
		//spin up part
		l_operator(&temp1[0],&Lysqrpsi[0],&temp2[0], &Lypsi[0], grid, ist, par, planfw, planbw, fftwpsi);
		//spin dn part
		l_operator(&temp1[0],&Lysqrpsi[ist->ngrid],&temp2[0], &Lypsi[ist->ngrid], grid, ist, par, planfw, planbw, fftwpsi);

		//Lzsqr parts
		//spin up part
		l_operator(&temp1[0],&temp2[0],&Lzsqrpsi[0], &Lzpsi[0], grid, ist, par, planfw, planbw, fftwpsi);
		//spin dn part
		l_operator(&temp1[0],&temp2[0],&Lzsqrpsi[ist->ngrid], &Lzpsi[ist->ngrid], grid, ist, par, planfw, planbw, fftwpsi);

		#pragma omp parallel for private(j, jstate_idx, jgrid, jgridup, jgriddn)
		for (j = 0; j < ist->n_holes; j++){
			zomplex *psi2 = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));;
			zomplex psi_j_tmp;
			zomplex lx_tmp, ly_tmp, lz_tmp, lsqr_tmp;
			lx_tmp.re = lx_tmp.im = 0.0;
			ly_tmp.re = ly_tmp.im = 0.0;
			lz_tmp.re = lz_tmp.im = 0.0;
			lsqr_tmp.re = lsqr_tmp.im = 0.0;
			jstate_idx = j * ist->nspinngrid * ist->complex_idx;

			memcpy(psi2, &psi_qp[jstate_idx], ist->nspinngrid*ist->complex_idx*sizeof(psi_qp[0]));
		
			for(jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
				psi_j_tmp.re = psi2[jgrid].re;
				psi_j_tmp.im = psi2[jgrid].im;
				
				// -<j|L_x|i>^*
				lx_tmp.re -= Lxpsi[jgrid].re*psi_j_tmp.re + Lxpsi[jgrid].im*psi_j_tmp.im;
				lx_tmp.im += Lxpsi[jgrid].im*psi_j_tmp.re - Lxpsi[jgrid].re*psi_j_tmp.im;
				// -<j|L_y|i>^*
				ly_tmp.re -= Lypsi[jgrid].re*psi_j_tmp.re + Lypsi[jgrid].im*psi_j_tmp.im;
				ly_tmp.im += Lypsi[jgrid].im*psi_j_tmp.re - Lypsi[jgrid].re*psi_j_tmp.im;
				// -<j|L_z|i>^*
				lz_tmp.re -= Lzpsi[jgrid].re*psi_j_tmp.re + Lzpsi[jgrid].im*psi_j_tmp.im;
				lz_tmp.im += Lzpsi[jgrid].im*psi_j_tmp.re - Lzpsi[jgrid].re*psi_j_tmp.im;

				// <j|L_x^2|i>^*
				lsqr_tmp.re += Lxsqrpsi[jgrid].re*psi_j_tmp.re + Lxsqrpsi[jgrid].im*psi_j_tmp.im;
				lsqr_tmp.im -= Lxsqrpsi[jgrid].im*psi_j_tmp.re - Lxsqrpsi[jgrid].re*psi_j_tmp.im;
				// <j|L_y^2|i>^*
				lsqr_tmp.re += Lysqrpsi[jgrid].re*psi_j_tmp.re + Lysqrpsi[jgrid].im*psi_j_tmp.im; 
				lsqr_tmp.im -= Lysqrpsi[jgrid].im*psi_j_tmp.re - Lysqrpsi[jgrid].re*psi_j_tmp.im;
				// <j|L_z^2|i>^*
				lsqr_tmp.re += Lzsqrpsi[jgrid].re*psi_j_tmp.re + Lzsqrpsi[jgrid].im*psi_j_tmp.im;
				lsqr_tmp.im -= Lzsqrpsi[jgrid].im*psi_j_tmp.re - Lzsqrpsi[jgrid].re*psi_j_tmp.im;
			}
			// normalize and store in allocated array
			l_mom[i*ist->n_holes + j].x_re = lx_tmp.re * grid->dv; l_mom[i*ist->n_holes + j].x_im = lx_tmp.im * grid->dv;
			l_mom[i*ist->n_holes + j].y_re = ly_tmp.re * grid->dv; l_mom[i*ist->n_holes + j].y_im = ly_tmp.im * grid->dv;
			l_mom[i*ist->n_holes + j].z_re = lz_tmp.re * grid->dv; l_mom[i*ist->n_holes + j].z_im = lz_tmp.im * grid->dv;
			l2_mom[i*ist->n_holes + j].re = lsqr_tmp.re * grid->dv; l2_mom[i*ist->n_holes + j].im = lsqr_tmp.im * grid->dv;
			
			//
			//  Compute LdotS, the dot product of orbital and spin angular momentum
			//
			zomplex psi_j_up_tmp, psi_j_dn_tmp;
			zomplex ldots_tmp;
			ldots_tmp.re = ldots_tmp.im = 0.00;

			for(jgrid = 0; jgrid < ist->ngrid; jgrid++){
				jgridup = jgrid;
				jgriddn = jgrid+ist->ngrid;
				psi_j_up_tmp.re = psi2[jgridup].re;
				psi_j_up_tmp.im = psi2[jgridup].im;
				psi_j_dn_tmp.re = psi2[jgriddn].re;
				psi_j_dn_tmp.im = psi2[jgriddn].im;
				
				//sx.lx
				//<j|r,dn> * <r,up|Lxi>
				ldots_tmp.re += psi_j_dn_tmp.re * Lxpsi[jgridup].re
										  + psi_j_dn_tmp.im * Lxpsi[jgridup].im;
				ldots_tmp.im += psi_j_dn_tmp.re * Lxpsi[jgridup].im
										  - psi_j_dn_tmp.im * Lxpsi[jgridup].re;
				//<j|r,up> * <r,dn|Lxi>
				ldots_tmp.re += psi_j_up_tmp.re * Lxpsi[jgriddn].re
										  + psi_j_up_tmp.im * Lxpsi[jgriddn].im;
				ldots_tmp.im += psi_j_up_tmp.re * Lxpsi[jgriddn].im
										  - psi_j_up_tmp.im * Lxpsi[jgriddn].re;

				//sy.ly
				//i*<j|r,dn> * <r,up|Lyi>
				ldots_tmp.im += psi_j_dn_tmp.re * Lypsi[jgridup].re
							  			  + psi_j_dn_tmp.im * Lypsi[jgridup].im;
				ldots_tmp.re -= psi_j_dn_tmp.re * Lypsi[jgridup].im
							  			  - psi_j_dn_tmp.im * Lypsi[jgridup].re;
				//-i*<j|r,up> * <r,dn|Lyi>
				ldots_tmp.im -= psi_j_up_tmp.re * Lypsi[jgriddn].re
							  			  + psi_j_up_tmp.im * Lypsi[jgriddn].im;
				ldots_tmp.re += psi_j_up_tmp.re * Lypsi[jgriddn].im
							  			  - psi_j_up_tmp.im * Lypsi[jgriddn].re;

  			  	//sz.lz
  			  	//<j|r,up> * <r,up|Lzi>
				ldots_tmp.re += psi_j_up_tmp.re * Lzpsi[jgridup].re
										  + psi_j_up_tmp.im * Lzpsi[jgridup].im;
				ldots_tmp.im += psi_j_up_tmp.re * Lzpsi[jgridup].im
										  - psi_j_up_tmp.im * Lzpsi[jgridup].re;
				//<j|r,dn> * <r,dn|Lzi>
				ldots_tmp.re -= psi_j_dn_tmp.re * Lzpsi[jgriddn].re
										  + psi_j_dn_tmp.im * Lzpsi[jgriddn].im;
				ldots_tmp.im -= psi_j_dn_tmp.re * Lzpsi[jgriddn].im
										  - psi_j_dn_tmp.im * Lzpsi[jgriddn].re;
			}
			//normalize and complex conjugate
			LdotS[i*ist->n_holes+j].re = 0.5 * ldots_tmp.re * grid->dv; 
			LdotS[i*ist->n_holes+j].im = -0.5 * ldots_tmp.im * grid->dv;

			fprintf(pfx,"%ld\t%ld\t%lf\t%lf\n",i,j,l_mom[i*ist->n_holes+j].x_re,l_mom[i*ist->n_holes+j].x_im);
			fprintf(pfy,"%ld\t%ld\t%lf\t%lf\n",i,j,l_mom[i*ist->n_holes+j].y_re,l_mom[i*ist->n_holes+j].y_im);
			fprintf(pfz,"%ld\t%ld\t%lf\t%lf\n",i,j,l_mom[i*ist->n_holes+j].z_re,l_mom[i*ist->n_holes+j].z_im);
			fprintf(pfsqr,"%ld\t%ld\t%lf\t%lf\n",i,j,l2_mom[i*ist->n_holes+j].re,l2_mom[i*ist->n_holes+j].im);
			fprintf(pfls,"%ld\t%ld\t%lf\t%lf\n",i,j,LdotS[i*ist->n_holes+j].re,LdotS[i*ist->n_holes+j].im); 
			if(i==j){
				printf("%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\n",i,j,l_mom[i*ist->n_holes+j].x_re,l_mom[i*ist->n_holes+j].y_re,l_mom[i*ist->n_holes+j].z_re,l2_mom[i*ist->n_holes+j].re,LdotS[i*ist->n_holes+j].re);
			}

			free(psi2);
		}
	}

	printf("Electron States:\n");
	printf("a       b         Lx.re           Ly.re           Lz.re          L^2.re           L.S.re\n");
	for (a = ist->lumo_idx;a<ist->n_elecs+ist->lumo_idx;a++){
		astate_idx = a * ist->nspinngrid * ist->complex_idx;
		memcpy(psi1, &psi_qp[astate_idx], ist->nspinngrid*ist->complex_idx*sizeof(psi_qp[0]));
		
		//spin up part
		l_operator(&Lxpsi[0], &Lypsi[0], &Lzpsi[0], &psi1[0], grid, ist, par, planfw, planbw, fftwpsi);
		//spin dn part
		l_operator(&Lxpsi[ist->ngrid],&Lypsi[ist->ngrid],&Lzpsi[ist->ngrid], &psi1[ist->ngrid], grid, ist, par, planfw, planbw, fftwpsi);
		
		//Lxsqr parts
		//spin up part
		l_operator(&Lxsqrpsi[0], &temp1[0], &temp2[0], &Lxpsi[0], grid, ist, par, planfw, planbw, fftwpsi);
		//spin dn part
		l_operator(&Lxsqrpsi[ist->ngrid],&temp1[0],&temp2[0], &Lxpsi[ist->ngrid], grid, ist, par, planfw, planbw, fftwpsi);

		//Lysqr parts
		//spin up part
		l_operator(&temp1[0], &Lysqrpsi[0], &temp2[0], &Lypsi[0], grid, ist, par, planfw, planbw, fftwpsi);
		//spin dn part
		l_operator(&temp1[0],&Lysqrpsi[ist->ngrid],&temp2[0], &Lypsi[ist->ngrid], grid, ist, par, planfw, planbw, fftwpsi);

		//Lzsqr parts
		//spin up part
		l_operator(&temp1[0], &temp2[0], &Lzsqrpsi[0], &Lzpsi[0], grid, ist, par, planfw, planbw, fftwpsi);
		//spin dn part
		l_operator(&temp1[0],&temp2[0],&Lzsqrpsi[ist->ngrid], &Lzpsi[ist->ngrid], grid, ist, par, planfw, planbw, fftwpsi);
		
		#pragma omp parallel for private(b, bstate_idx, jgrid, jgridup, jgriddn, index)
		for (b = ist->lumo_idx; b < ist->n_elecs+ist->lumo_idx; b++){
			index = sqr(ist->n_holes)+(a-ist->lumo_idx)*ist->n_elecs+(b-ist->lumo_idx);
			zomplex *psi2 = (zomplex*) calloc(ist->nspinngrid, sizeof(zomplex));;
			zomplex psi_b_tmp;
			zomplex lx_tmp, ly_tmp, lz_tmp, lsqr_tmp;
			lx_tmp.re = lx_tmp.im = 0.0;
			ly_tmp.re = ly_tmp.im = 0.0;
			lz_tmp.re = lz_tmp.im = 0.0;
			lsqr_tmp.re = lsqr_tmp.im = 0.0;
			bstate_idx = b * ist->nspinngrid * ist->complex_idx;

			memcpy(psi2, &psi_qp[bstate_idx], ist->nspinngrid*ist->complex_idx*sizeof(psi_qp[0]));

			for(jgrid = 0; jgrid < ist->nspinngrid; jgrid++){
				psi_b_tmp.re = psi2[jgrid].re;
				psi_b_tmp.im = psi2[jgrid].im;
				
				// <b|L_x|a>
				lx_tmp.re += Lxpsi[jgrid].re*psi_b_tmp.re + Lxpsi[jgrid].im*psi_b_tmp.im;
				lx_tmp.im += Lxpsi[jgrid].im*psi_b_tmp.re - Lxpsi[jgrid].re*psi_b_tmp.im;
				// <b|L_y|a>
				ly_tmp.re += Lypsi[jgrid].re*psi_b_tmp.re + Lypsi[jgrid].im*psi_b_tmp.im;
				ly_tmp.im += Lypsi[jgrid].im*psi_b_tmp.re - Lypsi[jgrid].re*psi_b_tmp.im;
				// <b|L_z|a>
				lz_tmp.re += Lzpsi[jgrid].re*psi_b_tmp.re + Lzpsi[jgrid].im*psi_b_tmp.im;
				lz_tmp.im += Lzpsi[jgrid].im*psi_b_tmp.re - Lzpsi[jgrid].re*psi_b_tmp.im;
				
				// <b|L_x^2|a>
				lsqr_tmp.re += Lxsqrpsi[jgrid].re*psi_b_tmp.re + Lxsqrpsi[jgrid].im*psi_b_tmp.im;
				lsqr_tmp.im += Lxsqrpsi[jgrid].im*psi_b_tmp.re - Lxsqrpsi[jgrid].re*psi_b_tmp.im;
				// <b|L_y^2|a>
				lsqr_tmp.re += Lysqrpsi[jgrid].re*psi_b_tmp.re + Lysqrpsi[jgrid].im*psi_b_tmp.im;
				lsqr_tmp.im += Lysqrpsi[jgrid].im*psi_b_tmp.re - Lysqrpsi[jgrid].re*psi_b_tmp.im;
				// <b|L_z^2|a>
				lsqr_tmp.re += Lzsqrpsi[jgrid].re*psi_b_tmp.re + Lzsqrpsi[jgrid].im*psi_b_tmp.im;
				lsqr_tmp.im += Lzsqrpsi[jgrid].im*psi_b_tmp.re - Lzsqrpsi[jgrid].re*psi_b_tmp.im;
			}
			// normalize and store in allocated array
			
			l_mom[index].x_re = lx_tmp.re * grid->dv; l_mom[index].x_im = lx_tmp.im * grid->dv;
			l_mom[index].y_re = ly_tmp.re * grid->dv; l_mom[index].y_im = ly_tmp.im * grid->dv;
			l_mom[index].z_re = lz_tmp.re * grid->dv; l_mom[index].z_im = lz_tmp.im * grid->dv;
			l2_mom[index].re = lsqr_tmp.re * grid->dv; l2_mom[index].im = lsqr_tmp.im * grid->dv;
			
			//
			// Compute LdotS for the electron states 
			//
			zomplex psi_b_up_tmp, psi_b_dn_tmp;
			zomplex ldots_tmp;
			ldots_tmp.re = ldots_tmp.im = 0.00;

			for(jgrid = 0; jgrid < ist->ngrid; jgrid++){
				jgridup = jgrid;
				jgriddn = jgrid+ist->ngrid;
				psi_b_up_tmp.re = psi2[jgridup].re;
				psi_b_up_tmp.im = psi2[jgridup].im;
				psi_b_dn_tmp.re = psi2[jgriddn].re;
				psi_b_dn_tmp.im = psi2[jgriddn].im;

				//s_mom.xlx
				//<b|r,dn> * <r,up|Lxa>
				ldots_tmp.re += psi_b_dn_tmp.re * Lxpsi[jgridup].re
							  + psi_b_dn_tmp.im * Lxpsi[jgridup].im;
				ldots_tmp.im += psi_b_dn_tmp.re * Lxpsi[jgridup].im
							  - psi_b_dn_tmp.im * Lxpsi[jgridup].re;
				//<b|r,up> * <r,dn|Lxa>
				ldots_tmp.re += psi_b_up_tmp.re * Lxpsi[jgriddn].re
							  + psi_b_up_tmp.im * Lxpsi[jgriddn].im;
				ldots_tmp.im += psi_b_up_tmp.re * Lxpsi[jgriddn].im
							  - psi_b_up_tmp.im * Lxpsi[jgriddn].re;

				//s_mom.yly
				//i*<j|r,dn> * <r,up|Lyi>
				ldots_tmp.im += psi_b_dn_tmp.re * Lypsi[jgridup].re
				  			  + psi_b_dn_tmp.im * Lypsi[jgridup].im;
				ldots_tmp.re -= psi_b_dn_tmp.re * Lypsi[jgridup].im
				  			  - psi_b_dn_tmp.im * Lypsi[jgridup].re;
				//-i*<j|r,up> * <r,dn|Lyi>
				ldots_tmp.im -= psi_b_up_tmp.re * Lypsi[jgriddn].re
				  			  + psi_b_up_tmp.im * Lypsi[jgriddn].im;
				ldots_tmp.re += psi_b_up_tmp.re * Lypsi[jgriddn].im
				  			  - psi_b_up_tmp.im * Lypsi[jgriddn].re;

  			  	//s_mom.zlz
  			  	//<j|r,up> * <r,up|Lzi>
				ldots_tmp.re += psi_b_up_tmp.re * Lzpsi[jgridup].re
							  + psi_b_up_tmp.im * Lzpsi[jgridup].im;
				ldots_tmp.im += psi_b_up_tmp.re * Lzpsi[jgridup].im
							  - psi_b_up_tmp.im * Lzpsi[jgridup].re;
				//<j|r,dn> * <r,dn|Lzi>
				ldots_tmp.re -= psi_b_dn_tmp.re * Lzpsi[jgriddn].re
							  + psi_b_dn_tmp.im * Lzpsi[jgriddn].im;
				ldots_tmp.im -= psi_b_dn_tmp.re * Lzpsi[jgriddn].im
							  - psi_b_dn_tmp.im * Lzpsi[jgriddn].re;
			}
			// normalize 
			LdotS[index].re = 0.5 * ldots_tmp.re * grid->dv;
			LdotS[index].im = 0.5 * ldots_tmp.im * grid->dv;

			fprintf(pfx,"%ld\t%ld\t%lf\t%lf\n",a,b,l_mom[index].x_re,l_mom[index].x_im);
			fprintf(pfy,"%ld\t%ld\t%lf\t%lf\n",a,b,l_mom[index].y_re,l_mom[index].y_im);
			fprintf(pfz,"%ld\t%ld\t%lf\t%lf\n",a,b,l_mom[index].z_re,l_mom[index].z_im);
			fprintf(pfsqr,"%ld\t%ld\t%lf\t%lf\n",a,b,l2_mom[index].re,l2_mom[index].im);
			fprintf(pfls,"%ld\t%ld\t%lf\t%lf\n",a,b,LdotS[index].re,LdotS[index].im);
			if(a==b){
				printf("%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\n",a,b,l_mom[index].x_re,l_mom[index].y_re,l_mom[index].z_re, l2_mom[index].re, LdotS[index].re);
			}

			free(psi2);
		}
	}
	
	fclose(pfx); 
	fclose(pfy); 
	fclose(pfz); 
	fclose(pfsqr); 
	fclose(pfls);
	
	free(psi1);
	free(Lxpsi); free(Lypsi); free(Lzpsi);
	free(Lxsqrpsi); free(Lysqrpsi); free(Lzsqrpsi);
	free(temp1); free(temp2);

	fftw_free(fftwpsi);
	fftw_destroy_plan(planfw);
	fftw_destroy_plan(planbw);
}


/************************************************************/
//Calculate the three vector components of the action of the L operator on the 
//spatial part of the grid (no spin part)
void l_operator(zomplex* Lxpsi, zomplex* Lypsi, zomplex* Lzpsi, zomplex* psi_qp, 
	grid_st *grid, index_st *ist, par_st *par, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex* fftwpsi){
	
	double *kx, *ky, *kz, *kindex;
	double density, x, y, z;
	long jx, jy, jz, jyz, jxyz;
	zomplex px, py, pz;

	
	//generate kx,ky,kz vectors
	if ((kx = (double *) malloc(grid->nx * sizeof(kx[0]))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for kx in angular.c\n"); exit(EXIT_FAILURE);
    }
	if ((ky = (double *) malloc(grid->nx * sizeof(ky[0]))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for ky in angular.c\n"); exit(EXIT_FAILURE);
    }
	if ((kz = (double *) malloc(grid->nx * sizeof(kz[0]))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for kz in angular.c\n"); exit(EXIT_FAILURE);
    }
	if ((kindex = (double *) malloc(3*grid->ngrid * sizeof(kindex[0]))) == NULL){
        fprintf(stderr, "ERROR: allocating memory for kindex in angular.c\n"); exit(EXIT_FAILURE);
    }

  	/***initializing the k vectors ***/
  	init_k_vecs(kindex, kx, ky, kz, grid, ist, par);

	// Take derivative along x direction
	p_operator("X", kindex, (zomplex*)psi_qp, Lxpsi, grid, ist, par, planfw, planbw, fftwpsi);
	// Take derivative along y direction
	p_operator("Y", kindex, (zomplex*)psi_qp, Lypsi, grid, ist, par, planfw, planbw, fftwpsi);
	// Take derivative along x direction
	p_operator("Z", kindex, (zomplex*)psi_qp, Lzpsi, grid, ist, par, planfw, planbw, fftwpsi);

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

	free(kx); free(ky); free(kz); free(kindex);

	return;
}

/*************************************************************************************************/

void init_k_vecs(double *kindex, double *kx, double *ky, double *kz, grid_st *grid, index_st *ist, par_st *par){

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