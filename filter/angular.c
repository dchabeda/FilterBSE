#include "fd.h"

/************************************************************/
void calc_angular_exp(double *psitot, grid_st *grid, int start, int stop, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel,
 	fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi){
	
	int i, ipsi; 
	long jx, jy, jz, jyz, jxyz, jxyz_real, jxyz_imag;
	double Jxexp, Jxsqr, Jyexp, Jysqr,Jzexp, Jzsqr, Jsqrexp, norm;
	zomplex *psi, *Jxpsi, *Jypsi, *Jzpsi;
	FILE* pf = fopen("angular.dat", "w");

	fprintf(pf, "###\tJx\t\tJx^2\t\tJy\t\tJy^2\t\tJz\t\tJz^2\t\tJ^2\t\tnorm\n");
	
	//allocate memoRy for the operated wavefucntions
	if ((Jxpsi  = (zomplex*) calloc(ist->nspinngrid,sizeof(zomplex)))==NULL) nerror("Jxpsi");
  	if ((Jypsi  = (zomplex*) calloc(ist->nspinngrid,sizeof(zomplex)))==NULL) nerror("Jypsi");
  	if ((Jzpsi  = (zomplex*) calloc(ist->nspinngrid,sizeof(zomplex)))==NULL) nerror("Jzpsi");
  	if ((psi = (zomplex*) calloc(ist->nspinngrid,sizeof(zomplex)))==NULL) nerror("psi");

	for(ipsi = start; ipsi < stop; ipsi++){	  	
	  	//psi = &psitot[ipsi*ist->nspinngrid];
	  	for (jz = 0; jz < grid->nz; jz++) { 
	  		for (jy = 0;jy<grid->ny;jy++){
	  			jyz = grid->nx * (grid->ny * jz + jy);
	  			for(jx = 0; jx<grid->nx;jx++){
					jxyz = jyz + jx;
					jxyz_real = ist->complex_idx * jxyz;
            		jxyz_imag = ist->complex_idx * jxyz + 1;

					psi[jxyz].re = psitot[ist->complex_idx*ipsi*ist->nspinngrid + jxyz_real];
					if (1 == flag->isComplex){
						psi[jxyz].im = psitot[ipsi*ist->nspinngrid+jxyz_imag];
						// Uncomment these lines if, for some reason, people start
						// using complex wavefunctions without spinors. I don't know why 
						// someone would use complex numbers without spinors because if there is
						// time reversal symmetry, the wavefunctions are guaranteed real,
						// which is at least 2x less expensive. But anyway, uncommenting these 
						// lines is how you would correctly handle that case.
						//if (1 == flag->useSpinors){
						psi[jxyz+ist->ngrid].re = psitot[ist->complex_idx*(ipsi*ist->nspinngrid + ist->ngrid) + jxyz_real];
						psi[jxyz+ist->ngrid].im = psitot[ist->complex_idx*(ipsi*ist->nspinngrid + ist->ngrid) + jxyz_imag];
						//}
					}
					
				}
			}		
		}

		
		low_pass_filter(&psi[0],grid,planfw,planbw,fftwpsi,ist,par,parallel);
		if (1 == flag->useSpinors){
			low_pass_filter(&psi[ist->ngrid],grid,planfw,planbw,fftwpsi,ist,par,parallel);
		}

		normalize(psi, ist->nspinngrid, ist, par, flag, parallel);

		char filename[20];
		double *rho = calloc(ist->ngrid, sizeof(double));
		
		for (long jgrid = 0; jgrid < ist->ngrid; jgrid++){
    		rho[jgrid] = (psi[jgrid].re);
		 }
		sprintf(filename, "smpsi%iUpRe.cube", ipsi);
	  	write_cube_file(rho, grid, filename);

		if (1 == flag->isComplex){
			for (long jgrid = 0; jgrid < ist->ngrid; jgrid++){
				rho[jgrid] = (psi[jgrid].im);
			}
			sprintf(filename, "smpsi%iUpIm.cube", ipsi);
			write_cube_file(rho, grid, filename);

			//if (1 == flag->useSpinors){
			for (long jgrid = 0; jgrid < ist->ngrid; jgrid++){
				rho[jgrid] = (psi[ist->ngrid+jgrid].re);
			}
			sprintf(filename, "smpsi%iDnRe.cube", ipsi);
			write_cube_file(rho, grid, filename);

			for (long jgrid = 0; jgrid < ist->ngrid; jgrid++){
				rho[jgrid] = (psi[ist->ngrid+jgrid].im);
			}
			sprintf(filename, "smpsi%iDnIm.cube", ipsi);
			write_cube_file(rho, grid, filename);
			//}
		}

		//compute the action of the J operator
	  	apply_J_op(Jxpsi, Jypsi, Jzpsi, psi, grid, planfw, planbw, fftwpsi, ist, par, parallel);

	  	//compute Ji expectation as <psi|Ji_psi> and J^2 as \sum_i <Ji_psi|Ji_psi>
	  	Jxexp = Jxsqr = Jyexp = Jysqr = Jzexp = Jzsqr = Jsqrexp = norm = 0.00;
	  	omp_set_dynamic(0);
	  	omp_set_num_threads(parallel->nthreads);
		#pragma omp parallel for reduction (+:Jxexp, Jxsqr, Jyexp, Jysqr,Jzexp, Jzsqr, Jsqrexp, norm)
		for (i = 0; i < ist->nspinngrid; i++){
			//<psi|Jxpsi>
			Jxexp += (psi[i].re * Jxpsi[i].re) + (psi[i].im * Jxpsi[i].im);
			Jxsqr += (Jxpsi[i].re * Jxpsi[i].re) + (Jxpsi[i].im * Jxpsi[i].im);

			//<psi|Jypsi>
			Jyexp += (psi[i].re * Jypsi[i].re) + (psi[i].im * Jypsi[i].im);
			Jysqr += (Jypsi[i].re * Jypsi[i].re) + (Jypsi[i].im * Jypsi[i].im);

			//<psi|Jzpsi>
			Jzexp += (psi[i].re * Jzpsi[i].re) + (psi[i].im * Jzpsi[i].im);	
			Jzsqr += (Jzpsi[i].re * Jzpsi[i].re) + (Jzpsi[i].im * Jzpsi[i].im);	

			Jsqrexp += (Jxpsi[i].re*Jxpsi[i].re) + (Jxpsi[i].im*Jxpsi[i].im) 
					+ (Jypsi[i].re*Jypsi[i].re) + (Jypsi[i].im*Jypsi[i].im)
					+ (Jzpsi[i].re*Jzpsi[i].re) + (Jzpsi[i].im*Jzpsi[i].im);
		
			norm += (psi[i].re*psi[i].re)+(psi[i].im*psi[i].im);
		}
		Jxexp*=par->dv;
		Jxsqr*=par->dv;
		Jyexp*=par->dv;
		Jysqr*=par->dv;
		Jzexp*=par->dv;
		Jzsqr*=par->dv;
		Jsqrexp*=par->dv;
		norm*=par->dv;


		fprintf(pf, "%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
					ipsi,Jxexp, Jxsqr, Jyexp, Jysqr, Jzexp, Jzsqr, Jsqrexp, norm);
	}

	fclose(pf);
	free(Jxpsi); free(Jypsi); free(Jzpsi);


}

/************************************************************/
//Calculate the three vector components of the action of the J=L+S operator
void apply_J_op( zomplex* Jxpsi, zomplex* Jypsi, zomplex* Jzpsi,zomplex* psi, 
	grid_st* grid, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi,index_st *ist, par_st *par, parallel_st *parallel){
	long jz,jy,jyz,jx,jxyz;
	double density,x,y,z;

	/*** Find electron center of mass ***/
	x = y = z = 0;
	omp_set_dynamic(0);
  	omp_set_num_threads(parallel->nthreads);
  	#pragma omp parallel for private (jz,jy,jyz,jx,jxyz) reduction(+:x,y,z)
  	for (jz = 0; jz < grid->nz; jz++) { 
  		for (jy = 0;jy<grid->ny;jy++){
  			jyz = grid->nx * (grid->ny * jz + jy);
  			for(jx = 0; jx<grid->nx;jx++){
				jxyz = jyz + jx;
				density = psi[jxyz].re *psi[jxyz].re + psi[jxyz].im * psi[jxyz].im 
				+ psi[jxyz+ist->ngrid].re *psi[jxyz+ist->ngrid].re + psi[jxyz+ist->ngrid].im * psi[jxyz+ist->ngrid].im;
				
				x += grid->x[jx] * density;
				y += grid->y[jy] * density;
				z += grid->z[jz] * density;
			}
		}
	}
	/*
	center.x = x*par->dv; center.y=y*par->dv;center.z=z*par->dv;
	printf("Center: <%f, %f, %f>\n", center.x, center.y , center.z);
	*/
	
	//spin up part
	apply_L_op(&Jxpsi[0],&Jypsi[0],&Jzpsi[0], &psi[0], grid,planfw,planbw, fftwpsi,ist, par, parallel);
	

	//spin dn part
	if (2 == ist->complex_idx){
		apply_L_op(&Jxpsi[ist->ngrid],&Jypsi[ist->ngrid],&Jzpsi[ist->ngrid], &psi[ist->ngrid], 
			grid,planfw,planbw, fftwpsi,ist, par, parallel);
	}

	//add the spin part which acts locally
	omp_set_dynamic(0);
  	omp_set_num_threads(parallel->nthreads);
	#pragma omp parallel for private (jz,jy,jyz,jx,jxyz)
  	for (jz = 0; jz < grid->nz; jz++) { 
  		for (jy = 0;jy<grid->ny;jy++){
  			jyz = grid->nx * (grid->ny * jz + jy);
  			for(jx = 0; jx<grid->nx;jx++){
				jxyz = jyz + jx;
				//Jx|r,up> = Lx|r,up> + 1/2 |r, dn>
				Jxpsi[jxyz].re += 0.5*psi[jxyz+ist->ngrid].re;
				Jxpsi[jxyz].im += 0.5*psi[jxyz+ist->ngrid].im;
				//Jx|r,dn> = Lx|r,dn> + 1/2 |r, up>
				Jxpsi[jxyz+ist->ngrid].re += 0.5*psi[jxyz].re;
				Jxpsi[jxyz+ist->ngrid].im += 0.5*psi[jxyz].im;

				//Jy|r,up> = Ly|r,up> - i/2 |r, dn>
				Jypsi[jxyz].re += 0.5* psi[jxyz+ist->ngrid].im;
				Jypsi[jxyz].im -= 0.5* psi[jxyz+ist->ngrid].re;
				//Jy|r,dn> = Ly|r,dn> + i/2 |r, up>
				Jypsi[jxyz+ist->ngrid].re -= 0.5* psi[jxyz].im;
				Jypsi[jxyz+ist->ngrid].im += 0.5* psi[jxyz].re;

				//Jz|r,up> = Lz|r,up> + 1/2 |r,up>
				Jzpsi[jxyz].re += 0.5*psi[jxyz].re;
				Jzpsi[jxyz].im += 0.5*psi[jxyz].im;
				//Jz|r,dn> = Lz|r,dn> - 1/2 |r,dn>
				Jzpsi[jxyz+ist->ngrid].re -= 0.5*psi[jxyz+ist->ngrid].re;
				Jzpsi[jxyz+ist->ngrid].im -= 0.5*psi[jxyz+ist->ngrid].im;

			}
		}
	}

}


/************************************************************/
//Calculate the three vector components of the action of the L operator on the 
//spatial part of the grid (no spin part)
void apply_L_op(zomplex* Lxpsi, zomplex* Lypsi, zomplex* Lzpsi, zomplex* psi, 
	grid_st* grid, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi,index_st *ist, par_st *par, parallel_st *parallel){
	
	double *kx, *ky, *kz, x, y, z;
	long jx,jy,jz, jyz, jxyz;

	xyz_st center;
	center.x = center.y = center.z = 0.00;


	/*** First use the fft to get the action of the p operator on each axis (-i d/dx) imag part from the fft definition***/


	//generate kx,ky,kz vectors
	if ((kx  = (double*)calloc(grid->nx,sizeof(double)))==NULL)nerror("kx");
  	if ((ky  = (double*)calloc(grid->ny,sizeof(double)))==NULL)nerror("ky");
  	if ((kz  = (double*)calloc(grid->nz,sizeof(double)))==NULL)nerror("kz");

  	/***initializing the k vectors ***/
  	//negative frequencies count backwards from end
  	//overall transform fw and bw aquires factor of nx*ny*nz so normalize here
  	for (kx[0] = 0.0, jx = 1; jx <= grid->nx / 2; jx++){
    	kx[grid->nx-jx] = -1.00 * (kx[jx] = (double)(jx) * grid->dkx * 
    		grid->nx_1 * grid->ny_1 * grid->nz_1);
  	}

  	for (ky[0] = 0.0, jy = 1; jy <= grid->ny / 2; jy++){
    	ky[grid->ny-jy] = -1.00 * (ky[jy] = (double)(jy) * grid->dky *
			grid->nx_1 * grid->ny_1 * grid->nz_1);
  	}

  	for (kz[0] = 0.0, jz = 1; jz <= grid->nz / 2; jz++){
    	kz[grid->nz-jz] = -1.00 * (kz[jz] = (double)(jz) * grid->dkz *
			grid->nx_1 * grid->ny_1 * grid->nz_1);
  	}


	// Copy psi to fftwpsi
  	memcpy(&fftwpsi[0], &psi[0], ist->ngrid*sizeof(fftwpsi[0]));
  	// FT from r-space to k-space
  	fftw_execute(planfw);
  	
  	omp_set_dynamic(0);
  	omp_set_num_threads(parallel->nthreads);
  	#pragma omp parallel for private (jz,jy,jyz,jx,jxyz)
  	for (jz = 0; jz < grid->nz; jz++) { 
  		for (jy = 0;jy<grid->ny;jy++){
  			jyz = grid->nx * (grid->ny * jz + jy);
  			for(jx = 0; jx<grid->nx;jx++){
				jxyz = jyz + jx;
				//multiply by kx to get partial along x-axis
				fftwpsi[jxyz][0] *= kx[jx];
				fftwpsi[jxyz][1] *= kx[jx];
			}
		}
	}

	// Inverse FT back to r-space
	fftw_execute(planbw);
	// Copy fftwpsi to psi to store Lx|psi> into |Lxpsi>
	memcpy(&Lxpsi[0], &fftwpsi[0], ist->ngrid*sizeof(Lxpsi[0]));
	

	// Copy psi to fftwpsi
  	memcpy(&fftwpsi[0], &psi[0], ist->ngrid*sizeof(fftwpsi[0]));

  	// FT from r-space to k-space
  	fftw_execute(planfw);
  	
  	omp_set_dynamic(0);
  	omp_set_num_threads(parallel->nthreads);
 	#pragma omp parallel for private (jz,jy,jyz,jx,jxyz)
  	for (jz = 0; jz < grid->nz; jz++) { 
  		for (jy = 0;jy<grid->ny;jy++){
  			jyz = grid->nx * (grid->ny * jz + jy);
  			for(jx = 0; jx<grid->nx;jx++){
				jxyz = jyz + jx;
				//multiply by ky to get partial along y-axis
				fftwpsi[jxyz][0] *= ky[jy];
				fftwpsi[jxyz][1] *= ky[jy];
			}
		}
	}

	// Inverse FT back to r-space
	fftw_execute(planbw);
  
	// Copy fftwpsi to psi to store Ly|psi> into |Lypsi>
	memcpy(&Lypsi[0], &fftwpsi[0], ist->ngrid*sizeof(Lypsi[0]));
	
	

	// Copy psi to fftwpsi
  	memcpy(&fftwpsi[0], &psi[0], ist->ngrid*sizeof(fftwpsi[0]));

  	// FT from r-space to k-space
  	fftw_execute(planfw);
  	omp_set_dynamic(0);
  	omp_set_num_threads(parallel->nthreads);
  	#pragma omp parallel for private (jz,jy,jyz,jx,jxyz)
  	for (jz = 0; jz < grid->nz; jz++) { 
  		for (jy = 0;jy<grid->ny;jy++){
  			jyz = grid->nx * (grid->ny * jz + jy);
  			for(jx = 0; jx<grid->nx;jx++){
				jxyz = jyz + jx;
				//multiply by kz to get partial along z-axis
				fftwpsi[jxyz][0] *= kz[jz];
				fftwpsi[jxyz][1] *= kz[jz];
			}
		}
	}

	// Inverse FT back to r-space
	fftw_execute(planbw);
  
	// Copy fftwpsi to psi to store Lz|psi> into |Lzpsi>
	memcpy(&Lzpsi[0], &fftwpsi[0], ist->ngrid*sizeof(Lzpsi[0]));

	free(kx); free(ky); free(kz);
	


	/*** Now do the cross product part at each grid point***/
	zomplex gradx, grady, gradz;
	omp_set_dynamic(0);
  	omp_set_num_threads(parallel->nthreads);
	#pragma omp parallel for private (jz,jy,jyz,jx,jxyz,gradx,grady,gradz,x,y,z)
  	for (jz = 0; jz < grid->nz; jz++) { 
  		z = grid->z[jz]-center.z;
  		for (jy = 0;jy<grid->ny;jy++){
  			jyz = grid->nx * (grid->ny * jz + jy);
  			y = grid->y[jy]- center.y;
  			for(jx = 0; jx<grid->nx;jx++){
  				x = grid->x[jx]-center.x;
				jxyz = jyz + jx;
				//copy over tmp varibales
				gradx.re = Lxpsi[jxyz].re; gradx.im = Lxpsi[jxyz].im;
				grady.re = Lypsi[jxyz].re; grady.im = Lypsi[jxyz].im;
				gradz.re = Lzpsi[jxyz].re; gradz.im = Lzpsi[jxyz].im;
				
				Lxpsi[jxyz].re = (y*gradz.re - z*grady.re);

				Lxpsi[jxyz].im = (y*gradz.im - z*grady.im);

				Lypsi[jxyz].re = (z*gradx.re - x*gradz.re);
		
				Lypsi[jxyz].im = (z*gradx.im - x*gradz.im);

				Lzpsi[jxyz].re = (x*grady.re - y*gradx.re);

				Lzpsi[jxyz].im = (x*grady.im - y*gradx.im);


			}
		}
	}


}
/************************************************************/
void low_pass_filter(zomplex* psi, grid_st *grid, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi,
	
	index_st *ist, par_st *par, parallel_st *parallel){
	double *kx,*ky,*kz;
	long jx,jy,jz,jyz,jxyz;
	//generate kx,ky,kz vectors
	if ((kx  = (double*)calloc(grid->nx,sizeof(double)))==NULL)nerror("kx");
  	if ((ky  = (double*)calloc(grid->ny,sizeof(double)))==NULL)nerror("ky");
  	if ((kz  = (double*)calloc(grid->nz,sizeof(double)))==NULL)nerror("kz");

  	/***initializing the k vectors ***/
  	//negative frequencies count backwards from end
  	//overall transform fw and bw aquires factor of nx*ny*nz so normalize here
  	for (kx[0] = 0.0, jx = 1; jx <= grid->nx / 2; jx++){
    	kx[grid->nx-jx] = -1.00 * (kx[jx] = (double)(jx) * grid->dkx );
  	}

  	for (ky[0] = 0.0, jy = 1; jy <= grid->ny / 2; jy++){
    	ky[grid->ny-jy] = -1.00 * (ky[jy] = (double)(jy) * grid->dky );
  	}

  	for (kz[0] = 0.0, jz = 1; jz <= grid->nz / 2; jz++){
    	kz[grid->nz-jz] = -1.00 * (kz[jz] = (double)(jz) * grid->dkz );
  	}


	// Copy psi to fftwpsi
  	memcpy(&fftwpsi[0], &psi[0], ist->ngrid*sizeof(fftwpsi[0]));
  	// FT from r-space to k-space
  	fftw_execute(planfw);
  	
  	omp_set_dynamic(0);
  	omp_set_num_threads(parallel->nthreads);
  	#pragma omp parallel for private (jz,jy,jyz,jx,jxyz)
  	for (jz = 0; jz < grid->nz; jz++) { 
  		for (jy = 0;jy<grid->ny;jy++){
  			jyz = grid->nx * (grid->ny * jz + jy);
  			for(jx = 0; jx<grid->nx;jx++){
				jxyz = jyz + jx;
				
				if(sqr(kx[jx])+sqr(ky[jy])+sqr(kz[jz]) < sqr(3.0*grid->dkz)){
					fftwpsi[jxyz][0] *=1.00;
					fftwpsi[jxyz][1] *=1.00;
				}
				else{
					fftwpsi[jxyz][0] =0;
					fftwpsi[jxyz][1] =0;
				}
			}
		}
	}

	// Inverse FT back to r-space
	fftw_execute(planbw);
	// Copy fftwpsi to psi to store Lx|psi> into |Lxpsi>
	memcpy(&psi[0], &fftwpsi[0], ist->ngrid*sizeof(psi[0]));
	
}

/************************************************************/
