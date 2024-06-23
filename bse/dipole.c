/****************************************************************************/

#include "fd.h"

/****************************************************************************/

void dipole(double *vx,double *vy,double *vz,zomplex *psi,zomplex *mux,zomplex *muy,zomplex *muz,double *eval,index_st ist,par_st par)
{
  FILE *pf1, *pf2, *pf3, *pf; 
  long a, i, jx, jy, jz, jgridup,jgriddn, jyz; 
  zomplex sumX, sumY, sumZ, tmp;
  double dz, dx, dy, ev, os;
  
  // Output will be written to these files
  pf = fopen("OS0.dat" , "w"); 
  pf1 = fopen("ux.dat", "w"); pf2 = fopen("uy.dat", "w"); pf3 = fopen("uz.dat", "w");

  // Make sure mux, muy and muz are zero to begin with
  for (i = 0; i < ist.total_lumo*ist.total_homo; i++) mux[i].re = muy[i].re = muz[i].re = mux[i].im = muy[i].im = muz[i].im = 0.0;

  // Main computional work of function performed here - must loop over all electron-hole (i-a) pairs
  // compute $$\vec{\mu} = <i|\vec{r}|a> = \int dr \sum_\simga \psi_{i}^{*}(r,\sigma) \vec{r} \psi_a*(r,\sigma)$$
  for (i = 0; i < ist.total_homo; i++, fprintf(pf1,"\n"),fprintf(pf2,"\n"),fprintf(pf3,"\n")){
    for (a = ist.lumo_idx; a < ist.lumo_idx+ist.total_lumo; a++) {
      sumX.re = sumY.re = sumZ.re = sumX.im = sumY.im = sumZ.im = 0.0;
      for (jz = 0; jz < ist.nz; jz++) {
        dz = vz[jz];
        for (jy = 0; jy < ist.ny; jy++) {
          dy = vy[jy];
          jyz = ist.nx * (ist.ny * jz + jy);
          for (jx = 0; jx < ist.nx; jx++) {
            dx = vx[jx];
            jgridup = jyz + jx;
            jgriddn = jgridup+ist.ngrid;
      	    tmp.re =  par.dv*(psi[i*ist.nspinngrid+jgridup].re * psi[a*ist.nspinngrid+jgridup].re + psi[i*ist.nspinngrid+jgridup].im * psi[a*ist.nspinngrid+jgridup].im
                             +psi[i*ist.nspinngrid+jgriddn].re * psi[a*ist.nspinngrid+jgriddn].re + psi[i*ist.nspinngrid+jgriddn].im * psi[a*ist.nspinngrid+jgriddn].im );
            tmp.im = par.dv * (-psi[i*ist.nspinngrid+jgridup].im * psi[a*ist.nspinngrid+jgridup].re + psi[i*ist.nspinngrid+jgridup].re * psi[a*ist.nspinngrid+jgridup].im
                              -psi[i*ist.nspinngrid+jgriddn].im * psi[a*ist.nspinngrid+jgriddn].re + psi[i*ist.nspinngrid+jgriddn].re * psi[a*ist.nspinngrid+jgriddn].im) ;
      	    

            sumX.re +=  dx*tmp.re; sumX.im +=  dx*tmp.im;
      	    sumY.re +=  dy*tmp.re; sumY.im +=  dy*tmp.im;
      	    sumZ.re +=  dz*tmp.re; sumZ.im +=  dz*tmp.im;
            //if(sqr(sumX.re)+sqr(sumX.im)>0){printf("a=%d i=%d",a,i);}
          }
        }
      }
      mux[i*ist.total_lumo+(a-ist.lumo_idx)].re = sumX.re;  mux[i*ist.total_lumo+(a-ist.lumo_idx)].im = sumX.im; 
      muy[i*ist.total_lumo+(a-ist.lumo_idx)].re = sumY.re;  muy[i*ist.total_lumo+(a-ist.lumo_idx)].im = sumY.im;
      muz[i*ist.total_lumo+(a-ist.lumo_idx)].re = sumZ.re;  muz[i*ist.total_lumo+(a-ist.lumo_idx)].im = sumZ.im;
      fprintf (pf1,"%ld %ld %g %g\n",i,a,mux[i*ist.total_lumo+(a-ist.lumo_idx)].re, mux[i*ist.total_lumo+(a-ist.lumo_idx)].im);
      fprintf (pf2,"%ld %ld %g %g\n",i,a,muy[i*ist.total_lumo+(a-ist.lumo_idx)].re, muy[i*ist.total_lumo+(a-ist.lumo_idx)].im);
      fprintf (pf3,"%ld %ld %g %g\n",i,a,muz[i*ist.total_lumo+(a-ist.lumo_idx)].re, muz[i*ist.total_lumo+(a-ist.lumo_idx)].im);

      os=((sqr(sumX.re)+sqr(sumX.im)) + (sqr(sumY.re)+sqr(sumY.im)) + (sqr(sumZ.re)+sqr(sumZ.im)));
      ev = eval[a] - eval[i];
      fprintf(pf,"%ld %ld %.8f %.12f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", i, a, sqrt(os), ev, (2.0/3.0)*ev*os,
	       mux[i*ist.total_lumo+(a-ist.lumo_idx)].re, mux[i*ist.total_lumo+(a-ist.lumo_idx)].im,
	       muy[i*ist.total_lumo+(a-ist.lumo_idx)].re, muy[i*ist.total_lumo+(a-ist.lumo_idx)].im,
	       muz[i*ist.total_lumo+(a-ist.lumo_idx)].re, muz[i*ist.total_lumo+(a-ist.lumo_idx)].im);
    }
  }
  fclose(pf); fclose(pf1); fclose(pf2); fclose(pf3);

  return;
}

/****************************************************************************/
// This function calculates the magnetic dipole matrix elements between the
// single-particle electron (a) and hole (i) states: <psi_a|m|psi_i>
// where m = -1/2*L = -1/2*r x p where x is the cross product.


void mag_dipole(double *vx, double *vy, double *vz, double *psi, double *mx, double *my, double *mz, 
  double *eval, fftw_plan_loc *planfw,fftw_plan_loc *planbw,fftw_complex *fftwpsi, index_st ist, par_st par)
{
  FILE *pf1, *pf2, *pf3, *pf; 
  long a, i, jx, jy, jz, jgrid, jyz; 
  double sumX, sumY, sumZ;
  double z, x, y, ev, ms, tmp;
  double *kx, *ky, *kz;
  double *kindex;
  double *psidx, *psidy, *psidz;
  zomplex *cpsi;
  
  // Output will be written to these files
  pf = fopen("M0.dat" , "w"); 
  pf1 = fopen("mx.dat", "w"); 
  pf2 = fopen("my.dat", "w"); 
  pf3 = fopen("mz.dat", "w");

  // Make sure mx, my and mz are zero to begin with
  for (i = 0; i < ist.total_lumo*ist.total_homo; i++) mx[i] = my[i] = mz[i] = 0.0;

  // Allocate memory 
  if ((kx = (double *) calloc(ist.nx, sizeof(double)))==NULL) nerror("kx");
  if ((ky = (double *) calloc(ist.ny, sizeof(double)))==NULL) nerror("ky");
  if ((kz = (double *) calloc(ist.nz, sizeof(double)))==NULL) nerror("kz");
  if ((kindex = (double *) calloc(3*ist.ngrid, sizeof(double)))==NULL) nerror("kindex");
  if ((psidx = (double *) calloc(ist.ngrid, sizeof(double)))==NULL) nerror("psidx");
  if ((psidy = (double *) calloc(ist.ngrid, sizeof(double)))==NULL) nerror("psidy");
  if ((psidz = (double *) calloc(ist.ngrid, sizeof(double)))==NULL) nerror("psidz");

  // Calculate kx 
  for (jx = 1; jx <= ist.nx / 2; jx++) {
    kx[jx] = (double)(jx) * par.dkx * ist.nx_1 * ist.ny_1 * ist.nz_1;
    kx[ist.nx-jx] = -kx[jx];
  }
  // Calculate ky 
  for (jy = 1; jy <= ist.ny / 2; jy++) {
    ky[jy] = (double)(jy) * par.dky * ist.nx_1 * ist.ny_1 * ist.nz_1;
    ky[ist.ny-jy] = -ky[jy];
  }
  // Calculate kz 
  for (jz = 1; jz <= ist.nz / 2; jz++) {
    kz[jz] = (double)(jz) * par.dkz * ist.nx_1 * ist.ny_1 * ist.nz_1;
    kz[ist.nz-jz] = -kz[jz];
  }  
  for (jz = 0; jz < ist.nz; jz++) {
    for (jy = 0; jy < ist.ny; jy++) {
      jyz = ist.nx * (ist.ny * jz + jy);
      for (jx = 0; jx < ist.nx; jx++) {   
        jgrid = jyz + jx;
        kindex[3*jgrid]   = kx[jx]; 
        kindex[3*jgrid+1] = ky[jy];
        kindex[3*jgrid+2] = kz[jz];
      } 
    }
  }

  // Allocate memory for Fourier transforms
  if ((cpsi = (zomplex *) calloc(ist.ngrid, sizeof(zomplex)))==NULL) nerror("cpsi");

  // Main computional work of function performed here - must loop over all electron-hole (i-a) pairs
  for (i = 0; i < ist.total_homo; i++, fprintf(pf1,"\n"),fprintf(pf2,"\n"),fprintf(pf3,"\n")) {
    // Fourier transform the hole wavefunctions
    for (jgrid = 0; jgrid < ist.ngrid; jgrid++) {
      cpsi[jgrid].re = psi[i*ist.ngrid+jgrid]; 
      cpsi[jgrid].im = 0.0;
    }
    memcpy(&fftwpsi[0], &cpsi[0], ist.ngrid*sizeof(fftwpsi[0]));
    fftw_execute(planfw[0]);
    memcpy(&cpsi[0], &fftwpsi[0], ist.ngrid*sizeof(cpsi[0])); // store the FT of the hole wavefunction in cpsi
    for(jgrid = 0; jgrid < ist.ngrid; jgrid++) {
      fftwpsi[jgrid][0] = cpsi[jgrid].im * kindex[3*jgrid];
      fftwpsi[jgrid][1] = -cpsi[jgrid].re * kindex[3*jgrid];
    }
    fftw_execute(planbw[0]);
    for(jgrid = 0; jgrid < ist.ngrid; jgrid++) {
      psidx[jgrid] = fftwpsi[jgrid][0]; // x-derivative of the hole wavefunction
    }

    for(jgrid = 0; jgrid < ist.ngrid; jgrid++) {
      fftwpsi[jgrid][0] = cpsi[jgrid].im * kindex[3*jgrid+1];
      fftwpsi[jgrid][1] = -cpsi[jgrid].re * kindex[3*jgrid+1];
    }
    fftw_execute(planbw[0]);
    for(jgrid = 0; jgrid < ist.ngrid; jgrid++) {
      psidy[jgrid] = fftwpsi[jgrid][0]; // y-derivative of the hole wavefunction
    }

    for(jgrid = 0; jgrid < ist.ngrid; jgrid++) {
      fftwpsi[jgrid][0] = cpsi[jgrid].im * kindex[3*jgrid+2];
      fftwpsi[jgrid][1] = -cpsi[jgrid].re * kindex[3*jgrid+2];
    }
    fftw_execute(planbw[0]);
    for(jgrid = 0; jgrid < ist.ngrid; jgrid++) {
      psidz[jgrid] = fftwpsi[jgrid][0]; // z-derivative of the hole wavefunction
    }

    for (a = ist.lumo_idx; a < ist.lumo_idx+ist.total_lumo; a++) {
      sumX = sumY = sumZ = 0.0;
      for (jz = 0; jz < ist.nz; jz++) {
        z = vz[jz];
        for (jy = 0; jy < ist.ny; jy++) {
          y = vy[jy];
          jyz = ist.nx * (ist.ny * jz + jy);
          for (jx = 0; jx < ist.nx; jx++) {
            x = vx[jx];
            jgrid = jyz + jx;
            tmp = psi[a*ist.ngrid+jgrid] * par.dv * 0.5;
            sumX += tmp * ( y * psidz[jgrid] - z * psidy[jgrid] );
            sumY += tmp * ( z * psidx[jgrid] - x * psidz[jgrid] );
            sumZ += tmp * ( x * psidy[jgrid] - y * psidx[jgrid] );
          }
        }
      }
      mx[i*ist.total_lumo+(a-ist.lumo_idx)] =  sumX;
      my[i*ist.total_lumo+(a-ist.lumo_idx)] =  sumY;
      mz[i*ist.total_lumo+(a-ist.lumo_idx)] =  sumZ;
      fprintf(pf1,"%ld %ld %g\n",i,a,sumX);
      fprintf(pf2,"%ld %ld %g\n",i,a,sumY);
      fprintf(pf3,"%ld %ld %g\n",i,a,sumZ);

      ms=(sqr(mx[i*ist.total_lumo+(a-ist.lumo_idx)]) + sqr(my[i*ist.total_lumo+(a-ist.lumo_idx)]) + sqr(mz[i*ist.total_lumo+(a-ist.lumo_idx)]));
      ev = eval[a] - eval[i];
      fprintf(pf,"%ld %ld %.8f %.12f %.8f %.8f %.8f %.8f\n", i, a, sqrt(ms), ev, (4.0/3.0)*ev*ms,
         sqr(mx[i*ist.total_lumo+(a-ist.lumo_idx)]),
         sqr(my[i*ist.total_lumo+(a-ist.lumo_idx)]),
         sqr(mz[i*ist.total_lumo+(a-ist.lumo_idx)]));
    }
  }
  fclose(pf); fclose(pf1); fclose(pf2); fclose(pf3);

  // Free dynamically allocated memory
  free(kx); free(ky); free(kz);
  free(kindex);
  free(cpsi);
  free(psidx); free(psidy); free(psidz);

  return;
}

/****************************************************************************/

void rotational_strength(double *rs, double *mux, double *muy, double *muz, double *mx, 
  double *my, double *mz, double *eval, index_st ist) {
  FILE *pf;
  int i, a, index;

  pf = fopen("rs0.dat", "w");
  for (i = 0; i < ist.total_homo; i++) {
    for (a = ist.lumo_idx; a < ist.lumo_idx+ist.total_lumo; a++) {
      index = i*ist.total_lumo + (a-ist.lumo_idx);
      rs[index] = mux[index]*mx[index] + muy[index]*my[index] + muz[index]*mz[index];
      fprintf(pf, "%d %d %.12f  %.16f\n", i, a, eval[a] - eval[i], rs[index]);
    }
  }
  fclose(pf);

  return;
}

/****************************************************************************/
