#include "fd.h"

/*****************************************************************************/

#define filt_func(x,dt)   (sqrt(dt / PIE) * exp(-sqr(x) * dt))
//#define filt_func(x,dt)   (0.04 / (PIE * (sqr(x) + sqr(0.04))))

void gen_newton_coeff(zomplex *an, double *samp, double *ene_targets, index_st *ist, par_st *par, parallel_st *parallel){
  /*******************************************************************
  * This function generates coefficients for a Newton interpolation  *
  * poolynomial approximation of Gaussian filter functions. A total  *
  * of ie series are generated, one for each of the energy targets   *
  * inputs:                                                          *
  *  [an] arr that will hold all ncheby*m_states_per_filter coeffs   *
  *  [samp] ncheby-long arr (zn) pf the Chebyshev support points     *
  *  [ist] ptr to counters, indices, and lengths                     *
  *  [par] ptr to par_st holding VBmin, VBmax... params              *
  *  [parallel] holds options for parallelization                    *
  * outputs: void                                                    *
  ********************************************************************/

  FILE *pf;
  zomplex *samploc;
  double scale, res = 1.0, Smin = -2.0, Smax = 2.0, x, rho, sumre, sumim;
  long i, j, ie;
  char str[20];
  
  scale = (Smax - Smin) / par->dE;

  //chebyshev_reordered(samp,Smin,Smax,ist->ncheby);
  samploc = (zomplex*)calloc(ist->ncheby, sizeof(zomplex));
  rho = samp_points_ashkenazy(samploc, Smin, Smax, ist->ncheby);
  for (j = 0; j < ist->ncheby; j++) {
    samp[j] = samploc[j].re;
  }

  pf = fopen("zn.dat", "w");
  for (j = 0; j < ist->ncheby; j++) {fprintf(pf, "%ld %lg\n", j, samp[j]);}
  fclose(pf);

  
  omp_set_dynamic(0);
  omp_set_num_threads(parallel->nthreads);
#pragma omp parallel for private(ie,x,j,i,res,sumre,sumim)
  for (ie = 0; ie < ist->m_states_per_filter; ie++){
    // open file to print Newton interpolation coefficients for the ie-th energy target
    sprintf(str, "coeff-%ld.dat", ie);
    pf = fopen(str , "w");
    // Before recursion, need to define term 0
    x = (samp[0] + 2.0) / scale + par->Vmin;
    an[ist->ncheby*ie+0].re = filt_func(x - ene_targets[ie], par->dt);
    an[ist->ncheby*ie+0].im = 0.0;
    fprintf (pf,"%d %g %g\n", 0, an[ist->ncheby*ie+0].re, an[ist->ncheby*ie+0].im);
    
    // term 1
    x = (samp[1] + 2.0) / scale + par->Vmin;
    an[ist->ncheby*ie+1].re = (filt_func(x - ene_targets[ie], par->dt) - an[ist->ncheby*ie+0].re) / (samp[1] - samp[0]);
    an[ist->ncheby*ie+1].im = (- an[ist->ncheby*ie+0].im) / (samp[1] - samp[0]);
    fprintf (pf,"%d %g %g\n", 1, an[ist->ncheby*ie+1].re, an[ist->ncheby*ie+1].im);
    
    // recursively define the remaining coefficients
    for (j = 2; j < ist->ncheby; j++) {
      for (i = 1, res = 1.0, sumre = sumim = 0.0; i < j; i++) {
        res *= (samp[j] - samp[i-1]);
        sumre += an[ist->ncheby*ie+i].re * res;
        sumim += an[ist->ncheby*ie+i].im * res;
      }
      res *= (samp[j] - samp[j-1]);
      x = (samp[j] + 2.0) / scale + par->Vmin;
      an[ist->ncheby*ie+j].re = (filt_func(x - ene_targets[ie], par->dt) - an[ist->ncheby*ie+0].re - sumre) / res;
      an[ist->ncheby*ie+j].im = (-an[ist->ncheby*ie+0].im - sumim) / res;
      fprintf (pf,"%ld %g %g\n", j, an[ist->ncheby*ie+j].re, an[ist->ncheby*ie+j].im);
    }

    fclose(pf);
  }
  
  check_function(an, samploc, ist, par, ene_targets[0]);
  free(samploc);

  return;
}

/**************************************************************************/

void chebyshev_reordered(double *point,double min,double max,long ncheby)
{
  long j,k,m;

  point[0] = 1;
  for (k = 1, m = ncheby/2; k < ncheby; k *= 2, m /= 2)
    for (j = 0; j < ncheby - k; j++)
      point[k+j] = point[j] + m;
  for (j = 0; j < ncheby; j++)
    point[j]= (max+min)/2.0 + (max-min)/2.0*cos(PIE*(point[j]-0.5)/(ncheby));
  return;
}

/**************************************************************************/

double samp_points_ashkenazy(zomplex *point,double min,double max,long ncheby)
{
  long j,k, nc3=32*ncheby, jrnd, kmax, imfrac;
  double dsre, dsim, fkmax, fk, minim, maxim, range = 0.0, *veca, del;
  zomplex *samp;

  samp = (zomplex*)calloc(nc3,sizeof(zomplex));  
  veca = (double*)calloc(nc3,sizeof(double));  

  minim = range * min;
  maxim = range * max;
  imfrac = nc3 / 8;
  dsre = (max-min)/(double)(nc3-1-imfrac);
  dsim = (maxim-minim)/(double)(imfrac-1);
  for (j = 0; j < nc3; j++) samp[j].re = samp[j].im = 0.0;
  for (j = 0; j < nc3-imfrac; j++) samp[j].re = min + (double)(j) * dsre;
  for (j = nc3-imfrac; j < nc3; j++) samp[j].im = minim + (double)(j-nc3+imfrac) * dsim;

  jrnd = 0;
  point[0] = samp[jrnd];

  for (k = 0; k < nc3; k++){
    del = sqr(samp[k].re-point[0].re) + sqr(samp[k].im-point[0].im);
    if (del < 1.0e-10) veca[k] = -1.0e30;
    else veca[k] = log(del);
  }
    
  for (j = 1; j < ncheby; j++){
    kmax = 0; fkmax = -2.0e30;
    for (k = 0; k < nc3; k++){
      if (veca[k] > fkmax){
	fkmax = veca[k];
	kmax = k;
      }
    }
    point[j] = samp[kmax];
    for (k = 0; k < nc3; k++){
      del = sqr(samp[k].re-point[j].re) + sqr(samp[k].im-point[j].im);
      if (del < 1.0e-10) veca[k] = -1.0e30;
      else veca[k] += log(del);
    }
    
  }
  for (fk = 1.0, j = 0; j < ncheby; j++) fk *= (sqr(point[j].re) + sqr(point[j].im));
  fk = sqrt(fk);
  fk = pow(fk,1.0/(double)(ncheby));
  /*for (j = 0; j < ncheby; j++) {
    point[j].re /= fk;
    point[j].im /= fk;
    }*/
  free(samp); free(veca);
  return (fk);
}

/**************************************************************************/

void check_function(zomplex *an, zomplex *samp, index_st *ist, par_st *par, double ene_target)
{
  FILE *pf; long j;
  zomplex f, x, xn, xm1, ctmp;
  double dx = 0.01, xunsc;

  pf = fopen("func.dat" , "w"); 
  for (x.im = 0.0, x.re = -2.0; x.re < 2.0; x.re += dx){
    xn.re = 1.0;
    xn.im = 0.0;
    for (f = an[0], j = 1; j < ist->ncheby; j++){
      xm1.re = x.re - samp[j-1].re;
      xm1.im = x.im - samp[j-1].im;
      cmul(xm1,xn,xn);
      cmul(an[j],xn,ctmp);
      
      f.re += ctmp.re;
      f.im += ctmp.im;
    }
    xunsc = ((x.re+2.0)* par->dE / 4.0 + par->Vmin);

    fprintf (pf,"%g %g %g\n",xunsc,f.re,filt_func(xunsc - ene_target, par->dt));
  }
  fclose(pf);
  return;
}

/*****************************************************************************/
