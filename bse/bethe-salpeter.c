/*****************************************************************************/

#include <float.h>
#include "fd.h"

/*****************************************************************************/

void bethe_salpeter(zomplex *bsmat, zomplex *direct, zomplex *exchange, double *h0mat, zomplex *psi, double *vz, zomplex *mux, zomplex *muy, zomplex * muz,
					double *mx, double *my, double *mz,zomplex *sx, zomplex *sy, zomplex *sz,zomplex *lx, zomplex *ly, zomplex *lz, zomplex* lsqr, zomplex* ls, long_st ist, par_st par)
{
  FILE *pf, *pf1, *pf2; 
  long a, b, i, j, k, l, ibs, jbs, jgamma, jgrid;
  double   *eval, *ev, os, msumx, msumy, msumz, mos, *pgrid;
  char str[100]; 
  zomplex sum, *mat,*h, *u, sumx, sumy, sumz;

  mat = (zomplex *) calloc(ist.ms2*ist.ms2, sizeof(zomplex));
  h = (zomplex *) calloc(ist.ms2*ist.ms2, sizeof(zomplex));
  u = (zomplex *) calloc(ist.ms2*ist.ms2, sizeof(zomplex));
  eval = (double *) calloc(ist.ms2, sizeof(double));

  printf("The number of electron-hole pairs in the exciton basis = %ld\n", ist.ms2);

#pragma omp parallel for private(i)
  for (i = 0; i < ist.ms2*ist.ms2; i++) {
    h[i].re = u[i].re = h0mat[i] - bsmat[i].re;
    h[i].im = u[i].im = -bsmat[i].im;
  }  
  diag((int)ist.ms2, ist.nthreads, u, eval);

  pf = fopen("HBSmatRE.dat", "w");
  for (i = 0; i < ist.ms2; i++, fprintf(pf,"\n")) {
    for (j = 0; j < ist.ms2; j++) {
      fprintf (pf,"%.*g \t", DBL_DIG, h[i*ist.ms2+j].re);
	}
  }
  fclose(pf);

  pf = fopen("HBSmatIM.dat", "w");
  for (i = 0; i < ist.ms2; i++, fprintf(pf,"\n")) {
    for (j = 0; j < ist.ms2; j++) {
      fprintf (pf,"%.*g \t", DBL_DIG,h[i*ist.ms2+j].im);
  }
  }
  fclose(pf);
  
  // Prints the coefficients for the 100 (or ist.ms2) lowest energy excitonic states
  long numExcStatesToPrint = 100;
  if (ist.ms2 < numExcStatesToPrint) numExcStatesToPrint = ist.ms2;
  pf = fopen("BSEcoeff.dat", "w");
  for (i = 0; i < ist.ms2; i++) {
    for (j = 0; j < numExcStatesToPrint; j++) {  
	  fprintf (pf,"{%.*g, %.*g}\t", DBL_DIG, u[i*ist.ms2+j].re, DBL_DIG, u[i*ist.ms2+j].im);
    } 
    fprintf (pf,"\n");	
  }
  fclose(pf);

  for (sum.re = sum.im = 0.0, i = 0; i < ist.ms2; i++) {
    for (j = 0; j < ist.ms2; j++) {
      sum.re += u[j*ist.ms2].re * (u[i*ist.ms2].re * h[j*ist.ms2+i].re - u[i*ist.ms2].im * h[j*ist.ms2+i].im)
          +  u[j*ist.ms2].im * (u[i*ist.ms2].re * h[j*ist.ms2+i].im + u[i*ist.ms2].im * h[j*ist.ms2+i].re);

      sum.im += u[j*ist.ms2].re * (u[i*ist.ms2].re * h[j*ist.ms2+i].im + u[i*ist.ms2].im * h[j*ist.ms2+i].re)
          -  u[j*ist.ms2].im * (u[i*ist.ms2].re * h[j*ist.ms2+i].re - u[i*ist.ms2].im * h[j*ist.ms2+i].im);
    }
  }
  printf("Ground state exciton has energy = %.10f (%.10f eV) (%.10f Imag)\n", sum.re, sum.re*AUTOEV, sum.im);


  //printf("eval[0]=%g\n",eval[0] );
  long *listibs = (long *) calloc(ist.ms2, sizeof(long));

  for (ibs = 0, a = ist.nlumo; a < ist.nlumo+ist.totallumo; a++) {
      for (i = 0; i < ist.totalhomo; i++, ibs++) {
          listibs[(a - ist.nlumo) * ist.totalhomo + i] = ibs;
          //printf("a:%ld i:%ld ibs:%ld\n",a,i,ibs);
      }
  }

//compute spins:
  FILE* spinpf =fopen("spins.dat", "w");  
  //printf("Spins:\n");
  zomplex spinx, spiny, spinz, spintot;
  zomplex tmpx, tmpy, tmpz, tmp, tmp2;
  long index, indexba, indexji,n;
  for (n=0;n<ist.ms2;n++){
    spinx.re = spinx.im = 0;
    spiny.re = spiny.im = 0;
    spinz.re = spinz.im = 0;
    spintot.re = spintot.im = 0;
    for (a = ist.nlumo;a<ist.nlumo+ist.totallumo;a++){
      for (i = 0; i < ist.totalhomo; i++) {
        ibs = listibs[(a - ist.nlumo) * ist.totalhomo + i];
        

        tmpx.re = tmpx.im = 0;
        tmpy.re = tmpy.im = 0;
        tmpz.re = tmpz.im = 0;
        //sum over b
        for (b = ist.nlumo;b<ist.nlumo+ist.totallumo;b++){
          jbs = listibs[(b - ist.nlumo) * ist.totalhomo + i];
          index = sqr(ist.totalhomo)+(a-ist.nlumo)*ist.totallumo+(b-ist.nlumo);
          
          //c_bi^* * <b| Sx |a>
          tmpx.re += u[jbs*ist.ms2+n].re*sx[index].re + u[jbs*ist.ms2+n].im*sx[index].im;
          tmpx.im += u[jbs*ist.ms2+n].re*sx[index].im - u[jbs*ist.ms2+n].im*sx[index].re;

          //c_bi^* * <b| Sy |a>
          tmpy.re += u[jbs*ist.ms2+n].re*sy[index].re + u[jbs*ist.ms2+n].im*sy[index].im;
          tmpy.im += u[jbs*ist.ms2+n].re*sy[index].im - u[jbs*ist.ms2+n].im*sy[index].re;

          //c_bi^* * <b| Sz |a>
          tmpz.re += u[jbs*ist.ms2+n].re*sz[index].re + u[jbs*ist.ms2+n].im*sz[index].im;
          tmpz.im += u[jbs*ist.ms2+n].re*sz[index].im - u[jbs*ist.ms2+n].im*sz[index].re;
        }
        
        //sum over j
        for (j = 0; j < ist.totalhomo; j++) {
          jbs = listibs[(a - ist.nlumo) * ist.totalhomo + j];
          index = i*ist.totalhomo+j;

          //c_aj^* * <j| Sx |i> 
          tmpx.re += u[jbs*ist.ms2+n].re*sx[index].re + u[jbs*ist.ms2+n].im*sx[index].im;
          tmpx.im += u[jbs*ist.ms2+n].re*sx[index].im - u[jbs*ist.ms2+n].im*sx[index].re;

          //c_aj^* * <j| Sy |i> 
          tmpy.re += u[jbs*ist.ms2+n].re*sy[index].re + u[jbs*ist.ms2+n].im*sy[index].im;
          tmpy.im += u[jbs*ist.ms2+n].re*sy[index].im - u[jbs*ist.ms2+n].im*sy[index].re;

          //c_aj^* * <j| Sz |i> 
          tmpz.re += u[jbs*ist.ms2+n].re*sz[index].re + u[jbs*ist.ms2+n].im*sz[index].im;
          tmpz.im += u[jbs*ist.ms2+n].re*sz[index].im - u[jbs*ist.ms2+n].im*sz[index].re;

        }

        //multiply by the c_ai coeff
        spinx.re += tmpx.re * u[ibs*ist.ms2+n].re - tmpx.im * u[ibs*ist.ms2+n].im;
        spinx.im += tmpx.im * u[ibs*ist.ms2+n].re + tmpx.re * u[ibs*ist.ms2+n].im;

        spiny.re += tmpy.re * u[ibs*ist.ms2+n].re - tmpy.im * u[ibs*ist.ms2+n].im;
        spiny.im += tmpy.im * u[ibs*ist.ms2+n].re + tmpy.re * u[ibs*ist.ms2+n].im;

        spinz.re += tmpz.re * u[ibs*ist.ms2+n].re - tmpz.im * u[ibs*ist.ms2+n].im;
        spinz.im += tmpz.im * u[ibs*ist.ms2+n].re + tmpz.re * u[ibs*ist.ms2+n].im;
        

        for (b = ist.nlumo;b<ist.nlumo+ist.totallumo;b++){
          for (j = 0; j < ist.totalhomo; j++) {
            
            indexba = sqr(ist.totalhomo)+(a-ist.nlumo)*ist.totallumo+(b-ist.nlumo);
            indexji = i*ist.totalhomo+j;

            jbs = listibs[(b - ist.nlumo) * ist.totalhomo + j];
            tmp.re=u[ibs*ist.ms2+n].re*u[jbs*ist.ms2+n].re
                  + u[ibs*ist.ms2+n].im*u[jbs*ist.ms2+n].im;

            tmp.im=u[ibs*ist.ms2+n].im*u[jbs*ist.ms2+n].re
                  - u[ibs*ist.ms2+n].re*u[jbs*ist.ms2+n].im;

            tmpx.re=sx[indexba].re*sx[indexji].re-sx[indexba].im*sx[indexji].im
             + sy[indexba].re*sy[indexji].re-sy[indexba].im*sy[indexji].im
             + sz[indexba].re*sz[indexji].re-sz[indexba].im*sz[indexji].im;

            tmpx.im=sx[indexba].im*sx[indexji].re+sx[indexba].re*sx[indexji].im
            + sy[indexba].im*sy[indexji].re+sy[indexba].re*sy[indexji].im
            + sz[indexba].im*sz[indexji].re+sz[indexba].re*sz[indexji].im;   


            spintot.re+=tmp.re*tmpx.re-tmp.im*tmpx.im;
            spintot.im+=tmp.im*tmpx.re+tmp.re*tmpx.im;


          }

        }



      }

    }
    fprintf(spinpf,"%ld\t%-10.5lf\t%-10.5lf\t%-10.5lf\t",n,spinx.re,spiny.re, spinz.re);
    fprintf(spinpf,"%-10.5lf\t (%-10.5lf)\n",1.5+2.0*spintot.re, 2.0*spintot.im);
  }
  fclose(spinpf);


//compute orbital momentum
  FILE* orbitpf =fopen("orbital.dat", "w");  
  //printf("Spins:\n");
  zomplex orbitx, orbity, orbitz, orbittot;
//  zomplex tmpx, tmpy, tmpz, tmp, tmp2;
//  long index, indexba, indexji,n;
  for (n=0;n<ist.ms2;n++){
    orbitx.re = orbitx.im = 0;
    orbity.re = orbity.im = 0;
    orbitz.re = orbitz.im = 0;
    orbittot.re = orbittot.im = 0;
    for (a = ist.nlumo;a<ist.nlumo+ist.totallumo;a++){
      for (i = 0; i < ist.totalhomo; i++) {
        ibs = listibs[(a - ist.nlumo) * ist.totalhomo + i];
        

        tmpx.re = tmpx.im = 0;
        tmpy.re = tmpy.im = 0;
        tmpz.re = tmpz.im = 0;
        //sum over b
        for (b = ist.nlumo;b<ist.nlumo+ist.totallumo;b++){
          jbs = listibs[(b - ist.nlumo) * ist.totalhomo + i];
          index = sqr(ist.totalhomo)+(a-ist.nlumo)*ist.totallumo+(b-ist.nlumo);
          
          //c_bi^* * <b| Lx |a>
          tmpx.re += u[jbs*ist.ms2+n].re*lx[index].re + u[jbs*ist.ms2+n].im*lx[index].im;
          tmpx.im += u[jbs*ist.ms2+n].re*lx[index].im - u[jbs*ist.ms2+n].im*lx[index].re;

          //c_bi^* * <b| Ly |a>
          tmpy.re += u[jbs*ist.ms2+n].re*ly[index].re + u[jbs*ist.ms2+n].im*ly[index].im;
          tmpy.im += u[jbs*ist.ms2+n].re*ly[index].im - u[jbs*ist.ms2+n].im*ly[index].re;

          //c_bi^* * <b| Lz |a>
          tmpz.re += u[jbs*ist.ms2+n].re*lz[index].re + u[jbs*ist.ms2+n].im*lz[index].im;
          tmpz.im += u[jbs*ist.ms2+n].re*lz[index].im - u[jbs*ist.ms2+n].im*lz[index].re;
        }
        
        //sum over j
        for (j = 0; j < ist.totalhomo; j++) {
          jbs = listibs[(a - ist.nlumo) * ist.totalhomo + j];
          index = i*ist.totalhomo+j;

          //c_aj^* * <j| Lx |i> 
          tmpx.re += u[jbs*ist.ms2+n].re*lx[index].re + u[jbs*ist.ms2+n].im*lx[index].im;
          tmpx.im += u[jbs*ist.ms2+n].re*lx[index].im - u[jbs*ist.ms2+n].im*lx[index].re;

          //c_aj^* * <j| Ly |i> 
          tmpy.re += u[jbs*ist.ms2+n].re*ly[index].re + u[jbs*ist.ms2+n].im*ly[index].im;
          tmpy.im += u[jbs*ist.ms2+n].re*ly[index].im - u[jbs*ist.ms2+n].im*ly[index].re;

          //c_aj^* * <j| Lz |i> 
          tmpz.re += u[jbs*ist.ms2+n].re*lz[index].re + u[jbs*ist.ms2+n].im*lz[index].im;
          tmpz.im += u[jbs*ist.ms2+n].re*lz[index].im - u[jbs*ist.ms2+n].im*lz[index].re;

        }

        //multiply by the c_ai coeff
        orbitx.re += tmpx.re * u[ibs*ist.ms2+n].re - tmpx.im * u[ibs*ist.ms2+n].im;
        orbitx.im += tmpx.im * u[ibs*ist.ms2+n].re + tmpx.re * u[ibs*ist.ms2+n].im;

        orbity.re += tmpy.re * u[ibs*ist.ms2+n].re - tmpy.im * u[ibs*ist.ms2+n].im;
        orbity.im += tmpy.im * u[ibs*ist.ms2+n].re + tmpy.re * u[ibs*ist.ms2+n].im;

        orbitz.re += tmpz.re * u[ibs*ist.ms2+n].re - tmpz.im * u[ibs*ist.ms2+n].im;
        orbitz.im += tmpz.im * u[ibs*ist.ms2+n].re + tmpz.re * u[ibs*ist.ms2+n].im;
        
        //Lsqr part
        for (b = ist.nlumo;b<ist.nlumo+ist.totallumo;b++){
          for (j = 0; j < ist.totalhomo; j++) {
            
            indexba = sqr(ist.totalhomo)+(a-ist.nlumo)*ist.totallumo+(b-ist.nlumo);
            indexji = i*ist.totalhomo+j;

            jbs = listibs[(b - ist.nlumo) * ist.totalhomo + j];
            
            //c_{ai}^n * (c_{bj}^n)^*
            tmp.re=u[ibs*ist.ms2+n].re*u[jbs*ist.ms2+n].re
                  + u[ibs*ist.ms2+n].im*u[jbs*ist.ms2+n].im;

            tmp.im=u[ibs*ist.ms2+n].im*u[jbs*ist.ms2+n].re
                  - u[ibs*ist.ms2+n].re*u[jbs*ist.ms2+n].im;


            tmpx.re=lx[indexba].re*lx[indexji].re-lx[indexba].im*lx[indexji].im
             + ly[indexba].re*ly[indexji].re-ly[indexba].im*ly[indexji].im
             + lz[indexba].re*lz[indexji].re-lz[indexba].im*lz[indexji].im;

            tmpx.im=lx[indexba].im*lx[indexji].re+lx[indexba].re*lx[indexji].im
            + ly[indexba].im*ly[indexji].re+ly[indexba].re*ly[indexji].im
            + lz[indexba].im*lz[indexji].re+lz[indexba].re*lz[indexji].im;   

            tmpx.re*=2.0; tmpx.im*=2.0;
            
            
            if (i==j){
              tmpx.re+=lsqr[indexba].re; tmpx.im+=lsqr[indexba].im;
            }

            if (a==b){
              tmpx.re+=lsqr[indexji].re; tmpx.im+=lsqr[indexji].im;
            }
            
            //printf("a:%ld b:%ld i:%ld j:%ld    tmpx: (%lf, %lf)\n", a,b,i,j,tmpx.re, tmpx.im);

            orbittot.re+=tmp.re*tmpx.re-tmp.im*tmpx.im;
            orbittot.im+=tmp.im*tmpx.re+tmp.re*tmpx.im;


          }

        }



      }

    }
    fprintf(orbitpf,"%ld\t%-10.5lf\t%-10.5lf\t%-10.5lf\t",n,orbitx.re,orbity.re, orbitz.re);
    fprintf(orbitpf,"%-10.5lf\t (%-10.5lf)\n",orbittot.re, orbittot.im);
  }
  fclose(orbitpf);




  //compute ls momentum
  FILE* lspf =fopen("couple.dat", "w");  
  //printf("Spins:\n");
  zomplex lstot;
//  zomplex tmpx, tmpy, tmpz, tmp, tmp2;
//  long index, indexba, indexji,n;
  for (n=0;n<ist.ms2;n++){
    lstot.re = lstot.im = 0;
    for (a = ist.nlumo;a<ist.nlumo+ist.totallumo;a++){
      for (i = 0; i < ist.totalhomo; i++) {
        ibs = listibs[(a - ist.nlumo) * ist.totalhomo + i];
        for (b = ist.nlumo;b<ist.nlumo+ist.totallumo;b++){
          for (j = 0; j < ist.totalhomo; j++) {
            
            tmpx.re = 0.0;tmpx.im = 0.0;

            indexba = sqr(ist.totalhomo)+(a-ist.nlumo)*ist.totallumo+(b-ist.nlumo);
            indexji = i*ist.totalhomo+j;

            jbs = listibs[(b - ist.nlumo) * ist.totalhomo + j];
            
            //c_{ai}^n * (c_{bj}^n)^*
            tmp.re=u[ibs*ist.ms2+n].re*u[jbs*ist.ms2+n].re
                  + u[ibs*ist.ms2+n].im*u[jbs*ist.ms2+n].im;

            tmp.im=u[ibs*ist.ms2+n].im*u[jbs*ist.ms2+n].re
                  - u[ibs*ist.ms2+n].re*u[jbs*ist.ms2+n].im;

            
            // <a|L|b>*<j|S|i>
            tmpx.re=lx[indexba].re*sx[indexji].re-lx[indexba].im*sx[indexji].im
                  + ly[indexba].re*sy[indexji].re-ly[indexba].im*sy[indexji].im
                  + lz[indexba].re*sz[indexji].re-lz[indexba].im*sz[indexji].im;

            tmpx.im=lx[indexba].im*sx[indexji].re+lx[indexba].re*sx[indexji].im
                  + ly[indexba].im*sy[indexji].re+ly[indexba].re*sy[indexji].im
                  + lz[indexba].im*sz[indexji].re+lz[indexba].re*sz[indexji].im;

            // <a|S|b>*<j|L|i>
            tmpx.re+=sx[indexba].re*lx[indexji].re-sx[indexba].im*lx[indexji].im
                   + sy[indexba].re*ly[indexji].re-sy[indexba].im*ly[indexji].im
                   + sz[indexba].re*lz[indexji].re-sz[indexba].im*lz[indexji].im;

            tmpx.im+=sx[indexba].im*lx[indexji].re+sx[indexba].re*lx[indexji].im
                   + sy[indexba].im*ly[indexji].re+sy[indexba].re*ly[indexji].im
                   + sz[indexba].im*lz[indexji].re+sz[indexba].re*lz[indexji].im;    

            
            
            //delta ij part
            if (i==j){
              tmpx.re+=ls[indexba].re; tmpx.im+=ls[indexba].im;
              
            }
            //delta ab part
            if (a==b){
              tmpx.re+=ls[indexji].re; tmpx.im+=ls[indexji].im;
            }
            
            //printf("a:%ld b:%ld i:%ld j:%ld    tmpx: (%lf, %lf)\n", a,b,i,j,tmpx.re, tmpx.im);

            lstot.re+=tmp.re*tmpx.re-tmp.im*tmpx.im;
            lstot.im+=tmp.im*tmpx.re+tmp.re*tmpx.im;


          }

        }



      }

    }
    fprintf(lspf,"%ld\t%-10.5lf\t (%-10.5lf)\n",n,lstot.re, lstot.im);
  }
  fclose(lspf);
/*
//compute direct and exchange parts
  FILE* enepartspf = fopen("eneparts.dat", "w");
  fprintf(enepartspf,"# n \t <H_dir> \t <H_exc> \t <H0>\n");
  double directE, exchangeE, E0;
  for (n=0;n<ist.ms2;n++){
    directE = exchangeE = E0= 0.0;

    for (a = ist.nlumo;a<ist.nlumo+ist.totallumo;a++){
      for (i = 0; i < ist.totalhomo; i++) {
        ibs = listibs[(a - ist.nlumo) * ist.totalhomo + i];
        for (b = ist.nlumo;b<ist.nlumo+ist.totallumo;b++){
          for (j = 0; j < ist.totalhomo; j++) {

            
            
            //(c_{ai}^n)^* * (c_{bj}^n)
            tmp.re=u[ibs*ist.ms2+n].re*u[jbs*ist.ms2+n].re
                  + u[ibs*ist.ms2+n].im*u[jbs*ist.ms2+n].im;

            tmp.im=u[ibs*ist.ms2+n].re*u[jbs*ist.ms2+n].im
                  - u[ibs*ist.ms2+n].im*u[jbs*ist.ms2+n].re;

            // (tmp)*H_{ai,bj}
            //directE+= direct[ibs * ist.ms2 + jbs].re * tmp.re 
            //        - direct[ibs * ist.ms2 + jbs].im * tmp.im;

            //exchangeE+= exchange[ibs * ist.ms2 + jbs].re * tmp.re 
            //          - exchange[ibs * ist.ms2 + jbs].im * tmp.im;

            //E0+= h0mat[ibs * ist.ms2 + jbs] * tmp.re;
            
            E0+=h[jbs+ist.ms2+ibs].re *tmp.re - h[jbs+ist.ms2+ibs].im *tmp.im;

          }
        }
      }
    }

    fprintf(enepartspf,"%ld %lg %lg %lg\n", n, directE, exchangeE, E0);
  }
  fclose(enepartspf);*/









  //compute $mat = h \cdot u$
#pragma omp parallel for private(l,j,k,sum)
  for (l = 0; l < ist.ms2; l++) {
    for (j = 0; j < ist.ms2; j++) {
      for (sum.re=sum.im = 0, k = 0; k < ist.ms2; k++) {
      	sum.re +=  h[l*ist.ms2+k].re * u[k*ist.ms2+j].re - h[l*ist.ms2+k].im * u[k*ist.ms2+j].im;
        sum.im +=  h[l*ist.ms2+k].im * u[k*ist.ms2+j].re + h[l*ist.ms2+k].re * u[k*ist.ms2+j].im;
      }
      mat[l*ist.ms2+j].re = sum.re;
      mat[l*ist.ms2+j].im = sum.im;
    }
  }


  //compute $u^\dagger \cdot mat = u^\dagger \cdot h \cdot u$
#pragma omp parallel for private(i,j,k,sum)
  for (i = 0; i < ist.ms2; i++) {
    for (j = 0; j < ist.ms2; j++) {
      for (sum.re=sum.im = 0, l = 0; l < ist.ms2; l++) {
      	sum.re +=   u[l*ist.ms2+i].re * mat[l*ist.ms2+j].re + u[l*ist.ms2+i].im * mat[l*ist.ms2+j].im;
        sum.im +=  -u[l*ist.ms2+i].im * mat[l*ist.ms2+j].re + u[l*ist.ms2+i].re * mat[l*ist.ms2+j].im;
      }
      h[i*ist.ms2+j].re = sum.re;
      h[i*ist.ms2+j].im = sum.im;
    }
  }

#pragma omp parallel for private(l,j,k,sum)
  for (l = 0; l < ist.ms2; l++) {
    for (j = 0; j < ist.ms2; j++) {
      for (sum.re=sum.im = 0, k = 0; k < ist.ms2; k++) {
        sum.re +=  h0mat[l*ist.ms2+k] * u[k*ist.ms2+j].re;
        sum.im +=  h0mat[l*ist.ms2+k] * u[k*ist.ms2+j].im;
      }
      mat[l*ist.ms2+j].re = sum.re;
      mat[l*ist.ms2+j].im = sum.im;
    }
  }

#pragma omp parallel for private(i,l,sum)
  for (i = 0; i < ist.ms2; i++) {
    for (j = 0; j < ist.ms2; j++) {
      for (sum.re=sum.im = 0, l = 0; l < ist.ms2; l++) {
        sum.re +=   u[l*ist.ms2+i].re * mat[l*ist.ms2+j].re + u[l*ist.ms2+i].im * mat[l*ist.ms2+j].im;
        sum.im +=  -u[l*ist.ms2+i].im * mat[l*ist.ms2+j].re + u[l*ist.ms2+i].re * mat[l*ist.ms2+j].im;

      }
      h0mat[i*ist.ms2+j] = sum.re;
    }
  }


#pragma omp parallel for private(l,j,k,sum)
  for (l = 0; l < ist.ms2; l++) {
    for (j = 0; j < ist.ms2; j++) {
      for (sum.re=sum.im = 0, k = 0; k < ist.ms2; k++) {
        sum.re +=  direct[l*ist.ms2+k].re * u[k*ist.ms2+j].re - direct[l*ist.ms2+k].im * u[k*ist.ms2+j].im;
        sum.im +=  direct[l*ist.ms2+k].im * u[k*ist.ms2+j].re + direct[l*ist.ms2+k].re * u[k*ist.ms2+j].im;
      }
      mat[l*ist.ms2+j].re = sum.re;
      mat[l*ist.ms2+j].im = sum.im;
    }
  }


  //compute $u^\dagger \cdot mat = u^\dagger \cdot h \cdot u$
#pragma omp parallel for private(i,j,k,sum)
  for (i = 0; i < ist.ms2; i++) {
    for (j = 0; j < ist.ms2; j++) {
      for (sum.re=sum.im = 0, l = 0; l < ist.ms2; l++) {
        sum.re +=   u[l*ist.ms2+i].re * mat[l*ist.ms2+j].re + u[l*ist.ms2+i].im * mat[l*ist.ms2+j].im;
        sum.im +=  -u[l*ist.ms2+i].im * mat[l*ist.ms2+j].re + u[l*ist.ms2+i].re * mat[l*ist.ms2+j].im;
      }
      direct[i*ist.ms2+j].re = sum.re;
      direct[i*ist.ms2+j].im = sum.im;
    }
  }

  #pragma omp parallel for private(l,j,k,sum)
  for (l = 0; l < ist.ms2; l++) {
    for (j = 0; j < ist.ms2; j++) {
      for (sum.re=sum.im = 0, k = 0; k < ist.ms2; k++) {
        sum.re +=  exchange[l*ist.ms2+k].re * u[k*ist.ms2+j].re - exchange[l*ist.ms2+k].im * u[k*ist.ms2+j].im;
        sum.im +=  exchange[l*ist.ms2+k].im * u[k*ist.ms2+j].re + exchange[l*ist.ms2+k].re * u[k*ist.ms2+j].im;
      }
      mat[l*ist.ms2+j].re = sum.re;
      mat[l*ist.ms2+j].im = sum.im;
    }
  }


  //compute $u^\dagger \cdot mat = u^\dagger \cdot h \cdot u$
#pragma omp parallel for private(i,j,k,sum)
  for (i = 0; i < ist.ms2; i++) {
    for (j = 0; j < ist.ms2; j++) {
      for (sum.re=sum.im = 0, l = 0; l < ist.ms2; l++) {
        sum.re +=   u[l*ist.ms2+i].re * mat[l*ist.ms2+j].re + u[l*ist.ms2+i].im * mat[l*ist.ms2+j].im;
        sum.im +=  -u[l*ist.ms2+i].im * mat[l*ist.ms2+j].re + u[l*ist.ms2+i].re * mat[l*ist.ms2+j].im;
      }
      exchange[i*ist.ms2+j].re = sum.re;
      exchange[i*ist.ms2+j].im = sum.im;
    }
  }





  // Print out the energies of the excitonic states
  pf = fopen("exciton.dat" , "w");
  fprintf(pf,"#n \t E_n \t <H> \t <H_dir> \t <H_exc> \t <H_0> \t E_B (eV)\n");
  for (i = 0; i < ist.ms2; i++) {  
    fprintf(pf,"%ld % .12f % .12f % .12f  % .12f  % .12f  % .12f\n", i, eval[i], h[i*ist.ms2+i].re, direct[i*ist.ms2+i].re, exchange[i*ist.ms2+i].re,
     h0mat[i*ist.ms2+i], (eval[i]-h0mat[i*ist.ms2+i])*AUTOEV);
  }
  fclose(pf);

  // Calculate and print the electric and magnetic dipole strengths and the rotational strength
  pf = fopen("OS.dat", "w");
  pf1 = fopen("M.dat", "w");
  pf2 = fopen("rs.dat", "w");
  for (jgamma = 0; jgamma < ist.ms2; jgamma++) {
    sumx.re = sumx.im= 0.0; sumy.re = sumy.im = 0.0; sumz.re =sumz.im = 0.0;
    //msumx = 0.0; msumy = 0.0; msumz = 0.0;
    for (ibs = 0, a = ist.nlumo; a < ist.nlumo+ist.totallumo; a++) {
      for (i = 0; i < ist.totalhomo; i++, ibs++) {
        sumx.re += u[ibs*ist.ms2 + jgamma].re * mux[i*ist.totallumo + (a - ist.nlumo)].re
                 - u[ibs*ist.ms2 + jgamma].im * mux[i*ist.totallumo + (a - ist.nlumo)].im;
        sumx.im += u[ibs*ist.ms2 + jgamma].re * mux[i*ist.totallumo + (a - ist.nlumo)].im
                 + u[ibs*ist.ms2 + jgamma].im * mux[i*ist.totallumo + (a - ist.nlumo)].re;
        
        sumy.re += u[ibs*ist.ms2 + jgamma].re * muy[i*ist.totallumo + (a - ist.nlumo)].re
                 - u[ibs*ist.ms2 + jgamma].im * muy[i*ist.totallumo + (a - ist.nlumo)].im;
        sumy.im += u[ibs*ist.ms2 + jgamma].re * muy[i*ist.totallumo + (a - ist.nlumo)].im
                 + u[ibs*ist.ms2 + jgamma].im * muy[i*ist.totallumo + (a - ist.nlumo)].re;
        
        sumz.re += u[ibs*ist.ms2 + jgamma].re * muz[i*ist.totallumo + (a - ist.nlumo)].re
                 - u[ibs*ist.ms2 + jgamma].im * muz[i*ist.totallumo + (a - ist.nlumo)].im;
        sumz.im += u[ibs*ist.ms2 + jgamma].re * muz[i*ist.totallumo + (a - ist.nlumo)].im
                 + u[ibs*ist.ms2 + jgamma].im * muz[i*ist.totallumo + (a - ist.nlumo)].re;
        //msumx += u[jgamma*ist.ms2 + ibs] * mx[i*ist.totallumo + (a - ist.nlumo)];
		    //msumy += u[jgamma*ist.ms2 + ibs] * my[i*ist.totallumo + (a - ist.nlumo)];
		    //msumz += u[jgamma*ist.ms2 + ibs] * mz[i*ist.totallumo + (a - ist.nlumo)];
      }
    } 
    os  = (sqr(sumx.re) +sqr(sumx.im) + sqr(sumy.re) + sqr(sumy.im) + sqr(sumz.re) +sqr(sumz.im));
    //mos = (msumx*msumx + msumy*msumy + msumz*msumz);
    fprintf(pf,  "%ld %.8f %.8f % .8f % .12f % .12f % .12f % .12f % .12f % .12f\n", jgamma, sqrt(os), eval[jgamma], (2.0/3.0)*eval[jgamma]*os, 
      				sumx.re, sumx.im, sumy.re, sumy.im, sumz.re, sumz.im);
    //fprintf(pf1, "%ld %.8f %.8f % .8f % .12f % .12f % .12f\n", jgamma, sqrt(mos), eval[jgamma], (4.0/3.0)*eval[jgamma]*mos,
		//			    msumx*msumx, msumy*msumy, msumz*msumz);
    //fprintf(pf2, "%ld %.8f % .16f\n", jgamma, eval[jgamma], (sumx*msumx + sumy*msumy + sumz*msumz));
  }
  fclose(pf); fclose(pf1); fclose(pf2);
  /*
  // Electron and hole carrier densities 
  pgrid = (double*)calloc(ist.ngrid,sizeof(double));
  for (jgamma = 0; jgamma < 10; jgamma++) {
	// Electron densities
    for (jgrid = 0; jgrid < ist.ngrid; jgrid++) { 
	  pgrid[jgrid] = 0.0;
	}
    for (ibs = 0, a = ist.nlumo; a < ist.nlumo+ist.totallumo; a++) {
      for (i = 0; i < ist.totalhomo; i++, ibs++) {
      	for (jbs = 0, b = ist.nlumo; b < ist.nlumo+ist.totallumo; b++) {
      	  for (j = 0; j < ist.totalhomo; j++, jbs++) {
      	    if (i == j) { 
      	      sum = u[jgamma*ist.ms2+ibs] * u[jgamma*ist.ms2+jbs];
            #pragma omp parallel for private(jgrid)
      	      for (jgrid = 0; jgrid < ist.ngrid; jgrid++)
            		pgrid[jgrid] += sum * psi[a*ist.ngrid+jgrid] * psi[b*ist.ngrid+jgrid];
      	    }
      	  }
      	}
      }
    }
    sprintf(str, "pcz%d-bs.dat", jgamma);
    print_pz_one(pgrid, vz, par, ist, str);
    sprintf(str, "pe-bs-%d.cub", jgamma);
    print_cube(pgrid, ist, par, str);
    
	// Hole densities
    for (jgrid = 0; jgrid < ist.ngrid; jgrid++) {
	  pgrid[jgrid] = 0.0;
    }
	for (ibs = 0, a = ist.nlumo; a < ist.nlumo+ist.totallumo; a++) {
      for (i = 0; i < ist.totalhomo; i++, ibs++) {
      	for (jbs = 0, b = ist.nlumo; b < ist.nlumo+ist.totallumo; b++) {
      	  for (j = 0; j < ist.totalhomo; j++, jbs++) {
      	    if (a == b) { 
      	      sum = u[jgamma*ist.ms2+ibs] * u[jgamma*ist.ms2+jbs];
            #pragma omp parallel for private(jgrid)
      	      for (jgrid = 0; jgrid < ist.ngrid; jgrid++)
            		pgrid[jgrid] += sum * psi[i*ist.ngrid+jgrid] * psi[j*ist.ngrid+jgrid];
      	    }
      	  }
      	}
      }
    }
    sprintf(str, "pvz%d-bs.dat", jgamma);
    print_pz_one(pgrid, vz, par, ist, str);
    sprintf(str, "ph-bs-%d.cub", jgamma);
    print_cube(pgrid, ist, par, str);
  }

  // Calculate pegrid and phgrid for the lowest 15 electron and hole state, respectively 
  // this is a noninteracting result -> test locations does not influence the result
  long nNonIntStates = 20;
  if (ist.totalhomo < nNonIntStates) nNonIntStates = ist.totalhomo;
  if (ist.totallumo < nNonIntStates) nNonIntStates = ist.totallumo;
  for (i = 0; i < nNonIntStates; i++) {
  #pragma omp parallel for private(jgrid)
    for (jgrid = 0; jgrid < ist.ngrid; jgrid++) { // electron first
      pgrid[jgrid] = par.dv *  psi[(ist.nlumo+i)*ist.ngrid + jgrid] * psi[(ist.nlumo+i)*ist.ngrid + jgrid];
    }
    norm_vector(pgrid, par.dv, ist.ngrid);
    sprintf(str, "pe-ni-%d.cub", i);
    print_cube(pgrid, ist, par, str);
    sprintf(str, "pcz%d-ni.dat", i); 
    print_pz_one(pgrid, vz, par, ist, str);
  #pragma omp parallel for private(jgrid)
    for (jgrid = 0; jgrid < ist.ngrid; jgrid++) { // hole second
      pgrid[jgrid] = par.dv *  psi[(ist.nlumo-1-i)*ist.ngrid + jgrid] * psi[(ist.nlumo-1-i)*ist.ngrid + jgrid];
    }
    norm_vector(pgrid, par.dv, ist.ngrid);
    sprintf(str, "ph-ni-%d.cub", i);
    print_cube(pgrid, ist, par, str);
    sprintf(str, "pvz%d-ni.dat", i); 
    print_pz_one(pgrid, vz, par, ist, str);
  }  

  // Calculate the pegrid and phgrid for the lowest excitonic state 
  // for fixed locations of the other carrier - an interacting result
  if (ist.printFPDensity) {
    print_fixed_qp_density(psi, u, vz, ist, par);
  } 
  
  free(pgrid);
  */
  free(u); free(h); free(eval); free(mat);

  return;
}

/*****************************************************************************/
