/*****************************************************************************/

#include <float.h>
#include "fd.h"

/*****************************************************************************/

void bethe_salpeter(zomplex *bsmat, zomplex *direct, zomplex *exchange, zomplex *bs_coeff, double *h0mat, double *xton_ene, zomplex *psi, 
					xyz_st *s_mom, xyz_st *l_mom, zomplex* l2_mom, zomplex* LdotS, grid_st *grid, index_st *ist, par_st *par)
{
  FILE *pf; 
  long a, b, i, j, k, l, ibs, jbs; 
  zomplex sum, *mat,*h;

  mat = (zomplex *) calloc(ist->n_xton*ist->n_xton, sizeof(zomplex));
  h = (zomplex *) calloc(ist->n_xton*ist->n_xton, sizeof(zomplex));

#pragma omp parallel for private(i)
  for (i = 0; i < ist->n_xton*ist->n_xton; i++) {
    h[i].re = bs_coeff[i].re = h0mat[i] - bsmat[i].re;
    h[i].im = bs_coeff[i].im = -bsmat[i].im;
  }  
  
  //
  // Diagonalize the BSE matrix to obtain the coefficients
  // But first, convert from row major to column major order
  
  diag((int)ist->n_xton, ist->nthreads,(double _Complex*) bs_coeff, xton_ene);
  zomplex *tmp_psi = calloc(ist->n_xton * ist->n_xton, sizeof(zomplex));
  memcpy(tmp_psi, bs_coeff, ist->n_xton * ist->n_xton * sizeof(zomplex));
  for (ibs = 0; ibs < ist->n_xton; ibs++){
    for (jbs = 0; jbs < ist->n_xton; jbs++){
      bs_coeff[jbs*ist->n_xton + ibs].re = tmp_psi[ibs*ist->n_xton + jbs].re;
      bs_coeff[jbs*ist->n_xton + ibs].im = -tmp_psi[ibs*ist->n_xton + jbs].im;
    }
  }
  free(tmp_psi);
  //
  //

  pf = fopen("HBSmatRE.dat", "w");
  for (i = 0; i < ist->n_xton; i++, fprintf(pf,"\n")) {
    for (j = 0; j < ist->n_xton; j++) {
      fprintf (pf,"%.*g \t", PR_LEN, h[i*ist->n_xton+j].re);
	}
  }
  fclose(pf);

  pf = fopen("HBSmatIM.dat", "w");
  for (i = 0; i < ist->n_xton; i++, fprintf(pf,"\n")) {
    for (j = 0; j < ist->n_xton; j++) {
      fprintf (pf,"%.*g \t", PR_LEN,h[i*ist->n_xton+j].im);
  }
  }
  fclose(pf);
  
  // Prints the coefficients for the 100 (or ist->n_xton) lowest energy excitonic states
  long numExcStatesToPrint = 100;
  if (ist->n_xton < numExcStatesToPrint) numExcStatesToPrint = ist->n_xton;

  pf = fopen("BSEcoeffRE.dat", "w");
  for (i = 0; i < ist->n_xton; i++, fprintf(pf,"\n")) {
    for (j = 0; j < ist->n_xton; j++) {  
	  fprintf (pf,"%.*g\t", PR_LEN, bs_coeff[i*ist->n_xton+j].re);
    }
  }
  fclose(pf);

  pf = fopen("BSEcoeffIM.dat", "w");
  for (i = 0; i < ist->n_xton; i++, fprintf(pf,"\n")) {
    for (j = 0; j < ist->n_xton; j++) {  
	  fprintf (pf,"%.*g\t", PR_LEN, bs_coeff[i*ist->n_xton+j].im);
    }
  }
  fclose(pf);

  sum.re = sum.im = 0.0;
  for (i = 0; i < ist->n_xton; i++) {
    for (j = 0; j < ist->n_xton; j++) {
      sum.re += bs_coeff[j*ist->n_xton].re * (bs_coeff[i*ist->n_xton].re * h[j*ist->n_xton+i].re - bs_coeff[i*ist->n_xton].im * h[j*ist->n_xton+i].im)
          +  bs_coeff[j*ist->n_xton].im * (bs_coeff[i*ist->n_xton].re * h[j*ist->n_xton+i].im + bs_coeff[i*ist->n_xton].im * h[j*ist->n_xton+i].re);

      sum.im += bs_coeff[j*ist->n_xton].re * (bs_coeff[i*ist->n_xton].re * h[j*ist->n_xton+i].im + bs_coeff[i*ist->n_xton].im * h[j*ist->n_xton+i].re)
          -  bs_coeff[j*ist->n_xton].im * (bs_coeff[i*ist->n_xton].re * h[j*ist->n_xton+i].re - bs_coeff[i*ist->n_xton].im * h[j*ist->n_xton+i].im);
    }
  }
  
  printf("\nGround state exciton has energy = %.5f a.u. | %.5f eV (%.5f Imag)\n", sum.re, sum.re*AUTOEV, sum.im);


  //printf("xton_ene[0]=%g\n",xton_ene[0] );
  long *listibs = (long *) calloc(ist->n_xton, sizeof(long));

  for (ibs = 0, a = ist->lumo_idx; a < ist->lumo_idx + ist->n_elecs; a++) {
      for (i = 0; i < ist->n_holes; i++, ibs++) {
          listibs[(a - ist->lumo_idx) * ist->n_holes + i] = ibs;
          //printf("a:%ld i:%ld ibs:%ld\n",a,i,ibs);
      }
  }

//compute spins:
  FILE* spinpf =fopen("spins.dat", "w");  
  //printf("Spins:\n");
  zomplex spinx, spiny, spinz, spintot;
  zomplex tmpx, tmpy, tmpz, tmp;
  long index, indexba, indexji,n;
  // printf("\nspins block\n\n"); fflush(0);
  for (n = 0; n < ist->n_xton; n++){
    spinx.re = spinx.im = 0;
    spiny.re = spiny.im = 0;
    spinz.re = spinz.im = 0;
    spintot.re = spintot.im = 0;
    for (a = ist->lumo_idx; a < ist->lumo_idx + ist->n_elecs; a++){
      for (i = 0; i < ist->n_holes; i++) {
        ibs = listibs[(a - ist->lumo_idx) * ist->n_holes + i];
        

        tmpx.re = tmpx.im = 0;
        tmpy.re = tmpy.im = 0;
        tmpz.re = tmpz.im = 0;
        //sum over b
        for (b = ist->lumo_idx; b < ist->lumo_idx + ist->n_elecs; b++){
          jbs = listibs[(b - ist->lumo_idx) * ist->n_holes + i];
          index = sqr(ist->n_holes) + (a - ist->lumo_idx) * ist->n_elecs + (b - ist->lumo_idx);
          
          //c_bi^* * <b| Sx |a>
          tmpx.re += bs_coeff[jbs*ist->n_xton+n].re*s_mom[index].x_re + bs_coeff[jbs*ist->n_xton+n].im*s_mom[index].x_im;
          tmpx.im += bs_coeff[jbs*ist->n_xton+n].re*s_mom[index].x_im - bs_coeff[jbs*ist->n_xton+n].im*s_mom[index].x_re;

          //c_bi^* * <b| Sy |a>
          tmpy.re += bs_coeff[jbs*ist->n_xton+n].re*s_mom[index].y_re + bs_coeff[jbs*ist->n_xton+n].im*s_mom[index].y_im;
          tmpy.im += bs_coeff[jbs*ist->n_xton+n].re*s_mom[index].y_im - bs_coeff[jbs*ist->n_xton+n].im*s_mom[index].y_re;

          //c_bi^* * <b| Sz |a>
          tmpz.re += bs_coeff[jbs*ist->n_xton+n].re*s_mom[index].z_re + bs_coeff[jbs*ist->n_xton+n].im*s_mom[index].z_im;
          tmpz.im += bs_coeff[jbs*ist->n_xton+n].re*s_mom[index].z_im - bs_coeff[jbs*ist->n_xton+n].im*s_mom[index].z_re;
        }
        
        //sum over j
        for (j = 0; j < ist->n_holes; j++) {
          jbs = listibs[(a - ist->lumo_idx) * ist->n_holes + j];
          index = i*ist->n_holes+j;

          //c_aj^* * <j| Sx |i> 
          tmpx.re += bs_coeff[jbs*ist->n_xton+n].re*s_mom[index].x_re + bs_coeff[jbs*ist->n_xton+n].im*s_mom[index].x_im;
          tmpx.im += bs_coeff[jbs*ist->n_xton+n].re*s_mom[index].x_im - bs_coeff[jbs*ist->n_xton+n].im*s_mom[index].x_re;

          //c_aj^* * <j| Sy |i> 
          tmpy.re += bs_coeff[jbs*ist->n_xton+n].re*s_mom[index].y_re + bs_coeff[jbs*ist->n_xton+n].im*s_mom[index].y_im;
          tmpy.im += bs_coeff[jbs*ist->n_xton+n].re*s_mom[index].y_im - bs_coeff[jbs*ist->n_xton+n].im*s_mom[index].y_re;

          //c_aj^* * <j| Sz |i> 
          tmpz.re += bs_coeff[jbs*ist->n_xton+n].re*s_mom[index].z_re + bs_coeff[jbs*ist->n_xton+n].im*s_mom[index].z_im;
          tmpz.im += bs_coeff[jbs*ist->n_xton+n].re*s_mom[index].z_im - bs_coeff[jbs*ist->n_xton+n].im*s_mom[index].z_re;

        }

        //multiply by the c_ai coeff
        spinx.re += tmpx.re * bs_coeff[ibs*ist->n_xton+n].re - tmpx.im * bs_coeff[ibs*ist->n_xton+n].im;
        spinx.im += tmpx.im * bs_coeff[ibs*ist->n_xton+n].re + tmpx.re * bs_coeff[ibs*ist->n_xton+n].im;

        spiny.re += tmpy.re * bs_coeff[ibs*ist->n_xton+n].re - tmpy.im * bs_coeff[ibs*ist->n_xton+n].im;
        spiny.im += tmpy.im * bs_coeff[ibs*ist->n_xton+n].re + tmpy.re * bs_coeff[ibs*ist->n_xton+n].im;

        spinz.re += tmpz.re * bs_coeff[ibs*ist->n_xton+n].re - tmpz.im * bs_coeff[ibs*ist->n_xton+n].im;
        spinz.im += tmpz.im * bs_coeff[ibs*ist->n_xton+n].re + tmpz.re * bs_coeff[ibs*ist->n_xton+n].im;
        

        for (b = ist->lumo_idx;b<ist->lumo_idx+ist->n_elecs;b++){
          for (j = 0; j < ist->n_holes; j++) {
            
            indexba = sqr(ist->n_holes)+(a-ist->lumo_idx)*ist->n_elecs+(b-ist->lumo_idx);
            indexji = i*ist->n_holes+j;

            jbs = listibs[(b - ist->lumo_idx) * ist->n_holes + j];
            tmp.re=bs_coeff[ibs*ist->n_xton+n].re*bs_coeff[jbs*ist->n_xton+n].re
                  + bs_coeff[ibs*ist->n_xton+n].im*bs_coeff[jbs*ist->n_xton+n].im;

            tmp.im=bs_coeff[ibs*ist->n_xton+n].im*bs_coeff[jbs*ist->n_xton+n].re
                  - bs_coeff[ibs*ist->n_xton+n].re*bs_coeff[jbs*ist->n_xton+n].im;

            tmpx.re=s_mom[indexba].x_re*s_mom[indexji].x_re-s_mom[indexba].x_im*s_mom[indexji].x_im
             + s_mom[indexba].y_re*s_mom[indexji].y_re-s_mom[indexba].y_im*s_mom[indexji].y_im
             + s_mom[indexba].z_re*s_mom[indexji].z_re-s_mom[indexba].z_im*s_mom[indexji].z_im;

            tmpx.im=s_mom[indexba].x_im*s_mom[indexji].x_re+s_mom[indexba].x_re*s_mom[indexji].x_im
            + s_mom[indexba].y_im*s_mom[indexji].y_re+s_mom[indexba].y_re*s_mom[indexji].y_im
            + s_mom[indexba].z_im*s_mom[indexji].z_re+s_mom[indexba].z_re*s_mom[indexji].z_im;   


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
  // printf("\norbital block\n\n"); fflush(0);
  for (n=0;n<ist->n_xton;n++){
    orbitx.re = orbitx.im = 0;
    orbity.re = orbity.im = 0;
    orbitz.re = orbitz.im = 0;
    orbittot.re = orbittot.im = 0;
    for (a = ist->lumo_idx;a<ist->lumo_idx+ist->n_elecs;a++){
      for (i = 0; i < ist->n_holes; i++) {
        ibs = listibs[(a - ist->lumo_idx) * ist->n_holes + i];
        

        tmpx.re = tmpx.im = 0;
        tmpy.re = tmpy.im = 0;
        tmpz.re = tmpz.im = 0;
        //sum over b
        for (b = ist->lumo_idx;b<ist->lumo_idx+ist->n_elecs;b++){
          jbs = listibs[(b - ist->lumo_idx) * ist->n_holes + i];
          index = sqr(ist->n_holes)+(a-ist->lumo_idx)*ist->n_elecs+(b-ist->lumo_idx);
          
          //c_bi^* * <b| Lx |a>
          tmpx.re += bs_coeff[jbs*ist->n_xton+n].re*l_mom[index].x_re + bs_coeff[jbs*ist->n_xton+n].im*l_mom[index].x_im;
          tmpx.im += bs_coeff[jbs*ist->n_xton+n].re*l_mom[index].x_im - bs_coeff[jbs*ist->n_xton+n].im*l_mom[index].x_re;

          //c_bi^* * <b| Ly |a>
          tmpy.re += bs_coeff[jbs*ist->n_xton+n].re*l_mom[index].y_re + bs_coeff[jbs*ist->n_xton+n].im*l_mom[index].y_im;
          tmpy.im += bs_coeff[jbs*ist->n_xton+n].re*l_mom[index].y_im - bs_coeff[jbs*ist->n_xton+n].im*l_mom[index].y_re;

          //c_bi^* * <b| Lz |a>
          tmpz.re += bs_coeff[jbs*ist->n_xton+n].re*l_mom[index].z_re + bs_coeff[jbs*ist->n_xton+n].im*l_mom[index].z_im;
          tmpz.im += bs_coeff[jbs*ist->n_xton+n].re*l_mom[index].z_im - bs_coeff[jbs*ist->n_xton+n].im*l_mom[index].z_re;
        }
        
        //sum over j
        for (j = 0; j < ist->n_holes; j++) {
          jbs = listibs[(a - ist->lumo_idx) * ist->n_holes + j];
          index = i*ist->n_holes+j;

          //c_aj^* * <j| Lx |i> 
          tmpx.re += bs_coeff[jbs*ist->n_xton+n].re*l_mom[index].x_re + bs_coeff[jbs*ist->n_xton+n].im*l_mom[index].x_im;
          tmpx.im += bs_coeff[jbs*ist->n_xton+n].re*l_mom[index].x_im - bs_coeff[jbs*ist->n_xton+n].im*l_mom[index].x_re;

          //c_aj^* * <j| Ly |i> 
          tmpy.re += bs_coeff[jbs*ist->n_xton+n].re*l_mom[index].y_re + bs_coeff[jbs*ist->n_xton+n].im*l_mom[index].y_im;
          tmpy.im += bs_coeff[jbs*ist->n_xton+n].re*l_mom[index].y_im - bs_coeff[jbs*ist->n_xton+n].im*l_mom[index].y_re;

          //c_aj^* * <j| Lz |i> 
          tmpz.re += bs_coeff[jbs*ist->n_xton+n].re*l_mom[index].z_re + bs_coeff[jbs*ist->n_xton+n].im*l_mom[index].z_im;
          tmpz.im += bs_coeff[jbs*ist->n_xton+n].re*l_mom[index].z_im - bs_coeff[jbs*ist->n_xton+n].im*l_mom[index].z_re;

        }

        //multiply by the c_ai coeff
        orbitx.re += tmpx.re * bs_coeff[ibs*ist->n_xton+n].re - tmpx.im * bs_coeff[ibs*ist->n_xton+n].im;
        orbitx.im += tmpx.im * bs_coeff[ibs*ist->n_xton+n].re + tmpx.re * bs_coeff[ibs*ist->n_xton+n].im;

        orbity.re += tmpy.re * bs_coeff[ibs*ist->n_xton+n].re - tmpy.im * bs_coeff[ibs*ist->n_xton+n].im;
        orbity.im += tmpy.im * bs_coeff[ibs*ist->n_xton+n].re + tmpy.re * bs_coeff[ibs*ist->n_xton+n].im;

        orbitz.re += tmpz.re * bs_coeff[ibs*ist->n_xton+n].re - tmpz.im * bs_coeff[ibs*ist->n_xton+n].im;
        orbitz.im += tmpz.im * bs_coeff[ibs*ist->n_xton+n].re + tmpz.re * bs_coeff[ibs*ist->n_xton+n].im;
        
        //Lsqr part
        for (b = ist->lumo_idx;b<ist->lumo_idx+ist->n_elecs;b++){
          for (j = 0; j < ist->n_holes; j++) {
            
            indexba = sqr(ist->n_holes)+(a-ist->lumo_idx)*ist->n_elecs+(b-ist->lumo_idx);
            indexji = i*ist->n_holes+j;

            jbs = listibs[(b - ist->lumo_idx) * ist->n_holes + j];
            
            //c_{ai}^n * (c_{bj}^n)^*
            tmp.re=bs_coeff[ibs*ist->n_xton+n].re*bs_coeff[jbs*ist->n_xton+n].re
                  + bs_coeff[ibs*ist->n_xton+n].im*bs_coeff[jbs*ist->n_xton+n].im;

            tmp.im=bs_coeff[ibs*ist->n_xton+n].im*bs_coeff[jbs*ist->n_xton+n].re
                  - bs_coeff[ibs*ist->n_xton+n].re*bs_coeff[jbs*ist->n_xton+n].im;


            tmpx.re=l_mom[indexba].x_re*l_mom[indexji].x_re-l_mom[indexba].x_im*l_mom[indexji].x_im
             + l_mom[indexba].y_re*l_mom[indexji].y_re-l_mom[indexba].y_im*l_mom[indexji].y_im
             + l_mom[indexba].z_re*l_mom[indexji].z_re-l_mom[indexba].z_im*l_mom[indexji].z_im;

            tmpx.im=l_mom[indexba].x_im*l_mom[indexji].x_re+l_mom[indexba].x_re*l_mom[indexji].x_im
            + l_mom[indexba].y_im*l_mom[indexji].y_re+l_mom[indexba].y_re*l_mom[indexji].y_im
            + l_mom[indexba].z_im*l_mom[indexji].z_re+l_mom[indexba].z_re*l_mom[indexji].z_im;   

            tmpx.re*=2.0; tmpx.im*=2.0;
            
            
            if (i==j){
              tmpx.re+=l2_mom[indexba].re; tmpx.im+=l2_mom[indexba].im;
            }

            if (a==b){
              tmpx.re+=l2_mom[indexji].re; tmpx.im+=l2_mom[indexji].im;
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
  // printf("\ncouple block\n\n"); fflush(0);
  for (n=0;n<ist->n_xton;n++){
    lstot.re = lstot.im = 0;
    for (a = ist->lumo_idx;a<ist->lumo_idx+ist->n_elecs;a++){
      for (i = 0; i < ist->n_holes; i++) {
        ibs = listibs[(a - ist->lumo_idx) * ist->n_holes + i];
        for (b = ist->lumo_idx;b<ist->lumo_idx+ist->n_elecs;b++){
          for (j = 0; j < ist->n_holes; j++) {
            
            tmpx.re = 0.0;tmpx.im = 0.0;

            indexba = sqr(ist->n_holes)+(a-ist->lumo_idx)*ist->n_elecs+(b-ist->lumo_idx);
            indexji = i*ist->n_holes+j;

            jbs = listibs[(b - ist->lumo_idx) * ist->n_holes + j];
            
            //c_{ai}^n * (c_{bj}^n)^*
            tmp.re=bs_coeff[ibs*ist->n_xton+n].re*bs_coeff[jbs*ist->n_xton+n].re
                  + bs_coeff[ibs*ist->n_xton+n].im*bs_coeff[jbs*ist->n_xton+n].im;

            tmp.im=bs_coeff[ibs*ist->n_xton+n].im*bs_coeff[jbs*ist->n_xton+n].re
                  - bs_coeff[ibs*ist->n_xton+n].re*bs_coeff[jbs*ist->n_xton+n].im;

            
            // <a|L|b>*<j|S|i>
            tmpx.re=l_mom[indexba].x_re*s_mom[indexji].x_re-l_mom[indexba].x_im*s_mom[indexji].x_im
                  + l_mom[indexba].y_re*s_mom[indexji].y_re-l_mom[indexba].y_im*s_mom[indexji].y_im
                  + l_mom[indexba].z_re*s_mom[indexji].z_re-l_mom[indexba].z_im*s_mom[indexji].z_im;

            tmpx.im=l_mom[indexba].x_im*s_mom[indexji].x_re+l_mom[indexba].x_re*s_mom[indexji].x_im
                  + l_mom[indexba].y_im*s_mom[indexji].y_re+l_mom[indexba].y_re*s_mom[indexji].y_im
                  + l_mom[indexba].z_im*s_mom[indexji].z_re+l_mom[indexba].z_re*s_mom[indexji].z_im;

            // <a|S|b>*<j|L|i>
            tmpx.re+=s_mom[indexba].x_re*l_mom[indexji].x_re-s_mom[indexba].x_im*l_mom[indexji].x_im
                   + s_mom[indexba].y_re*l_mom[indexji].y_re-s_mom[indexba].y_im*l_mom[indexji].y_im
                   + s_mom[indexba].z_re*l_mom[indexji].z_re-s_mom[indexba].z_im*l_mom[indexji].z_im;

            tmpx.im+=s_mom[indexba].x_im*l_mom[indexji].x_re+s_mom[indexba].x_re*l_mom[indexji].x_im
                   + s_mom[indexba].y_im*l_mom[indexji].y_re+s_mom[indexba].y_re*l_mom[indexji].y_im
                   + s_mom[indexba].z_im*l_mom[indexji].z_re+s_mom[indexba].z_re*l_mom[indexji].z_im;    

            
            
            //delta ij part
            if (i==j){
              tmpx.re+=LdotS[indexba].re; tmpx.im+=LdotS[indexba].im;
              
            }
            //delta ab part
            if (a==b){
              tmpx.re+=LdotS[indexji].re; tmpx.im+=LdotS[indexji].im;
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



  // printf("\nOMP block 1\n\n"); fflush(0);
  //compute $mat = h \cdot u$
#pragma omp parallel for private(l,j,k,sum)
  for (l = 0; l < ist->n_xton; l++) {
    for (j = 0; j < ist->n_xton; j++) {
      sum.re = sum.im = 0.0;
      for (k = 0; k < ist->n_xton; k++) {
      	sum.re +=  h[l*ist->n_xton+k].re * bs_coeff[k*ist->n_xton+j].re - h[l*ist->n_xton+k].im * bs_coeff[k*ist->n_xton+j].im;
        sum.im +=  h[l*ist->n_xton+k].im * bs_coeff[k*ist->n_xton+j].re + h[l*ist->n_xton+k].re * bs_coeff[k*ist->n_xton+j].im;
      }
      mat[l*ist->n_xton+j].re = sum.re;
      mat[l*ist->n_xton+j].im = sum.im;
    }
  }

  // printf("\nOMP block 2\n\n"); fflush(0);
  //compute $u^\dagger \cdot mat = u^\dagger \cdot h \cdot u$
#pragma omp parallel for private(i,j,l,sum)
  for (i = 0; i < ist->n_xton; i++) {
    for (j = 0; j < ist->n_xton; j++) {
      sum.re = sum.im = 0.0;
      for (l = 0; l < ist->n_xton; l++) {
      	sum.re +=   bs_coeff[l*ist->n_xton+i].re * mat[l*ist->n_xton+j].re + bs_coeff[l*ist->n_xton+i].im * mat[l*ist->n_xton+j].im;
        sum.im +=  -bs_coeff[l*ist->n_xton+i].im * mat[l*ist->n_xton+j].re + bs_coeff[l*ist->n_xton+i].re * mat[l*ist->n_xton+j].im;
      }
      h[i*ist->n_xton+j].re = sum.re;
      h[i*ist->n_xton+j].im = sum.im;
    }
  }
  // printf("\nOMP block 3\n\n"); fflush(0);
#pragma omp parallel for private(l,j,k,sum)
  for (l = 0; l < ist->n_xton; l++) {
    for (j = 0; j < ist->n_xton; j++) {
      for (sum.re=sum.im = 0, k = 0; k < ist->n_xton; k++) {
        sum.re +=  h0mat[l*ist->n_xton+k] * bs_coeff[k*ist->n_xton+j].re;
        sum.im +=  h0mat[l*ist->n_xton+k] * bs_coeff[k*ist->n_xton+j].im;
      }
      mat[l*ist->n_xton+j].re = sum.re;
      mat[l*ist->n_xton+j].im = sum.im;
    }
  }
  // printf("\nOMP block 4\n\n"); fflush(0);
#pragma omp parallel for private(i,j,l,sum)
  for (i = 0; i < ist->n_xton; i++) {
    for (j = 0; j < ist->n_xton; j++) {
      for (sum.re=sum.im = 0, l = 0; l < ist->n_xton; l++) {
        sum.re +=   bs_coeff[l*ist->n_xton+i].re * mat[l*ist->n_xton+j].re + bs_coeff[l*ist->n_xton+i].im * mat[l*ist->n_xton+j].im;
        sum.im +=  -bs_coeff[l*ist->n_xton+i].im * mat[l*ist->n_xton+j].re + bs_coeff[l*ist->n_xton+i].re * mat[l*ist->n_xton+j].im;

      }
      h0mat[i*ist->n_xton+j] = sum.re;
    }
  }

  // printf("\nOMP block 5\n\n"); fflush(0);
#pragma omp parallel for private(l,j,k,sum)
  for (l = 0; l < ist->n_xton; l++) {
    for (j = 0; j < ist->n_xton; j++) {
      sum.re = sum.im = 0.0;
      for (k = 0; k < ist->n_xton; k++) {
        sum.re +=  direct[l*ist->n_xton+k].re * bs_coeff[k*ist->n_xton+j].re - direct[l*ist->n_xton+k].im * bs_coeff[k*ist->n_xton+j].im;
        sum.im +=  direct[l*ist->n_xton+k].im * bs_coeff[k*ist->n_xton+j].re + direct[l*ist->n_xton+k].re * bs_coeff[k*ist->n_xton+j].im;
      }
      mat[l*ist->n_xton+j].re = sum.re;
      mat[l*ist->n_xton+j].im = sum.im;
    }
  }

  // printf("\nOMP block 6\n\n"); fflush(0);
  //compute $u^\dagger \cdot mat = u^\dagger \cdot h \cdot u$
#pragma omp parallel for private(i,j,l,sum)
  for (i = 0; i < ist->n_xton; i++) {
    for (j = 0; j < ist->n_xton; j++) {
      sum.re = sum.im = 0.0;
      for (l = 0; l < ist->n_xton; l++) {
        sum.re +=   bs_coeff[l*ist->n_xton+i].re * mat[l*ist->n_xton+j].re + bs_coeff[l*ist->n_xton+i].im * mat[l*ist->n_xton+j].im;
        sum.im +=  -bs_coeff[l*ist->n_xton+i].im * mat[l*ist->n_xton+j].re + bs_coeff[l*ist->n_xton+i].re * mat[l*ist->n_xton+j].im;
      }
      direct[i*ist->n_xton+j].re = sum.re;
      direct[i*ist->n_xton+j].im = sum.im;
    }
  }
  // printf("\nOMP block 7\n\n"); fflush(0);
  #pragma omp parallel for private(l,j,k,sum)
  for (l = 0; l < ist->n_xton; l++) {
    for (j = 0; j < ist->n_xton; j++) {
      for (sum.re=sum.im = 0, k = 0; k < ist->n_xton; k++) {
        sum.re +=  exchange[l*ist->n_xton+k].re * bs_coeff[k*ist->n_xton+j].re - exchange[l*ist->n_xton+k].im * bs_coeff[k*ist->n_xton+j].im;
        sum.im +=  exchange[l*ist->n_xton+k].im * bs_coeff[k*ist->n_xton+j].re + exchange[l*ist->n_xton+k].re * bs_coeff[k*ist->n_xton+j].im;
      }
      mat[l*ist->n_xton+j].re = sum.re;
      mat[l*ist->n_xton+j].im = sum.im;
    }
  }

  // printf("\nOMP block 8\n\n"); fflush(0);
  //compute $u^\dagger \cdot mat = u^\dagger \cdot h \cdot u$
#pragma omp parallel for private(i,j,l,sum)
  for (i = 0; i < ist->n_xton; i++) {
    for (j = 0; j < ist->n_xton; j++) {
      for (sum.re=sum.im = 0, l = 0; l < ist->n_xton; l++) {
        sum.re +=   bs_coeff[l*ist->n_xton+i].re * mat[l*ist->n_xton+j].re + bs_coeff[l*ist->n_xton+i].im * mat[l*ist->n_xton+j].im;
        sum.im +=  -bs_coeff[l*ist->n_xton+i].im * mat[l*ist->n_xton+j].re + bs_coeff[l*ist->n_xton+i].re * mat[l*ist->n_xton+j].im;
      }
      exchange[i*ist->n_xton+j].re = sum.re;
      exchange[i*ist->n_xton+j].im = sum.im;
    }
  }




  // printf("\nexciton\n\n"); fflush(0);
  // Print out the energies of the excitonic states
  pf = fopen("exciton.dat" , "w");
  fprintf(pf,"#n \t E_n \t <H> \t <H_dir> \t <H_exc> \t <H_0> \t E_B (eV)\n");
  for (i = 0; i < ist->n_xton; i++) {  
    fprintf(pf,"%ld % .12f % .12f % .12f  % .12f  % .12f  % .12f\n", i, xton_ene[i], h[i*ist->n_xton+i].re, direct[i*ist->n_xton+i].re, exchange[i*ist->n_xton+i].re,
     h0mat[i*ist->n_xton+i], (xton_ene[i]-h0mat[i*ist->n_xton+i])*AUTOEV);
  }
  fclose(pf);

  
  /*
  // Electron and hole carrier densities 
  pgrid = (double*)calloc(ist->ngrid,sizeof(double));
  for (jgamma = 0; jgamma < 10; jgamma++) {
	// Electron densities
    for (jgrid = 0; jgrid < ist->ngrid; jgrid++) { 
	  pgrid[jgrid] = 0.0;
	}
    for (ibs = 0, a = ist->lumo_idx; a < ist->lumo_idx+ist->n_elecs; a++) {
      for (i = 0; i < ist->n_holes; i++, ibs++) {
      	for (jbs = 0, b = ist->lumo_idx; b < ist->lumo_idx+ist->n_elecs; b++) {
      	  for (j = 0; j < ist->n_holes; j++, jbs++) {
      	    if (i == j) { 
      	      sum = bs_coeff[jgamma*ist->n_xton+ibs] * bs_coeff[jgamma*ist->n_xton+jbs];
            #pragma omp parallel for private(jgrid)
      	      for (jgrid = 0; jgrid < ist->ngrid; jgrid++)
            		pgrid[jgrid] += sum * psi[a*ist->ngrid+jgrid] * psi[b*ist->ngrid+jgrid];
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
    for (jgrid = 0; jgrid < ist->ngrid; jgrid++) {
	  pgrid[jgrid] = 0.0;
    }
	for (ibs = 0, a = ist->lumo_idx; a < ist->lumo_idx+ist->n_elecs; a++) {
      for (i = 0; i < ist->n_holes; i++, ibs++) {
      	for (jbs = 0, b = ist->lumo_idx; b < ist->lumo_idx+ist->n_elecs; b++) {
      	  for (j = 0; j < ist->n_holes; j++, jbs++) {
      	    if (a == b) { 
      	      sum = bs_coeff[jgamma*ist->n_xton+ibs] * bs_coeff[jgamma*ist->n_xton+jbs];
            #pragma omp parallel for private(jgrid)
      	      for (jgrid = 0; jgrid < ist->ngrid; jgrid++)
            		pgrid[jgrid] += sum * psi[i*ist->ngrid+jgrid] * psi[j*ist->ngrid+jgrid];
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
  if (ist->n_holes < nNonIntStates) nNonIntStates = ist->n_holes;
  if (ist->n_elecs < nNonIntStates) nNonIntStates = ist->n_elecs;
  for (i = 0; i < nNonIntStates; i++) {
  #pragma omp parallel for private(jgrid)
    for (jgrid = 0; jgrid < ist->ngrid; jgrid++) { // electron first
      pgrid[jgrid] = par->dv *  psi[(ist->lumo_idx+i)*ist->ngrid + jgrid] * psi[(ist->lumo_idx+i)*ist->ngrid + jgrid];
    }
    norm_vector(pgrid, par->dv, ist->ngrid);
    sprintf(str, "pe-ni-%d.cub", i);
    print_cube(pgrid, ist, par, str);
    sprintf(str, "pcz%d-ni.dat", i); 
    print_pz_one(pgrid, vz, par, ist, str);
  #pragma omp parallel for private(jgrid)
    for (jgrid = 0; jgrid < ist->ngrid; jgrid++) { // hole second
      pgrid[jgrid] = par->dv *  psi[(ist->lumo_idx-1-i)*ist->ngrid + jgrid] * psi[(ist->lumo_idx-1-i)*ist->ngrid + jgrid];
    }
    norm_vector(pgrid, par->dv, ist->ngrid);
    sprintf(str, "ph-ni-%d.cub", i);
    print_cube(pgrid, ist, par, str);
    sprintf(str, "pvz%d-ni.dat", i); 
    print_pz_one(pgrid, vz, par, ist, str);
  }  

  // Calculate the pegrid and phgrid for the lowest excitonic state 
  // for fixed locations of the other carrier - an interacting result
  if (ist->printFPDensity) {
    print_fixed_qp_density(psi, u, vz, ist, par);
  } 
  
  free(pgrid);
  */
  free(h); free(mat);
  free(listibs);
  
  return;
}

/*****************************************************************************/
