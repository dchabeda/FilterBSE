/*****************************************************************************/

#include "bethe-salpeter.h"

/*****************************************************************************/

void bethe_salpeter(
  double complex*  direct, 
  double complex*  exchange,
  double complex*  bsmat, 
  double complex*  bs_coeff, 
  double*          h0mat, 
  double*          xton_ene, 
	xyz_st*          s_mom, 
  xyz_st*          l_mom, 
  double complex*  l2_mom, 
  double complex*  ldots, 
  grid_st*         grid, 
  index_st*        ist, 
  par_st*          par,
  flag_st*         flag,
  parallel_st*     parallel
  ){

  /************************************************************/
	/*******************  DECLARE VARIABLES   *******************/
	/************************************************************/

  FILE*           pf;

  unsigned long   a;
  unsigned long   b;
  unsigned long   i;
  unsigned long   j;
  unsigned long   k;
  unsigned long   l;
  unsigned long   ibs;
  unsigned long   jbs;
  unsigned long   idx;

  const long      n_xton = ist->n_xton;

  double complex* mat;
  double complex* h;
  double complex  sum;

  /************************************************************/
	/********************   BUILD BSE "H"    ********************/
	/************************************************************/

  ALLOCATE(&mat, n_xton * n_xton, "mat in bethe-salpeter");
  ALLOCATE(&h,   n_xton * n_xton, "h in bethe-salpeter");
  
  #pragma omp parallel for private(i)
  for (i = 0; i < n_xton*n_xton; i++) {
    h[i] = bs_coeff[i] = h0mat[i] - bsmat[i];
  }

  // Print out the full BSE mat
  pf = fopen("HBSmatRE.dat", "w");
  for (i = 0; i < n_xton; i++, fprintf(pf,"\n")) {
    for (j = 0; j < n_xton; j++) {
      fprintf (pf,"%.*g \t", PR_LEN, creal(h[i*n_xton+j]));
	}
  }
  fclose(pf);

  pf = fopen("HBSmatIM.dat", "w");
  for (i = 0; i < n_xton; i++, fprintf(pf,"\n")) {
    for (j = 0; j < n_xton; j++) {
      fprintf (pf,"%.*g \t", PR_LEN, cimag(h[i*n_xton+j]));
  }
  }
  fclose(pf);

  /************************************************************/
	/********************    DIAG BSE MAT    ********************/
	/************************************************************/

  // Diagonalize the BSE matrix to obtain the coefficients
  
  diag((int)n_xton, ist->nthreads, bs_coeff, xton_ene);

  // Convert from lapack column major to row major order
  
  double complex *tmp_mat;
  ALLOCATE(&tmp_mat, n_xton * n_xton, "tmp_mat");

  memcpy(tmp_mat, bs_coeff, n_xton * n_xton * sizeof(double complex));

  for (ibs = 0; ibs < n_xton; ibs++){
    for (jbs = 0; jbs < n_xton; jbs++){
      bs_coeff[jbs*n_xton + ibs] = conj(tmp_mat[ibs*n_xton + jbs]);
    }
  }
  free(tmp_mat);
  
  
  // Prints the coefficients for the first 100 (or n_xton) lowest energy excitonic states
  long numExcStatesToPrint = 100;
  if (n_xton < numExcStatesToPrint) numExcStatesToPrint = n_xton;

  pf = fopen("BSEcoeffRE.dat", "w");
  for (i = 0; i < n_xton; i++, fprintf(pf,"\n")) {
    for (j = 0; j < n_xton; j++) {  
	  fprintf (pf,"%.*g\t", PR_LEN, creal(bs_coeff[i*n_xton+j]));
    }
  }
  fclose(pf);

  pf = fopen("BSEcoeffIM.dat", "w");
  for (i = 0; i < n_xton; i++, fprintf(pf,"\n")) {
    for (j = 0; j < n_xton; j++) {  
	  fprintf (pf,"%.*g\t", PR_LEN, cimag(bs_coeff[i*n_xton+j]));
    }
  }
  fclose(pf);


  /************************************************************/
	/******************    GROUND XTON ENE    *******************/
	/************************************************************/

  sum = 0.0 + 0.0*I;
  double complex tmp;

  for (i = 0; i < n_xton; i++) {
    for (j = 0; j < n_xton; j++) {
      tmp = bs_coeff[i * n_xton] * h[j*n_xton + i];
      sum += conjmul(bs_coeff[j * n_xton], tmp);
    }
  }
  
  printf("\nGround state exciton has energy = %.5f a.u. | %.5f eV (%.5f Imag)\n", creal(sum), creal(sum)*AUTOEV, cimag(sum));


  //printf("xton_ene[0]=%g\n",xton_ene[0] );
//   long *listibs = (long *) calloc(n_xton, sizeof(long));

//   for (ibs = 0, a = ist->lumo_idx; a < ist->lumo_idx + ist->n_elecs; a++) {
//       for (i = 0; i < ist->n_holes; i++, ibs++) {
//           listibs[(a - ist->lumo_idx) * ist->n_holes + i] = ibs;
//           //printf("a:%ld i:%ld ibs:%ld\n",a,i,ibs);
//       }
//   }

//   //compute spins:
//   FILE* spinpf =fopen("spins.dat", "w");  
//   //printf("Spins:\n");
//   double complex spinx, spiny, spinz, spintot;
//   double complex tmpx, tmpy, tmpz, tmp;
//   long index, indexba, indexji,n;
//   // printf("\nspins block\n\n"); fflush(0);
//   for (n = 0; n < n_xton; n++){
//     spinx.re = spinx.im = 0;
//     spiny.re = spiny.im = 0;
//     spinz.re = spinz.im = 0;
//     spintot.re = spintot.im = 0;
//     for (a = ist->lumo_idx; a < ist->lumo_idx + ist->n_elecs; a++){
//       for (i = 0; i < ist->n_holes; i++) {
//         ibs = listibs[(a - ist->lumo_idx) * ist->n_holes + i];
        

//         tmpx.re = tmpx.im = 0;
//         tmpy.re = tmpy.im = 0;
//         tmpz.re = tmpz.im = 0;
//         //sum over b
//         for (b = ist->lumo_idx; b < ist->lumo_idx + ist->n_elecs; b++){
//           jbs = listibs[(b - ist->lumo_idx) * ist->n_holes + i];
//           index = sqr(ist->n_holes) + (a - ist->lumo_idx) * ist->n_elecs + (b - ist->lumo_idx);
          
//           //c_bi^* * <b| Sx |a>
//           tmpx.re += bs_coeff[jbs*n_xton+n].re*s_mom[index].x_re + bs_coeff[jbs*n_xton+n].im*s_mom[index].x_im;
//           tmpx.im += bs_coeff[jbs*n_xton+n].re*s_mom[index].x_im - bs_coeff[jbs*n_xton+n].im*s_mom[index].x_re;

//           //c_bi^* * <b| Sy |a>
//           tmpy.re += bs_coeff[jbs*n_xton+n].re*s_mom[index].y_re + bs_coeff[jbs*n_xton+n].im*s_mom[index].y_im;
//           tmpy.im += bs_coeff[jbs*n_xton+n].re*s_mom[index].y_im - bs_coeff[jbs*n_xton+n].im*s_mom[index].y_re;

//           //c_bi^* * <b| Sz |a>
//           tmpz.re += bs_coeff[jbs*n_xton+n].re*s_mom[index].z_re + bs_coeff[jbs*n_xton+n].im*s_mom[index].z_im;
//           tmpz.im += bs_coeff[jbs*n_xton+n].re*s_mom[index].z_im - bs_coeff[jbs*n_xton+n].im*s_mom[index].z_re;
//         }
        
//         //sum over j
//         for (j = 0; j < ist->n_holes; j++) {
//           jbs = listibs[(a - ist->lumo_idx) * ist->n_holes + j];
//           index = i*ist->n_holes+j;

//           //c_aj^* * <j| Sx |i> 
//           tmpx.re += bs_coeff[jbs*n_xton+n].re*s_mom[index].x_re + bs_coeff[jbs*n_xton+n].im*s_mom[index].x_im;
//           tmpx.im += bs_coeff[jbs*n_xton+n].re*s_mom[index].x_im - bs_coeff[jbs*n_xton+n].im*s_mom[index].x_re;

//           //c_aj^* * <j| Sy |i> 
//           tmpy.re += bs_coeff[jbs*n_xton+n].re*s_mom[index].y_re + bs_coeff[jbs*n_xton+n].im*s_mom[index].y_im;
//           tmpy.im += bs_coeff[jbs*n_xton+n].re*s_mom[index].y_im - bs_coeff[jbs*n_xton+n].im*s_mom[index].y_re;

//           //c_aj^* * <j| Sz |i> 
//           tmpz.re += bs_coeff[jbs*n_xton+n].re*s_mom[index].z_re + bs_coeff[jbs*n_xton+n].im*s_mom[index].z_im;
//           tmpz.im += bs_coeff[jbs*n_xton+n].re*s_mom[index].z_im - bs_coeff[jbs*n_xton+n].im*s_mom[index].z_re;

//         }

//         //multiply by the c_ai coeff
//         spinx.re += tmpx.re * bs_coeff[ibs*n_xton+n].re - tmpx.im * bs_coeff[ibs*n_xton+n].im;
//         spinx.im += tmpx.im * bs_coeff[ibs*n_xton+n].re + tmpx.re * bs_coeff[ibs*n_xton+n].im;

//         spiny.re += tmpy.re * bs_coeff[ibs*n_xton+n].re - tmpy.im * bs_coeff[ibs*n_xton+n].im;
//         spiny.im += tmpy.im * bs_coeff[ibs*n_xton+n].re + tmpy.re * bs_coeff[ibs*n_xton+n].im;

//         spinz.re += tmpz.re * bs_coeff[ibs*n_xton+n].re - tmpz.im * bs_coeff[ibs*n_xton+n].im;
//         spinz.im += tmpz.im * bs_coeff[ibs*n_xton+n].re + tmpz.re * bs_coeff[ibs*n_xton+n].im;
        

//         for (b = ist->lumo_idx;b<ist->lumo_idx+ist->n_elecs;b++){
//           for (j = 0; j < ist->n_holes; j++) {
            
//             indexba = sqr(ist->n_holes)+(a-ist->lumo_idx)*ist->n_elecs+(b-ist->lumo_idx);
//             indexji = i*ist->n_holes+j;

//             jbs = listibs[(b - ist->lumo_idx) * ist->n_holes + j];
//             tmp.re=bs_coeff[ibs*n_xton+n].re*bs_coeff[jbs*n_xton+n].re
//                   + bs_coeff[ibs*n_xton+n].im*bs_coeff[jbs*n_xton+n].im;

//             tmp.im=bs_coeff[ibs*n_xton+n].im*bs_coeff[jbs*n_xton+n].re
//                   - bs_coeff[ibs*n_xton+n].re*bs_coeff[jbs*n_xton+n].im;

//             tmpx.re=s_mom[indexba].x_re*s_mom[indexji].x_re-s_mom[indexba].x_im*s_mom[indexji].x_im
//              + s_mom[indexba].y_re*s_mom[indexji].y_re-s_mom[indexba].y_im*s_mom[indexji].y_im
//              + s_mom[indexba].z_re*s_mom[indexji].z_re-s_mom[indexba].z_im*s_mom[indexji].z_im;

//             tmpx.im=s_mom[indexba].x_im*s_mom[indexji].x_re+s_mom[indexba].x_re*s_mom[indexji].x_im
//             + s_mom[indexba].y_im*s_mom[indexji].y_re+s_mom[indexba].y_re*s_mom[indexji].y_im
//             + s_mom[indexba].z_im*s_mom[indexji].z_re+s_mom[indexba].z_re*s_mom[indexji].z_im;   


//             spintot.re+=tmp.re*tmpx.re-tmp.im*tmpx.im;
//             spintot.im+=tmp.im*tmpx.re+tmp.re*tmpx.im;


//           }

//         }



//       }

//     }
//     fprintf(spinpf,"%ld\t%-10.5lf\t%-10.5lf\t%-10.5lf\t",n,spinx.re,spiny.re, spinz.re);
//     fprintf(spinpf,"%-10.5lf\t (%-10.5lf)\n",1.5+2.0*spintot.re, 2.0*spintot.im);
//   }
//   fclose(spinpf);


// //compute orbital momentum
//   FILE* orbitpf =fopen("orbital.dat", "w");  
//   //printf("Spins:\n");
//   double complex orbitx, orbity, orbitz, orbittot;
// //  double complex tmpx, tmpy, tmpz, tmp, tmp2;
// //  long index, indexba, indexji,n;
//   // printf("\norbital block\n\n"); fflush(0);
//   for (n=0;n<n_xton;n++){
//     orbitx.re = orbitx.im = 0;
//     orbity.re = orbity.im = 0;
//     orbitz.re = orbitz.im = 0;
//     orbittot.re = orbittot.im = 0;
//     for (a = ist->lumo_idx;a<ist->lumo_idx+ist->n_elecs;a++){
//       for (i = 0; i < ist->n_holes; i++) {
//         ibs = listibs[(a - ist->lumo_idx) * ist->n_holes + i];
        

//         tmpx.re = tmpx.im = 0;
//         tmpy.re = tmpy.im = 0;
//         tmpz.re = tmpz.im = 0;
//         //sum over b
//         for (b = ist->lumo_idx;b<ist->lumo_idx+ist->n_elecs;b++){
//           jbs = listibs[(b - ist->lumo_idx) * ist->n_holes + i];
//           index = sqr(ist->n_holes)+(a-ist->lumo_idx)*ist->n_elecs+(b-ist->lumo_idx);
          
//           //c_bi^* * <b| Lx |a>
//           tmpx.re += bs_coeff[jbs*n_xton+n].re*l_mom[index].x_re + bs_coeff[jbs*n_xton+n].im*l_mom[index].x_im;
//           tmpx.im += bs_coeff[jbs*n_xton+n].re*l_mom[index].x_im - bs_coeff[jbs*n_xton+n].im*l_mom[index].x_re;

//           //c_bi^* * <b| Ly |a>
//           tmpy.re += bs_coeff[jbs*n_xton+n].re*l_mom[index].y_re + bs_coeff[jbs*n_xton+n].im*l_mom[index].y_im;
//           tmpy.im += bs_coeff[jbs*n_xton+n].re*l_mom[index].y_im - bs_coeff[jbs*n_xton+n].im*l_mom[index].y_re;

//           //c_bi^* * <b| Lz |a>
//           tmpz.re += bs_coeff[jbs*n_xton+n].re*l_mom[index].z_re + bs_coeff[jbs*n_xton+n].im*l_mom[index].z_im;
//           tmpz.im += bs_coeff[jbs*n_xton+n].re*l_mom[index].z_im - bs_coeff[jbs*n_xton+n].im*l_mom[index].z_re;
//         }
        
//         //sum over j
//         for (j = 0; j < ist->n_holes; j++) {
//           jbs = listibs[(a - ist->lumo_idx) * ist->n_holes + j];
//           index = i*ist->n_holes+j;

//           //c_aj^* * <j| Lx |i> 
//           tmpx.re += bs_coeff[jbs*n_xton+n].re*l_mom[index].x_re + bs_coeff[jbs*n_xton+n].im*l_mom[index].x_im;
//           tmpx.im += bs_coeff[jbs*n_xton+n].re*l_mom[index].x_im - bs_coeff[jbs*n_xton+n].im*l_mom[index].x_re;

//           //c_aj^* * <j| Ly |i> 
//           tmpy.re += bs_coeff[jbs*n_xton+n].re*l_mom[index].y_re + bs_coeff[jbs*n_xton+n].im*l_mom[index].y_im;
//           tmpy.im += bs_coeff[jbs*n_xton+n].re*l_mom[index].y_im - bs_coeff[jbs*n_xton+n].im*l_mom[index].y_re;

//           //c_aj^* * <j| Lz |i> 
//           tmpz.re += bs_coeff[jbs*n_xton+n].re*l_mom[index].z_re + bs_coeff[jbs*n_xton+n].im*l_mom[index].z_im;
//           tmpz.im += bs_coeff[jbs*n_xton+n].re*l_mom[index].z_im - bs_coeff[jbs*n_xton+n].im*l_mom[index].z_re;

//         }

//         //multiply by the c_ai coeff
//         orbitx.re += tmpx.re * bs_coeff[ibs*n_xton+n].re - tmpx.im * bs_coeff[ibs*n_xton+n].im;
//         orbitx.im += tmpx.im * bs_coeff[ibs*n_xton+n].re + tmpx.re * bs_coeff[ibs*n_xton+n].im;

//         orbity.re += tmpy.re * bs_coeff[ibs*n_xton+n].re - tmpy.im * bs_coeff[ibs*n_xton+n].im;
//         orbity.im += tmpy.im * bs_coeff[ibs*n_xton+n].re + tmpy.re * bs_coeff[ibs*n_xton+n].im;

//         orbitz.re += tmpz.re * bs_coeff[ibs*n_xton+n].re - tmpz.im * bs_coeff[ibs*n_xton+n].im;
//         orbitz.im += tmpz.im * bs_coeff[ibs*n_xton+n].re + tmpz.re * bs_coeff[ibs*n_xton+n].im;
        
//         //Lsqr part
//         for (b = ist->lumo_idx;b<ist->lumo_idx+ist->n_elecs;b++){
//           for (j = 0; j < ist->n_holes; j++) {
            
//             indexba = sqr(ist->n_holes)+(a-ist->lumo_idx)*ist->n_elecs+(b-ist->lumo_idx);
//             indexji = i*ist->n_holes+j;

//             jbs = listibs[(b - ist->lumo_idx) * ist->n_holes + j];
            
//             //c_{ai}^n * (c_{bj}^n)^*
//             tmp.re=bs_coeff[ibs*n_xton+n].re*bs_coeff[jbs*n_xton+n].re
//                   + bs_coeff[ibs*n_xton+n].im*bs_coeff[jbs*n_xton+n].im;

//             tmp.im=bs_coeff[ibs*n_xton+n].im*bs_coeff[jbs*n_xton+n].re
//                   - bs_coeff[ibs*n_xton+n].re*bs_coeff[jbs*n_xton+n].im;


//             tmpx.re=l_mom[indexba].x_re*l_mom[indexji].x_re-l_mom[indexba].x_im*l_mom[indexji].x_im
//              + l_mom[indexba].y_re*l_mom[indexji].y_re-l_mom[indexba].y_im*l_mom[indexji].y_im
//              + l_mom[indexba].z_re*l_mom[indexji].z_re-l_mom[indexba].z_im*l_mom[indexji].z_im;

//             tmpx.im=l_mom[indexba].x_im*l_mom[indexji].x_re+l_mom[indexba].x_re*l_mom[indexji].x_im
//             + l_mom[indexba].y_im*l_mom[indexji].y_re+l_mom[indexba].y_re*l_mom[indexji].y_im
//             + l_mom[indexba].z_im*l_mom[indexji].z_re+l_mom[indexba].z_re*l_mom[indexji].z_im;   

//             tmpx.re*=2.0; tmpx.im*=2.0;
            
            
//             if (i==j){
//               tmpx.re+=l2_mom[indexba].re; tmpx.im+=l2_mom[indexba].im;
//             }

//             if (a==b){
//               tmpx.re+=l2_mom[indexji].re; tmpx.im+=l2_mom[indexji].im;
//             }
            
//             //printf("a:%ld b:%ld i:%ld j:%ld    tmpx: (%lf, %lf)\n", a,b,i,j,tmpx.re, tmpx.im);

//             orbittot.re+=tmp.re*tmpx.re-tmp.im*tmpx.im;
//             orbittot.im+=tmp.im*tmpx.re+tmp.re*tmpx.im;


//           }

//         }



//       }

//     }
//     fprintf(orbitpf,"%ld\t%-10.5lf\t%-10.5lf\t%-10.5lf\t",n,orbitx.re,orbity.re, orbitz.re);
//     fprintf(orbitpf,"%-10.5lf\t (%-10.5lf)\n",orbittot.re, orbittot.im);
//   }
//   fclose(orbitpf);




//   //compute ls momentum
//   FILE* lspf =fopen("couple.dat", "w");  
//   //printf("Spins:\n");
//   double complex lstot;
// //  double complex tmpx, tmpy, tmpz, tmp, tmp2;
// //  long index, indexba, indexji,n;
//   // printf("\ncouple block\n\n"); fflush(0);
//   for (n=0;n<n_xton;n++){
//     lstot.re = lstot.im = 0;
//     for (a = ist->lumo_idx;a<ist->lumo_idx+ist->n_elecs;a++){
//       for (i = 0; i < ist->n_holes; i++) {
//         ibs = listibs[(a - ist->lumo_idx) * ist->n_holes + i];
//         for (b = ist->lumo_idx;b<ist->lumo_idx+ist->n_elecs;b++){
//           for (j = 0; j < ist->n_holes; j++) {
            
//             tmpx.re = 0.0;tmpx.im = 0.0;

//             indexba = sqr(ist->n_holes)+(a-ist->lumo_idx)*ist->n_elecs+(b-ist->lumo_idx);
//             indexji = i*ist->n_holes+j;

//             jbs = listibs[(b - ist->lumo_idx) * ist->n_holes + j];
            
//             //c_{ai}^n * (c_{bj}^n)^*
//             tmp.re=bs_coeff[ibs*n_xton+n].re*bs_coeff[jbs*n_xton+n].re
//                   + bs_coeff[ibs*n_xton+n].im*bs_coeff[jbs*n_xton+n].im;

//             tmp.im=bs_coeff[ibs*n_xton+n].im*bs_coeff[jbs*n_xton+n].re
//                   - bs_coeff[ibs*n_xton+n].re*bs_coeff[jbs*n_xton+n].im;

            
//             // <a|L|b>*<j|S|i>
//             tmpx.re=l_mom[indexba].x_re*s_mom[indexji].x_re-l_mom[indexba].x_im*s_mom[indexji].x_im
//                   + l_mom[indexba].y_re*s_mom[indexji].y_re-l_mom[indexba].y_im*s_mom[indexji].y_im
//                   + l_mom[indexba].z_re*s_mom[indexji].z_re-l_mom[indexba].z_im*s_mom[indexji].z_im;

//             tmpx.im=l_mom[indexba].x_im*s_mom[indexji].x_re+l_mom[indexba].x_re*s_mom[indexji].x_im
//                   + l_mom[indexba].y_im*s_mom[indexji].y_re+l_mom[indexba].y_re*s_mom[indexji].y_im
//                   + l_mom[indexba].z_im*s_mom[indexji].z_re+l_mom[indexba].z_re*s_mom[indexji].z_im;

//             // <a|S|b>*<j|L|i>
//             tmpx.re+=s_mom[indexba].x_re*l_mom[indexji].x_re-s_mom[indexba].x_im*l_mom[indexji].x_im
//                    + s_mom[indexba].y_re*l_mom[indexji].y_re-s_mom[indexba].y_im*l_mom[indexji].y_im
//                    + s_mom[indexba].z_re*l_mom[indexji].z_re-s_mom[indexba].z_im*l_mom[indexji].z_im;

//             tmpx.im+=s_mom[indexba].x_im*l_mom[indexji].x_re+s_mom[indexba].x_re*l_mom[indexji].x_im
//                    + s_mom[indexba].y_im*l_mom[indexji].y_re+s_mom[indexba].y_re*l_mom[indexji].y_im
//                    + s_mom[indexba].z_im*l_mom[indexji].z_re+s_mom[indexba].z_re*l_mom[indexji].z_im;    

            
            
//             //delta ij part
//             if (i==j){
//               tmpx.re+=LdotS[indexba].re; tmpx.im+=LdotS[indexba].im;
              
//             }
//             //delta ab part
//             if (a==b){
//               tmpx.re+=LdotS[indexji].re; tmpx.im+=LdotS[indexji].im;
//             }
            
//             //printf("a:%ld b:%ld i:%ld j:%ld    tmpx: (%lf, %lf)\n", a,b,i,j,tmpx.re, tmpx.im);

//             lstot.re+=tmp.re*tmpx.re-tmp.im*tmpx.im;
//             lstot.im+=tmp.im*tmpx.re+tmp.re*tmpx.im;


//           }

//         }



//       }

//     }
//     fprintf(lspf,"%ld\t%-10.5lf\t (%-10.5lf)\n",n,lstot.re, lstot.im);
//   }
//   fclose(lspf);



//   // printf("\nOMP block 1\n\n"); fflush(0);
//   //compute $mat = h \cdot u$
// #pragma omp parallel for private(l,j,k,sum)
//   for (l = 0; l < n_xton; l++) {
//     for (j = 0; j < n_xton; j++) {
//       sum.re = sum.im = 0.0;
//       for (k = 0; k < n_xton; k++) {
//       	sum.re +=  h[l*n_xton+k].re * bs_coeff[k*n_xton+j].re - h[l*n_xton+k].im * bs_coeff[k*n_xton+j].im;
//         sum.im +=  h[l*n_xton+k].im * bs_coeff[k*n_xton+j].re + h[l*n_xton+k].re * bs_coeff[k*n_xton+j].im;
//       }
//       mat[l*n_xton+j].re = sum.re;
//       mat[l*n_xton+j].im = sum.im;
//     }
//   }

//   // printf("\nOMP block 2\n\n"); fflush(0);
//   //compute $u^\dagger \cdot mat = u^\dagger \cdot h \cdot u$
// #pragma omp parallel for private(i,j,l,sum)
//   for (i = 0; i < n_xton; i++) {
//     for (j = 0; j < n_xton; j++) {
//       sum.re = sum.im = 0.0;
//       for (l = 0; l < n_xton; l++) {
//       	sum.re +=   bs_coeff[l*n_xton+i].re * mat[l*n_xton+j].re + bs_coeff[l*n_xton+i].im * mat[l*n_xton+j].im;
//         sum.im +=  -bs_coeff[l*n_xton+i].im * mat[l*n_xton+j].re + bs_coeff[l*n_xton+i].re * mat[l*n_xton+j].im;
//       }
//       h[i*n_xton+j].re = sum.re;
//       h[i*n_xton+j].im = sum.im;
//     }
//   }
//   // printf("\nOMP block 3\n\n"); fflush(0);
// #pragma omp parallel for private(l,j,k,sum)
//   for (l = 0; l < n_xton; l++) {
//     for (j = 0; j < n_xton; j++) {
//       for (sum.re=sum.im = 0, k = 0; k < n_xton; k++) {
//         sum.re +=  h0mat[l*n_xton+k] * bs_coeff[k*n_xton+j].re;
//         sum.im +=  h0mat[l*n_xton+k] * bs_coeff[k*n_xton+j].im;
//       }
//       mat[l*n_xton+j].re = sum.re;
//       mat[l*n_xton+j].im = sum.im;
//     }
//   }
//   // printf("\nOMP block 4\n\n"); fflush(0);
// #pragma omp parallel for private(i,j,l,sum)
//   for (i = 0; i < n_xton; i++) {
//     for (j = 0; j < n_xton; j++) {
//       for (sum.re=sum.im = 0, l = 0; l < n_xton; l++) {
//         sum.re +=   bs_coeff[l*n_xton+i].re * mat[l*n_xton+j].re + bs_coeff[l*n_xton+i].im * mat[l*n_xton+j].im;
//         sum.im +=  -bs_coeff[l*n_xton+i].im * mat[l*n_xton+j].re + bs_coeff[l*n_xton+i].re * mat[l*n_xton+j].im;

//       }
//       h0mat[i*n_xton+j] = sum.re;
//     }
//   }

//   // printf("\nOMP block 5\n\n"); fflush(0);
// #pragma omp parallel for private(l,j,k,sum)
//   for (l = 0; l < n_xton; l++) {
//     for (j = 0; j < n_xton; j++) {
//       sum.re = sum.im = 0.0;
//       for (k = 0; k < n_xton; k++) {
//         sum.re +=  direct[l*n_xton+k].re * bs_coeff[k*n_xton+j].re - direct[l*n_xton+k].im * bs_coeff[k*n_xton+j].im;
//         sum.im +=  direct[l*n_xton+k].im * bs_coeff[k*n_xton+j].re + direct[l*n_xton+k].re * bs_coeff[k*n_xton+j].im;
//       }
//       mat[l*n_xton+j].re = sum.re;
//       mat[l*n_xton+j].im = sum.im;
//     }
//   }

//   // printf("\nOMP block 6\n\n"); fflush(0);
//   //compute $u^\dagger \cdot mat = u^\dagger \cdot h \cdot u$
// #pragma omp parallel for private(i,j,l,sum)
//   for (i = 0; i < n_xton; i++) {
//     for (j = 0; j < n_xton; j++) {
//       sum.re = sum.im = 0.0;
//       for (l = 0; l < n_xton; l++) {
//         sum.re +=   bs_coeff[l*n_xton+i].re * mat[l*n_xton+j].re + bs_coeff[l*n_xton+i].im * mat[l*n_xton+j].im;
//         sum.im +=  -bs_coeff[l*n_xton+i].im * mat[l*n_xton+j].re + bs_coeff[l*n_xton+i].re * mat[l*n_xton+j].im;
//       }
//       direct[i*n_xton+j].re = sum.re;
//       direct[i*n_xton+j].im = sum.im;
//     }
//   }
//   // printf("\nOMP block 7\n\n"); fflush(0);
//   #pragma omp parallel for private(l,j,k,sum)
//   for (l = 0; l < n_xton; l++) {
//     for (j = 0; j < n_xton; j++) {
//       for (sum.re=sum.im = 0, k = 0; k < n_xton; k++) {
//         sum.re +=  exchange[l*n_xton+k].re * bs_coeff[k*n_xton+j].re - exchange[l*n_xton+k].im * bs_coeff[k*n_xton+j].im;
//         sum.im +=  exchange[l*n_xton+k].im * bs_coeff[k*n_xton+j].re + exchange[l*n_xton+k].re * bs_coeff[k*n_xton+j].im;
//       }
//       mat[l*n_xton+j].re = sum.re;
//       mat[l*n_xton+j].im = sum.im;
//     }
//   }

//   // printf("\nOMP block 8\n\n"); fflush(0);
//   //compute $u^\dagger \cdot mat = u^\dagger \cdot h \cdot u$
// #pragma omp parallel for private(i,j,l,sum)
//   for (i = 0; i < n_xton; i++) {
//     for (j = 0; j < n_xton; j++) {
//       for (sum.re=sum.im = 0, l = 0; l < n_xton; l++) {
//         sum.re +=   bs_coeff[l*n_xton+i].re * mat[l*n_xton+j].re + bs_coeff[l*n_xton+i].im * mat[l*n_xton+j].im;
//         sum.im +=  -bs_coeff[l*n_xton+i].im * mat[l*n_xton+j].re + bs_coeff[l*n_xton+i].re * mat[l*n_xton+j].im;
//       }
//       exchange[i*n_xton+j].re = sum.re;
//       exchange[i*n_xton+j].im = sum.im;
//     }
//   }




//   // printf("\nexciton\n\n"); fflush(0);
//   // Print out the energies of the excitonic states
//   pf = fopen("exciton.dat" , "w");
//   fprintf(pf,"#n \t E_n \t <H> \t <H_dir> \t <H_exc> \t <H_0> \t E_B (eV)\n");
//   for (i = 0; i < n_xton; i++) {  
//     fprintf(pf,"%ld % .12f % .12f % .12f  % .12f  % .12f  % .12f\n", i, xton_ene[i], h[i*n_xton+i].re, direct[i*n_xton+i].re, exchange[i*n_xton+i].re,
//      h0mat[i*n_xton+i], (xton_ene[i]-h0mat[i*n_xton+i])*AUTOEV);
//   }
//   fclose(pf);

  
  free(h); free(mat);
  // free(listibs);
  
  return;
}

/*****************************************************************************/
