#include "fd.h"
#include <float.h>

/***************************************************************************************/

void calc_eh_kernel_cplx(zomplex       *psi_qp, 
                        zomplex       *pot_bare,
                        zomplex       *pot_screened,
                        zomplex       *pot_htree,
                        zomplex       *bsmat,
                        zomplex       *direct,
                        zomplex       *exchange,
                        double        *h0mat,
                        double        *eval,
                        index_st      *ist,
                        par_st        *par,
                        flag_st       *flag,
                        fftw_plan_loc *planfw,
                        fftw_plan_loc *planbw,
                        fftw_complex  *fftwpsi
                        ){
    
    FILE   *pf;  
    long   i, j, a, b, ibs, jbs;
    int    tid; //ispin; 
    long   *listibs;
    double ene, ene1, ene2;
    zomplex *rho;
    
    
    rho = (zomplex *) malloc(ist->ngrid * ist->nthreads * sizeof(zomplex));
    listibs = (long *) malloc(ist->n_xton * sizeof(long));

    for (ibs = 0, a = ist->lumo_idx; a < ist->lumo_idx+ist->n_elecs; a++) {
        for (i = 0; i < ist->n_holes; i++, ibs++) {
            listibs[(a - ist->lumo_idx) * ist->n_holes + i] = ibs;
            // printf("a:%ld i:%ld ibs:%ld\n",a,i,ibs);
        }
    }

    omp_set_dynamic(0);
    omp_set_num_threads(ist->nthreads);

    printf("Computing screened direct matrix, K^d_(ab,ji)\n");
    pf = fopen("direct.dat" , "w");
    /*** vabji direct ***/
    // The two electron integrals have a trivial 4-fold permutation symmetry if complex
    // 8-fold if real! We avoid performing additional computational work by recognizing
    // [ij|ab] = [ab|ij] = [ji|ba]^* = [ba|ji]^*
    // using Chemist's notation from Szabo and Ostlund
    // This loop indexing scheme effectively avoids trivial extra computation
    FILE *pf1, *pf2;
    pf1 = fopen("rho_htree.dat", "w");
    pf2 = fopen("pot_htree.dat", "w");

    #pragma omp parallel for private(ibs,jbs,ene1,ene2,ene,tid,a,b,i,j)
    for (a = ist->lumo_idx; a < ist->lumo_idx+ist->n_elecs; a++) {
        long jgrid, ispingrid, jgrid_imag, jgrid_real, jgriddn_imag, jgriddn_real, jgridup_imag, jgridup_real;
        long i_up_real, i_up_imag, i_dn_real, i_dn_imag, j_up_real, j_up_imag, j_dn_real, j_dn_imag;
        long a_up_real, a_up_imag, a_dn_real, a_dn_imag, b_up_real, b_up_imag, b_dn_real, b_dn_imag;
        int ispin;
        //loop over electron states b, b <= a for symmetry in (ab) permutation
        for (b = ist->lumo_idx; b < ist->lumo_idx+ist->n_elecs; b++) {
            tid = omp_get_thread_num();
            
            //get joint density \rho_{ab}(r) = \sum_{\sigma} psi_{a}^{*}(r,\sigma) psi_{b}(r,\sigma)
            for (jgrid = 0; jgrid < ist->ngrid; jgrid++) {
                // // If the wavefunctions are scalar, only one of these indices is necessary (jgridup_real)
                // // If psi_qp are spinors, then all four indices are necessary
                // jgridup_real = jgrid * ist->complex_idx;
                // a_up_real = a*ist->nspinngrid*ist->complex_idx+jgridup_real;
                // b_up_real = b*ist->nspinngrid*ist->complex_idx+jgridup_real; 
                
                // // compute joint density for real, scalar wavefuntions
                // rho[tid*ist->ngrid+jgrid].re = psi_qp[a_up_real] * psi_qp[b_up_real];
                // // rho[tid*ist->ngrid+jgrid].im = 0.0;

                // // If using spinors, add the down component of the integral (and the remaining imaginary piece from the up spin)
                // if (1 == flag->useSpinors){
                //     // get indices for computing matrix elem with spinor
                //     jgridup_imag = jgrid * ist->complex_idx + 1;
                //     jgriddn_real = jgridup_real + ist->ngrid*ist->complex_idx;
                //     jgriddn_imag = jgridup_imag + ist->ngrid*ist->complex_idx;
                //     a_up_imag = a*ist->nspinngrid*ist->complex_idx+jgridup_imag; b_up_imag = b*ist->nspinngrid*ist->complex_idx+jgridup_imag;
                //     a_dn_real = a*ist->nspinngrid*ist->complex_idx+jgriddn_real; b_dn_real = b*ist->nspinngrid*ist->complex_idx+jgriddn_real;
                //     a_dn_imag = a*ist->nspinngrid*ist->complex_idx+jgriddn_imag; b_dn_imag = b*ist->nspinngrid*ist->complex_idx+jgriddn_imag;

                //     rho[tid*ist->ngrid+jgrid].re += psi_qp[a_up_imag]*psi_qp[b_up_imag] + psi_qp[a_dn_real]*psi_qp[b_dn_real] + psi_qp[a_dn_imag]*psi_qp[b_dn_imag];
                //     rho[tid*ist->ngrid+jgrid].im = psi_qp[a_up_real]*psi_qp[b_up_imag] + psi_qp[a_dn_real]*psi_qp[b_dn_imag]
                //                                         - psi_qp[a_up_imag]*psi_qp[b_up_real] - psi_qp[a_dn_imag]*psi_qp[b_dn_real];
                //}
                for (ispin = 0; ispin < 2; ispin++){
                    ispingrid=jgrid+ist->ngrid*ispin;
    	            rho[tid*ist->ngrid+jgrid].re += psi_qp[a*ist->nspinngrid+ispingrid].re * psi_qp[b*ist->nspinngrid+ispingrid].re
                                                         + psi_qp[a*ist->nspinngrid+ispingrid].im * psi_qp[b*ist->nspinngrid+ispingrid].im;
    	            rho[tid*ist->ngrid+jgrid].im += psi_qp[a*ist->nspinngrid+ispingrid].re * psi_qp[b*ist->nspinngrid+ispingrid].im
                                                         - psi_qp[a*ist->nspinngrid+ispingrid].im * psi_qp[b*ist->nspinngrid+ispingrid].re;
                }
                if (a == 4 && b == 5) {
                    fprintf(pf1, "%ld %.10g %.10g\n", jgrid, rho[tid*ist->ngrid+jgrid].re, rho[tid*ist->ngrid+jgrid].im);
                }
            }
            
            //this should populate the pot_htree array with h_d(r) = \int W(r,r') \rho_{ab}(r') d^3r' via fourier transform
            hartree(&rho[tid*ist->ngrid], pot_screened, &pot_htree[tid*ist->ngrid], ist, planfw[tid], planbw[tid], &fftwpsi[tid*ist->ngrid]);            
            
            if (a == 4 && b == 5) {
                for (jgrid = 0; jgrid < ist->ngrid; jgrid++){
                    fprintf(pf2, "%ld %.10g %.10g\n", jgrid, pot_htree[tid*ist->ngrid+jgrid].re, pot_htree[tid*ist->ngrid+jgrid].im);
                }    
                
            }
            //loop over hole states i
            for (i = 0; i < ist->n_holes; i++) {
	            
                //loop over hole states j
                for (j = 0; j < ist->n_holes; j++) {
                    zomplex sum1, sum2, tmp;
	                // printf("(%ld %ld|%ld %ld)\n", a, b, i, j);
                    //get pair state excitation energy
                    ene1 = eval[a] - eval[i];
	                ene2 = eval[b] - eval[j];
	                ene = ene1 - ene2;
                    
                    //integrate the effective potential to get K^d_{ai,bj}=\int h_d(r) \sum_\sigma psi_{i}(r,\sigma) psi_{j}^{*}(r,\sigma) d^3r
                    sum1.re = sum1.im = tmp.re = tmp.im = 0.0;
					for (jgrid = 0; jgrid < ist->ngrid; jgrid++){
                        // jgridup_real = jgrid * ist->complex_idx;
                        // i_up_real = i*ist->nspinngrid*ist->complex_idx + jgridup_real;
                        // j_up_real = j*ist->nspinngrid*ist->complex_idx + jgridup_real;
                            
				        // tmp.re = psi_qp[i_up_real] * psi_qp[j_up_real];
                        
                        // // If using spinors, add the down component of the integral (and the remaining imaginary piece from the up spin)
                        // if (1 == flag->useSpinors){
                        //     // get indices for computing matrix elem with spinor
                        //     jgridup_imag = jgrid * ist->complex_idx + 1;
                        //     jgriddn_real = jgridup_real + ist->ngrid*ist->complex_idx;
                        //     jgriddn_imag = jgridup_imag + ist->ngrid*ist->complex_idx;
                        //     i_up_imag = i*ist->nspinngrid*ist->complex_idx+jgridup_imag; j_up_imag = j*ist->nspinngrid*ist->complex_idx+jgridup_imag;
                        //     i_dn_real = i*ist->nspinngrid*ist->complex_idx+jgriddn_real; j_dn_real = j*ist->nspinngrid*ist->complex_idx+jgriddn_real;
                        //     i_dn_imag = i*ist->nspinngrid*ist->complex_idx+jgriddn_imag; j_dn_imag = j*ist->nspinngrid*ist->complex_idx+jgriddn_imag;

                        //     tmp.re = psi_qp[i_up_real] * psi_qp[j_up_real] + psi_qp[i_up_imag]*psi_qp[j_up_imag] 
                        //              + psi_qp[i_dn_real]*psi_qp[j_dn_real] + psi_qp[i_dn_imag]*psi_qp[j_dn_imag];
                        //     tmp.im = (psi_qp[j_up_real]*psi_qp[i_up_imag] + psi_qp[j_dn_real]*psi_qp[i_dn_imag]
                        //              - psi_qp[j_up_imag]*psi_qp[i_up_real] - psi_qp[j_dn_imag]*psi_qp[i_dn_real]);
                        // }
                        for (ispin = 0; ispin<2; ispin++){
                            ispingrid=jgrid+ist->ngrid*ispin;
				            tmp.re = (psi_qp[j*ist->nspinngrid + ispingrid].re * psi_qp[i*ist->nspinngrid + ispingrid].re
                                     +psi_qp[j*ist->nspinngrid + ispingrid].im * psi_qp[i*ist->nspinngrid + ispingrid].im);
                            tmp.im = (psi_qp[j*ist->nspinngrid + ispingrid].re * psi_qp[i*ist->nspinngrid + ispingrid].im
                                     -psi_qp[j*ist->nspinngrid + ispingrid].im * psi_qp[i*ist->nspinngrid + ispingrid].re);
                            
                        
                            sum1.re += pot_htree[tid*ist->ngrid + jgrid].re * tmp.re 
                                    -  pot_htree[tid*ist->ngrid + jgrid].im * tmp.im;
                            
                            sum1.im += pot_htree[tid*ist->ngrid + jgrid].re * tmp.im
                                    +  pot_htree[tid*ist->ngrid + jgrid].im * tmp.re;
                        }
                    }
	                sum1.re *= par->dv;
                    sum1.im *= par->dv;

                    //get the matrix indicies for {ai,bj} and set bsmat
	                ibs = listibs[(a - ist->lumo_idx)*ist->n_holes + i];
	                jbs = listibs[(b - ist->lumo_idx)*ist->n_holes + j];
                    //printf("Index bsmas: ibs:%ld jbs:%ld index:%ld\n",ibs, jbs,ibs * ist->n_xton + jbs);
	                bsmat[ibs * ist->n_xton + jbs].re = sum1.re;
                    bsmat[ibs * ist->n_xton + jbs].im = sum1.im;
                    
                    direct[ibs * ist->n_xton + jbs].re = sum1.re;
                    direct[ibs * ist->n_xton + jbs].im = sum1.im;
	                
                    //if diagonal put the energy difference in the h0mat
                    if (ibs == jbs) 
                        h0mat[ibs*ist->n_xton+jbs] = eval[a] - eval[i];
	                else 
                        h0mat[ibs*ist->n_xton+jbs] = 0.0;
	                
                    fprintf(pf,"%ld %ld %ld %ld %ld %ld %.10f %.10f %.10f %.10f\n",a,i,b,j,
		                     listibs[(a-ist->lumo_idx)*ist->n_holes+i],
		                     listibs[(b-ist->lumo_idx)*ist->n_holes+j],
		                     ene1, ene2, sum1.re,sum1.im);
	                //fflush(0);
	            }
            }
        }
    }
    //fprintf(pf,"#############################################\n");
    fclose(pf1); fclose(pf2);
    fclose(pf);
    /*** vjbai exchange ***/
    printf("Computing bare exchange matrix, K^d_(jb,ai)\n");
    pf = fopen("exchange.dat" , "w");
    //loop over electron states a
    #pragma omp parallel for private(ibs,jbs,ene1,ene2,ene,tid,a,b,i,j)
    for (a = ist->lumo_idx; a < ist->lumo_idx+ist->n_elecs; a++) {
        long jgrid, ispingrid, jgrid_imag, jgrid_real, jgriddn_imag, jgriddn_real, jgridup_imag, jgridup_real;
        long i_up_real, i_up_imag, i_dn_real, i_dn_imag, j_up_real, j_up_imag, j_dn_real, j_dn_imag;
        long a_up_real, a_up_imag, a_dn_real, a_dn_imag, b_up_real, b_up_imag, b_dn_real, b_dn_imag;
        int ispin;
        //loop over hole states i
        for (i = 0; i < ist->n_holes; i++) {
            tid = omp_get_thread_num();	
            ene1 = eval[a] - eval[i];
            
            //get joint density \rho_{ai}(r) = \sum_{\sigma} psi_{a}^*(r,\sigma) psi_{i}(r,\sigma)
            for (jgrid = 0; jgrid < ist->ngrid; jgrid++) {
                // If the wavefunctions are scalar, only one of these indices is necessary (jgridup_real)
                // If psi_qp are spinors, then all four indices are necessary
                // jgridup_real = jgrid * ist->complex_idx;
                // a_up_real = a*ist->nspinngrid*ist->complex_idx+jgridup_real;
                // i_up_real = i*ist->nspinngrid*ist->complex_idx+jgridup_real; 
                
                // // compute joint density for real, scalar wavefuntions
                // rho[tid*ist->ngrid+jgrid].re = psi_qp[a_up_real] * psi_qp[i_up_real];
                // // rho[tid*ist->ngrid+jgrid].im = 0.0;

                // // If using spinors, add the down component of the integral (and the remaining imaginary piece from the up spin)
                // if (1 == flag->useSpinors){
                //     // get indices for computing matrix elem with spinor
                //     jgridup_imag = jgrid * ist->complex_idx + 1;
                //     jgriddn_real = jgridup_real + ist->ngrid*ist->complex_idx;
                //     jgriddn_imag = jgridup_imag + ist->ngrid*ist->complex_idx;
                //     a_up_imag = a*ist->nspinngrid*ist->complex_idx+jgridup_imag; i_up_imag = i*ist->nspinngrid*ist->complex_idx+jgridup_imag;
                //     a_dn_real = a*ist->nspinngrid*ist->complex_idx+jgriddn_real; i_dn_real = i*ist->nspinngrid*ist->complex_idx+jgriddn_real;
                //     a_dn_imag = a*ist->nspinngrid*ist->complex_idx+jgriddn_imag; i_dn_imag = i*ist->nspinngrid*ist->complex_idx+jgriddn_imag;

                //     rho[tid*ist->ngrid+jgrid].re += psi_qp[a_up_imag]*psi_qp[i_up_imag] + psi_qp[a_dn_real]*psi_qp[i_dn_real] + psi_qp[a_dn_imag]*psi_qp[i_dn_imag];
                //     rho[tid*ist->ngrid+jgrid].im = psi_qp[i_up_real]*psi_qp[i_up_imag] + psi_qp[i_dn_real]*psi_qp[i_dn_imag]
                //                                         - psi_qp[a_up_imag]*psi_qp[i_up_real] - psi_qp[a_dn_imag]*psi_qp[i_dn_real];
                // }
                for (ispin = 0; ispin<2; ispin++){
                    ispingrid=jgrid+ist->ngrid*ispin;
                    rho[tid*ist->ngrid+jgrid].re += psi_qp[a*ist->nspinngrid+ispingrid].re * psi_qp[i*ist->nspinngrid+ispingrid].re
                                                 + psi_qp[a*ist->nspinngrid+ispingrid].im * psi_qp[i*ist->nspinngrid+ispingrid].im;
                    rho[tid*ist->ngrid+jgrid].im += psi_qp[a*ist->nspinngrid+ispingrid].re * psi_qp[i*ist->nspinngrid+ispingrid].im
                                                 - psi_qp[a*ist->nspinngrid+ispingrid].im * psi_qp[i*ist->nspinngrid+ispingrid].re;
                } 
            }
            //this should populate the pot_htree array with h_x(r) = \int v(r,r') \rho_{ai}(r') d^3r' but i don't understand how
            hartree(&rho[tid*ist->ngrid], pot_bare, &pot_htree[tid*ist->ngrid], ist, planfw[tid], planbw[tid], &fftwpsi[tid*ist->ngrid]);
            
            //loop over electron states b
            for (b = ist->lumo_idx; b < ist->lumo_idx+ist->n_elecs; b++) {

                //loop over hole states j
	            for (j = 0; j < ist->n_holes; j++) {
                    zomplex sum1, sum2, tmp;
	                // printf("(%ld %ld|%ld %ld)\n", a, b, i, j);
                    
	                ene2 = eval[b] - eval[j];
	                ene = ene1 - ene2;
	                //
                    //integrate the effective potential to get K^x_{ai,bj}=\int h_x(r) \sum_\sigma psi_{b}(r,\sigma) psi_{j}^{*}(r,\sigma) d^3r
                    sum2.re = sum2.im = tmp.re = tmp.im = 0.0;
                    for (jgrid = 0; jgrid < ist->ngrid; jgrid++){
                        // jgridup_real = jgrid * ist->complex_idx;
                        // j_up_real = j*ist->nspinngrid*ist->complex_idx + jgridup_real;
                        // b_up_real = b*ist->nspinngrid*ist->complex_idx + jgridup_real;
                            
				        // tmp.re = psi_qp[j_up_real] * psi_qp[b_up_real];
                        
                        // // If using spinors, add the down component of the integral (and the remaining imaginary piece from the up spin)
                        // if (1 == flag->useSpinors){
                        //     // get indices for computing matrix elem with spinor
                        //     jgridup_imag = jgrid * ist->complex_idx + 1;
                        //     jgriddn_real = jgridup_real + ist->ngrid*ist->complex_idx;
                        //     jgriddn_imag = jgridup_imag + ist->ngrid*ist->complex_idx;
                        //     j_up_imag = j*ist->nspinngrid*ist->complex_idx+jgridup_imag; b_up_imag = b*ist->nspinngrid*ist->complex_idx+jgridup_imag;
                        //     j_dn_real = j*ist->nspinngrid*ist->complex_idx+jgriddn_real; b_dn_real = b*ist->nspinngrid*ist->complex_idx+jgriddn_real;
                        //     j_dn_imag = j*ist->nspinngrid*ist->complex_idx+jgriddn_imag; b_dn_imag = b*ist->nspinngrid*ist->complex_idx+jgriddn_imag;

                        //     tmp.re += psi_qp[j_up_imag]*psi_qp[b_up_imag] + psi_qp[j_dn_real]*psi_qp[b_dn_real] + psi_qp[j_dn_imag]*psi_qp[b_dn_imag];
                        //     tmp.im = (psi_qp[j_up_real]*psi_qp[b_up_imag] + psi_qp[j_dn_real]*psi_qp[b_dn_imag]
                        //              - psi_qp[j_up_imag]*psi_qp[b_up_real] - psi_qp[j_dn_imag]*psi_qp[b_dn_real]);
                        // }
                                     
                        // sum2.re += pot_htree[tid*ist->ngrid + jgrid].re * tmp.re 
                        //         - pot_htree[tid*ist->ngrid + jgrid].im * tmp.im;
                        
                        // sum2.im += pot_htree[tid*ist->ngrid + jgrid].re * tmp.im
                        //         +  pot_htree[tid*ist->ngrid + jgrid].im * tmp.re;
                        for (ispin = 0; ispin<2; ispin++){
                            ispingrid=jgrid+ist->ngrid*ispin;
                            tmp.re = (psi_qp[j*ist->nspinngrid + ispingrid].re * psi_qp[b*ist->nspinngrid + ispingrid].re
                                     +psi_qp[j*ist->nspinngrid + ispingrid].im * psi_qp[b*ist->nspinngrid + ispingrid].im);
                            tmp.im = (psi_qp[j*ist->nspinngrid + ispingrid].re * psi_qp[b*ist->nspinngrid + ispingrid].im
                                     -psi_qp[j*ist->nspinngrid + ispingrid].im * psi_qp[b*ist->nspinngrid + ispingrid].re);

                            sum2.re += pot_htree[tid*ist->ngrid + jgrid].re * tmp.re
                                    -  pot_htree[tid*ist->ngrid + jgrid].im * tmp.im;
                            
                            sum2.im += pot_htree[tid*ist->ngrid + jgrid].re * tmp.im
                                    +  pot_htree[tid*ist->ngrid + jgrid].im * tmp.re;                                        
                        }
                    }
                    
                    sum2.re *= par->dv;
                    sum2.im *= par->dv;

                    ibs = listibs[(a-ist->lumo_idx)*ist->n_holes+i];
	                jbs = listibs[(b-ist->lumo_idx)*ist->n_holes+j];
                    //NOTE: scalar version has a 2 as to calc for bright only. Don't want for full matrix. Took out 11/10 --DW
					bsmat[ibs*ist->n_xton+jbs].re -=  sum2.re;
                    bsmat[ibs*ist->n_xton+jbs].im -=  sum2.im;

                    exchange[ibs*ist->n_xton+jbs].re = -1.0* sum2.re;
                    exchange[ibs*ist->n_xton+jbs].im = -1.0* sum2.im;

                    if(ibs==jbs){bsmat[ibs*ist->n_xton+jbs].im = 0.0;}
	                fprintf(pf,"%ld %ld %ld %ld %ld %ld %.10f %.10f %.10f %.10f\n",
                            a,i,b,j,ibs,jbs,
                            ene1, ene2, sum2.re, sum2.im);
	                // fprintf(pf,"%ld %ld %ld %ld %ld %ld %.*g %.*g %.*g %.*g\n",
                    //         a,i,b,j,ibs,jbs,
                    //         DBL_DIG, ene1, DBL_DIG, ene2, DBL_DIG, sum2.re, DBL_DIG, sum2.im);
	                // fflush(0);
	            }
            }
        }
    }
    
	fclose(pf);
    
	free(rho); free(listibs);
	
    return;
}

/***************************************************************************************/
