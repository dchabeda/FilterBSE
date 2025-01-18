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
    long jgrid, ispingrid;
    int ispin;
    int    tid; //ispin; 
    long   *listibs;
    double ene, ene1, ene2;
    zomplex *rho;
    
    rho = (zomplex *) calloc(ist->ngrid * ist->nthreads, sizeof(zomplex));
    listibs = (long *) malloc(ist->n_xton * sizeof(long));

    for (ibs = 0, a = ist->lumo_idx; a < ist->lumo_idx+ist->n_elecs; a++) {
        for (i = 0; i < ist->n_holes; i++, ibs++) {
            listibs[(a - ist->lumo_idx) * ist->n_holes + i] = ibs;
            
        }
    }

    omp_set_dynamic(0);
    omp_set_num_threads(ist->nthreads);

    printf("\nComputing screened direct matrix, K^d_(ab,ij)\n");
    
    /*** vabji direct ***/
    // The two electron integrals have a trivial 4-fold permutation symmetry if complex
    // 8-fold if real! We avoid performing additional computational work by recognizing
    // [ij|ab] = [ab|ij] = [ji|ba]^* = [ba|ji]^*
    // using Chemist's notation from Szabo and Ostlund
    // This loop indexing scheme effectively avoids trivial extra computation
    
    // #pragma omp parallel for private(ibs,jbs,ene1,ene2,ene,tid,a,b,i,j)
    for (a = ist->lumo_idx; a < ist->lumo_idx+ist->n_elecs; a++) {
        
        // long jgrid, ispingrid, jgrid_imag, jgrid_real, jgriddn_imag, jgriddn_real, jgridup_imag, jgridup_real;
        // long i_up_real, i_up_imag, i_dn_real, i_dn_imag, j_up_real, j_up_imag, j_dn_real, j_dn_imag;
        // long a_up_real, a_up_imag, a_dn_real, a_dn_imag, b_up_real, b_up_imag, b_dn_real, b_dn_imag;
        // int ispin;
        //loop over electron states b, b <= a for symmetry in (ab) permutation
        for (b = ist->lumo_idx; b < ist->lumo_idx+ist->n_elecs; b++) {
            tid = omp_get_thread_num();
            
            //get joint density \rho_{ab}(r) = \sum_{\sigma} psi_{a}^{*}(r,\sigma) psi_{b}(r,\sigma)
            for (jgrid = 0; jgrid < ist->ngrid; jgrid++) {
                rho[tid*ist->ngrid+jgrid].re = 0.0;
                rho[tid*ist->ngrid+jgrid].im = 0.0;
                for (ispin = 0; ispin < ist->nspin; ispin++){
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
                    ispingrid = jgrid + ist->ngrid*ispin;
    	            rho[tid*ist->ngrid+jgrid].re += psi_qp[a*ist->nspinngrid+ispingrid].re * psi_qp[b*ist->nspinngrid+ispingrid].re
                                                         + psi_qp[a*ist->nspinngrid+ispingrid].im * psi_qp[b*ist->nspinngrid+ispingrid].im;
    	            rho[tid*ist->ngrid+jgrid].im += psi_qp[a*ist->nspinngrid+ispingrid].re * psi_qp[b*ist->nspinngrid+ispingrid].im
                                                         - psi_qp[a*ist->nspinngrid+ispingrid].im * psi_qp[b*ist->nspinngrid+ispingrid].re;
                }
            }
            
            //this should populate the pot_htree array with h_d(r) = \int W(r,r') \rho_{ab}(r') d^3r' via fourier transform
            hartree(&rho[tid*ist->ngrid], pot_screened, &pot_htree[tid*ist->ngrid], ist, planfw[tid], planbw[tid], &fftwpsi[tid*ist->ngrid]);            
            
            //loop over hole states i
            for (i = 0; i < ist->n_holes; i++) {
	            
                //loop over hole states j
                for (j = 0; j < ist->n_holes; j++) {
                    zomplex sum1, sum2, tmp;
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
                        for (ispin = 0; ispin < ist->nspin; ispin++){
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
                    
	                bsmat[ibs * ist->n_xton + jbs].re = sum1.re;
                    bsmat[ibs * ist->n_xton + jbs].im = sum1.im;
                    
                    direct[ibs * ist->n_xton + jbs].re = sum1.re;
                    direct[ibs * ist->n_xton + jbs].im = sum1.im;
	                
                    //if diagonal put the energy difference in the h0mat
                    if (ibs == jbs) 
                        h0mat[ibs*ist->n_xton+jbs] = eval[a] - eval[i];
	                else 
                        h0mat[ibs*ist->n_xton+jbs] = 0.0;
	                
                    
	            }
            }
        }
    }

    pf = fopen("direct.dat" , "w");
    for (ibs = 0; ibs < ist->n_xton; ibs++){
        for (jbs = 0; jbs < ist->n_xton; jbs++){
            fprintf(pf,"%ld %ld %.10f %.10f\n", ibs, jbs, \
            direct[ibs * ist->n_xton + jbs].re, direct[ibs * ist->n_xton + jbs].re);
        }
    }
    fclose(pf);
    printf("  Done computing direct mat\n");

    
    /*** vjbai exchange ***/
    printf("\nComputing bare exchange matrix, K^x_(ai,bj)\n");
    
    //loop over electron states a
    for (a = ist->lumo_idx; a < ist->lumo_idx+ist->n_elecs; a++) {
        // long jgrid, ispingrid, jgrid_imag, jgrid_real, jgriddn_imag, jgriddn_real, jgridup_imag, jgridup_real;
        // long i_up_real, i_up_imag, i_dn_real, i_dn_imag, j_up_real, j_up_imag, j_dn_real, j_dn_imag;
        // long a_up_real, a_up_imag, a_dn_real, a_dn_imag, b_up_real, b_up_imag, b_dn_real, b_dn_imag;
        // int ispin;
        //loop over hole states i
        #pragma omp parallel for private(ibs,jbs,ene1,ene2,ene,tid,b,i,j)
        for (i = 0; i < ist->n_holes; i++) {
            tid = omp_get_thread_num();	
            ene1 = eval[a] - eval[i];
            
            //get joint density \rho_{ai}(r) = \sum_{\sigma} psi_{a}^*(r,\sigma) psi_{i}(r,\sigma)
            for (jgrid = 0; jgrid < ist->ngrid; jgrid++) {
                rho[tid*ist->ngrid+jgrid].re = 0.0;
                rho[tid*ist->ngrid+jgrid].im = 0.0;
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
                for (ispin = 0; ispin < 2; ispin++){
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

                    if(ibs==jbs){
                        bsmat[ibs*ist->n_xton+jbs].im = 0.0;
                    }
	            }
            }
        }
    }

    pf = fopen("exchange.dat" , "w");
    for (ibs = 0; ibs < ist->n_xton; ibs++){
        for (jbs = 0; jbs < ist->n_xton; jbs++){
            fprintf(pf,"%ld %ld %.10f %.10f\n", ibs, jbs, \
            exchange[ibs * ist->n_xton + jbs].re, exchange[ibs * ist->n_xton + jbs].im);
        }
    }
	fclose(pf);
    printf("  Done computing exchange mat\n");
    

    FILE *ppsi;
    ppsi = fopen("bsRE.dat", "w");
    for (i = 0; i < ist->n_xton; i++, fprintf(ppsi,"\n"))
        for (j = 0; j < ist->n_xton; j++)
            fprintf(ppsi,"%.*g ", PR_LEN, bsmat[i*ist->n_xton+j].re);
    fclose(ppsi);

    ppsi = fopen("bsIM.dat", "w");
    for (i = 0; i < ist->n_xton; i++, fprintf(ppsi,"\n"))
        for (j = 0; j < ist->n_xton; j++)
            fprintf(ppsi,"%.*g ", PR_LEN, bsmat[i*ist->n_xton+j].im);
    fclose(ppsi);

    ppsi = fopen("h0.dat", "w");
    for (i = 0; i < ist->n_xton; i++, fprintf(ppsi,"\n"))
        for (j = 0; j < ist->n_xton; j++)
             fprintf(ppsi,"%.*g ", PR_LEN, h0mat[i*ist->n_xton+j]);
    fclose(ppsi);

    free(rho); free(listibs);
	
    return;
}

/***************************************************************************************/
