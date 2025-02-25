#include "fd.h"
#include <float.h>
#include <mpi.h>

/***************************************************************************************/

void calc_eh_kernel_cplx(
	zomplex       *psi_qp, 
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
	fftw_complex  *fftwpsi,
	parallel_st *parallel
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

	
	/*** vabji direct ***/
	// We avoid performing additional computational work by recognizing that
	// the two electron integrals have a 2-fold permutation symmetry if complex
	// [ij|ab] = [ji|ba]^*
	// 4-fold if real! 
	// [ij|ab] = [ji|ab] = [ji|ba] = [ij|ba]
	// using Chemist's notation from Szabo and Ostlund
	// This loop indexing scheme effectively avoids trivial extra computation
	
	if (parallel->mpi_rank == 0)
	{

	// Allocate storage for cached potentials (assuming n_elecs is manageable)
	zomplex **pot_cache = malloc(ist->n_elecs * ist->n_elecs * sizeof(zomplex *));
	for (i = 0; i < sqr(ist->n_elecs); i++) {
		pot_cache[i] = calloc(ist->ngrid, sizeof(zomplex));
	}
	// Cache for the integrals
	zomplex Kd_cache[ist->n_elecs][ist->n_elecs][ist->n_holes][ist->n_holes];
	for (a = 0; a < ist->n_elecs; a++){
	for (b = 0; b < ist->n_elecs; b++){
	for (i = 0; i < ist->n_holes; i++){
	for (j = 0; j < ist->n_holes; j++){
		Kd_cache[a][b][i][j].re = 0.0;
		Kd_cache[a][b][i][j].im = 0.0;
	}
	}
	}
	}

	printf("\nComputing screened direct matrix, K^d_(ab,ij) on rank %d\n", parallel->mpi_rank);
	pf = fopen("direct.dat" , "w");

	for (a = ist->lumo_idx; a < ist->lumo_idx+ist->n_elecs; a++) {
		//loop over electron states b, b <= a for symmetry in (ab) permutation
		#pragma omp parallel for private(ibs,jbs,ene1,ene2,ene,tid,b,i,j)
		for (b = ist->lumo_idx; b < ist->lumo_idx+ist->n_elecs; b++) {
			tid = omp_get_thread_num();

			// Check if pot_htree[a][b] was already computed
			if (pot_cache[(a - ist->lumo_idx)*ist->n_elecs + (b - ist->lumo_idx)][0].re == 0.0 &&
					pot_cache[(a - ist->lumo_idx)*ist->n_elecs + (b - ist->lumo_idx)][0].im == 0.0) 
			{
			
				//get joint density \rho_{ab}(r) = \sum_{\sigma} psi_{a}^{*}(r,\sigma) psi_{b}(r,\sigma)
				for (jgrid = 0; jgrid < ist->ngrid; jgrid++) {
					rho[tid*ist->ngrid+jgrid].re = 0.0;
					rho[tid*ist->ngrid+jgrid].im = 0.0;

					for (ispin = 0; ispin < ist->nspin; ispin++){
						ispingrid = jgrid + ist->ngrid*ispin;

						rho[tid*ist->ngrid+jgrid].re += \
							psi_qp[a*ist->nspinngrid+ispingrid].re * psi_qp[b*ist->nspinngrid+ispingrid].re
							+ psi_qp[a*ist->nspinngrid+ispingrid].im * psi_qp[b*ist->nspinngrid+ispingrid].im;

						rho[tid*ist->ngrid+jgrid].im += \
							psi_qp[a*ist->nspinngrid+ispingrid].re * psi_qp[b*ist->nspinngrid+ispingrid].im
							- psi_qp[a*ist->nspinngrid+ispingrid].im * psi_qp[b*ist->nspinngrid+ispingrid].re;
					}
				}
			
				//this should populate the pot_htree array with h_d(r) = \int W(r,r') \rho_{ab}(r') d^3r' via fourier transform
				hartree(&rho[tid*ist->ngrid], pot_screened, &pot_htree[tid*ist->ngrid], ist, planfw[tid], planbw[tid], &fftwpsi[tid*ist->ngrid]);   
				
				// Cache the result for ab
				memcpy(pot_cache[(a - ist->lumo_idx)*ist->n_elecs + (b - ist->lumo_idx)],
					&pot_htree[tid * ist->ngrid], ist->ngrid * sizeof(zomplex));
				// Cache the result for (ba)*
				for (jgrid = 0; jgrid < ist->ngrid; jgrid++){
					pot_cache[(b - ist->lumo_idx)*ist->n_elecs + (a - ist->lumo_idx)][jgrid].re = pot_htree[tid*ist->ngrid + jgrid].re;
					pot_cache[(b - ist->lumo_idx)*ist->n_elecs + (a - ist->lumo_idx)][jgrid].re = - pot_htree[tid*ist->ngrid + jgrid].im;
				}
			} else {
				// Reuse the cached value
				// printf("Reusing cached pot_htree for a = %ld b = %ld\n", (a - ist->lumo_idx), (b - ist->lumo_idx) );
				memcpy(&pot_htree[tid * ist->ngrid],
					pot_cache[(a - ist->lumo_idx)*ist->n_elecs + (b - ist->lumo_idx)], ist->ngrid * sizeof(zomplex));
			}

			//loop over hole states i
			for (i = 0; i < ist->n_holes; i++) {
				//loop over hole states j
				for (j = 0; j < ist->n_holes; j++) {
					//get the matrix indicies for {ai,bj}
					ibs = listibs[(a - ist->lumo_idx)*ist->n_holes + i];
					jbs = listibs[(b - ist->lumo_idx)*ist->n_holes + j];

					zomplex sum1, tmp;
					//get pair state excitation energy
					ene1 = eval[a] - eval[i];
					ene2 = eval[b] - eval[j];
					ene = ene1 - ene2;

					//if diagonal put the energy difference in the h0mat
					if (ibs == jbs) {
            h0mat[ibs*ist->n_xton+jbs] = eval[a] - eval[i];
          }
					else {
            h0mat[ibs*ist->n_xton+jbs] = 0.0;
          }

					// If an equivalent integral was computed before, reuse it
					if ( (Kd_cache[(a-ist->lumo_idx)][b-ist->lumo_idx][i][j].re != 0.0) && (Kd_cache[a-ist->lumo_idx][b-ist->lumo_idx][i][j].im != 0.0) ) {
						// printf("Kd_cache for a = %ld b = %ld i = %ld j = %ld: %lg %lg\n", a, b, i, j, Kd_cache[a-ist->lumo_idx][b-ist->lumo_idx][i][j].re, Kd_cache[a-ist->lumo_idx][b-ist->lumo_idx][i][j].im);
						
            bsmat[ibs * ist->n_xton + jbs].re += Kd_cache[a-ist->lumo_idx][b-ist->lumo_idx][i][j].re;
            bsmat[ibs * ist->n_xton + jbs].im += Kd_cache[a-ist->lumo_idx][b-ist->lumo_idx][i][j].im;
						direct[ibs * ist->n_xton + jbs] = Kd_cache[a-ist->lumo_idx][b-ist->lumo_idx][i][j];
            continue;
        	}
					
						
					//integrate the effective potential to get K^d_{ai,bj}=\int h_d(r) \sum_\sigma psi_{i}(r,\sigma) psi_{j}^{*}(r,\sigma) d^3r
					sum1.re = sum1.im = tmp.re = tmp.im = 0.0;
					for (jgrid = 0; jgrid < ist->ngrid; jgrid++){
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
						
					bsmat[ibs * ist->n_xton + jbs].re = sum1.re;
					bsmat[ibs * ist->n_xton + jbs].im = sum1.im;
					
					direct[ibs * ist->n_xton + jbs].re = sum1.re;
					direct[ibs * ist->n_xton + jbs].im = sum1.im;
				
					// Cache this integral to avoid recomputing
					Kd_cache[a-ist->lumo_idx][b-ist->lumo_idx][i][j].re = sum1.re;
					Kd_cache[a-ist->lumo_idx][b-ist->lumo_idx][i][j].im = sum1.im;
					// Use the two-fold symmetry
					Kd_cache[b-ist->lumo_idx][a-ist->lumo_idx][j][i].re = sum1.re;
					Kd_cache[b-ist->lumo_idx][a-ist->lumo_idx][j][i].im = - sum1.im;
				}
			}
		}

		for (long aux = a; aux < a+1; aux++){
			for (b = 0; b < ist->n_elecs; b++){
			for (i = 0; i < ist->n_holes; i++){
			for (j = 0; j < ist->n_holes; j++){
				ibs = listibs[(aux - ist->lumo_idx)*ist->n_holes + i];
				jbs = listibs[b*ist->n_holes + j];
				fprintf(pf,"%ld %ld %ld %ld %ld %ld %.10f %.10f\n", aux, b, i, j, ibs, jbs, \
				direct[ibs * ist->n_xton + jbs].re, direct[ibs * ist->n_xton + jbs].re);
			}
			}
			}
		}
		fflush(0);
	}

	fclose(pf);
	printf("  Done computing direct mat\n"); 
	fflush(0);

	for (i = 0; i < sqr(ist->n_elecs) ; i++) {
		free(pot_cache[i]);
	}
	free(pot_cache);

	} // end of mpi rank 0
	

	/*******************************************************************/
	/*******************************************************************/
	/*******************************************************************/

	if ( (parallel->mpi_rank == 1 ) || (parallel->mpi_size == 1) ) //
	{
	// Cache for the integrals
	zomplex Kx_cache[ist->n_elecs][ist->n_holes][ist->n_elecs][ist->n_holes];
	for (a = 0; a < ist->n_elecs; a++){
	for (i = 0; i < ist->n_holes; i++){
	for (b = 0; b < ist->n_elecs; b++){
	for (j = 0; j < ist->n_holes; j++){
		Kx_cache[a][i][b][j].re = 0.0;
		Kx_cache[a][i][b][j].im = 0.0;
	}
	}
	}
	}

	/*** vjbai exchange ***/
	printf("\nComputing bare exchange matrix, K^x_(ai,bj) on rank %d\n", parallel->mpi_rank);
	pf = fopen("exchange.dat" , "w");

	//loop over electron states a
	for (a = ist->lumo_idx; a < ist->lumo_idx+ist->n_elecs; a++) {
		//loop over hole states i
		#pragma omp parallel for private(ibs,jbs,ene1,ene2,ene,tid,b,i,j)
		for (i = 0; i < ist->n_holes; i++) {
				tid = omp_get_thread_num();	
				ene1 = eval[a] - eval[i];
				
				
				//get joint density \rho_{ai}(r) = \sum_{\sigma} psi_{a}^*(r,\sigma) psi_{i}(r,\sigma)
				for (jgrid = 0; jgrid < ist->ngrid; jgrid++) {
					rho[tid*ist->ngrid+jgrid].re = 0.0;
					rho[tid*ist->ngrid+jgrid].im = 0.0;
					
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
					ibs = listibs[(a-ist->lumo_idx) * ist->n_holes + i];
					jbs = listibs[(b-ist->lumo_idx) * ist->n_holes + j];
					
					if(ibs==jbs){
						bsmat[ibs*ist->n_xton+jbs].im = 0.0;
					}

					// If an equivalent integral was computed before, reuse it
					if ( (Kx_cache[(a-ist->lumo_idx)][i][b-ist->lumo_idx][j].re != 0.0) && (Kx_cache[a-ist->lumo_idx][i][b-ist->lumo_idx][j].im != 0.0) ) {
						// printf("Kx_cache for a = %ld b = %ld i = %ld j = %ld: %lg %lg\n", a, b, i, j, Kx_cache[a-ist->lumo_idx][i][b-ist->lumo_idx][j].re, Kx_cache[a-ist->lumo_idx][i][b-ist->lumo_idx][j].im);
						
            exchange[ibs * ist->n_xton + jbs] = Kx_cache[a-ist->lumo_idx][i][b-ist->lumo_idx][j];

            bsmat[ibs * ist->n_xton + jbs].re += Kx_cache[a-ist->lumo_idx][i][b-ist->lumo_idx][j].re;
            if(ibs==jbs){
              bsmat[ibs*ist->n_xton+jbs].im = 0.0;
            } else{
            bsmat[ibs * ist->n_xton + jbs].im += Kx_cache[a-ist->lumo_idx][i][b-ist->lumo_idx][j].im;
            }
            continue;
        	}

					zomplex sum2, tmp;
					// printf("(%ld %ld|%ld %ld)\n", a, b, i, j);
						
					ene2 = eval[b] - eval[j];
					ene = ene1 - ene2;
					//
					//integrate the effective potential to get K^x_{ai,bj}=\int h_x(r) \sum_\sigma psi_{b}(r,\sigma) psi_{j}^{*}(r,\sigma) d^3r
					sum2.re = sum2.im = tmp.re = tmp.im = 0.0;
					for (jgrid = 0; jgrid < ist->ngrid; jgrid++){
						for (ispin = 0; ispin < 2; ispin++){
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
					
					//NOTE: scalar version has a 2 as to calc for bright only. Don't want for full matrix. Took out 11/10 --DW
					bsmat[ibs*ist->n_xton+jbs].re -=  sum2.re;
					bsmat[ibs*ist->n_xton+jbs].im -=  sum2.im;

					exchange[ibs*ist->n_xton+jbs].re = -1.0* sum2.re;
					exchange[ibs*ist->n_xton+jbs].im = -1.0* sum2.im;
					
					// Cache this integral to avoid recomputing
					Kx_cache[a-ist->lumo_idx][i][b-ist->lumo_idx][j].re = -1.0 * sum2.re;
					Kx_cache[a-ist->lumo_idx][i][b-ist->lumo_idx][j].im = -1.0 * sum2.im;
					// Utilize the two-fold symmetry
					Kx_cache[b-ist->lumo_idx][j][a-ist->lumo_idx][i].re = - 1.0 * sum2.re;
					Kx_cache[b-ist->lumo_idx][j][a-ist->lumo_idx][i].im = sum2.im;
				}
			}
		}
		
		for (long aux = a; aux < a+1; aux++){
			for (b = 0; b < ist->n_elecs; b++){
			for (i = 0; i < ist->n_holes; i++){
			for (j = 0; j < ist->n_holes; j++){
				ibs = listibs[(aux-ist->lumo_idx)*ist->n_holes+i];
				jbs = listibs[b*ist->n_holes+j];
			
				fprintf(pf,"%ld %ld %.10f %.10f\n", ibs, jbs, \
					exchange[ibs * ist->n_xton + jbs].re, exchange[ibs * ist->n_xton + jbs].im);
			}
			}
			}
		}
		fflush(stdout);
	}

	fclose(pf);
  printf("  Done computing exchange mat\n"); 
	fflush(0);

	} // close mpi rank 2

	// Synchronize all ranks here
	MPI_Barrier(MPI_COMM_WORLD);
  
  // If multiple ranks were used to compute the kernel
  // Send data from rank 1 to rank 0
  if (parallel->mpi_size > 1){
	if (parallel->mpi_rank == 1) {
		MPI_Send(exchange, 2*sqr(ist->n_xton), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	  }
	  if (parallel->mpi_rank == 0) {
		MPI_Recv(exchange, 2*sqr(ist->n_xton), MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  }
  }
  
	
	if ( (parallel->mpi_rank == 0) && (parallel->mpi_size != 1) ){
		long xton_idx;
		for (a = 0; a < ist->n_elecs; a++){
			for (i = 0; i < ist->n_holes; i++){
      for (b = 0; b < ist->n_elecs; b++){
      for (j = 0; j < ist->n_holes; j++){
        ibs = listibs[a*ist->n_holes + i];
        jbs = listibs[b*ist->n_holes + j];
			  xton_idx = ibs * ist->n_xton + jbs;

				bsmat[xton_idx].re += exchange[xton_idx].re;
        if(ibs==jbs){
          bsmat[ibs*ist->n_xton+jbs].im = 0.0;
        } else{
				  bsmat[xton_idx].im += exchange[xton_idx].im;
        }
      }
      }
			}
		}
  }

  if (parallel->mpi_rank == 0){
		FILE *ppsi;
		ppsi = fopen("bsRE.dat", "w");
		for (i = 0; i < ist->n_xton; i++, fprintf(ppsi,"\n")){
			for (j = 0; j < ist->n_xton; j++){
				fprintf(ppsi,"%.*g ", PR_LEN, bsmat[i*ist->n_xton+j].re);
			}
		}
		fclose(ppsi);

		ppsi = fopen("bsIM.dat", "w");
		for (i = 0; i < ist->n_xton; i++, fprintf(ppsi,"\n")){
			for (j = 0; j < ist->n_xton; j++){
				fprintf(ppsi,"%.*g ", PR_LEN, bsmat[i*ist->n_xton+j].im);
			}
		}
		fclose(ppsi);

		ppsi = fopen("h0.dat", "w");
		for (i = 0; i < ist->n_xton; i++, fprintf(ppsi,"\n")){
			for (j = 0; j < ist->n_xton; j++){
				fprintf(ppsi,"%.*g ", PR_LEN, h0mat[i*ist->n_xton+j]);
			}
		}
		fclose(ppsi);
	}

	free(rho); free(listibs);

	return;
}

/***************************************************************************************/