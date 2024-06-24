/*****************************************************************************/

#include "fd.h"

/*****************************************************************************/

void get_qp_basis_indices(double *eig_vals, double *sigma_E, long **eval_hole_idxs, long **eval_elec_idxs, index_st *ist, par_st *par, flag_st *flag){
  //this is where we set which is eleectron and what is hole

  FILE *pf;
  long i, cntr, eval_homo_idx, eval_lumo_idx, old_n_holes;
  double deltaE;
  
  // Allocate memory for eval_hole_idxs etc if it is NULL
  if (NULL == *eval_hole_idxs){
    printf("\tAllocating new memory for eval_hole_idxs\n");
    *eval_hole_idxs = malloc(500 * sizeof(long));
  }
  if (NULL == *eval_elec_idxs){
    printf("\tAllocating new memory for eval_hole_idxs\n");
    *eval_elec_idxs = malloc(500 * sizeof(long));
  }

  // Set the quasiparticle basis indices to zero
  ist->homo_idx = ist->lumo_idx = ist->n_holes = ist->n_elecs = 0;

  cntr = 0; // We will reorder eig_vals and sigma_E to only contain the converged eigenstates
  // Get the indices of the HOMO, LUMO, and also the number of eigstates in VB/CB 
  for (i = 0; i < ist->mn_states_tot; i++){
    // Find HOMO and total number of VB states
    if (sigma_E[i] < par->sigma_E_cut && eig_vals[i] < par->fermi_E){
      ist->n_holes++; // increment the counter for number of hole states
      ist->homo_idx = cntr; // homo_idx in RAM is the value of cntr for which condition is met
      eval_homo_idx = i; // the value of homo_idx from the eval.dat file
      
      eig_vals[cntr] = eig_vals[i]; // reorder the eig_vals
      sigma_E[cntr] = sigma_E[i]; // reorder the sigma_E
      (*eval_hole_idxs)[cntr] = i; // add this index to the array holding all indices from eval.dat
      cntr++;
    }
    // Find LUMO and total number of CB states
    if (sigma_E[i] < par->sigma_E_cut && eig_vals[i] > par->fermi_E){
      if (0 == ist->lumo_idx){
        eval_lumo_idx = i;
        ist->lumo_idx = ist->homo_idx + 1; // Get first value of i for which condition is met
      }
      ist->n_elecs++; 
      
      eig_vals[cntr] = eig_vals[i];
      sigma_E[cntr] = sigma_E[i];
      (*eval_elec_idxs)[(long) (cntr - ist->n_holes)] = i;
      cntr++;
    }
  }
  old_n_holes = ist->n_holes;

  printf("\n\tTotal # of filtered hole eigenstates = %ld\n", ist->n_holes);
  printf("\tTotal # of filtered electron eigenstates = %ld\n", ist->n_elecs);
  printf("\tThe filter eval.dat index of the HOMO state = %ld  LUMO state = %ld\n", eval_homo_idx, eval_lumo_idx);
  
  // Check how the quasiparticle basis compares to the desired energy range
  // First, in the VB
  deltaE = eig_vals[ist->homo_idx] - eig_vals[ist->homo_idx - ist->n_holes + 1];
  if (deltaE < par->delta_E_hole){
    printf("\n\tUnconstrained energy span of holes, %lg a.u. < desired span = %lg a.u.\n", deltaE, par->delta_E_hole);
    printf("\tIncrease size of VB basis states to reach desired result\n");
  } else {
    printf("\n\tUnconstrained energy span of holes %lg a.u. > desired span = %lg a.u.\n", deltaE, par->delta_E_hole);
    printf("\tFewer VB basis states would reach the desired result\n");
  }
  // Then in the conduction band
  deltaE = eig_vals[ist->n_elecs + ist->lumo_idx - 1] - eig_vals[ist->lumo_idx];
  if (deltaE < par->delta_E_hole){
    printf("\tUnconstrained energy span of elecs, %lg a.u. < desired span = %lg a.u.\n", deltaE, par->delta_E_hole);
    printf("\tIncrease size of CB basis states to reach desired result\n");
  } else {
    printf("\tUnconstrained energy span of elecs, %lg a.u. > desired span = %lg a.u.\n", deltaE, par->delta_E_hole);
    printf("\tFewer CB basis states would reach the desired result\n");
  }

  // If the max number of electron or hole states was set, then constrain the basis
  // We will also modify the eig_vals and sigma_E array again so they only contain the
  // energies of the selected qp basis states
  
  if ((-1 != ist->max_hole_states) && (ist->n_holes > ist->max_hole_states)){
    
    ist->n_holes = ist->max_hole_states;
    printf("\n\tConstraining hole basis states to maxHoleStates\n\t  new ist->n_holes = %ld\n", ist->n_holes);
    
    // Reorder eig_vals and sigma_E to only contain eigenstates
    cntr = 0;
    for (i = 0; i < ist->n_holes; i++){
      eig_vals[cntr] = eig_vals[ist->homo_idx - ist->n_holes + i + 1];
      sigma_E[cntr] = sigma_E[ist->homo_idx - ist->n_holes + i + 1];
      (*eval_hole_idxs)[cntr] = (*eval_hole_idxs)[(long) (ist->homo_idx - ist->n_holes + i + 1)];
      cntr++;
    }
    // check that the counter has the same value as ist->n_holes
    if ((cntr - ist->n_holes) != 0){
      printf("ERROR: something went wrong reordering eig_vals in get_qp_basis_indices\n");
      exit(EXIT_FAILURE);
    }

    // Determine the new energy span
    deltaE = eig_vals[ist->homo_idx] - eig_vals[ist->homo_idx - ist->n_holes + 1];
    if (deltaE < par->delta_E_hole){
      printf("\tConstrained energy span of holes, %lg a.u. < desired span = %lg a.u.\n", deltaE, par->delta_E_hole);
      printf("\tIncrease size of VB basis states to reach desired result\n");
    } else {
      printf("\n\tConstrained energy span of holes %lg a.u. > desired span = %lg a.u.\n", deltaE, par->delta_E_hole);
      printf("\tFewer VB basis states would reach the desired result\n");
    }
    ist->homo_idx = ist->homo_idx - (old_n_holes - ist->n_holes);

  }
  // If the user does not specify max_elec or max_hole_states, then use all the quasiparticle states for BSE
  if ((-1 != ist->max_elec_states) && (ist->n_elecs > ist->max_elec_states) ){
    
    ist->n_elecs = ist->max_elec_states;
    ist->lumo_idx = ist->homo_idx + 1;
    printf("\tConstraining elec basis states to maxElecStates\n\t  new ist->n_elecs = %ld\n", ist->n_elecs);
    
    // Reorder eig_vals and sigma_E to only contain eigenstates
    cntr = ist->homo_idx + 1;
    for (i = ist->n_holes; i < ist->n_holes + ist->n_elecs; i++){
      eig_vals[cntr] = eig_vals[i];
      sigma_E[cntr] = sigma_E[i];
      // (*eval_elec_idxs)[(long) cntr - ist->n_holes] = (*eval_elec_idxs)[i - ist->n_holes];
      // we have to subtract by n_holes because we start the iterator from n_holes + 1
      cntr++;
    }
    // for (i = 0; i < ist->max_elec_states; i++){
    //   (*eval_elec_idxs)[(long) cntr - ist->n_holes] = (*eval_elec_idxs)[i];
    // }
 
    // Determine the new energy span
    deltaE = eig_vals[ist->n_elecs + ist->lumo_idx - 1] - eig_vals[ist->lumo_idx];
    if (deltaE < par->delta_E_elec){
      printf("\tConstrained energy span of elecs, %lg a.u. < desired span = %lg a.u.\n", deltaE, par->delta_E_hole);
      printf("\tIncrease size of CB basis states to reach desired result\n");
    } 
    else {
      printf("\tConstrained energy span of elecs %lg a.u. > desired span = %lg a.u.\n", deltaE, par->delta_E_hole);
      printf("\tFewer CB basis states would reach the desired result\n");
    }

  } 
  // If the user does not specify max_elec or max_hole_states, then use all the quasiparticle states for BSE
  if ((-1 == ist->max_elec_states) && (-1 == ist->max_elec_states)) {
    printf("\n\tAll filtered quasiparticle states will be used for exciton basis\n");
  }

  ist->n_qp = ist->n_holes + ist->n_elecs;

  // Print QP basis info
  printf("\n\tSelected # of filtered h+ qp basis states = %ld\n", ist->n_holes);
  printf("\tSelected # of filtered e- qp basis states = %ld\n", ist->n_elecs);
  printf("\tTotal number of quasiparticle states, n_qp = %ld\n", ist->n_qp);
  printf("\tThe BSEeval.par index of the HOMO state = %ld  LUMO state = %ld\n", ist->homo_idx, ist->lumo_idx);
  printf("\tThe HOMO energy = % .6g a.u. % .5f eV\n", eig_vals[ist->homo_idx], eig_vals[ist->homo_idx]*AUTOEV);
  printf("\tThe LUMO energy = % .6g a.u. % .5f eV\n", eig_vals[ist->lumo_idx], eig_vals[ist->lumo_idx]*AUTOEV);
  printf("\tFundamental gap = %.6g a.u. %.5f eV\n", eig_vals[ist->lumo_idx]-eig_vals[ist->homo_idx], (eig_vals[ist->lumo_idx]-eig_vals[ist->homo_idx])*AUTOEV);

  // Print BSEeval.par
  pf = fopen("BSEeval.par", "w");
  for (i = 0; i < ist->n_qp; i++){
    fprintf(pf, "% .12f %lg\n", eig_vals[i], sigma_E[i]);
  }
  fclose(pf);

  return;
}

/*****************************************************************************/

// void init(double *potl, double *vx, double *vy, double *vz, double *ksqr, double *rx, double *ry, double *rz, par_st *par, index_st *ist){
//   FILE *pf; 
//   long ntmp, jx, jy, jz, jyz, jxyz, ie, ntot, jp, *npot, nn, flags=0;
//   double del, mx, my, mz, xd, yd, zd, dx, dy, dz, *ksqrx, *ksqry, *ksqrz;
//   double *vr, *potatom, *dr;

//   printf("Final Box Dimensions: xd = %g yd = %g zd = %g\n",par->xmax,par->ymax,par->zmax);

//   if ((xd < yd) && (xd < zd))  {
//     par->minr = xd;
//     par->boxl = (double)(ist->nx) * par->dx;
//   }
//   else if ((yd < xd) && (yd < zd))  {
//     par->minr = yd;
//     par->boxl = (double)(ist->ny) * par->dy;
//   }
//   else {
//     par->minr = zd;
//     par->boxl = (double)(ist->nz) * par->dz;
//   }    
  
//   par->gamma = 7.0 / (2.0 * par->minr);
//   par->gamma2 = sqr(par->gamma);

//   printf("Gamma = %.6f\n", par->gamma);
  
//   return;
// }


// /****************************************************************************/

// void init_pot(zomplex *potq,zomplex *potqx, grid_st *grid, par_st *par, index_st *ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
// {
//   long jx, jy, jz, jyz, jxyz, sx, sy, sz;
//   double dr, dr2, x2, y2, z2, *kx2, *ky2, *kz2, alpha, cosa, sina;
//   double ex_1, ey_1, ez_1, sqrtexeyez_1, sqrtaveps, boxl2 = sqr(par->boxl), sqrk0;
//   double gammaeps;
//   zomplex *potr, *potrx, tmp;

//   ex_1 = 1.0 / par->epsX;
//   ey_1 = 1.0 / par->epsY;
//   ez_1 = 1.0 / par->epsZ;
//   sqrtexeyez_1 = 1.0 / sqrt(par->epsX * par->epsY * par->epsZ);
//   sqrtaveps = sqrt((par->epsX + par->epsY + par->epsZ) / 3.0);
//   gammaeps = par->gamma * sqrtaveps;

//   /*** no yukawa screening for the exchange ***/
//   sqrk0 = par->gamma2 * sqr(sqrtaveps);

//   if ((kx2  = (double*)calloc(ist->nx,sizeof(double)))==NULL)nerror("kx2");
//   if ((ky2  = (double*)calloc(ist->ny,sizeof(double)))==NULL)nerror("ky2");
//   if ((kz2  = (double*)calloc(ist->nz,sizeof(double)))==NULL)nerror("kz2");
//   if ((potr = (zomplex*)calloc(ist->ngrid,sizeof(zomplex)))==NULL)nerror("potr");
//   if ((potrx = (zomplex*)calloc(ist->ngrid,sizeof(zomplex)))==NULL)nerror("potrx");
  
//   for (jxyz = 0; jxyz < ist->ngrid; jxyz++) potr[jxyz].re = potr[jxyz].im = 0.0;
//   for (jxyz = 0; jxyz < ist->ngrid; jxyz++) potrx[jxyz].re = potrx[jxyz].im = 0.0;

//   for (jz = 0; jz < ist->nz; jz++) {
//     z2 =sqr(vz[jz]);
//     for (jy = 0; jy < ist->ny; jy++) {
//       y2 = sqr(vy[jy]);
//       jyz = ist->nx * (ist->ny * jz + jy);
//       for (jx = 0; jx < ist->nx; jx++) {
//       	x2 = sqr(vx[jx]);
//       	jxyz = jyz + jx;
//       	dr2 = (x2 + y2 + z2);
//       	if (dr2 < boxl2) {
//       	  dr = sqrt(x2 + y2 + z2);
//       	  potr[jxyz].re = screenedcoulomb(dr, par->gamma);
//       	  dr = sqrt(ex_1 * x2 + ey_1 * y2 + ez_1 * z2);
//       	  potrx[jxyz].re = sqrtexeyez_1 * screenedcoulomb(dr, gammaeps);
//       	}
//       }
//     }
//   }

//   for (jxyz = 0; jxyz < ist->ngrid; jxyz++)
//     potq[jxyz].re = potq[jxyz].im = potqx[jxyz].re = potqx[jxyz].im = 0.0;
  
//   memcpy(&fftwpsi[0], &potr[0], ist->ngrid*sizeof(fftwpsi[0]));
//   fftw_execute(planfw);
//   memcpy(&potq[0], &fftwpsi[0], ist->ngrid*sizeof(potq[0]));

//   memcpy(&fftwpsi[0], &potrx[0], ist->ngrid*sizeof(fftwpsi[0]));
//   fftw_execute(planfw);
//   memcpy(&potqx[0], &fftwpsi[0], ist->ngrid*sizeof(potqx[0]));
  
//   for (kx2[0] = 0.0, jx = 1; jx <= ist->nx / 2; jx++)
//     kx2[jx] = (kx2[ist->nx-jx] = sqr((double)(jx) * par->dkx));
//   for (ky2[0] = 0.0, jy = 1; jy <= ist->ny / 2; jy++)
//     ky2[jy] = (ky2[ist->ny-jy] = sqr((double)(jy) * par->dky));
//   for (kz2[0] = 0.0, jz = 1; jz <= ist->nz / 2; jz++)
//     kz2[jz] = (kz2[ist->nz-jz] = sqr((double)(jz) * par->dkz));

//   for (sz = 1.0, jz = 0; jz < ist->nz; jz++, sz = -sz) {
//     z2 = kz2[jz];
//     for (sy = 1.0, jy = 0; jy < ist->ny; jy++, sy = -sy) {
//       y2 = ky2[jy];
//       jyz = ist->nx * (ist->ny * jz + jy);
//       for (sx = 1.0, jx = 0; jx < ist->nx; jx++, sx = -sx) {
//       	x2 = kx2[jx];
//       	jxyz = jyz + jx;
//       	alpha = PIE * (double)(jx + jy + jz + ist->ngrid / 2);
//       	cosa = cos(alpha);
//       	sina = sin(alpha);

//       	/*** hartree term ***/
//       	tmp.re = potq[jxyz].re;
//       	tmp.im = potq[jxyz].im;
      	
//       	potq[jxyz].re = (tmp.re * cosa - tmp.im * sina) * par->dv;
//       	potq[jxyz].im = (tmp.re * sina + tmp.im * cosa) * par->dv;
//       	potq[jxyz].re += (FOURPI / (x2 + y2 + z2 + par->gamma2));
//       	potq[jxyz].re *= ist->ngrid_1;
//       	potq[jxyz].im *= ist->ngrid_1;

//       	/*** screened exchange term ***/
//       	tmp.re = potqx[jxyz].re;
//       	tmp.im = potqx[jxyz].im;

//       	potqx[jxyz].re = (tmp.re * cosa - tmp.im * sina) * par->dv;
//       	potqx[jxyz].im = (tmp.re * sina + tmp.im * cosa) * par->dv;
//       	//potqx[jxyz].re += FOURPI / (par->epsX * x2 + par->epsY * y2 + par->epsZ * z2 + sqrk0);
//       	potqx[jxyz].re += (FOURPI * (1.0 - exp(-0.25* (par->epsX * x2 + par->epsY * y2 + par->epsZ * z2) / sqrk0)) / (par->epsX * x2 + par->epsY * y2 + par->epsZ * z2 + EPSR));
//         //printf("denominator = % .12f\n", par->epsX * x2 + par->epsY * y2 + par->epsZ * z2);
// 		    potqx[jxyz].re *= ist->ngrid_1;
//       	potqx[jxyz].im *= ist->ngrid_1;
//       }
//     }
//   }

//   free(potr);  free(potrx); free(kx2); free(ky2); free(kz2);

//   return;
// }
  
// /****************************************************************************/

// void init_pot_old(double *vx,double *vy,double *vz,zomplex *potq,par_st par,index_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
// {
//   long jx, jy, jz, jyz, jxyz, sx, sy, sz;
//   double dr, x2, y2, z2, *kx2, *ky2, *kz2, alpha, cosa, sina;
//   double boxl = (double)(ist->nx) * par->dx;
//   zomplex *potr, tmp;

//   if ((kx2  = (double*)calloc(ist->nx,sizeof(double)))==NULL)nerror("kx2");
//   if ((ky2  = (double*)calloc(ist->ny,sizeof(double)))==NULL)nerror("ky2");
//   if ((kz2  = (double*)calloc(ist->nz,sizeof(double)))==NULL)nerror("kz2");
//   if ((potr = (zomplex*)calloc(ist->ngrid,sizeof(zomplex)))==NULL)nerror("potr");
  
//   for (jx = 0; jx < ist->ngrid; jx++) potr[jx].re = potr[jx].im = 0.0;
//   for (jz = 0; jz < ist->nz; jz++){
//     z2 = sqr(vz[jz]);
//     for (jy = 0; jy < ist->ny; jy++){
//       y2 = sqr(vy[jy]);
//       jyz = ist->nx * (ist->ny * jz + jy);
//       for (jx = 0; jx < ist->nx; jx++) {
//       	x2 = sqr(vx[jx]); 
//       	jxyz = jyz + jx;
//       	dr = sqrt(x2 + y2 + z2);
//       	if (dr < boxl) potr[jxyz].re = screenedcoulomb(dr, par->gamma);
//       }
//     }
//   }

  
//   memcpy(&fftwpsi[0],&potr[0],ist->ngrid*sizeof(fftwpsi[0]));
//   fftw_execute(planfw);
//   memcpy(&potq[0],&fftwpsi[0],ist->ngrid*sizeof(potq[0]));
//   /*fftwnd_one(planfw,potr,potq);*/

//   /*for (jx = 0; jx < ist->nx; jx++)
//     kx2[jx] = sqr(par->kxmin + (double)(jx) * par->dkx);
//   for (jy = 0; jy < ist->ny; jy++)
//     ky2[jy] = sqr(par->kymin + (double)(jy) * par->dky);
//   for (jz = 0; jz < ist->nz; jz++)
//   kz2[jz] = sqr(par->kzmin + (double)(jz) * par->dkz);*/
  
//   for (kx2[0] = 0.0, jx = 1; jx <= ist->nx / 2; jx++)
//     kx2[jx] = (kx2[ist->nx-jx] = sqr((double)(jx) * par->dkx));
//   for (ky2[0] = 0.0, jy = 1; jy <= ist->ny / 2; jy++)
//     ky2[jy] = (ky2[ist->ny-jy] = sqr((double)(jy) * par->dky));
//   for (kz2[0] = 0.0, jz = 1; jz <= ist->nz / 2; jz++)
//     kz2[jz] = (kz2[ist->nz-jz] = sqr((double)(jz) * par->dkz));

//   for (sz = 1.0, jz = 0; jz < ist->nz; jz++, sz = -sz){
//     z2 = kz2[jz];
//     for (sy = 1.0, jy = 0; jy < ist->ny; jy++, sy = -sy){
//       y2 = ky2[jy];
//       jyz = ist->nx * (ist->ny * jz + jy);
//       for (sx = 1.0, jx = 0; jx < ist->nx; jx++, sx = -sx){
// 	x2 = kx2[jx];
// 	jxyz = jyz + jx;
// 	alpha = PIE * (double)(jx + jy + jz + ist->ngrid / 2);
// 	cosa = cos(alpha);
// 	sina = sin(alpha);
	
// 	/*potq[jxyz].re *= (double)(sx * sy * sz) * par->dv;
// 	  potq[jxyz].im *= (double)(sx * sy * sz) * par->dv;*/

// 	tmp.re = potq[jxyz].re;
// 	tmp.im = potq[jxyz].im;

// 	potq[jxyz].re = (tmp.re * cosa - tmp.im * sina) * par->dv;
// 	potq[jxyz].im = (tmp.re * sina + tmp.im * cosa) * par->dv;
	
// 	potq[jxyz].re += (FOURPI / (x2 + y2 + z2 + par->gamma2));
    
// 	potq[jxyz].re *= ist->ngrid_1;
// 	potq[jxyz].im *= ist->ngrid_1;

//       }
//     }
//   }
//   free(potr);  free(kx2); free(ky2); free(kz2);
//   return;
// }

// /************************************************************************/

// #define SQRTPI       (sqrt(3.14159265358979323846))

// double screenedcoulomb(double dr, double gamma)
// {
//   //if (dr < EPSR) return (gamma);
//   if (dr < EPSR) return (2.0*gamma/SQRTPI);
//   return (erf(gamma * dr) / dr);
//   //return ((1.0 - exp(-gamma * dr)) / dr);
// }

// /************************************************************************/

// void init_psi(zomplex *psi,double *vx,double *vy,double *vz,index_st ist,par_st par,long *idum)
// {
//   long jx, jy, jz, jzy, jxyz;
//   long tidum = (*idum);

//   for (jz = 0; jz < ist->nz; jz++) for (jy = 0; jy < ist->ny; jy++){
//     for (jzy = ist->nx * (ist->ny * jz + jy), jx = 0; jx < ist->nx; jx++){
//       jxyz = jzy + jx;
//       psi[jxyz].re = (-1.0 + 2.0 * ran_nrc(&tidum));
//       psi[jxyz].im = 0.0;
//     }
//   }
//   normalize_zomplex(psi, par->dv, ist->ngrid);
//   (*idum) = tidum;
//   return;
// }

/************************************************************************/
