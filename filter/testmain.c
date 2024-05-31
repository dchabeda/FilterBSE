
#include "fd.h"
#include <unistd.h>

/*****************************************************************************/

int main(int argc, char *argv[])
{
  FILE *ppsi; zomplex *psitot;
  par_st par;   index_st  ist; atom_info *atm;
  double  *gridX, *gridY, *gridZ;
  double *Rx, *Ry, *Rz;
  double *dr, *vr, *potForAtom, dx,dy,dz; 
  long nStatesInEval, i, ieof,flags=0;
  long jx,jy,jz;
  fftw_plan_loc planfw, planbw; fftw_complex *fftwpsi;
  time_t currentTime = time(NULL);


  printf("This calculation began at: %s", ctime(&currentTime)); 
  fflush(stdout);

  /*** read initial setup from input.par ***/
  init_size(argc, argv, &par, &ist);

  /*** allocating memoRy ***/
  /*** the positions of the atom_infos in the x, y, and z directions ***/
  if ((Rx = (double *) calloc(ist.natoms, sizeof(double))) == NULL) nerror("Rx");
  if ((Ry = (double *) calloc(ist.natoms, sizeof(double))) == NULL) nerror("Ry");
  if ((Rz = (double *) calloc(ist.natoms, sizeof(double))) == NULL) nerror("Rz");
  if ((atm = (atom_info *) calloc(ist.natoms, sizeof(atom_info))) == NULL) nerror("atm");
  init_conf(Rx, Ry, Rz, atm, &par, &ist);

    /*** the grid in the x, y, and z directions ***/
  if ((gridX = (double *) calloc(ist.nx, sizeof(double))) == NULL) nerror("gridX");
  if ((gridY = (double *) calloc(ist.ny, sizeof(double))) == NULL) nerror("gridY");
  if ((gridZ = (double *) calloc(ist.nz, sizeof(double))) == NULL) nerror("gridZ");
  
  par.dv = par.dx * par.dy * par.dz;
  /***initializing the potential vector  ***/
  for (jx = 0, dx = par.xmin; jx < ist.nx; jx++, dx += par.dx) gridX[jx] = dx;
  for (jy = 0, dy = par.ymin; jy < ist.ny; jy++, dy += par.dy) gridY[jy] = dy;
  for (jz = 0, dz = par.zmin; jz < ist.nz; jz++, dz += par.dz) gridZ[jz] = dz;


/*** initialization for the fast Fourier transform ***/
  fftwpsi = fftw_malloc(sizeof (fftw_complex )*ist.ngrid);
  planfw = fftw_plan_dft_3d(ist.nz, ist.ny, ist.nx, fftwpsi, fftwpsi, FFTW_FORWARD, flags);
  planbw = fftw_plan_dft_3d(ist.nz, ist.ny, ist.nx, fftwpsi, fftwpsi, FFTW_BACKWARD, flags);
  

  /*** read in all wavefunctions in the energy range  ***/
    // Get the number of states in eval.par
  nStatesInEval = 0;
  if ( access("eval.dat", F_OK) != -1 ) {
    ppsi = fopen("eval.dat" , "r");
    for ( i = ieof = 0; ieof != EOF; i++) {
      ieof = fscanf(ppsi, "%*ld %*lg %*lg");
    }
    nStatesInEval = i-1;
    fclose(ppsi);
  }
  else {
    printf("\n\nNo eval.dat file detected in current working directoRy - the program is exiting!!!\n\n");
    exit(EXIT_FAILURE);
  }

  ist.mn_states_tot = nStatesInEval;

  if ((psitot = (zomplex *) calloc(ist.nspinngrid*ist.mn_states_tot, sizeof(zomplex))) == NULL) nerror("psitot");
  printf("Reading %d states from psi.dat from filter run\n", ist.mn_states_tot);
  if ( access("psi.dat", F_OK) != -1 ){
    ppsi = fopen("psi.dat", "r");
    fread(psitot, sizeof(zomplex), ist.nspinngrid*ist.mn_states_tot,ppsi);
  }
  else{
    printf("\n\nNo psi.dat file detected in cwd - exiting!!!\n\n");
    nerror("readPsi");
  }

  /*** calculate the spin and angular momentum statistics ***/
  calc_angular_exp(psitot, gridX, gridY, gridZ, planfw, planbw, fftwpsi, ist, par, 300, 333);
  
  
  double* rho = calloc(ist.ngrid,sizeof(double));
  char filename[20];
  for (int j = 300; j<333; j++){ 
  
  for (long i = 0;i<ist.ngrid; i++){
    rho[i] = (psitot[j*ist.nspinngrid+i].re);
  }
  sprintf(filename, "psi%iUpRe.cube", j);
  write_cube_file(rho, par,ist, filename);

  for (long i = 0;i<ist.ngrid; i++){
    rho[i] = (psitot[j*ist.nspinngrid+i].im);
  }
  sprintf(filename, "psi%iUpIm.cube", j);
  write_cube_file(rho, par,ist, filename);

  for (long i = 0;i<ist.ngrid; i++){
    rho[i] = (psitot[j*ist.nspinngrid+ist.ngrid+i].re);
  }
  sprintf(filename, "psi%iDnRe.cube", j);
  write_cube_file(rho, par,ist, filename);

  for (long i = 0;i<ist.ngrid; i++){
    rho[i] = (psitot[j*ist.nspinngrid+ist.ngrid+i].im);
  }
  sprintf(filename, "psi%iDnIm.cube", j);
  write_cube_file(rho, par,ist, filename);
  }

  free(Rx); free(Ry); free(Rz); free(atm);
  free(gridX); free(gridY); free(gridZ); free(psitot);
}
