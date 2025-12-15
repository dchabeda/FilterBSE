/***************************************************************************************/

#include "hartree.h"

/***************************************************************************************/

void hartree(
  double complex* rho, 
  double complex* potq, 
  double complex* poth, 
  index_st*       ist, 
  fftw_plan_loc   planfw, 
  fftw_plan_loc   planbw, 
  fftw_complex*   fftwpsi
  ){
    
  long j;

  memcpy(&fftwpsi[0],&rho[0],ist->ngrid*sizeof(fftwpsi[0]));
  fftw_execute(planfw);

  for (j = 0; j < ist->ngrid; j++) {
    fftwpsi[j] *= potq[j];
  }

  fftw_execute(planbw);
  memcpy(&poth[0], &fftwpsi[0], ist->ngrid*sizeof(fftwpsi[0]));
  

  return;
}

/***************************************************************************************/
