#include "fd.h"

/*****************************************************************************/

double interpolate(double r,double dr,double *vr,double *vr_LR,double *pot,double *pot_LR,long potFileLen,long n,long j, int scale_LR, double scale_LR_par, double strain_factor){

  double a, b;
  double a_LR, b_LR;
  double val_at_gridpt;
  long i;

  i = (long)(r / dr);
  if (i > (n - 2)) return (0.0);

  a = strain_factor * (pot[j*potFileLen+i+1] - pot[j*potFileLen+i]) / (vr[j*potFileLen+i+1] - vr[j*potFileLen+i]);
  b = strain_factor * pot[j*potFileLen+i] - vr[j*potFileLen+i] * a;
  
  val_at_gridpt = a * r + b;
  
  if (1 == scale_LR){
    a_LR = (scale_LR_par * (pot_LR[j*potFileLen+i+1] - pot_LR[j*potFileLen+i])) / (vr_LR[j*potFileLen+i+1] - vr_LR[j*potFileLen+i]);
    b_LR = (scale_LR_par * pot_LR[j*potFileLen+i]) - vr_LR[j*potFileLen+i] * a_LR;

    val_at_gridpt += ( (a_LR * r + b_LR) );
  }

  return val_at_gridpt;
}

/*****************************************************************************/
