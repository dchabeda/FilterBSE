#include "fd.h"

/*****************************************************************************/

double interpolate(double r,double dr,double *vr,double *vr_LR,double *pot,double *pot_LR,long pot_file_len,long n,long j, int scale_LR, double scale_LR_par, double strain_factor, int is_LR){

  double a, b;
  double a_LR, b_LR;
  double val_at_gridpt;
  long i;

  i = (long)(r / dr);

  // This line causing the potential to be set to zero at some distance cutoff from an atomic center.
  // This is incorrect behavior for the long range potentials. 
  if (0 == is_LR){
    if (i > (n - 2)) return (0.0);
  }
  // Daniel C: patch fix just to run some jobs for a presentation. The distance in psuedopotential files only goes to 200 Bohr
  if (r > 199.0){
    return 0.0;
  }
  
  a = strain_factor * (pot[j*pot_file_len+i+1] - pot[j*pot_file_len+i]) / (vr[j*pot_file_len+i+1] - vr[j*pot_file_len+i]);
  b = strain_factor * pot[j*pot_file_len+i] - vr[j*pot_file_len+i] * a;
  
  val_at_gridpt = a * r + b;
  
  if (1 == scale_LR){
    a_LR = (scale_LR_par * (pot_LR[j*pot_file_len+i+1] - pot_LR[j*pot_file_len+i])) / (vr_LR[j*pot_file_len+i+1] - vr_LR[j*pot_file_len+i]);
    b_LR = (scale_LR_par * pot_LR[j*pot_file_len+i]) - vr_LR[j*pot_file_len+i] * a_LR;

    val_at_gridpt += ( (a_LR * r + b_LR) );
  }

  return val_at_gridpt;
}

/*****************************************************************************/
