/*****************************************************************************/

#include "fd.h"

/*****************************************************************************/


double norma(zomplex *psi, double dr,long ngrid)
{
  long i;
  double norm = 0.0;

  for (i = 0; i < ngrid; i++) norm += (sqr(psi[i].re) + sqr(psi[i].im));
  return (sqrt(norm * dr));
}

/*****************************************************************************/

double normalize_zomplex(zomplex *psi, double dr, long ngrid) {
	long i;
	double N = norma(psi, dr, ngrid);

	for (i = 0; i < ngrid; i++) {
		psi[i].re /= N;
		psi[i].im /= N;
	}

	return N;
}

/*****************************************************************************/

void normalize_all(double *psi,double dr,long n_qp,long ngrid)
{
  long i, ie, ieg;
  double nrm;

  for (ie = 0; ie < n_qp; ie++){
    for (ieg = ie * ngrid, nrm = 0.0, i = 0; i < ngrid; i++)
      nrm += sqr(psi[ieg+i]);
    if (nrm != 0.0){
      nrm = 1.0 / sqrt(nrm * dr);
//      printf ("%ld norm = %g\n",ie,nrm); fflush(0);
    }
    for (ieg = ie * ngrid, i = 0; i < ngrid; i++)
      psi[ieg+i] *= nrm;
  }
  return;
}

/*****************************************************************************/

void norm_vector(double *vector, double dV, long length) {
	long i;
	double norm = 0.0;

	for (i = 0; i < length; i++) norm += vector[i];
	if (norm > 1e-10) {
		norm *= dV;
		for (i = 0; i < length; i++) {
			vector[i] /= norm;
		}
	}

	return;
}

/*****************************************************************************/
