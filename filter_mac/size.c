/*****************************************************************************/
// File contains functions that involve calculating the total size of the system 

#include "fd.h"

/*****************************************************************************/

double get_dot_ligand_size(double *Rx, double *Ry, double *Rz, long n_atoms) {
  /*******************************************************************
  * This function calculates the max distance between two NC atoms   *
  * inputs:                                                          *
  *  [Rx] array of coordinates along the x direction                 *
  *  [Ry] array of coordinates along the y direction                 *
  *  [Rz] array of coordinates along the z direction                 *
  *  [n_atoms] number of atoms in geometry                           *
  * outputs: [double] distance between atoms in a.u.                 *
  ********************************************************************/

  long i, j;
  double dr2, dis2;
  
  for (dr2 = 0.0, i = 0; i < n_atoms-1; i++) for (j = i+1; j < n_atoms; j++) {
    dis2 = sqr(Rx[i]-Rx[j]) + sqr(Ry[i]-Ry[j]) + sqr(Rz[i]-Rz[j]);
    if (dis2 > dr2) dr2 = dis2;
  }

  return(sqrt(dr2));
}

/***************************************************************************/

double get_dot_ligand_size_z(double *R, long n_atoms){
  /*******************************************************************
  * This function calculates the size of the NC along one dimension  *
  * inputs:                                                          *
  *  [R] array of coordinates along one dimension                    *
  *  [n_atoms] number of atoms in geometry                           *
  * outputs: [double] distance between atoms in a.u.                 *
  ********************************************************************/

  long i, j;
  double dr2, dis2;
  
  for (dr2 = 0.0, i = 0; i < n_atoms-1; i++) for (j = i+1; j < n_atoms; j++) {
      dis2 = sqr(R[i]-R[j]);
      if (dis2 > dr2) dr2 = dis2;
    }

  return(sqrt(dr2));
}

/*****************************************************************************/ 
