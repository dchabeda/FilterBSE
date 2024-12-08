#include "pot_coupling.h"

/****************************************************************************/
// 

void write_cube_file(double *rho, grid_st *grid, char *fileName) {
  /*****************************************************************
  * This function prints out cube files to visualize wavefunctions *
  * inputs: [double *rho] is a pointer to an ngrid long array      *
  *         [grid_st *grid] is a pointer to the grid struct        *
  *         [char *fileName] is a pointer to the output file name  *
  * outputs: void                                                  *
  ******************************************************************/ 

  FILE *pf, *pConfFile;
  
  long jgrid, iX, iY, iZ, iYZ, natoms, atomType;
  double x, y, z;
  char line[80], atomSymbol[10];

  if( access("conf.dat", F_OK) == -1 ){
        printf("ERROR: no conf.dat file exists in directory\n");
        fprintf(stderr, "ERROR: no conf.dat file exists in directory\n");
        exit(EXIT_FAILURE);
  } else{
    pConfFile = fopen("conf.dat", "r");
  }
  fscanf(pConfFile, "%ld", &natoms);

  pf = fopen(fileName, "w");
  
  fprintf(pf, "CUBE FILE\n");
  fprintf(pf, "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n");
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", natoms, grid->xmin, grid->ymin, grid->zmin);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", grid->nx, grid->dx, 0.0, 0.0);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", grid->ny, 0.0, grid->dy, 0.0);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", grid->nz, 0.0, 0.0, grid->dz);
  fgets(line, 80, pConfFile); 
  while(fgets(line, 80, pConfFile) != NULL) {
    sscanf(line, "%2s %lf %lf %lf", (char*)&atomSymbol, &x, &y, &z);
    
    //TODO: make sane atom_info handling.... 
    if (! strcmp(atomSymbol, "Cd")) { 
      atomType = 48;
    }
    else if (! strcmp(atomSymbol, "S")) {  
      atomType = 16;
    }
    else if (! strcmp(atomSymbol, "Se")) { 
      atomType = 34;
    }
    else if (! strcmp(atomSymbol, "Zn")) {
      atomType = 30;
    }
    else if (! strcmp(atomSymbol, "Te")) {
      atomType = 52;
    }
  	else if (! strcmp(atomSymbol, "C")) {
  	  atomType = 6;
  	}
    else if (! strcmp(atomSymbol, "Si")) {
	  atomType = 14;
    }
    else if (! strcmp(atomSymbol, "Cs")) {
    atomType = 55;
    }
    else if (! strcmp(atomSymbol, "Pb")) {
    atomType = 82;
    }
    else if (! strcmp(atomSymbol, "I")) {
    atomType = 53;
    }
    else { 
      atomType = 1; 
    }
    fprintf(pf, "%5li%12.6f%12.6f%12.6f%12.6f\n", atomType, 0.0, x, y, z);
  }
  for (iX = 0; iX < grid->nx; iX++) {
    for (iY = 0; iY < grid->ny; iY++) {
      for (iZ = 0; iZ < grid->nz; iZ++) {
        iYZ = grid->nx * (grid->ny * iZ + iY);
        jgrid = iYZ + iX;
        fprintf(pf, "%12.5f ", rho[jgrid]);
        if (iZ % 6 == 5) {
          fprintf(pf, "\n");
        }
      }
      fprintf(pf, "\n");
    }
  }
  fclose(pConfFile);
  fclose(pf);

  return;
}

void write_separation(FILE *pf, char *top_bttm) {
  /*****************************************************************
  * This function prints asterisk separation lines in stdout       *
  * inputs:                                                        *
  *   [FILE *pf] pointer to output file stream                     *
  *   [char *top_bttm] pointer to char for top/bottom formatting   *
  * outputs: void                                                  *
  ******************************************************************/

  char *top_key; top_key = malloc(2*sizeof(top_key[0]));
  char *bttm_key; bttm_key = malloc(2*sizeof(bttm_key[0]));
  
  strcpy(top_key, "T");
  strcpy(bttm_key, "B");

  if ( 0 == strcmp(top_bttm, (const char *) top_key) ){
    fprintf(pf, "\n\n******************************************************************************\n");
  } else if ( 0 == strcmp(top_bttm, (const char *) bttm_key) ){
    fprintf(pf, "\n******************************************************************************\n");
  } else {
    fprintf(stderr, "Invalid string supplied to writeSeparation. Exiting!\n");
    exit(EXIT_FAILURE);
  }

  free(top_key); free(bttm_key);
  
  return;
}
