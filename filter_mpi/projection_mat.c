#include "fd.h"

int countlines(const char *filename);
void allocate_memory(void **ptr, size_t length, size_t type_size, char* message);
void read_eval(FILE*, double*, double*, long*, double, double, int, int*, long*);
void read_psi(FILE*, zomplex*, long*, long, long, long, long);
void calc_U_proj(zomplex*, zomplex*, zomplex*, double, long, long, long, long);
void write_cube(double *rho, char *fileName);

int main(int argc, char *argv[]){
	if (argc!=10) {
    printf("Usage: projection_mat.x eval1.dat psi1.dat eval2.dat psi2.dat ngrid nhomo nlumo deps dv\n"); 
    exit(EXIT_FAILURE);
  }

  // ----------------------------------------
  // Declare variables
  FILE *peval1;
  FILE *peval2;
  FILE *ppsi1;
  FILE *ppsi2;

  zomplex *psi1;
  zomplex *psi2;

  double *eval1, *deval1;
  double *eval2, *deval2;
  double evalloc, deloc;
	double fermiEnergy = -0.16;

  long i, j;
  long a, b;
  long k;
  long ieof;
  long ihomo1;
  long ilumo1;
  long ihomo2;
  long ilumo2;
  long *eval1_idxs;
  long *eval2_idxs;

	long ngrid = atol(argv[5]);
	long nhomo = atol(argv[6]);
	long nlumo = atol(argv[7]);
	double deps = atof(argv[8]);
  double dv = atof(argv[9]);

  long nspinngrid = 2*ngrid;
  long n_qp = nhomo + nlumo;

	int n_lines_1, n_lines_2;
	int n_states_1, n_states_2;
  
  // ----------------------------------------
	// Count the number of lines in the eval.dat files
	n_lines_1 = countlines(argv[1]);
	n_lines_2 = countlines(argv[3]);

  // ----------------------------------------
  // Print job setup info
	printf("\nRunning projection_mat.c with the following inputs:\n");
	printf("eval1 file: %s, length: %d lines\n", argv[1], n_lines_1);
	printf("psi1 file: %s\n", argv[2]);
	printf("eval2 file: %s, length: %d lines\n", argv[3], n_lines_2);
	printf("psi2 file: %s\n", argv[4]);
	printf("ngrid: %ld\n", ngrid);
	printf("nspinngrid: %ld\n", nspinngrid);
	printf("nhomo = %ld; nlumo = %ld\n", nhomo, nlumo);
	printf("deps = %g\n", deps);
  printf("dv = %lg\n", dv);
	printf("nhomo+nlumo = %ld\n", nhomo+nlumo);

  // ----------------------------------------
	// allocate memory for the energies and variances
  ALLOCATE(&eval1, n_lines_1, "eval1");
  ALLOCATE(&eval2, n_lines_2, "eval2");
  ALLOCATE(&deval1, n_lines_1, "deval1");
  ALLOCATE(&deval2, n_lines_2, "deval2");
  ALLOCATE(&eval1_idxs, n_lines_1, "eval1_idxs");
  ALLOCATE(&eval2_idxs, n_lines_2, "eval2_idxs");
	ALLOCATE(&psi1, n_qp * nspinngrid, "psi1");
  ALLOCATE(&psi2, n_qp * nspinngrid, "psi2");

	// ----------------------------------------
	// Read the eval files into arrays
  peval1 = fopen(argv[1], "r");
	read_eval(peval1, eval1, deval1, eval1_idxs, fermiEnergy, deps, n_lines_1, &n_states_1, &ihomo1);
	printf("%d eval1 energies successfully acquired\n", n_states_1);
  fclose(peval1);

  peval2 = fopen(argv[3], "r");
  read_eval(peval2, eval2, deval2, eval2_idxs, fermiEnergy, deps, n_lines_2, &n_states_2, &ihomo2);
	printf("%d eval2 energies successfully acquired\n", n_states_2);
  fclose(peval2);

  // ----------------------------------------
  // Read in the wavefunctions
  ppsi1 = fopen(argv[2], "r");
  read_psi(ppsi1, psi1, eval1_idxs, nspinngrid, nhomo, nlumo, ihomo1);
  printf("Psi1 states successfully loaded\n");
  fclose(ppsi1);

  ppsi2 = fopen(argv[4], "r");
  read_psi(ppsi2, psi2, eval2_idxs, nspinngrid, nhomo, nlumo, ihomo2);
  printf("Psi2 states successfully loaded\n");
	fclose(ppsi2);

  // ----------------------------------------
  // For debugging, print out the wavefunctions
  double* rho;
  ALLOCATE(&rho, ngrid, "rho");
  char str[40];

  for (i = 0; i < n_qp; i++) {
    for (k = 0; k < ngrid; k++) {
      rho[k] = sqr(psi1[i*nspinngrid + k].re) + sqr(psi1[i*nspinngrid + k].im);
    }
    sprintf(str, "psi1_%ld.cube", i);
    write_cube(rho, str);
  }

  for (i = 0; i < n_qp; i++) {
    for (k = 0; k < ngrid; k++) {
      rho[k] = sqr(psi2[i*nspinngrid + k].re) + sqr(psi2[i*nspinngrid + k].im);
    }
    sprintf(str, "psi2_%ld.cube", i);
    write_cube(rho, str);
  }
  

  // ----------------------------------------
  // Compute and print the unitary transformation matrix
  printf("Computing unitary transformation matrix elements\n");
  FILE* pU;
  zomplex* U;

  ALLOCATE(&U, n_qp*n_qp*sizeof(zomplex), "U");

  calc_U_proj(U, psi1, psi2, dv, nspinngrid, 0, nhomo, n_qp);

  pU = fopen("U_ij.dat", "w");
  for (i = 0; i < nhomo; i++) {
    for (j = 0; j < nhomo; j++) {
      fprintf(pU, "%ld %ld %lg %lg\n", i, j, U[i*n_qp + j].re, U[i*n_qp + j].im);
    }
  }
  fclose(pU);

  calc_U_proj(U, psi1, psi2, dv, nspinngrid, nhomo, nhomo+nlumo, n_qp);

  pU = fopen("U_ab.dat", "w");
  for (a = nhomo; a < nhomo + nlumo; a++) {
    for (b = nhomo; b < nhomo + nlumo; b++) {
      fprintf(pU, "%ld %ld %lg %lg\n", a, b, U[a*n_qp + b].re, U[a*n_qp + b].im);
    }
  }
  fclose(pU);

  // Confirm that the matrix is unitary by computing UU+ and U+U.
  
  zomplex* UU_dagger;
  zomplex* U_daggerU;
  ALLOCATE(&UU_dagger, n_qp * n_qp, "UU_dagger");
  ALLOCATE(&U_daggerU, n_qp * n_qp, "U_daggerU");

  for (i = 0; i < n_qp; i++) {
    for (a = 0; a < n_qp; a++) {
      // UU+ = sum_a U_ia * U^*_ai
      UU_dagger[i*n_qp+a].re += U[i*n_qp+a].re * U[a*n_qp+i].re + U[i*n_qp+a].im * U[a*n_qp+i].im;
      UU_dagger[i*n_qp+a].im += U[i*n_qp+a].im * U[a*n_qp+i].re - U[i*n_qp+a].re * U[a*n_qp+i].im;

      // U+U = U^*_ai * U_ia
      U_daggerU[i*n_qp+a].re += U[a*n_qp+i].re * U[i*n_qp+a].re + U[a*n_qp+i].im * U[i*n_qp+a].im;
      U_daggerU[i*n_qp+a].im += U[a*n_qp+i].re * U[i*n_qp+a].im - U[a*n_qp+i].im * U[i*n_qp+a].re;
    }
  }

  FILE* pUUd;
  FILE* pUdU;
  pUUd = fopen("UU+.dat", "w");
  pUdU = fopen("U+U.dat", "w");

  for (i = 0; i < n_qp; i++) {
    for (a = 0; a < n_qp; a++) {
      fprintf(pUUd, "%ld %ld %lg %lg\n", i, a, UU_dagger[i*n_qp+a].re, UU_dagger[i*n_qp+a].im);
      fprintf(pUdU, "%ld %ld %lg %lg\n", i, a, U_daggerU[i*n_qp+a].re, U_daggerU[i*n_qp+a].im);
    }
  }
  fclose(pUUd);
  fclose(pUdU);

  printf("Done with projection_mat.x!\n");

  // ----------------------------------------
  // Test projection onto new basis

  zomplex* psi2_proj;
  zomplex u_ia;
  long state_i, state_a;
  ALLOCATE(&psi2_proj, n_qp * nspinngrid, "psi2_proj");
  FILE *pf = fopen("u_ia.dat", "w");
  // psi2proj_a = sum_i U_ia psi1_i
  for (a = 0; a < n_qp; a++) {
    state_a = a * nspinngrid;
    for (i = 0; i < n_qp; i++) {
      state_i = i * nspinngrid;
      u_ia.re = U[i*n_qp + a].re;
      u_ia.im = U[i*n_qp + a].im;
      
      for (k = 0; k < nspinngrid; k++) {
        psi2_proj[state_a + k].re += u_ia.re * psi1[state_i + k].re - u_ia.im * psi1[state_i + k].im;
        psi2_proj[state_a + k].im += u_ia.re * psi1[state_i + k].im + u_ia.im * psi1[state_i + k].re;
      }
    }
  }
  fclose(pf);

  // Print out psi2_proj as a cube

  for (a = 0; a < n_qp; a++) {
    for (k = 0; k < ngrid; k++) {
      rho[k] = sqr(psi2_proj[a*nspinngrid + k].re) + sqr(psi2_proj[a*nspinngrid + k].im);
    }
    sprintf(str, "psi2proj_%ld.cube", a);
    write_cube(rho, str);
  }
  free(rho);

  free(eval1);
  free(eval2);
  free(deval1);
  free(deval2);
  free(eval1_idxs);
  free(eval2_idxs);
  free(psi1);
  free(psi2);
  free(U);
}

int countlines(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Could not open file");
        return -1;
    }

    int count = 0;
    char ch;

    // Read through the file character by character
    while ((ch = fgetc(file)) != EOF) {
        if (ch == '\n') {
            count++;
        }
    }

    fclose(file);
    return count;
}

void allocate_memory(void **ptr, size_t length, size_t type_size, char* message) {
    *ptr = calloc(length, type_size);
    if (*ptr == NULL) {
        fprintf(stderr, "Memory allocation failed: %s\n", message);
        exit(EXIT_FAILURE);
    }
    return;
}

void read_eval(
  FILE*    peval,
  double*  eig_vals,
  double*  sigma_E,
  long*    eval_idxs,
  double   fermiE,
  double   secut,
  int      n_lines,
  int*     n_states,
  long*    ihomo
  ){

  int i;
  int cntr;
  int tmp;

  long homo_idx;

  double eval;
  double deval;

  if (peval == NULL){
		printf("Error opening eval file!\n");
		exit(EXIT_FAILURE);
	}

  cntr = 0;
  // Read eval file into arrays
  for (i = 0; i < n_lines; i++){
    fscanf(peval, "%d %lg %lg", &tmp, &eval, &deval);
    // If the sigma_E val > secut, skip
    if ( (deval < secut) ){
      eig_vals[cntr] = eval;
      sigma_E[cntr] = deval;
      eval_idxs[cntr] = i;
      if (eval < fermiE){
        homo_idx = cntr;
      }
      cntr++;
    }
  }

  *n_states = cntr;
  *ihomo = homo_idx;

  return;
}

void read_psi(
  FILE*   ppsi,
  zomplex* psi,
  long*   eval_idxs,
  long    ngrid,
  long    nhomo,
  long    nlumo,
  long    ihomo  
  ) {

  int cntr = 0;
  int i;
  int state_idx;
  long file_offset;
  long ilumo = ihomo + 1;
  

  if (ppsi == NULL){
    printf("Error opening psi file!\n");
		exit(EXIT_FAILURE);
  }
  // ----------------------------------------
  // Read in nhomo hole states
  for (i = ihomo; i > ihomo - nhomo; i--){
    state_idx = eval_idxs[i];
    if ((state_idx == 0) && (i != 0)) {
      printf("Error: overrun converged hole eigenvalues!\n");
      exit(EXIT_FAILURE);
    }

    file_offset = state_idx * ngrid * sizeof(zomplex);
    fseek(ppsi, file_offset, SEEK_SET);
    fread(&psi[cntr*ngrid], sizeof(zomplex), ngrid, ppsi);
    cntr++;
  }

  // ----------------------------------------
  // Read in nlumo electron states
  for (i = ilumo; i < ilumo + nlumo; i++) {
    state_idx = eval_idxs[i];
    if ((state_idx == 0) && (i != 0)) {
      printf("Error: overrun converged elec eigenvalues!\n");
      exit(EXIT_FAILURE);
    }
    file_offset = state_idx * ngrid * sizeof(zomplex);
    fseek(ppsi, file_offset, SEEK_SET);
    fread(&psi[cntr*ngrid], sizeof(zomplex), ngrid, ppsi);
    cntr++;
  }

  return;

}

void calc_U_proj(
  zomplex* U,
  zomplex* psi1,
  zomplex* psi2,
  double dv,
  long nspinngrid,
  long start,
  long end,
  long tot_dim
  ) {

  long i, j, k;
  zomplex u_ij;

  for (i = start; i < end; i++) {
    for (j = start; j < end; j++) {
      u_ij.re = u_ij.im = 0;
      for (k = 0; k < nspinngrid; k++){
        u_ij.re += 
          psi1[i*nspinngrid + k].re * psi2[j*nspinngrid + k].re 
          + psi1[i*nspinngrid + k].im * psi2[j*nspinngrid + k].im;
        
        u_ij.im += 
          psi1[i*nspinngrid + k].re * psi2[j*nspinngrid + k].im 
          - psi1[i*nspinngrid + k].im * psi2[j*nspinngrid + k].re;
      }
      U[i * tot_dim + j].re = u_ij.re * dv;
      U[i * tot_dim + j].im = u_ij.im * dv;
    }
  }

  return;
}

void write_cube(double *rho, char *fileName) {
  /*****************************************************************
  * This function prints out cube files to visualize wavefunctions *
  * inputs: [double *rho] is a pointer to an ngrid long array      *
  *         [grid_st *grid] is a pointer to the grid struct        *
  *         [char *fileName] is a pointer to the output file name  *
  * outputs: void                                                  *
  ******************************************************************/ 

  long nx = 72;
  long ny = 72;
  long nz = 72;
  double dx = 0.5;
  double dy = 0.5;
  double dz = 0.5;
  double xmin = -18.0;
  double ymin = -18.0;
  double zmin = -18.0;


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
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", natoms, xmin, ymin, zmin);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", nx, dx, 0.0, 0.0);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", ny, 0.0, dy, 0.0);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", nz, 0.0, 0.0, dz);
  
  fgets(line, 80, pConfFile); 
  while(fgets(line, 80, pConfFile) != NULL) {
    sscanf(line, "%2s %lf %lf %lf", (char*)&atomSymbol, &x, &y, &z);
    
    if (! strcmp(atomSymbol, "Cd")) atomType = 48;
    else if (! strcmp(atomSymbol, "Se")) atomType = 34;
    else if (! strcmp(atomSymbol, "Cs")) atomType = 55;
    else if (! strcmp(atomSymbol, "Pb")) atomType = 82;
    else if (! strcmp(atomSymbol, "I")) atomType = 53;
    else { atomType = 1;}

    fprintf(pf, "%5li%12.6f%12.6f%12.6f%12.6f\n", atomType, 0.0, x, y, z);
  }
  for (iX = 0; iX < nx; iX++) {
    for (iY = 0; iY < ny; iY++) {
      for (iZ = 0; iZ < nz; iZ++) {
        iYZ = nx * (ny * iZ + iY);
        jgrid = iYZ + iX;
        fprintf(pf, "%.5g ", rho[jgrid]);
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