/*****************************************************************************/
// Main file for cube printing utility.
#include "fd.h"

/*****************************************************************************/
int countlines(char *filename);


/*****************************************************************************/
int main(int argc, char *argv[])
{
  FILE *ppsi, *pf;  
  zomplex *psi, *coeff_bs;
  // custom structs 
  grid_st grid; par_st par; index_st ist; parallel_st parallel; flag_st flag;
  xyz_st *R; atom_info *atom;
  // double arrays
  double *rho; 
  // long int arrays and counters
  long i, a, jms, ibs, n;
  int numExcStatesToPrint = 100;
  int j, start, end, n_qp, ms2, nhomo, nlumo, totalhomo, totallumo;
  ist.atom_types = malloc(N_MAX_ATOM_TYPES*sizeof(ist.atom_types[0]));
  time_t currentTime = time(NULL);


  //command line input parsing
  if (argc!=9){
    printf("Usage: carrier_dens. xton_start xton_end n_qp n_ms2 nhomo totalhomo nlumo totallumo\n");
    exit(EXIT_FAILURE);
  }

  start = atoi(argv[1]);
  end = atoi(argv[2]);
  n_qp = atoi(argv[3]);
  ms2 = atoi(argv[4]);
  nhomo = atoi(argv[5]);
  totalhomo = atoi(argv[6]);
  nlumo = atoi(argv[7]);
  totallumo = atoi(argv[8]);

  if (start > end){
    printf("Invaid start (%d), end(%d): start > end\n", start,end);
    exit(EXIT_FAILURE);
  }
  if (start < 0){
    printf("Invaid start (%d): start < 0\n", start);
    exit(EXIT_FAILURE);
  }

  printf("This calculation began at: %s", ctime(&currentTime)); 
  fflush(stdout);

  coeff_bs = (zomplex *) calloc(ms2*ms2, sizeof(zomplex));
  long *listibs = (long *) calloc(ms2, sizeof(long));

  for (ibs = 0, a = nlumo; a < nlumo+totallumo; a++) {
      for (i = 0; i < totalhomo; i++, ibs++) {
          listibs[(a - nlumo) * totalhomo + i] = ibs;
          printf("a:%ld i:%ld ibs:%ld\n",a,i,ibs);
      }
  }

  /*** read initial setup from input.par ***/
  printf("\nReading job specifications from input.par:\n");
  read_input(&flag, &grid, &ist, &par, &parallel);

  /*** allocating memory ***/
  // the positions of the atoms in the x, y, and z directions 
  if ((R = (xyz_st *) calloc(ist.natoms, sizeof(xyz_st))) == NULL) nerror("R");
  // the atom specific information 
  if ((atom = (atom_info *) calloc(ist.natoms, sizeof(atom_info))) == NULL) nerror("atom");
  
  /*** read the nanocrystal configuration ***/
  printf("\nReading atomic configuration from conf.par:\n"); fflush(0);
  read_conf(R, atom, &ist, &par, &flag);
  
  /*** initialize the grid ***/
  printf("\nInitializing the grid parameters:\n"); fflush(0);
  init_grid_params(&grid, R, &ist, &par);

  if ((rho = (double *)calloc(ist.ngrid, sizeof(double)))==NULL) nerror("rho");


  //count number of states found
  jms = countlines("eval.dat");
  printf("%ld total states in psi.dat\n", jms); fflush(0);
  
  //allocate memory for all psi
  if ((psi = (zomplex *) calloc(n_qp * ist.nspinngrid, sizeof(zomplex))) == NULL) nerror("psi");


  //read psi from file
  ppsi = fopen("BSEpsi.par" , "r");
  fread (&psi[0],sizeof(zomplex),n_qp*ist.nspinngrid,ppsi);
  /*char filename[20];
  for (j = 0; j <= n_qp; j++){ 
    printf("Printing state %d from BSEpsi.par\n", j);

    for (i = 0;i<ist.ngrid; i++){
      rho[i] = sqr(psi[j*ist.nspinngrid + i].re)+sqr(psi[j*ist.nspinngrid + i].im);
              
    }
    sprintf(filename, "rhoUp%i.cube", j);
    write_cube_file(rho, &grid, filename);
    
    // for (i = 0;i<ist.ngrid; i++){
    //   rho[i]=sqr(psi[ist.ngrid+i].re)+sqr(psi[ist.ngrid+i].im);
    // }
    // sprintf(filename, "rhoDn%i.cube", j);
    // write_cube_file(rho, &grid, filename);

    // for (i = 0;i<ist.ngrid; i++){
    //   rho[i]= sqr(psi[i].re)+sqr(psi[i].im)+
    //           sqr(psi[ist.ngrid+i].re)+sqr(psi[ist.ngrid+i].im);
    // }
    // sprintf(filename, "rhoTot%i.cube", j);
    // write_cube_file(rho, &grid, filename);

  }
  fclose(ppsi); */

  
  if (ms2 < numExcStatesToPrint) numExcStatesToPrint = ms2;
  pf = fopen("BSEcoeff.dat", "r");
  for (i = 0; i < ms2; i++) {
    printf("Reading exciton %ld of the BSE coeffs\n", i); fflush(0);
    for (j = 0; j < numExcStatesToPrint; j++) {  
	  fscanf (pf,"%lg, %lg\n", &coeff_bs[i*ms2+j].re, &coeff_bs[i*ms2+j].im);
      if (i == 0) printf ("%lg %lg\n", coeff_bs[i*ms2+j].re, coeff_bs[i*ms2+j].im);
    } 
    fscanf (pf,"");	
  }
  fclose(pf);

  printf("Successfully read the BSE coeffs\n");
  printf("ist.ngrid = %ld\n", ist.ngrid);
//   for (ibs = 0, a = nlumo; a < nlumo + totallumo; a++){
//     for (i = 0; i < nhomo; i++, ibs++){
//         printf("xton %d i = %ld a = %ld C_ai = %lg %lg\n", 0, i, a, coeff_bs[ibs*ms2 + 0].re, coeff_bs[ibs*ms2 + 0].im);
//     }
//   }
  

  long num_fpd_states = 6; // the number of excitons that we will print fixed-point densities for.
  // Print the electron carrier density
  // |Psi_n(r_e)|^2 = sum_{iabs} C_{ia}^n* C_{ib} phi_{a,\sigma}^* * phi_{b,\sigma}
  // where sigma is the spin variable
  #pragma omp parallel for private(n)
  for (n = start; n < end; n++){
    long i, a, b, ia_idx, ib_idx, igrid, ispin, ispingrid;
    double *rho;
    char str[100];
    zomplex *elec_fpd, * elec_dens;

    elec_dens = (zomplex*) calloc(ist.ngrid, sizeof(elec_dens[0]));
    elec_fpd = (zomplex*) calloc(ist.ngrid, sizeof(elec_fpd[0]));
    rho = (double *) calloc(ist.ngrid, sizeof(rho[0]));

    // double rhoa = 0.0;
    // for (igrid=0; igrid < ist.ngrid; igrid++){
    //   elec_dens[igrid].re = elec_dens[igrid].im = 0.0;
    //   rhoa += sqr(psi[26*ist.nspinngrid+igrid].re) + sqr(psi[26*ist.nspinngrid+igrid].im) + sqr(psi[26*ist.nspinngrid+ist.ngrid+igrid].re) + sqr(psi[26*ist.nspinngrid+ist.ngrid+igrid].im);
    // }
    // rhoa *= par.dv;
    // //printf("Norm of elec qp state %ld = %lg\n", 26, rhoa);

    for (i = 0; i <= nhomo; i++){
      // Loop over electron states a
      for (a = nlumo; a < nlumo+totallumo; a++) {
        //loop over electron states b
        for (b = nlumo; b < nlumo+totallumo; b++) {
          ia_idx = listibs[(a - nlumo) * totalhomo + i];
          ib_idx = listibs[(b - nlumo) * totalhomo + i];
          if (a==b) printf("i: %ld a: %ld b: %ld C.re = %lg C.im = %lg\n", i, a, b, coeff_bs[ia_idx*ms2 + n].re * coeff_bs[ib_idx*ms2 + n].re + coeff_bs[ia_idx*ms2 + n].im * coeff_bs[ib_idx*ms2 + n].im, coeff_bs[ia_idx*ms2 + n].re * coeff_bs[ib_idx*ms2 + n].im - coeff_bs[ia_idx*ms2 + n].im * coeff_bs[ib_idx*ms2 + n].re); 
          //printf("i: %ld a: %ld ia_idx: %ld C_ia^n: %lg %lg\n", i, a, ia_idx, coeff_bs[ia_idx*ms2 + n].re, coeff_bs[ia_idx*ms2 + n].im);
          //printf("  i: %ld b: %ld ib_idx: %ld C_ib^n: %lg %lg\n", i, b, ib_idx, coeff_bs[ib_idx*ms2 + n].re, coeff_bs[ib_idx*ms2 + n].im);
          //get joint density \rho_{ab}(r) = \sum_{\sigma} psi_{a}^{*}(r,\sigma) psi_{b}(r,\sigma)
          for (igrid = 0; igrid < ist.ngrid; igrid++) {
            for (ispin = 0; ispin < 2; ispin++){
              ispingrid=ist.ngrid*ispin + igrid;
              elec_dens[igrid].re += psi[a*ist.nspinngrid+ispingrid].re * psi[b*ist.nspinngrid+ispingrid].re
                  + psi[a*ist.nspinngrid+ispingrid].im * psi[b*ist.nspinngrid+ispingrid].im;
              elec_dens[igrid].im += psi[a*ist.nspinngrid+ispingrid].re * psi[b*ist.nspinngrid+ispingrid].im
                  - psi[a*ist.nspinngrid+ispingrid].im * psi[b*ist.nspinngrid+ispingrid].re;
            }
            // Multiply the density by the Bethe-Salpether coefficients for these excitations
            // C_{ia}^n = coeff_bs[ia_idx*ms2 + n];   C_{ib}^n = coeff_bs[ib_idx*ms2 + n]
            elec_fpd[igrid].re += (coeff_bs[ia_idx*ms2 + n].re * coeff_bs[ib_idx*ms2 + n].re + coeff_bs[ia_idx*ms2 + n].im * coeff_bs[ib_idx*ms2 + n].im ) * elec_dens[igrid].re\
                      - (coeff_bs[ia_idx*ms2 + n].re * coeff_bs[ib_idx*ms2 + n].im - coeff_bs[ia_idx*ms2 + n].im * coeff_bs[ib_idx*ms2 + n].re) * elec_dens[igrid].im;
                     
            elec_fpd[igrid].im *= (coeff_bs[ia_idx*ms2 + n].re * coeff_bs[ib_idx*ms2 + n].re + coeff_bs[ia_idx*ms2 + n].im * coeff_bs[ib_idx*ms2 + n].im)* elec_dens[igrid].im\
                      + (coeff_bs[ia_idx*ms2 + n].re*coeff_bs[ib_idx*ms2 + n].im - coeff_bs[ia_idx*ms2 + n].im * coeff_bs[ib_idx*ms2 + n].re) * elec_dens[igrid].re;
            //if (0 == (igrid % 1000000)) printf("\telec_fpd[%ld] = %lg %lg\n", igrid, elec_fpd[igrid].re, elec_fpd[igrid].im ); 
          }
        }
      }
    }
    //printf("\nAbout to print xton density for state %ld\n", n);
    // Print the electron carrier density
    //
    // The carrier density is complex. Take the square magnitude in order to obtain
    // a single real-valued density at all points in space
    double sum = 0.0;
    for (igrid = 0; igrid < ist.ngrid; igrid++){
      rho[igrid] = elec_fpd[igrid].re;
      sum += elec_fpd[igrid].re + elec_fpd[igrid].im;
    }
    sum *= par.dv;
    printf("Norm of elec carrier dens of xton %ld = %lg\n", n, sum);
    
    // Write the cube file
    sprintf(str, "xton-%ld-elec-fpd.cube", n);
    write_cube_file(rho, &grid, str);

    free(elec_fpd); free(elec_dens); free(rho);
  }

  // Print the hole carrier density
  // |Psi_n(r_h)|^2 = sum_{ija,sigma} C_{ia}^n* C_{ja} phi_{i,\sigma}^* * phi_{j,\sigma}
  // where sigma is the spin variable
  #pragma omp parallel for private(n)
  for (n = start; n < end; n++){
    long i, j, a, ia_idx, ja_idx, igrid, ispin, ispingrid, dn_grid;
    double *rho_up, *rho_dn;
    double c_ia_norm = 0.0;
    char str[50], str_dn[50];
    zomplex *hole_fpd, *hole_dens;

    hole_fpd = (zomplex*) calloc(ist.ngrid, sizeof(hole_fpd[0]));
    hole_dens = (zomplex*) calloc(ist.ngrid, sizeof(hole_dens[0]));
    rho_up = (double *) calloc(ist.ngrid, sizeof(rho_up[0]));
    rho_dn = (double *) calloc(ist.ngrid, sizeof(rho_dn[0]));
    for (igrid=0; igrid < ist.ngrid; igrid++){
      hole_dens[igrid].re = hole_dens[igrid].im = 0.0;
    }
    for (i = 0; i <= nhomo; i++){
      sprintf(str, "i-state-%ld-up.cube", i);
      sprintf(str_dn, "i-state-%ld-dn.cube", i);
      for (igrid = 0; igrid < ist.ngrid; igrid++){
        dn_grid = ist.ngrid + igrid;
        rho_up[igrid] = sqr(psi[i*ist.nspinngrid + igrid].re) + sqr(psi[i*ist.nspinngrid + igrid].re);
        rho_dn[igrid] = sqr(psi[i*ist.nspinngrid + dn_grid].re) + sqr(psi[i*ist.nspinngrid + dn_grid].re);
      }
      write_cube_file(rho_up, &grid, str);
      write_cube_file(rho_dn, &grid, str_dn);
      // Loop over electron states a
      for (j = 0; j <= nhomo; j++) {
        //printf("i: %ld a: %ld ia_idx: %ld C_ia^n: %lg %lg\n", i, a, ia_idx, coeff_bs[ia_idx*ms2 + n].re, coeff_bs[ia_idx*ms2 + n].im);
        //loop over electron states b
        for (a = nlumo; a < nlumo+totallumo; a++) {
          ia_idx = listibs[(a - nlumo) * totalhomo + i];
          ja_idx = listibs[(a - nlumo) * totalhomo + j];
          if (i == j) printf("i: %ld j: %ld a: %ld C.re = %lg C.im = %lg\n", i, j, a, coeff_bs[ia_idx*ms2 + n].re * coeff_bs[ja_idx*ms2 + n].re + coeff_bs[ia_idx*ms2 + n].im * coeff_bs[ja_idx*ms2 + n].im, coeff_bs[ia_idx*ms2 + n].re * coeff_bs[ja_idx*ms2 + n].im - coeff_bs[ia_idx*ms2 + n].im * coeff_bs[ja_idx*ms2 + n].re);
          //printf("  j: %ld a: %ld ja_idx: %ld C_ja^n: %lg %lg\n", j, a, ja_idx, coeff_bs[ja_idx*ms2 + n].re, coeff_bs[ja_idx*ms2 + n].im); fflush(0);
          //get joint density \rho_{ab}(r) = \sum_{\sigma} psi_{a}^{*}(r,\sigma) psi_{b}(r,\sigma)
          for (igrid = 0; igrid < ist.ngrid; igrid++) {
            for (ispin = 0; ispin < 2; ispin++){
              ispingrid=ist.ngrid*ispin + igrid;
              hole_dens[igrid].re += psi[i*ist.nspinngrid+ispingrid].re * psi[j*ist.nspinngrid+ispingrid].re
                  + psi[i*ist.nspinngrid+ispingrid].im * psi[j*ist.nspinngrid+ispingrid].im;
              hole_dens[igrid].im += psi[i*ist.nspinngrid+ispingrid].re * psi[j*ist.nspinngrid+ispingrid].im
                  - psi[i*ist.nspinngrid+ispingrid].im * psi[j*ist.nspinngrid+ispingrid].re;
            }
            // Multiply the density by the Bethe-Salpether coefficients for these excitations
            // C_{ia}^n = coeff_bs[ia_idx*ms2 + n];   C_{ib}^n = coeff_bs[ib_idx*ms2 + n]
            hole_fpd[igrid].re += (coeff_bs[ia_idx*ms2 + n].re * coeff_bs[ja_idx*ms2 + n].re + coeff_bs[ia_idx*ms2 + n].im * coeff_bs[ja_idx*ms2 + n].im ) * hole_dens[igrid].re\
                      - (coeff_bs[ia_idx*ms2 + n].re * coeff_bs[ja_idx*ms2 + n].im - coeff_bs[ia_idx*ms2 + n].im * coeff_bs[ja_idx*ms2 + n].re) * hole_dens[igrid].im;
            
            hole_fpd[igrid].im *= (coeff_bs[ia_idx*ms2 + n].re * coeff_bs[ja_idx*ms2 + n].re + coeff_bs[ia_idx*ms2 + n].im * coeff_bs[ja_idx*ms2 + n].im)* hole_dens[igrid].im\
                      + (coeff_bs[ia_idx*ms2 + n].re*coeff_bs[ja_idx*ms2 + n].im - coeff_bs[ia_idx*ms2 + n].im * coeff_bs[ja_idx*ms2 + n].re) * hole_dens[igrid].re;
            //if (0 == (igrid % 100000)) printf("\thole_fpd[%ld] = %lg %lg\n", igrid, hole_fpd[igrid].re, hole_fpd[igrid].im ); fflush(0);
          }
          c_ia_norm += sqr(coeff_bs[ia_idx*ms2 + n].re) + sqr(coeff_bs[ia_idx*ms2 + n].im);
        }
      }
    }
    
    // Print the electron carrier density
    //

    // The carrier density is complex. Take the square magnitude in order to obtain
    // a single real-valued density at all points in space
    double sum = 0.0;
    for (igrid = 0; igrid < ist.ngrid; igrid++){
      rho[igrid] = hole_fpd[igrid].re;
      sum += hole_fpd[igrid].re + hole_fpd[igrid].im;
    }
    sum *= par.dv;
    printf("Norm of hole carrier dens of xton %ld = %lg\n", n, sum);
    printf("Norm of C_ia coeffs for xton %ld = %lg\n", n, c_ia_norm);
    // Write the cube file
    sprintf(str, "xton-%ld-hole-fpd.cube", n);
    write_cube_file(rho, &grid, str);

    free(hole_fpd); free(hole_dens); free(rho);
  }

  return 0;


}


/*****************************************************************************/
int countlines(char *filename){
  FILE* fp = fopen(filename,"r");
  int lines = 0;
  int ch;
  while(1){
    ch = fgetc(fp);
    if (feof(fp)){ break; }
    if(ch == '\n')
    {
      lines++;
    }
  }
  fclose(fp);
  return lines;
}


