/*****************************************************************************/
// Main file for cube printing utility.
#include "fd.h"

/*****************************************************************************/
int countlines(char *filename);


/*****************************************************************************/
int main(int argc, char *argv[])
{
  FILE *ppsi;  zomplex *psi;
  // custom structs 
  grid_st grid; par_st par; index_st ist; parallel_st parallel; flag_st flag;
  xyz_st *R; atom_info *atom;
  // double arrays
  double *rho; 
  // long int arrays and counters
  long i, jms;
  int j, start, end;
  ist.atom_types = malloc(N_MAX_ATOM_TYPES*sizeof(ist.atom_types[0]));
  time_t currentTime = time(NULL);


  //command line input parsing
  if (argc!=3){
    printf("Usage: makecube start end");
    exit(EXIT_FAILURE);
  }

  start = atoi(argv[1]);
  end = atoi(argv[2]);

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
  
  //allocate memory for psi
  if ((psi = (zomplex *) calloc(ist.nspinngrid, sizeof(zomplex))) == NULL) nerror("psi");


  //read psi from file
	ppsi = fopen("psi.dat" , "r");

	
  char filename[20];
  for (j = start; j <= end; j++){ 
    printf("Reading state %d from psi.dat\n", j);

    if(fseek(ppsi,j*ist.nspinngrid*sizeof(zomplex),SEEK_SET)!=0){
      printf("Error reading from psi.dat!\n"); exit(EXIT_FAILURE);
    }
    fread (&psi[0],sizeof(zomplex),ist.nspinngrid,ppsi);

    for (i = 0;i<ist.ngrid; i++){
      rho[i] = sqr(psi[i].re)+sqr(psi[i].im);
              
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
  fclose(ppsi);  

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


