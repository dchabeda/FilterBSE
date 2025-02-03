/*****************************************************************************/
// Main file for cube printing utility.
#include "fd.h"

/*****************************************************************************/
int countlines(char *filename);


/*****************************************************************************/
int main(int argc, char *argv[])
{
  FILE *ppsi;  double *psi;
  // custom structs 
  grid_st grid; par_st par; index_st ist; parallel_st parallel; flag_st flag;
  xyz_st *R; atom_info *atom;
  // double arrays
  double *rho; 
  // long int arrays and counters
  long i, i_re, i_im, jms;
  long j, start, end;
  long long offset;
  ist.atom_types = malloc(N_MAX_ATOM_TYPES*sizeof(ist.atom_types[0]));
  time_t currentTime = time(NULL);


  //command line input parsing
  if (argc!=3){
    printf("Usage: makecube start end\n");
    exit(EXIT_FAILURE);
  }

  start = atoi(argv[1]);
  end = atoi(argv[2]);

  if (start > end){
    printf("Invaid start (%ld), end(%ld): start > end\n", start, end);
    exit(EXIT_FAILURE);
  }
  if (start < 0){
    printf("Invaid start (%ld): start < 0\n", start);
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
  char *file_name; file_name = malloc(9*sizeof(file_name[0]));
  strcpy(file_name, "conf.par");
  read_conf(file_name, R, atom, &ist, &par, &flag);
  
  /*** initialize the grid ***/
  printf("\nInitializing the grid parameters:\n"); fflush(0);
  init_grid_params(&grid, R, &ist, &par, &flag);

  if ((rho = (double *)calloc(ist.ngrid, sizeof(double)))==NULL) nerror("rho");


  //count number of states found
  //jms = countlines("eval.dat");
  //printf("%ld total states in psi.dat\n", jms); fflush(0);
  
  //allocate memory for psi
  if ((psi = (double *) calloc(ist.complex_idx * ist.nspinngrid, sizeof(double))) == NULL) nerror("psi");


  //read psi from file
	ppsi = fopen("psi.dat" , "r");

	
  char filename[20];
  // long file_ptr_loc;
  for (j = start; j <= end; j++){ 
    printf("Reading state %ld from psi.dat\n", j);
    offset = j*ist.complex_idx*ist.nspinngrid*sizeof(double);
    // printf("File offset = %lld\n", offset);
    if(fseek(ppsi, offset, SEEK_SET) != 0){
      printf("Error reading from psi.dat!\n"); exit(EXIT_FAILURE);
    }
   
    // file_ptr_loc = ftell(ppsi);
    // printf("The file pointer is at %ld\n", file_ptr_loc);
    if (fread(&psi[0], sizeof(double), ist.complex_idx*ist.nspinngrid, ppsi) == 0){
      printf("Error in fread reading from psi.dat!\n"); exit(EXIT_FAILURE);
    }
    
    // file_ptr_loc = ftell(ppsi);
    // printf("After reading the file pointer is at %ld\n", file_ptr_loc);
    
    if (flag.isComplex == 0){
      for (i = 0; i < ist.ngrid; i++){
        rho[i] = sqr(psi[i]);      
      }
    }
    else {
      for (i = 0; i < ist.ngrid; i++){
        i_re = i * ist.complex_idx;
        i_im = i * ist.complex_idx + 1;
        rho[i]+=sqr(psi[i_re]) + sqr(psi[i_im]);
      }
    }
    sprintf(filename, "rhoUp%ld.cube", j);
    write_cube_file(rho, &grid, filename);

    if (flag.useSpinors == 1){
      for (i = 0; i < ist.ngrid; i++){
        i_re = i * ist.complex_idx;
        i_im = i * ist.complex_idx + 1;
        rho[i] += sqr(psi[ist.ngrid + i_re]) + sqr(psi[ist.ngrid + i_im]);
      }
      sprintf(filename, "rhoDn%ld.cube", j);
      write_cube_file(rho, &grid, filename);
    }
  }
  fclose(ppsi);  

  printf("Done with makecube.x\n\n");
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


