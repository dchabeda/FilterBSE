/*****************************************************************************/
// Main file for cube printing utility.
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

/*****************************************************************************/
int main(int argc, char *argv[]){
  
  FILE *pf_in;
  FILE *pf_out;

  double *psi;

  int j;
  int start, end;
  int is_cplx;
  
  long nspinngrid;
  long offset;
  
  char fileName[50];

  //command line input parsing
  if (argc!=6){
    printf("Usage: get_n_states start end nspinngrid cmplx filename");
    exit(EXIT_FAILURE);
  }

  start = atoi(argv[1]);
  end = atoi(argv[2]);
  nspinngrid = atol(argv[3]);
  is_cplx = atoi(argv[4]);
  strcpy(fileName, argv[5]);

  printf("start = %d\n", start);
  printf("end = %d\n", end);
  printf("nspinngrid = %ld\n", nspinngrid);
  printf("is_cplx = %d\n", is_cplx);
  printf("fileName = %s\n\n", fileName);

  if (start > end){
    printf("Invaid start (%d), end(%d): start > end\n", start,end);
    exit(EXIT_FAILURE);
  }
  if (start < 0){
    printf("Invaid start (%d): start < 0\n", start);
    exit(EXIT_FAILURE);
  }

  fflush(stdout);

  //allocate memory for psi
  if ((psi = (double *) calloc(is_cplx * nspinngrid, sizeof(double))) == NULL){
    printf("ERROR allocating psi\n"); exit(EXIT_FAILURE);
  }

  //read psi from file
	pf_in = fopen(fileName , "r");
  
  sprintf(fileName, "psi_%d-%d.dat", start, end);
  pf_out = fopen(fileName, "w");
  
  for (j = start; j <= end; j++){
    printf("Reading state %d\n", j);
    
    offset = j * is_cplx * nspinngrid * sizeof(double);

    if(fseek(pf_in, offset, SEEK_SET) != 0){
      printf("Error reading from psi.dat!\n"); exit(EXIT_FAILURE);
    }

    fread (&psi[0], sizeof(double), is_cplx * nspinngrid, pf_in);
    fwrite(&psi[0], sizeof(double), is_cplx * nspinngrid, pf_out);
  }

  fclose(pf_in);
  fclose(pf_out);  

  printf("Done with get_n_states.x\n");


  return 0;
}