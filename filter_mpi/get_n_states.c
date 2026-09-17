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

  char fileName[1024];
  char evalInName[1024];
  char evalOutName[1024];

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

  // Generate a clipped eval.dat holding only the selected window of states.
  // eval.dat is a plain-text companion to psi.dat: one line per state,
  // "index eigenvalue variance" (see write_eval_dat in write.c). It lives in
  // the same directory as the input psi file, so derive its path from argv[5].
  {
    FILE *pf_eval_in;
    FILE *pf_eval_out;
    char *slash;
    long idx;
    double eval_loc, sigma_loc;
    long written;

    strcpy(evalInName, argv[5]);
    slash = strrchr(evalInName, '/');
    if (slash != NULL){
      // keep the directory portion (including trailing '/'), append eval.dat
      strcpy(slash + 1, "eval.dat");
    }
    else {
      strcpy(evalInName, "eval.dat");
    }

    pf_eval_in = fopen(evalInName, "r");
    if (pf_eval_in == NULL){
      printf("Warning: could not open %s; skipping clipped eval.dat\n", evalInName);
    }
    else {
      sprintf(evalOutName, "eval_%d-%d.dat", start, end);
      pf_eval_out = fopen(evalOutName, "w");
      if (pf_eval_out == NULL){
        printf("ERROR opening %s for writing\n", evalOutName);
        exit(EXIT_FAILURE);
      }

      idx = 0;
      written = 0;
      // Read every line; copy only those in [start, end], renumbering the
      // index from 0 so line i matches state i of the clipped psi file.
      while (fscanf(pf_eval_in, "%ld %lg %lg", &idx, &eval_loc, &sigma_loc) == 3){
        if (idx >= start && idx <= end){
          fprintf(pf_eval_out, "%ld %.16lg %lg\n", written, eval_loc, sigma_loc);
          written++;
        }
      }

      fclose(pf_eval_in);
      fclose(pf_eval_out);
      printf("Wrote %s (%ld states)\n", evalOutName, written);
    }
  }

  printf("Done with get_n_states.x\n");


  return 0;
}