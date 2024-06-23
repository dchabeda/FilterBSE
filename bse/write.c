/****************************************************************************/
/* This file does the printing for the program */

#include "fd.h"

/****************************************************************************/
// prints out current time to stdout 
 
void write_current_time(FILE *pf) {
  time_t startTime;

  startTime = time(NULL);
  fprintf(pf,"%s", ctime(&startTime));

  return;
}

/****************************************************************************/

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
    fprintf(stderr, "Invalid string supplied to write_separation. Exiting!\n");
    exit(EXIT_FAILURE);
  }

  free(top_key); free(bttm_key);
  
  return;
}

/****************************************************************************/
