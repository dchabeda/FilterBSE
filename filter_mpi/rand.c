#include "fd.h"

/*****************************************************************************/

#define IM1 2147483563
#define IM2 2147483399
#define AM  (1.0 / (double)IM1)
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1 + IMM1 / NTAB)
#define REPS 1.2e-7
#define RNMX (1.0 - REPS)

/*****************************************************************************/

double ran_nrc(long *rand_seed)
{
  long    j;
  long   k;
  static long randSeed2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  double temp;
  
  if (*rand_seed <= 0){  /* Initialize */
    if (-(*rand_seed) < 1) *rand_seed = 1;    /* Be sure to prevent rand_seed = 0 */
    else *rand_seed = -(*rand_seed);
    randSeed2 = (*rand_seed);

    for (j = NTAB+7; j >= 0; j--){/* Load the shuffle table after 8 warm-ups */
      k = (*rand_seed) / IQ1;
      *rand_seed = IA1 * (*rand_seed - k * IQ1) - k * IR1;
      if (*rand_seed < 0) *rand_seed += IM1;
      if (j < NTAB) iv[j] = *rand_seed;
    }
    iy = iv[0];
  }

  k = (*rand_seed) / IQ1;      /* Start Here when not initializing */
  *rand_seed = IA1 * (*rand_seed - k * IQ1) - k * IR1; /* Compute rand_seed = IA1*rand_seed % IM1*/
  if (*rand_seed < 0) *rand_seed += IM1;               /* without overflow by Scharge  */

  k = randSeed2 / IQ2;                           /* method.                      */
  randSeed2 = IA2 * (randSeed2 - k * IQ2) - k * IR2; /* Compute randSeed2=IA2*rand_seed % IM2 */
  if (randSeed2 < 0) randSeed2 += IM2;               /* likewise.                    */
  
  j = iy / NDIV;          /* Will be in the range 0..NTAB-1 */
  iy = iv[j] - randSeed2;     /* Here rand_seed is shuffled, rand_seed and randSeed2 are */
  iv[j] = *rand_seed;          /* combined to generate output               */

  if (iy < 1) iy += IMM1;
  if ((temp = AM * iy) > RNMX) return RNMX; /* Because users don't expect */
  else return temp;                         /* endpoint value             */
}

/*****************************************************************************/

void Randomize()
{
  double tt = (double)time(0);
  long   rn = (long)tt % 1000;
  long    i;
  
  for (i = 0; i < rn; i++) random();
  return;
}

/*****************************************************************************/
