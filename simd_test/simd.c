#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <immintrin.h>

#define N 1000000000  // Large enough for performance testing
#define f32_VEC_LEN 8

void init_arrays(float *a, float *b, float *c) {
  #pragma omp parallel for
  for (long i = 0; i < N; i++) {
    a[i] = 1.0f * (i % 100);
    b[i] = 2.0f * (i % 100);
    c[i] = 0.0f;
  }
}

void init_rand_arrays(float *a, float *b) {
  #pragma omp parallel for
  for (long i = 0; i < N; i++) {
    a[i] = (2.0f * (float)rand()/RAND_MAX) - 1.0f;
    b[i] = (2.0f * (float)rand()/RAND_MAX) - 1.0f;
  }
}

void multiply_arrays(float *restrict a, float *restrict b, float *restrict c) {
  for (long i = 0; i < N; i++) {
    c[i] = a[i] * b[i];
  }
}

void multiply_arrays_simd(float *restrict a, float *restrict b, float *restrict c) {
  #pragma omp simd aligned(a, b, c: 32)
  for (long i = 0; i < N; i+=8) {
    c[i] = a[i] * b[i];
  }
}

void multiply_arrays_omp(float *restrict a, float *restrict b, float *restrict c) {
  #pragma omp parallel for aligned(a, b, c: 32)
  for (long i = 0; i < N; i++) {
    c[i] = a[i] * b[i];
  }
}

void complex_dot(float *a, float *b, float *Cr, float *Ci){
  long i;
  float Ar;
  float Ai;
  float Br;
  float Bi;
  float sum_re; 
  float sum_im;

  sum_re = 0.0f;
  sum_im = 0.0f;

  for (i = 0; i < N; i += 2){
    // Real and imaginary components of A
    Ar = a[i];
    Ai = a[i + 1];

    // Real and imaginary components of B
    Br = b[i];
    Bi = - b[i + 1];

    // Complex dot product
    sum_re += Ar * Br - Ai * Bi;
    sum_im += Ar * Bi + Ai * Br;
  }

  *Cr = sum_re;
  *Ci = sum_im;

  return;
}

void complex_dot_simd(float *A, float *B, float *Cr, float *Ci){

  long j;
  __m256 sum_re; 
  __m256 sum_im;

  const int n = N / f32_VEC_LEN;

  // Initialize the real and imaginary components of the summation for the a.b* product
  // Sets all elements of the 256 bit vector to 0
  sum_re = _mm256_set1_ps(0.0);
  sum_im = _mm256_set1_ps(0.0);

  // Create a vector with 1.0 in the real component and -1 in the imag for conjugating
  const __m256 conj = _mm256_set_ps(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0);

  // Alias the float arrays a and b as 256-bit packed single vectors (__m256)
  __m256* a = (__m256*)A;
  __m256* b = (__m256*)B;

  for (j = 0; j < n; j++){
    __m256 cr    = _mm256_mul_ps(a[j], b[j]);    // |a3i*b3i|a3r*b3r|a2i*b2i|a2r*b2r|a1i*b1i|a1r*b1r|a0i*b0i|a0r*b0r|
    __m256 bConj = _mm256_mul_ps(b[j], conj); // Conjugate b |-b3i|b3r|-b2i|b2r|-b1i|b1r|-b0i|b0r|

    // Swap the real and imaginary components of the conjugated b vector with an
    // in-lane [2, 3, 0, 1] swap. Binary control byte is read in blocks of two: 10(2), 11(3), 00(0), 01(1)
    __m256 bFlip = _mm256_permute_ps(bConj, 0b10110001); // |b3r|-b3i|b2r|-b2i|b1r|-b1i|b0r|-b0i| 
    __m256 ci    = _mm256_mul_ps(a[j], bFlip);           // |a3i*b3r|-a3r*b3i|a2i*b2r|-a2r*b2i|a1i*b1r|-a1r*b1i|a0i*b0r|-a0r*b0i|
    sum_re       = _mm256_add_ps(sum_re, cr);
    sum_im       = _mm256_add_ps(sum_im, ci);
  }

  // With fused multiply-add
  // for (j = 0; j < n; j++){
  //   sum_re = _mm256_fmadd_ps(a[j], b[j], sum_re);
  //   __m256 bConj = _mm256_mul_ps(b[j], conj); // Conjugate b |-b3i|b3r|-b2i|b2r|-b1i|b1r|-b0i|b0r|
  //   __m256 bFlip = _mm256_permute_ps(bConj, 0b10110001);
  //   sum_im = _mm256_fmadd_ps(a[j], bFlip, sum_im);
  // }

  // Use horizontal vector additions to add up every element in the vector
  // Real part
  sum_re = _mm256_hadd_ps(sum_re, sum_re);
  sum_re = _mm256_hadd_ps(sum_re, sum_re);
  __m256 sumr_flip = _mm256_permute2f128_ps(sum_re, sum_re, 1);
  sum_re = _mm256_add_ps(sum_re, sumr_flip);
  // Imaginary part
  sum_im = _mm256_hadd_ps(sum_im, sum_im);
  sum_im = _mm256_hadd_ps(sum_im, sum_im);
  __m256 sumi_flip = _mm256_permute2f128_ps(sum_im, sum_im, 1);
  sum_im = _mm256_add_ps(sum_im, sumi_flip);

  *Cr = sum_re[0];
  *Ci = sum_im[0];

  return;
}



int main(int arc, char* argv[]) {
  float *a = (float*) aligned_alloc(64, N * sizeof(float));
  float *b = (float*) aligned_alloc(64, N * sizeof(float));
  float *c = (float*) aligned_alloc(64, N * sizeof(float));

  init_arrays(a, b, c);

  int nthreads = atoi(argv[1]);

  printf("Number of OMP threads set to %d\n", nthreads);
  omp_set_num_threads(nthreads);
  

  double start = omp_get_wtime();
  multiply_arrays(a, b, c);
  double end = omp_get_wtime();
  printf("No parallel or simd: Time taken: %f seconds\n", end - start);
  printf("Select elems of c:\n");
  printf("c[0] = %f\tc[1230] = %f\tc[480247] = %f\n", c[0], c[1230], c[480247]);

  init_arrays(a, b, c);

  start = omp_get_wtime();
  multiply_arrays_omp(a, b, c);
  end = omp_get_wtime();

  printf("Parallel, no simd: Time taken: %f seconds\n", end - start);
  printf("Select elems of c:\n");
  printf("c[0] = %f\tc[1230] = %f\tc[480247] = %f\n", c[0], c[1230], c[480247]);

  init_arrays(a, b, c);

  start = omp_get_wtime();
  multiply_arrays_simd(a, b, c);
  end = omp_get_wtime();

  printf("Parallel + simd: Time taken: %f seconds\n", end - start);
  printf("Select elems of c:\n");
  printf("c[0] = %f\tc[1230] = %f\tc[480247] = %f\n", c[0], c[1230], c[480247]);


  /*************************************************************/
  /*************************************************************/
  // Complex dot product
  /*************************************************************/
  /*************************************************************/
  float Cr;
  float Ci;

  init_rand_arrays(a, b);

  start = omp_get_wtime();
  complex_dot(a, b, &Cr, &Ci);
  end = omp_get_wtime();
  printf("Complex dot: Time taken: %f seconds\n", end - start);
  printf("Cr = %f Ci = %f\n", Cr, Ci);

  start = omp_get_wtime();
  complex_dot_simd(a, b, &Cr, &Ci);
  end = omp_get_wtime();
  printf("Complex dot + simd: Time taken: %f seconds\n", end - start);
  printf("Cr = %f Ci = %f\n", Cr, Ci);
  


  free(a);
  free(b);
  free(c);

  return 0;
}
