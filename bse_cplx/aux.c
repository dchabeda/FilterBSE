#include "aux.h"

void print_progress_bar(int cur, int tot){
  // print the filtering progress to the output file
  int barWidth = 16; // Width of the progress bar
  float percent = (float)cur / tot * 100;
  int pos = barWidth * cur / tot;

  printf("\t  [");
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) printf("#");
    else printf(" ");
  }
  printf("] %3.0f%% | %s\n", percent, get_time());
  fflush(stdout);

  return;
}

/*****************************************************************************/

char* get_time(void){
  time_t current_time;
  char* c_time_string;

  current_time = time(NULL);
  c_time_string = ctime(&current_time);
  
  if (current_time == ((time_t)-1)) {
    fprintf(stderr, "Error: time() failed\n");
    return "UNKNOWN_TIME";
  }

  if (c_time_string == NULL) {
    fprintf(stderr, "Error: ctime() failed\n");
    return "UNKNOWN_TIME";
  }

  return c_time_string;
}

/*****************************************************************************/

void allocate_memory(void **ptr, size_t length, size_t type_size, char* message) {
  *ptr = calloc(length, type_size);
  if (*ptr == NULL) {
    fprintf(stderr, "Memory allocation failed: %s\n", message);
    exit(EXIT_FAILURE);
  }
  return;
}


/*****************************************************************************/
void allocate_aligned_memory(void **ptr, size_t length, size_t type_size, char* message) {
  // Align mem to 32 byte boundaries for AVX2 SIMD operations
  *ptr = aligned_alloc(BYTE_BOUNDARY, length * type_size); 
  if (*ptr == NULL) {
      fprintf(stderr, "Memory allocation failed: %s\n", message);
      exit(EXIT_FAILURE);
  }
  memset(*ptr, 0, length * type_size);  // Initialize all elements to 0
  return;
}

int get_dynamic_process_workload(long trip_count){
  int num_threads;
  int chunk_size;
  
  // Get the number of threads available
  #pragma omp parallel
  {
      #pragma omp master
      num_threads = omp_get_num_threads();
  }
  
  // Adjust chunk size based on problem size and thread count
  if (trip_count < 128) {
    // For very small workloads, use sequential execution
    omp_set_num_threads(1);
    chunk_size = trip_count;
  } else if (trip_count < 4 * num_threads) {
      // For small workloads, use smaller number of threads
      int adjusted_threads = (trip_count + 3) / 4; // Ceiling division
      omp_set_num_threads(adjusted_threads);
      chunk_size = 1; // Use smallest chunk size for best load balancing
  } else if (trip_count < 16 * num_threads) {
      // For medium workloads, use all threads but with small chunks
      chunk_size = 1;
  } else {
      // For larger workloads, calculate optimal chunk size
      // Aim for ~8-16 chunks per thread for good load balancing with lower overhead
      int target_chunks_per_thread = 12;
      chunk_size = trip_count / (num_threads * target_chunks_per_thread);
      
      // Set reasonable bounds
      if (chunk_size < 4) chunk_size = 4;
      if (chunk_size > 64) chunk_size = 64;
  }

  return chunk_size;
}

