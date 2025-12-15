#include "fd.h"

char *get_time();
char *format_duration(double elapsed_seconds);
void allocate_memory(void **ptr, size_t length, size_t type_size, char *message);
void allocate_aligned_memory(void **ptr, size_t length, size_t type_size, char *message);
void print_progress_bar(int cur, int tot);
int get_dynamic_process_workload(long trip_count);
void real_to_complex(const double *in, double complex *out, size_t n);