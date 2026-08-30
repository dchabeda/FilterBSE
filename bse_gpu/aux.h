#include "fd.h"

void allocate_memory(void **ptr, size_t length, size_t type_size, char* message);
void allocate_aligned_memory(void **ptr, size_t length, size_t type_size, char* message);
void print_progress_bar(int cur, int tot);
char* get_time();

