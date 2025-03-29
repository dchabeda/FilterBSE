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
