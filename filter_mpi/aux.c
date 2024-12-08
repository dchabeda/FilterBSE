#include "fd.h"

/*****************************************************************************/
char* format_duration(double elapsed_seconds) {
    // Calculate hours, minutes, and seconds
    int hours = (int)(elapsed_seconds / 3600);
    int minutes = (int)((elapsed_seconds - (hours * 3600)) / 60);
    int seconds = (int)(elapsed_seconds - (hours * 3600) - (minutes * 60));
    
    // Allocate memory for the string
    char* result = (char*)malloc(16 * sizeof(char)); // "hh:mm:ss" + null terminator

    // Format the string
    if (result != NULL) {
        snprintf(result, 16, "%02d:%02d:%02d", hours, minutes, seconds);
    }

    return result;
}

/*****************************************************************************/
// void mpi_print(const char *message) {
//     int rank;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     if (rank == 0) {
//         printf("%s", message);
//     }
// }