#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

#define BUFFER_SIZE (64 * 1024 * 1024)  // 64MB buffer for efficient file copying

void concatenate_files(int jr, int r_s, int r_e, const char *output_file) {
    FILE *out_fp = fopen(output_file, "wb");
    if (!out_fp) {
        fprintf(stderr, "Error: Cannot open output file %s: %s\n", output_file, strerror(errno));
        exit(EXIT_FAILURE);
    }

    // Determine reference file size
    char ref_file[256];
    snprintf(ref_file, sizeof(ref_file), "psi-filt-0-%d.dat", r_s);

    struct stat ref_stat;
    if (stat(ref_file, &ref_stat) != 0) {
        fprintf(stderr, "Error: Reference file %s does not exist.\n", ref_file);
        fclose(out_fp);
        exit(EXIT_FAILURE);
    }
    off_t ref_size = ref_stat.st_size;

    char file_name[256];
    char buffer[BUFFER_SIZE];
    int count = 0;
    FILE *concat_log = fopen("concat_psi_filt.out", "w");

    if (!concat_log) {
        fprintf(stderr, "Error: Cannot open log file.\n");
        fclose(out_fp);
        exit(EXIT_FAILURE);
    }

    // Iterate over files and concatenate if they match the reference size
    for (int i = 0; i < jr; i++) {
        for (int j = r_s; j < r_e; j++) {
            snprintf(file_name, sizeof(file_name), "psi-filt-%d-%d.dat", i, j);

            struct stat file_stat;
            if (stat(file_name, &file_stat) != 0) {
                fprintf(stderr, "Warning: %s does not exist!\n", file_name);
                continue;
            }

            if (file_stat.st_size != ref_size) {
                fprintf(stderr, "Warning: %s has a different size (%ld bytes) than reference (%ld bytes)!\n",
                        file_name, file_stat.st_size, ref_size);
                continue;
            }

            FILE *in_fp = fopen(file_name, "rb");
            if (!in_fp) {
                fprintf(stderr, "Error: Cannot open %s: %s\n", file_name, strerror(errno));
                continue;
            }

            size_t bytes_read;
            while ((bytes_read = fread(buffer, 1, BUFFER_SIZE, in_fp)) > 0) {
                fwrite(buffer, 1, bytes_read, out_fp);
            }

            fclose(in_fp);
            remove(file_name); // Delete successfully processed file

            printf("Concatenated %s\n", file_name);
            fprintf(concat_log, "Concatenated %s\n", file_name);
            fflush(concat_log);

            count++;
        }
    }

    fprintf(concat_log, "%d files were successfully concatenated into %s\n", count, output_file);
    printf("%d files were successfully concatenated into %s\n", count, output_file);

    fclose(out_fp);
    fclose(concat_log);
}

int main(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <j_per_rank> <rnk_start> <rank_end> <output_file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    int jr = atoi(argv[1]);
    int r_s = atoi(argv[2]);
    int r_e = atoi(argv[3]);
    const char *output_file = argv[4];

    concatenate_files(jr, r_s, r_e, output_file);
    return EXIT_SUCCESS;
}
