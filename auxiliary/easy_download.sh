#!/bin/bash

# Define paths
SCRATCH_DIR="$SCRATCH/filter_bse"  # Change this if needed
TARGET_DIR="$SCRATCH_DIR/easy_download"

# Ensure easy_download directory exists
mkdir -p "$TARGET_DIR"

# Loop over job directories matching the naming pattern
for job_dir in "$SCRATCH_DIR"/{r,o,c}_*x*x*_CsPbI3; do
    # Extract just the directory name
    job_name=$(basename "$job_dir")
    
    # Create the corresponding directory in easy_download
    mkdir -p "$TARGET_DIR/$job_name"

    # Loop over "bse" and "filter" directories inside the job directory
    for subdir in "bse" "filter"; do
        src="$job_dir/$subdir"
        dest="$TARGET_DIR/$job_name/$subdir"

        # Check if the subdirectory exists before copying
        if [[ -d "$src" ]]; then
            mkdir -p "$dest"

            # Copy only the required files
            cp "$src"/run*dat "src"/exci* "src"/input* "$src"/conf* "$src"/OS* "$src"/*eval* "$dest" 2>/dev/null
        fi
    done
done

echo "Lite versions of bse and filter directories created in $TARGET_DIR."
