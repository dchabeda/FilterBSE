#!/bin/bash

# Check if the user provided the maximum value for the loop
if [[ -z "$1" ]]; then
    echo "Usage: $0 <max_value>"
    exit 1
fi

max_value=$1

# Get the size of the first file as a reference
ref_size=$(stat --format="%s" psi-filt-0.dat 2>/dev/null)

# Check if the reference file exists
if [[ -z "$ref_size" ]]; then
    echo "Reference file psi-filt-0.dat does not exist."
    exit 1
fi

# Loop through all the files up to the max_value and compare their sizes
for i in $(seq 0 $max_value); do
    file="psi-filt-$i.dat"

    # Check if the file exists
    if [[ -f "$file" ]]; then
        size=$(stat --format="%s" "$file")
        if [[ "$size" -ne "$ref_size" ]]; then
            echo "File $file has a different size: $size bytes (expected $ref_size bytes)"
        fi
    fi
done

echo "Size check complete."
