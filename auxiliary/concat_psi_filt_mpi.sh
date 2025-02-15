#!/bin/bash

# Check if both arguments (nstates and output_file) are provided
if [ $# -ne 3 ]; then
    echo "Usage: $0 <j_per_rank> <n_rank> <output_file>"
    exit 1
fi

# Get the number of states and the output file name
jr=$1
nr=$2
output_file=$3

# Get the size of the first file as a reference
ref_file="psi-filt-0-0.dat"
ref_size=$(stat --format="%s" "$ref_file" 2>/dev/null)

# Check if the reference file exists
if [[ -z "$ref_size" ]]; then
    echo "Reference file $ref_file does not exist."
    exit 1
fi

# Initialize or clear the output file
> $output_file

# Initialize a counter for successfully concatenated files
count=0

# Loop over all the files and concatenate them if their sizes match
for ((i=0; i<jr; i++)); do
for ((j=0; j<nr; j++)); do
    file="psi-filt-${i}-${j}.dat"
    
    # Check if the file exists
    if [[ -f "$file" ]]; then
        size=$(stat --format="%s" "$file")
        
        # Only concatenate if the file size matches the reference
        if [[ "$size" -eq "$ref_size" ]]; then
            cat "$file" >> "$output_file"
            echo "Concatenated $file"
            echo "Concatenated $file" > concat_psi_filt.out
            ((count++))  # Increment the counter for each successfully concatenated file
        else
            echo "Warning: $file has a different size ($size bytes) than the reference ($ref_size bytes)!"
        fi
    else
        echo "Warning: $file does not exist!" 
    fi
done
done
# Print the number of successfully concatenated files
echo "$count files were successfully concatenated into $output_file" > concat_psi_filt.out
echo "$count files were successfully concatenated into $output_file"
