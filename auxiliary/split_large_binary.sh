#!/bin/bash

# User-defined parameters
MAX_FILE_LEN_GB=220  # Max file size per split (adjust as needed)
SEGMENT_SIZE_BYTES=$1  # Replace with your segment size
FILE="psi-filt.dat"

# Convert max file length to bytes
MAX_FILE_LEN_BYTES=$(( MAX_FILE_LEN_GB * 1024**3 ))

# Get total file size
FILE_SIZE_BYTES=$(stat --format="%s" "$FILE")

# Calculate the max number of full segments that fit in each split
SEGMENTS_PER_SPLIT=$(( MAX_FILE_LEN_BYTES / SEGMENT_SIZE_BYTES ))
ACTUAL_SPLIT_SIZE=$(( SEGMENTS_PER_SPLIT * SEGMENT_SIZE_BYTES ))

# Calculate the number of full splits
NUM_SEGMENTS_TOTAL=$(( FILE_SIZE_BYTES / SEGMENT_SIZE_BYTES ))
NUM_SPLITS=$(( NUM_SEGMENTS_TOTAL / SEGMENTS_PER_SPLIT ))

# If there's a remainder, add one more split
REMAINDER_SEGMENTS=$(( NUM_SEGMENTS_TOTAL % SEGMENTS_PER_SPLIT ))
if [[ $REMAINDER_SEGMENTS -ne 0 ]]; then
    ((NUM_SPLITS++))
fi

echo "Splitting '$FILE' ($FILE_SIZE_BYTES bytes) into $NUM_SPLITS parts, each containing $SEGMENTS_PER_SPLIT segments (~$ACTUAL_SPLIT_SIZE bytes per split)."

# Perform the splitting
for ((i=0; i<NUM_SPLITS; i++)); do
    #DIRNAME="ortho_holes$((i+1))"
    mkdir -p "$DIRNAME"

    SKIP_SEGMENTS=$(( i * SEGMENTS_PER_SPLIT ))
    COUNT_SEGMENTS=$SEGMENTS_PER_SPLIT

    # Adjust for the last segment
    if [[ $i -eq $((NUM_SPLITS-1)) && $REMAINDER_SEGMENTS -ne 0 ]]; then
        COUNT_SEGMENTS=$REMAINDER_SEGMENTS
    fi

    SPLIT_SIZE=$(( COUNT_SEGMENTS * SEGMENT_SIZE_BYTES ))

    echo "Writing $SPLIT_SIZE bytes to psi-filt$i.dat (Skipping $SKIP_SEGMENTS segments)..."

    dd if="$FILE" bs=$SEGMENT_SIZE_BYTES skip=$SKIP_SEGMENTS count=$COUNT_SEGMENTS of="psi-filt$i.dat" status=progress &
done

wait  # Wait for all background `dd` processes to finish
echo "Splitting complete!"
