#!/bin/bash

# Check if both input files are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <FILE_HIGH> <FILE_LOW> > <OUT>" >&2
    exit 1
fi

FILE_HIGH="$1"
FILE_LOW="$2"

# --- Header (Line 1: Sample IDs) Processing ---

# Get just the first line of sample IDs from both files
samples_high=$(head -n 1 "$FILE_HIGH")
samples_low=$(head -n 1 "$FILE_LOW")

# Combine the sample ID lines and prepend a tab for the empty ID header spot.
# This prints the header line to STDOUT.
echo -e "\t${samples_high}\t${samples_low}"

# --- Data (Line 2 onwards: ID + Data) Processing ---

# Extract data starting from the second line of both files.

# The first file keeps its ID column and all data columns.
data_high_pipe=$(tail -n +2 "$FILE_HIGH")

# The second file has its first column (ID) removed using 'cut -f 2-'.
data_low_stripped_pipe=$(tail -n +2 "$FILE_LOW" | cut -f 2-)

# Paste them together column-wise using process substitution <()
# and print the result to STDOUT.
paste <(echo "$data_high_pipe") <(echo "$data_low_stripped_pipe")
