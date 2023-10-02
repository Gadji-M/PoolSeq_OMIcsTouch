#!/bin/bash

# Default options
samples_dir=""

# Function to show usage information
usage() {
    echo "Usage: $(basename "$0") -i SAMPLES_DIR"
    echo "Options:"
    echo "  -i SAMPLES_DIR  Directory containing sample subdirectories"
    exit 1
}

# Function to process a FASTA file
process_fasta() {
    local fasta_file="$1"
    local sample_name="$2"

    # Check if the file exists and has a .fa extension
    if [ -f "$fasta_file" ] && [[ "$fasta_file" == *.fa ]]; then
        # Run TransDecoder.LongOrfs on the input FASTA file
        TransDecoder.LongOrfs -t "$fasta_file" -m 100

        echo "TransDecoder.LongOrfs completed for sample: $sample_name (FASTA file: $fasta_file)"
    fi
}

# Parse command-line options
while getopts "i:" opt; do
    case $opt in
        i)
            samples_dir="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
    esac
done

# Check if the required option is provided
if [ -z "$samples_dir" ]; then
    echo "Error: Missing required option."
    usage
fi

# Loop through each sample subdirectory
for sample_dir in "$samples_dir"/*; do
    if [ -d "$sample_dir" ]; then
        sample_name=$(basename "$sample_dir")  # Extract the sample name from the directory name

        # Loop through each file in the sample subdirectory
        for fasta_file in "$sample_dir"/*.fa; do
            process_fasta "$fasta_file" "$sample_name"
        done
    fi
done
