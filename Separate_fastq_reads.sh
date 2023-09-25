#!/bin/bash

# Default values
input_dir=""
output_dir=""
threads=5

# Function to print script usage
print_usage() {
    echo "Usage: $0 -i <input_directory> -o <output_directory> [-t <threads>]"
    echo "Options:"
    echo "  -i <input_directory>: Directory containing FASTQ files"
    echo "  -o <output_directory>: Directory to store paired-end FASTQ files"
    echo "  -t <threads>: Number of threads for parallel processing (default: 1)"
    exit 1
}

# Parse command-line options
while getopts "i:o:t:" opt; do
    case "$opt" in
        i) input_dir="$OPTARG" ;;
        o) output_dir="$OPTARG" ;;
        t) threads="$OPTARG" ;;
        *) print_usage ;;
    esac
done

# Check if required options are provided
if [ -z "$input_dir" ] || [ -z "$output_dir" ]; then
    echo "Error: Input and output directories are required."
    print_usage
fi

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Use find to locate FASTQ files in the input directory
find "$input_dir" -maxdepth 1 -type f -name "*.fastq" | while read -r fastq_file; do
    filename=$(basename "$fastq_file")
    base_filename="${filename%.*}"

    # Set output file paths for the paired-end FASTQ files
    output_file1="$output_dir/${base_filename}_R1.fastq"
    output_file2="$output_dir/${base_filename}_R2.fastq"

    # Use seqtk to split the FASTQ file into paired-end files
    seqtk seq -1 "$fastq_file" > "$output_file1"
    seqtk seq -2 "$fastq_file" > "$output_file2"

    echo "Splitting completed for $filename."
done

echo "All FASTQ files processed and split into paired-end files."


#AUTHOR: GADJI_MAHAMAT
