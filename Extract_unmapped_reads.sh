#!/bin/bash

# Default values
input_dir=""
output_dir=""
threads=5

# Function to print script usage
print_usage() {
    echo "Usage: $0 -i <input_directory> -o <output_directory> [-t <threads>]"
    echo "Options:"
    echo "  -i <input_directory>: Directory containing BAM files"
    echo "  -o <output_directory>: Directory to store output FastQ files"
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

# Use find to locate BAM files in the input directory (excluding subdirectories)
find "$input_dir" -maxdepth 1 -type f -name "*.bam" | while read -r input_bam; do
    # Extract the filename (without extension)
    filename=$(basename "$input_bam")
    filename="${filename%.*}"

    # Set the output FastQ file path
    output_fastq="$output_dir/$filename.fastq"

    # Use samtools to extract unmapped reads in parallel
    samtools view -@ "$threads" -b -f 4 "$input_bam" | samtools fastq -@ "$threads" - > "$output_fastq"

    echo "Unmapped reads extraction completed for $filename."
done

echo "All BAM files processed."
