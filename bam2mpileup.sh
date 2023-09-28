#!/bin/bash

# Default values
reference_genome=""
output_dir=""
output_prefix=""

# Function to show usage information
usage() {
    echo "Usage: $(basename "$0") -r REFERENCE_GENOME [-o OUTPUT_DIR] [-p OUTPUT_PREFIX] BAM_FILE [BAM_FILE ...]"
    echo "Options:"
    echo "  -r REFERENCE_GENOME   Path to the reference genome file (required)"
    echo "  -o OUTPUT_DIR         Output directory for MPileup files (default: $output_dir)"
    echo "  -p OUTPUT_PREFIX      Prefix for the combined MPileup file (default: $output_prefix)"
    exit 1
}

# Parse command-line options
while getopts "r:o:p:" opt; do
    case $opt in
        r)
            reference_genome="$OPTARG"
            ;;
        o)
            output_dir="$OPTARG"
            ;;
        p)
            output_prefix="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
    esac
done

# Shift the options and get the remaining arguments (BAM files)
shift $((OPTIND-1))
bam_files=("$@")

# Check if the reference genome file exists
if [ ! -e "$reference_genome" ]; then
    echo "Reference genome file '$reference_genome' does not exist."
    exit 1
fi

# Check if there are any BAM files specified
if [ ${#bam_files[@]} -eq 0 ]; then
    echo "No BAM files specified."
    usage
fi

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Generate the MPileup file for all specified BAM files
samtools mpileup -B -Q 0 -f "$reference_genome" "${bam_files[@]}" > "$output_dir/$output_prefix.mpileup"

echo "MPileup generation complete. Combined MPileup saved as '$output_dir/$output_prefix.mpileup'."
