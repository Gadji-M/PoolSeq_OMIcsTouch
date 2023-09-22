#!/bin/bash

# Default values
bam_dir=""
output_dir="coverage_results"

# Function to print script usage
usage() {
  echo "Usage: $0 -b <bam_dir> [-o <output_dir>]"
  exit 1
}

# Parse command-line options
while getopts ":b:o:" opt; do
  case $opt in
    b)
      bam_dir="$OPTARG"
      ;;
    o)
      output_dir="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      usage
      ;;
  esac
done

# Check if the BAM directory is provided
if [ -z "$bam_dir" ]; then
  echo "BAM directory (-b) must be specified."
  usage
fi

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# List of BAM files in the directory
bam_files=("$bam_dir"/*.bam)

# Iterate over the BAM files and calculate coverage
for bam_file in "${bam_files[@]}"; do
    # Extract the sample name from the BAM file
    sample_name=$(basename "$bam_file" .bam)

    # Run samtools coverage on the BAM file and save the output to a coverage file
    coverage_file="${output_dir}/${sample_name}_coverage.txt"
    samtools coverage "$bam_file" > "$coverage_file"

    echo "Coverage calculated for $bam_file. Output file: $coverage_file"
done

echo "Coverage calculation complete!"
