#!/bin/bash

# Default values
input_dir=""
output_dir=""
threads=60

# Parse command-line arguments
while getopts "i:o:t:" opt; do
  case $opt in
    i) input_dir="$OPTARG" ;;
    o) output_dir="$OPTARG" ;;
    t) threads="$OPTARG" ;;
    \?) echo "Usage: $0 -i input_dir -o output_dir -t threads" >&2
        exit 1 ;;
  esac
done

# Check if input and output directories are specified
if [ -z "$input_dir" ] || [ -z "$output_dir" ]; then
  echo "Input and output directories must be specified. Use -i and -o options." >&2
  exit 1
fi

# Create FastQC output directory if it doesn't exist
mkdir -p "$output_dir"

# Run FastQC on all FASTQ files in the input directory
find "$input_dir" \( -name "*.fq.gz" -o -name "*.fastq.gz" \) -print0 | xargs -0 -P "$threads" -I {} fastqc -o "$output_dir" {}

# Run MultiQC on the FastQC output
multiqc -p --fn_as_s_name -i QC_final_Report -o "$output_dir" "$output_dir"
