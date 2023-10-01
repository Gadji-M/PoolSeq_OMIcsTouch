#!/bin/bash

# Default values for parameters
input_dir=""
output_base_dir=""
threads=12

# Function to display usage instructions
usage() {
  echo "Usage: $0 -i <input_directory> -o <output_directory> [-t <threads>]"
  exit 1
}

# Parse command-line options
while getopts "i:o:t:" opt; do
  case "$opt" in
    i) input_dir="$OPTARG" ;;
    o) output_base_dir="$OPTARG" ;;
    t) threads="$OPTARG" ;;
    ?) usage ;;
  esac
done

# Check if required options are provided
if [ -z "$input_dir" ] || [ -z "$output_base_dir" ]; then
  echo "Error: Input and output directories are required."
  usage
fi

# Ensure the input directory exists
if [ ! -d "$input_dir" ]; then
  echo "Error: Input directory does not exist."
  usage
fi

# Create the output base directory if it doesn't exist
mkdir -p "$output_base_dir"

# List all pairs of forward and reverse reads in the input directory
samples=($(find "$input_dir" -type f -name "*_R1.fastq" -o -name "*_R1.fastq.gz"))

# Specify the parameters for megahit
megahit_params="-t $threads"

# Iterate over each pair of forward and reverse reads
for forward_read in "${samples[@]}"; do
  # Extract the sample name without the _R1 extension
  sample_name=$(basename -- "$forward_read")
  sample_name_noext="${sample_name%%_R1.fastq}"
  sample_name_noext="${sample_name_noext%%_R1.fastq.gz}"

  # Define the path to the corresponding reverse read
  reverse_read="$input_dir/${sample_name_noext}_R2.fastq"
  reverse_read_gz="$input_dir/${sample_name_noext}_R2.fastq.gz"

  # Create an output directory for the current sample
  output_dir="${output_base_dir}/${sample_name_noext}_output"
  mkdir -p "$output_dir"

  # Check if the reverse read exists as .fastq or .fastq.gz
  if [ -e "$reverse_read" ]; then
    megahit -1 "$forward_read" -2 "$reverse_read" -o "$output_dir" $megahit_params
  elif [ -e "$reverse_read_gz" ]; then
    megahit -1 "$forward_read" -2 "$reverse_read_gz" -o "$output_dir" $megahit_params
  else
    echo "Error: Reverse read not found for $sample_name_noext."
    continue
  fi

  echo "Finished processing $sample_name_noext"
done

echo "All samples processed"
