#!/bin/bash

# Default values for parameters
input_file=""
output_dir=""
min_count=""
min_coverage=""
max_coverage=""
min_covered_fraction=""
window_sizes=""
step_sizes=""
pool_size=""
popoolation2_path=""


# Function to display usage instructions
usage() {
  echo "Usage: $0 [-i input_file] [-o output_dir] [-c min_count] [-C min_coverage] [-f max_coverage] [-w window_sizes] [-s step_sizes] [-p pool_size] [-P popoolation2_path] [-F min_covered_fraction]"
  echo "Options:"
  echo "  -i input_file               Path to the input file."
  echo "  -o output_dir               Path to the output directory."
  echo "  -c min_count                Minimum count for SNP calling."
  echo "  -C min_coverage             Minimum coverage for SNP calling."
  echo "  -f max_coverage             Maximum coverage for SNP calling."
  echo "  -w window_sizes             Comma-separated list of window sizes."
  echo "  -s step_sizes               Comma-separated list of step sizes."
  echo "  -p pool_size                Pool size for SNP calling."
  echo "  -P popoolation2_path        Path to the Popoolation2 directory."
  echo "  -F min_covered_fraction     Minimum covered fraction for SNP calling."
  exit 1
}

# Parse command-line options
while getopts "i:o:c:C:f:w:s:p:P:F:" opt; do
  case "$opt" in
    i) input_file="$OPTARG";;
    o) output_dir="$OPTARG";;
    c) min_count="$OPTARG";;
    C) min_coverage="$OPTARG";;
    f) max_coverage="$OPTARG";;
    w) window_sizes="$OPTARG";;
    s) step_sizes="$OPTARG";;
    p) pool_size="$OPTARG";;
    P) popoolation2_path="$OPTARG";;
    F) min_covered_fraction="$OPTARG";;
    \?) echo "Usage: $0 [-i input_file] [-o output_dir] [-c min_count] [-C min_coverage] [-f max_coverage] [-w window_sizes] [-s step_sizes] [-p pool_size] [-P popoolation2_path] [-F min_covered_fraction]" >&2
        exit 1;;
  esac
done

# Validate mandatory options
if [ -z "$input_file" ] || [ -z "$output_dir" ] || [ -z "$window_sizes" ] || [ -z "$step_sizes" ]; then
  echo "Error: Missing required options. Use -h for help."
  exit 1
fi

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Convert comma-separated strings to arrays
IFS=',' read -ra window_sizes <<< "$window_sizes"
IFS=',' read -ra step_sizes <<< "$step_sizes"

# Loop through window sizes and step sizes
for ((i = 0; i < ${#window_sizes[@]}; i++)); do
  window_size="${window_sizes[i]}"
  step_size="${step_sizes[i]}"

  # Define output file name based on window size and step size
  output_file="$output_dir/fst_result_${window_size}_${step_size}.txt"

  # Run the fst-sliding.pl script for pairwise FST calculations
  perl "${popoolation2_path}/fst-sliding.pl" \
    --input "$input_file" \
    --output "$output_file" \
    --min-count "$min_count" \
    --min-coverage "$min_coverage" \
    --max-coverage "$max_coverage" \
    --min-covered-fraction "$min_covered_fraction" \
    --window-size "$window_size" \
    --step-size "$step_size" \
    --pool-size "$pool_size" \
    --suppress-noninformative

  echo "Pairwise FST values for window size $window_size and step size $step_size calculated and stored in $output_file"
done
