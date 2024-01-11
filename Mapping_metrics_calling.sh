#!/bin/bash

# Default values
input_dir=""
reference_sequence=""
output_dir=""
picard_path=""  # Specify the path to picard.jar

# Function to print script usage
usage() {
  echo "Usage: $0 -i <input_dir> [-R <reference_sequence>] [-o <output_dir>] [-p <picard_path>]"
  exit 1
}

# Parse command-line options
while getopts ":i:R:o:p:" opt; do
  case $opt in
    i)
      input_dir="$OPTARG"
      ;;
    R)
      reference_sequence="$OPTARG"
      ;;
    o)
      output_dir="$OPTARG"
      ;;
    p)
      picard_path="$OPTARG"
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

# Check if the input directory is provided
if [ -z "$input_dir" ]; then
  echo "Input directory (-i) must be specified."
  usage
fi

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over the BAM files in the input directory and calculate alignment summaries
for bam_file in "$input_dir"/*.bam; do
    # Extract the sample name from the BAM file
    sample_name=$(basename "$bam_file" .bam)

    # Run Picard's CollectAlignmentSummaryMetrics on each BAM file
    output_file="${output_dir}/${sample_name}_alignment_summary_metrics.txt"
    java -jar "$picard_path" CollectAlignmentSummaryMetrics \
      R="$reference_sequence" \
      I="$bam_file" \
      O="$output_file"


echo "Mapping metrics calculated for $bam_file. Output file: $output_file"
done

echo "Mapping metrics calculation complete!"
