#!/bin/bash

# Default values
bam_directory=""
output_directory=""
picard_path=""
reference_genome=""

# Function to print script usage
usage() {
  echo "Usage: $0 -b <bam_directory> -o <output_directory> -p <picard_path> [-r <reference_genome>]"
  exit 1
}

# Parse command-line options
while getopts ":b:o:p:r:" opt; do
  case $opt in
    b)
      bam_directory="$OPTARG"
      ;;
    o)
      output_directory="$OPTARG"
      ;;
    p)
      picard_path="$OPTARG"
      ;;
    r)
      reference_genome="$OPTARG"
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

# Check if required parameters are provided
if [ -z "$bam_directory" ] || [ -z "$output_directory" ] || [ -z "$picard_path" ]; then
  echo "BAM directory (-b), output directory (-o), and Picard path (-p) must be specified."
  usage
fi

# Create the output directory if it doesn't exist
mkdir -p "$output_directory"

# Iterate over BAM files in the directory
for bam_file in "$bam_directory"/*.bam; do
    # Get the filename without the extension
    filename=$(basename "$bam_file" .bam)

    # Output file name for the metrics
    output_file="${output_directory}/${filename}_alignment_metrics.txt"

    # Run Picard CollectAlignmentSummaryMetrics
    java -jar "${picard_path}/picard.jar" CollectAlignmentSummaryMetrics \
    --COLLECT_ALIGNMENT_INFORMATION true \
        -R "$reference_genome" \
        -I "$bam_file" \
        -O "$output_file" \
        --VALIDATION_STRINGENCY SILENT

    echo "Alignment metrics saved for $filename: $output_file"
done
