#!/bin/bash

# Default values
input_dir=""
output_dir=""
output_metrics=""
picard_jar=""

# Function to print script usage
usage() {
  echo "Usage: $0 -i <input_dir> -o <output_dir> -m <output_metrics> [-p <picard_jar>]"
  exit 1
}

# Parse command-line options
while getopts ":i:o:m:p:" opt; do
  case $opt in
    i)
      input_dir="$OPTARG"
      ;;
    o)
      output_dir="$OPTARG"
      ;;
    m)
      output_metrics="$OPTARG"
      ;;
    p)
      picard_jar="$OPTARG"
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
if [ -z "$input_dir" ] || [ -z "$output_dir" ] || [ -z "$output_metrics" ]; then
  echo "Input directory (-i), output directory (-o), and output metrics (-m) must be specified."
  usage
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over the BAM files in the input directory
for bam_file in "${input_dir}"/*.bam; do
    # Extract the file name without extension
    filename=$(basename "$bam_file" .bam)
    
    # Define the sorted BAM file path
    sorted_bam="${output_dir}/${filename}_sorted.bam"
    
    # Define the marked duplicates BAM file path
    marked_dup_bam="${output_dir}/${filename}_marked_duplicates.bam"
    
    # Sort the BAM file using Picard
    java -jar "$picard_jar" SortSam \
       -I "$bam_file" \
        -O "$sorted_bam" \
        -SO coordinate \
        --VALIDATION_STRINGENCY SILENT
    
    # Mark duplicates using Picard
    java -jar "$picard_jar" MarkDuplicates \
        -I "$sorted_bam" \
        -O "$marked_dup_bam" \
        -M "${output_dir}/${filename}_${output_metrics}" \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY SILENT
done
