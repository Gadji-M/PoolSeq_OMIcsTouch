#!/bin/bash

while getopts "r:i:o:" opt; do
  case $opt in
    r)
      reference_genome="$OPTARG"
      ;;
    i)
      IFS=',' read -ra bam_files <<< "$OPTARG"
      ;;
    o)
      output_dir="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

# Create output directory if it doesn't exit

mkdir -p "$output_dir"

if [ -z "$reference_genome" ] || [ -z "$bam_files" ] || [ -z "$output_dir" ]; then
  echo "Usage: $0 -r <reference_genome> -i <input_bam_files> -o <output_dir>"
  exit 1
fi

# Loop through the BAM files and call variants for each
for bam_file in "${bam_files[@]}"; do
    sample_name=$(basename "$bam_file" .bam)
    output_vcf="$output_dir/$sample_name.vcf"
    echo "Processing $bam_file..."

    # Call variants using bcftools
    bcftools mpileup -Ou -f "$reference_genome" "$bam_file" | bcftools call -O b -m -v --threads 50 -o "$output_vcf"

    echo "Variants called for $bam_file. Output saved to $output_vcf"
done


###############
####AUTHOR#####
####GADJI_M####
###############
