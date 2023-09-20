#!/bin/bash

# Default values for parameters
threads=30
reference_genome=""
output_dir=""

# Function to print script usage
usage() {
  echo "Usage: $0 -b <base_dir> [-t <threads>] [-r <reference_genome>] [-o <output_dir>]"
  exit 1
}

# Parse command-line options
while getopts ":b:t:r:o:" opt; do
  case $opt in
    b)
      base_dir="$OPTARG"
      ;;
    t)
      threads="$OPTARG"
      ;;
    r)
      reference_genome="$OPTARG"
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

# Check if base_dir is provided
if [ -z "$base_dir" ]; then
  echo "Error: -b <base_dir> is a required parameter."
  usage
fi

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over the sample folders
for sample_dir in "${base_dir}"/*; do
  # Extract the sample name from the sample directory name
  sample=$(basename "${sample_dir}")

  # Determine the number of forward and reverse reads for the sample
  num_forward_reads=$(ls "${sample_dir}"/*_1.fq.gz 2>/dev/null | wc -l)
  num_reverse_reads=$(ls "${sample_dir}"/*_2.fq.gz 2>/dev/null | wc -l)

  # Check if there is exactly one pair of forward and reverse reads
  if [ "$num_forward_reads" -eq 1 ] && [ "$num_reverse_reads" -eq 1 ]; then
    # BWA-MEM Alignment command for a single pair
    bwa mem -t "$threads" \
      "$reference_genome" \
      "${sample_dir}"/*_1.fq.gz \
      "${sample_dir}"/*_2.fq.gz \
      | samtools view -@ 30 -bS - > "${output_dir}/${sample}.paired.bam" 2>"${output_dir}/${sample}.bwa-mem.err"

    # Create a flag file to indicate completion
    touch "${output_dir}/${sample}.paired.bam.done"
  elif [ "$num_forward_reads" -ge 2 ] && [ "$num_reverse_reads" -ge 2 ]; then
    # Combine all forward reads into one file
    cat "${sample_dir}"/*_1.fq.gz > "${output_dir}/${sample}_combined_1.fq.gz"

    # Combine all reverse reads into one file
    cat "${sample_dir}"/*_2.fq.gz > "${output_dir}/${sample}_combined_2.fq.gz"

    # BWA-MEM Alignment command for combined reads
    bwa mem -t "$threads" -R "@RG\tID:${sample}\tLB:${sample}\tSM:${sample}\tPL:ILLUMINA" \
      "$reference_genome" \
      "${output_dir}/${sample}_combined_1.fq.gz" \
      "${output_dir}/${sample}_combined_2.fq.gz" \
      | samtools view -@ 3n0 -bS - > "${output_dir}/${sample}.paired.bam" 2>"${output_dir}/${sample}.bwa-mem.err"

    # Create a flag file to indicate completion
    touch "${output_dir}/${sample}.paired.bam.done"

    # Clean up combined read files
    rm "${output_dir}/${sample}_combined_1.fq.gz" "${output_dir}/${sample}_combined_2.fq.gz"
  else
    echo "Skipping sample $sample: Incorrect number of reads."
  fi
done
