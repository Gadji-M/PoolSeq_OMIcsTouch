#!/bin/bash

# Default values for parameters
threads=60
reference_genome="/mnt/46TB/Ghost/Space_Gadji/PoolSeq_2022/crosses/reference/fasta/VectorBase-61_AfunestusFUMOZ_Genome.fasta"
output_dir="/home/mahamat_g/gadjipoolseq/Raw_fastq/bam_files"

# Function to print script usage
usage() {
  echo "Usage: $0 [-t <threads>] [-r <reference_genome>] [-o <output_dir>]"
  exit 1
}

# Parse command-line options
while getopts ":t:r:o:" opt; do
  case $opt in
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

# Shift the parsed options so that "$@" only contains the sample directories
shift $((OPTIND-1))

# Specify the base directory where the sample folders are located
base_dir="/home/mahamat_g/gadjipoolseq/Raw_fastq"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over the sample folders
for sample_dir in "${base_dir}"/*; do
  # Extract the sample name from the sample directory name
  sample=$(basename "${sample_dir}")

  # Determine the number of reads for the sample
  num_reads=$(ls "${sample_dir}"/*_1.fq.gz | wc -l)

  # BWA-MEM Alignment command
  bwa mem -t "$threads" -R "@RG\tID:${sample}\tLB:${sample}\tSM:${sample}\tPL:ILLUMINA" \
    "$reference_genome" \
    <(cat "${sample_dir}"/*_1.fq.gz) \
    <(cat "${sample_dir}"/*_2.fq.gz) \
    | samtools view -@ 50 -bS - > "${output_dir}/${sample}.paired.bam" 2>"${output_dir}/${sample}.bwa-mem.err"

  # Create a flag file to indicate completion
  touch "${output_dir}/${sample}.paired.bam.done"
done

