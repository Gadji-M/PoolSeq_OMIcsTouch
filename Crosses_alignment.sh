#!/bin/bash

# Specify the base directory where the sample folders are located
base_dir="/home/mahamat_g/gadjipoolseq/Raw_fastq"

# Specify the output directory for BAM files and error logs
output_dir="/home/mahamat_g/gadjipoolseq/Raw_fastq/bam_files"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over the sample folders
for sample_dir in "${base_dir}"/*; do
  # Extract the sample name from the sample directory name
  sample=$(basename "${sample_dir}")

  # Determine the number of reads for the sample
  num_reads=$(ls "${sample_dir}"/*_1.fq.gz | wc -l)

  # BWA-MEM Alignment command
  bwa mem -t 60 -R "@RG\tID:${sample}\tLB:${sample}\tSM:${sample}\tPL:ILLUMINA" \
    /mnt/46TB/Ghost/Space_Gadji/PoolSeq_2022/crosses/reference/fasta/VectorBase-61_AfunestusFUMOZ_Genome.fasta \
    <(cat "${sample_dir}"/*_1.fq.gz) \
    <(cat "${sample_dir}"/*_2.fq.gz) \
    | samtools view -@ 50 -bS - > "${output_dir}/${sample}.paired.bam" 2>"${output_dir}/${sample}.bwa-mem.err"

  # Create a flag file to indicate completion
  touch "${output_dir}/${sample}.paired.bam.done"
done

