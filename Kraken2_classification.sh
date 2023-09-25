#!/bin/bash

# Default values
fastq_dir=""
kraken_db=""
output_dir=""
threads=10

# Variables to store paths for unclassified and classified sequences
unclassified_dir=""
classified_dir=""

# Function to print script usage
print_usage() {
    echo "Usage: $0 -i <fastq_directory> -d <kraken_database> -o <output_directory> -u <unclassified_dir> -c <classified_dir> [-t <threads>]"
    echo "Options:"
    echo "  -i <fastq_directory>: Directory containing FastQ files"
    echo "  -d <kraken_database>: Kraken 2 database directory"
    echo "  -o <output_directory>: Directory to store Kraken reports"
    echo "  -u <unclassified_dir>: Directory to store unclassified sequences"
    echo "  -c <classified_dir>: Directory to store classified sequences"
    echo "  -t <threads>: Number of threads for parallel processing (default: 1)"
    exit 1
}

# Parse command-line options
while getopts "i:d:o:u:c:t:" opt; do
    case "$opt" in
        i) fastq_dir="$OPTARG" ;;
        d) kraken_db="$OPTARG" ;;
        o) output_dir="$OPTARG" ;;
        u) unclassified_dir="$OPTARG" ;;
        c) classified_dir="$OPTARG" ;;
        t) threads="$OPTARG" ;;
        *) print_usage ;;
    esac
done

# Check if required options are provided
if [ -z "$fastq_dir" ] || [ -z "$kraken_db" ] || [ -z "$output_dir" ] || [ -z "$unclassified_dir" ] || [ -z "$classified_dir" ]; then
    echo "Error: FastQ directory, Kraken database, output directory, unclassified directory, and classified directory are required."
    print_usage
fi

# Create the output directories if they don't exist
mkdir -p "$output_dir" "$unclassified_dir" "$classified_dir"

# Iterate over the FastQ files in the directory
find "$fastq_dir" -maxdepth 1 -type f -name "*.fastq" | while read -r fastq_file; do
    # Extract the sample name from the FastQ file name
    sample=$(basename "${fastq_file}" .fastq)

    # Run Kraken 2 on the FastQ file
    kraken2 --db "${kraken_db}" --report "${output_dir}/${sample}.kraken2.report" --output "${output_dir}/${sample}.kraken2.out" --unclassified-out "${unclassified_dir}/${sample}.unclassified.fastq" --classified-out "${classified_dir}/${sample}.classified.fastq" --use-names --threads "$threads" "${fastq_file}"
done

echo "Kraken analysis completed."
