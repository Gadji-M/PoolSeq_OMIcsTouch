#!/bin/bash

# Default values for options
input_directory=""
output_directory=""
threads=1

while getopts ":i:o:t:" opt; do
  case $opt in
    i)
      input_directory="$OPTARG"
      ;;
    o)
      output_directory="$OPTARG"
      ;;
    t)
      threads="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Check for missing required options
if [ -z "$input_directory" ] || [ -z "$output_directory" ]; then
  echo "Usage: $0 -i input_directory -o output_directory [-t threads]"
  exit 1
fi

# Check if the output directory exists, and create it if not
if [ ! -d "$output_directory" ]; then
  mkdir -p "$output_directory"
fi

# Function to convert a single FASTQ file to FASTA format
convert_fastq_to_fasta() {
  input_fastq="$1"
  output_directory="$2"

  # Extract the filename without the path and extension
  file_name=$(basename "$input_fastq")
  file_name_no_ext="${file_name%.fastq}"

  # Define the corresponding output FASTA file
  output_fasta="$output_directory/$file_name_no_ext.fasta"

  # Convert the current FASTQ file to FASTA format
  awk 'NR%4==1 {print ">" substr($0,2)} NR%4==2 {print}' "$input_fastq" > "$output_fasta"

  echo "Converted $input_fastq to $output_fasta"
}

# Iterate through each FASTQ file in the input directory
for input_fastq in "$input_directory"/*.fastq; do
  # Check if there are any matching files
  if [ -f "$input_fastq" ]; then
    # Convert the current FASTQ file to FASTA format in parallel using specified threads
    convert_fastq_to_fasta "$input_fastq" "$output_directory" &
    ((counter++))

    # Limit the number of parallel processes based on the specified thread count
    if [ "$counter" -ge "$threads" ]; then
      wait
      counter=0
    fi
  fi
done

# Wait for any remaining background processes to finish
wait

echo "FASTQ to FASTA conversion complete."
