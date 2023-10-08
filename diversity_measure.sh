#!/bin/bash

# Default values for options
grenedalf_path=""
bam_dir=""
window_sizes=""
window_strides=""
output_base_dir=""
sam_min_map_qual=""
sam_min_base_qual=""
multi_file_locus_set=""
filter_sample_min_count=""
filter_sample_max_count=""
filter_sample_min_coverage=""
filter_sample_max_coverage=""
measure=""
separator_char=""
na_entry=""
popoolation-corrected-tajimas-d=""
threads=""
log_file=""
window_type=""

# Parse command-line options
while getopts "p:d:w:s:o:M:S:Q:q:L:C:c:F:f:m:z:n:t:v:l:T:" opt; do
  case "$opt" in
    p) grenedalf_path="$OPTARG";;
    d) bam_dir="$OPTARG";;
    w) window_sizes="$OPTARG";;
    s) window_strides="$OPTARG";;
    o) output_base_dir="$OPTARG";;
    M) sam_min_map_qual="$OPTARG";;
    S) sam_min_base_qual="$OPTARG";;
    Q) multi_file_locus_set="$OPTARG";;
    q) filter_sample_min_count="$OPTARG";;
    L) filter_sample_max_count="$OPTARG";;
    C) filter_sample_min_coverage="$OPTARG";;
    c) filter_sample_max_coverage="$OPTARG";;
    F) measure="$OPTARG";;
    f) separator_char="$OPTARG";;
    m) na_entry="$OPTARG";;
    z) popoolation-corrected-tajimas-d="$OPTARG";;
    n) threads="$OPTARG";;
    t) log_file="$OPTARG";;
    T) window_type="$OPTARG";;
    \?) echo "Usage: $0 -p <grenedalf_path> -d <bam_directory> -w <window_sizes> -s <window_strides> -o <output_base_dir> -M <sam_min_map_qual> -S <sam_min_base_qual> -Q <multi_file_locus_set> -q <filter_sample_min_count> -L <filter_sample_max_count> -C <filter_sample_min_coverage> -c <filter_sample_max_coverage> -F <measure> -f <separator_char> -m <na_entry> ->
        exit 1;;
  esac
done

# Ensure the path to grenedalf and the BAM directory are specified
if [ -z "$grenedalf_path" ] || [ -z "$bam_dir" ]; then
  echo "Error: Path to grenedalf (-p) and BAM directory (-d) are required. Use -h for help."
  exit 1
fi

# Loop through all .bam files in the specified directory
for bam_file in "$bam_dir"/*.bam; do
  if [ -e "$bam_file" ]; then
    # Extract the BAM file name without extension
    bam_name=$(basename -- "$bam_file")
    bam_name_no_extension="${bam_name%.*}"

    # Create a subdirectory for each BAM file
    output_dir="$output_base_dir/$bam_name_no_extension"
    mkdir -p "$output_dir"

    IFS=' ' read -ra window_sizes <<< "$window_sizes"
    IFS=' ' read -ra window_strides <<< "$window_strides"

    # Loop through window sizes and strides
    for window_size in "${window_sizes[@]}"; do
      for window_stride in "${window_strides[@]}"; do
        # Create a subdirectory for each window size and stride combination
        window_dir="$output_dir/window_${window_size}_stride_${window_stride}"
        mkdir -p "$window_dir"

        # Construct the full path to the 'grendalf' binary
        grendalf_cmd="$grenedalf_path/grenedalf"

        # Run the grendalf diversity command with specified options and path
        cmd="$grendalf_cmd diversity --sam-path $bam_file --window-sliding-width $window_size --window-sliding-stride $window_stride --window-type $window_type --pool-sizes 50 --sam-min-map-qual $sam_min_map_qual --sam-min-base-qual $sam_min_base_qual --multi-file-locus-set $multi_file_locus_set --filter-sample-min-count $filter_sample_min_count --filter-sample>

        # Execute the command
        echo "Running command: $cmd"
        $cmd
      done
    done
  fi
done

#AUTHOR#
####GADJI_M###
