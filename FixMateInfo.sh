#!/usr/bin/env bash

# Set the default values for options
PICARD_JAR=""
BAM_DIR=""
OUTPUT_DIR=""
VALIDATION_STRINGENCY="SILENT"
SORT_ORDER="coordinate"
CREATE_INDEX="false"

# Function to display usage instructions
usage() {
  echo "Usage: $0 -p <picard_jar_path> -i <input_directory> -o <output_directory> [-v <validation_stringency>] [-s <sort_order>] [-c <create_index>]"

  echo "Options:"

  echo "  -p <picard_jar_path>   represents the path to the picard file which is the java program to use for the JOB;"

  echo "  -i <input_directory>    Specifies the directory containing the input BAM files to process;"

  echo "  -o <output_directory>   Specifies the directory where the output BAM files with fixed mate information will be stored;"

  echo "  -v <validation_stringency>  Allows you to specify the validation stringency level for Picard's FixMateInformation tool. Setting stringency to SILENT can improve performance when processing BAM files in which variable-length data (read, qualities, tags) do not otherwise need to be decoded. The default is set to "SILENT," but you can change it to other levels like "LENIENT" or "STRICT" as needed;"

  echo "  -s <sort_order>         Specifies the sorting order for the output BAM files. The default is set to "coordinate," which is a common sorting order for BAM files;"

  echo "  -c <create_index>       Specifies whether to create an index for the output BAM files. The default is set to "false." If you want to create an index, set it to "true"."



  echo "
        ########AUTHOR: GADJI_M########
        ###############################"
  exit 1
}



# Parse command-line options
while getopts "p:i:o:v:s:c:" opt; do
  case $opt in
    p) PICARD_JAR="$OPTARG" ;;
    i) BAM_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    v) VALIDATION_STRINGENCY="$OPTARG" ;;
    s) SORT_ORDER="$OPTARG" ;;
    c) CREATE_INDEX="$OPTARG" ;;
    \?) usage ;;
  esac
done

# Check if the help option is provided
if [ "$1" = "--help" ]; then
  usage
fi


# Check if required options are provided
if [ -z "$PICARD_JAR" ] || [ -z "$BAM_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
  usage
fi

# Process each BAM file in the directory
for BAM_FILE in "$BAM_DIR"/*.bam; do
  if [ -f "$BAM_FILE" ]; then
    # Get the base name of the BAM file without extension
    BASE_NAME=$(basename "$BAM_FILE" .bam)

    # Generate the output file name for fixed mate information
    OUTPUT_BAM="$OUTPUT_DIR/${BASE_NAME}.paired_marked_duplicates.fixed_mate.bam"

    # Run Picard to fix mate information
    java -jar "$PICARD_JAR" FixMateInformation -I "$BAM_FILE" -O "$OUTPUT_BAM" -MC true \
      --VALIDATION_STRINGENCY "$VALIDATION_STRINGENCY" --VERBOSITY INFO -SO "$SORT_ORDER" --CREATE_INDEX "$CREATE_INDEX"

    # Print a message indicating the completion of the operation
    echo "Fixed mate information for $BAM_FILE and saved to $OUTPUT_BAM"
  fi
done
