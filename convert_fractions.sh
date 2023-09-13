#!/bin/bash

# Initialize variables with default values
input_file=""
output_file=""
columns_to_convert=""

# Parse command-line arguments
while getopts "i:o:c:" opt; do
    case $opt in
        i) input_file="$OPTARG";;
        o) output_file="$OPTARG";;
        c) columns_to_convert="$OPTARG";;
        \?) echo "Usage: $0 -i input_file -o output_file -c columns_to_convert"
            exit 1;;
    esac
done

# Check if required arguments are provided
if [ -z "$input_file" ] || [ -z "$output_file" ] || [ -z "$columns_to_convert" ]; then
    echo "Usage: $0 -i input_file -o output_file -c columns_to_convert"
    exit 1
fi

# Process the input file and generate the output
awk -F '\t' -v OFS='\t' -v columns="$columns_to_convert" 'BEGIN {
    split(columns, cols, " ");
    for (i in cols) {
        col = cols[i];
        column_numbers[col] = 1;
    }
}
{
    for (col in column_numbers) {
        if ($col ~ "/") {
            split($col, frac, "/");
            if (frac[2] != 0) {
                $col = (frac[1] == frac[2]) ? 1 : ((frac[1] == 0) ? 0 : (frac[1] / frac[2]));
                $col = sprintf("%.2f", $col); # Limit to two decimal places
            } else {
                $col = 0;
            }
        }
    }
} 1' "$input_file" > "$output_file"

echo "Conversion completed. Results saved in $output_file."
