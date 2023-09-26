#!/usr/bin/env python3
import argparse

# Function to parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract sequences from a FASTA file based on sequence IDs.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    parser.add_argument("-l", "--id-list", required=True, help="File containing a list of sequence IDs to extract")
    return parser.parse_args()

# Main function
def main():
    args = parse_arguments()
    id_set = set()

    # Read the list of sequence IDs from the file
    with open(args.id_list, "r") as id_file:
        for line in id_file:
            line = line.strip()
            if line:
                id_set.add(line)

    # Open the input and output files
    with open(args.input, "r") as input_file, open(args.output, "w") as output_file:
        sequence = ""
        write_sequence = False
        current_sequence_id = ""

        for line in input_file:
            if line.startswith(">"):
                current_sequence_id = line[1:].strip()
                # Check if the current sequence ID is in the set of IDs to extract
                if current_sequence_id in id_set:
                    write_sequence = True
                else:
                    write_sequence = False
                # Write the sequence identifier to the output
                output_file.write(line)
            elif write_sequence:
                # Write the sequence lines to the output
                output_file.write(line)

if __name__ == "__main__":
    main()

#AUTHOR: GADJI_M#
