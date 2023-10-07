#!/usr/bin/env python3
import argparse
import os

# Function to validate mandatory options
def validate_mandatory_options(args):
    if not args.input_file or not args.output_dir or not args.window_sizes or not args.step_sizes:
        parser.error("Error: Missing required options. Use -h for help.")

# Parse command-line options
parser = argparse.ArgumentParser(description="Pairwise FST calculations using Popoolation2")
parser.add_argument("-i", "--input_file", required=True, help="Input file")
parser.add_argument("-o", "--output_dir", required=True, help="Output directory")
parser.add_argument("-c", "--min_count", type=int, help="Minimum count of alleles required")
parser.add_argument("-C", "--min_coverage", type=int, help="Minimum coverage required")
parser.add_argument("-f", "--max_coverage", type=str, help="Maximum coverage as a percentage")
parser.add_argument("-w", "--window_sizes", required=True, help="Comma-separated list of window sizes")
parser.add_argument("-s", "--step_sizes", required=True, help="Comma-separated list of step sizes")
parser.add_argument("-p", "--pool_size", type=int, help="Pool size")
parser.add_argument("-P", "--popoolation2_path", default="", help="Path to Popoolation2")
parser.add_argument("-F", "--min_covered_fraction", type=float, help="Minimum covered fraction")
parser.add_argument("-n", "--suppress_noninformative", action="store_true", help="Suppress noninformative")

args = parser.parse_args()

# Validate mandatory options
validate_mandatory_options(args)

# Create the output directory if it doesn't exist
os.makedirs(args.output_dir, exist_ok=True)

# Convert comma-separated strings to lists
window_sizes = args.window_sizes.split(',')
step_sizes = args.step_sizes.split(',')

# Loop through window sizes and step sizes
for window_size in window_sizes:
    for step_size in step_sizes:

        # Define output file name based on window size and step size
        output_file = os.path.join(args.output_dir, f"fst_result_w{window_size}_s{step_size}.txt")

        # Run the fst-sliding.pl script for pairwise FST calculations
        cmd = [
            "perl",
            os.path.join(args.popoolation2_path, "fst-sliding.pl"),
            "--input", args.input_file,
            "--output", output_file,
            "--min-count", str(args.min_count) if args.min_count is not None else "",
            "--min-coverage", str(args.min_coverage) if args.min_coverage is not None else "",
            "--max-coverage", args.max_coverage if args.max_coverage else "",
            "--min-covered-fraction", str(args.min_covered_fraction) if args.min_covered_fraction is not None else "",
            "--window-size", window_size,
            "--step-size", step_size,
            "--pool-size", str(args.pool_size) if args.pool_size is not None else "",
            "--suppress-noninformative" if args.suppress_noninformative else "",
        ]

        os.system(" ".join(cmd))

      #python Fst_sliding_windows_equi.py -i path/to/input_file -o path/to/output_dir -c 2 -C 10 -f "5%" -w "500,1000,2000" -s "250,500,1000" -p 40 -P "path/to/popoolation2_1201/" -F 0.2

        print(f"Pairwise FST values for window size {window_size} and step size {step_size} calculated and stored in {output_file}")
