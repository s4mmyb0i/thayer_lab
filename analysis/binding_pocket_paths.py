"""
Samvit Prem Singhal | spsinghal@wesleyan.edu
A script that extracts only the paths that hit the binding pocket
"""

import csv
import os
import sys

if len(sys.argv) != 2:
    print("Error!\nUsage: python3 binding_pocket_paths.py <construct>")
    sys.exit(1)

CONSTRUCT = sys.argv[1]

data_folder = f"../../data/{CONSTRUCT}/random_walk"
output_folder = f"../../data/{CONSTRUCT}/filtered_walks"

def filter_termination_code(input_file, output_file, termination_code=-1):
    # Open the input file and prepare to read
    with open(input_file, mode='r') as infile:
        reader = infile.readlines()
        
        # Extract header from the first line
        header = reader[0].strip().split("; ")
        
        # Filter rows where the last element in the second column is the termination code
        filtered_rows = []
        for line in reader[1:]:
            row = line.strip().split("; ")
            full_path = row[1].split(', ')
            full_path = list(map(int, full_path))
            if full_path[-1] == termination_code:
                filtered_rows.append(row)

    # Write the filtered rows to the output file
    with open(output_file, mode='w', newline='') as outfile:
        # Manually write header and rows to avoid quotes
        outfile.write("; ".join(header) + "\n")
        for row in filtered_rows:
            outfile.write("; ".join(row) + "\n")

# Process each CSV file in the data folder
for filename in os.listdir(data_folder):
    if filename.endswith(".csv"):
        input_file = os.path.join(data_folder, filename)
        output_file = os.path.join(output_folder, f"{os.path.splitext(filename)[0]}_filtered_walks.csv")
        filter_termination_code(input_file, output_file)
        print(f"Processed file {filename} and extracted paths")