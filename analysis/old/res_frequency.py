"""
Samvit Prem Singhal | spsinghal@wesleyan.edu
A script that plots and describes the residues that were in the paths that hit the binding pocket residues.

Input: Directory of .csv files containing walks
Output: .txt file, .svg file, and .csv file with residue frequencies
"""

import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import sys

if len(sys.argv) != 2:
    print("Usage: python3 res_frequency.py <construct>")
    sys.exit(1)

CONSTRUCT = sys.argv[1]
# Step 1: Specify the directory containing the CSV files
INPUT_DIR = f'../../data/{CONSTRUCT}/random_walk'
OUTPUT_DIR = f'../../analysis/{CONSTRUCT}'

# Create the output directory if it does not exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Step 2: Initialize variables to store all residues and their frequencies
residue_frequencies = {}

# Step 3: Iterate over all CSV files in the specified directory
for filename in os.listdir(INPUT_DIR):
    if filename.endswith('.csv'):
        csv_file_path = os.path.join(INPUT_DIR, filename)

        # Load the CSV file, using the correct delimiter ';'
        data = np.genfromtxt(csv_file_path, delimiter=';', dtype=str, skip_header=0)

        # Filter paths that hit the binding pocket (ending with -1)
        paths = data[:, 1]  # Assuming paths are in the second column
        binding_pocket_paths = [path for path in paths if path.strip().endswith('-1')]

        # Count frequencies of residues in the filtered paths
        for path in binding_pocket_paths:
            # Clean up path string to extract residues
            path = path.strip()  # Remove any leading or trailing whitespace
            residues = path.split(', ')  # Split by comma and space to get residues

            for residue in residues:
                residue = residue.strip()  # Clean up residue (if needed)
                if int(residue) < 0:  # Ignore negative residues
                    continue
                if residue in residue_frequencies:
                    residue_frequencies[residue] += 1
                else:
                    residue_frequencies[residue] = 1

# Step 4: Sort residues and counts by residue number
residues = list(residue_frequencies.keys())
counts = [residue_frequencies[residue] for residue in residues]

# Convert residues to integers for sorting
residues = sorted(residues, key=lambda x: int(x))

# Step 5: Create a CSV file for the residue frequencies
output_csv_path = os.path.join(OUTPUT_DIR, 'res_frequencies.csv')
with open(output_csv_path, mode='w', newline='') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(['Residue', 'Frequency'])  # Write header
    for residue, count in zip(residues, counts):
        writer.writerow([residue, count])

# Step 6: Plot the frequencies of all residues in ascending order
plt.figure(figsize=(10, 6))
plt.bar(residues, counts)
plt.xlabel('Residue')
plt.ylabel('Frequency')
plt.title(f'Frequency of Residues Hitting the Binding Pocket (Sorted) for {CONSTRUCT}")')
# plt.ylim(0, 50000)
plt.xticks(rotation=45)
plt.xticks(residues[::2])  # Label every other residue
plt.tight_layout()

# Step 7: Save the plot as an SVG file
output_svg_path = os.path.join(OUTPUT_DIR, 'res_frequencies.svg')
plt.savefig(output_svg_path)
plt.close()

print(f"Script executed successfully. Output files generated in '{OUTPUT_DIR}': residue_frequencies.csv, residue_frequencies.svg")