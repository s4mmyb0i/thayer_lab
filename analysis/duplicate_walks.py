"""
Samvit Prem Singhal | spsinghal@wesleyan.edu
A script that counts the number of duplicate random walks for each construct
"""

import pandas as pd
import os
import sys

if len(sys.argv) != 2:
    print("Error!\nUsage: python3 duplicate_walk.py <construct>")
    sys.exit(1)

CONSTRUCT = sys.argv[1]

# Folder containing the CSV files for each construct
data_folder = f"../../data/{CONSTRUCT}/random_walk"
output_folder = "./"

def count_duplicate_walks(folder_path):
    total_duplicates = 0
    total_walks = 0
    
    # Loop through each CSV file in the construct's random_walk folder
    for filename in os.listdir(folder_path):
        if filename.endswith(".csv"):
            file_path = os.path.join(folder_path, filename)
            print(f"file path: {file_path}\n")
            
            # Load the CSV file
            df = pd.read_csv(file_path, delimiter='; ', engine='python')
            
            # Count total random walks in the file
            total_walks += len(df)
            
            # Identify duplicated full paths (excluding the first occurrence)
            duplicates_in_file = df[df.duplicated(subset="Full Path", keep=False)]
            num_duplicates = duplicates_in_file.duplicated(subset="Full Path", keep='first').sum()
            
            # Add the duplicates from this file to the total count
            total_duplicates += num_duplicates
    
    print(f"total dup: {total_duplicates}\ntotal walk: {total_walks}")
    return total_duplicates, total_walks

# # Ensure output folder exists
# if not os.path.exists(output_folder):
#     os.makedirs(output_folder)

# Count duplicated random walks for the current construct
num_duplicates, total_walks = count_duplicate_walks(data_folder)

print(f"total walkss: {total_walks}")
# Prepare the output text
output_text = (f"Total random walks: {total_walks}\nDuplicated random walks: {num_duplicates}\nPercentage duplicated: {100 * num_duplicates / total_walks}%\n")

# Print results to console
print(output_text)

# Save results to a text file in the output folder
output_file_path = os.path.join(output_folder, "duplicate_walks.txt")
with open(output_file_path, "w") as f:
    f.write(output_text)
