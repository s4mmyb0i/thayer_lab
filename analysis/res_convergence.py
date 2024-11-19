"""
Samvit Prem Singhal | spsinghal@wesleyan.edu
A script that graphs the normalized residue hitting times for all simulations for each construct

Input:  .csv files
Output: .svg plot
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

if len(sys.argv) != 2:
    print("Error!\nUsage: python3 res_convergence.py <construct>")
    sys.exit(1)

CONSTRUCT = sys.argv[1]

# Folder containing the CSV files for each construct
data_folder = f"../../data/{CONSTRUCT}/random_walk"
output_folder = "./"

# Loop through each construct
for filename in os.listdir(data_folder):
    if filename.endswith(".csv"):
        file_path = os.path.join(data_folder, filename)
        print(f"file path: {file_path}\n")
        
        # Load the CSV file
        df = pd.read_csv(file_path, delimiter='; ', engine='python')
        
        





for construct in constructs:
    plt.figure(figsize=(10, 6))  # Create a new figure for each construct

    # Loop through each simulation count for the current construct
    for sim in simulations:
        # Construct the file path for the current construct and simulation count
        file_path = f'project_{sim}/analysis/{construct}/res_frequencies.csv'

        # Load data using numpy
        data = np.genfromtxt(file_path, delimiter=',', skip_header=1)
        residues = data[:, 0]  # First column for Residue
        frequencies = data[:, 1]  # Second column for Frequency

        # Normalize frequencies by the number of simulations
        normalized_frequencies = frequencies / sim

        print(f"SIM: {sim}")
        
        # Check if any normalized frequency is greater than 0
        for i, freq in enumerate(normalized_frequencies):
            if freq < 1:
                print(f"Frequency for residue {residues[i]} is {freq:.4f} (original {frequencies[i]})")

        # Plot each normalized dataset with a label
        plt.plot(residues, normalized_frequencies, label=f'{sim} Simulations')

    # Add titles, labels, and legend
    plt.title(f'Normalized Hitting Times Across Simulations for {construct} Construct')
    plt.xlabel('Residue')
    plt.ylabel('Normalized Frequency (Proportion)')
    plt.legend()  # Add a legend to distinguish between the different simulations

    # Save the plot as an SVG file in the specified output directory
    output_file = os.path.join(output_dir, f'{construct}_residue_comparison.svg')
    plt.savefig(output_file, format='svg')
    plt.close()  # Close the figure to prevent overlap when creating multiple plots

print(f"Normalized graphs saved as SVG files in the directory: {output_dir}")










# """
# Samvit Prem Singhal | spsinghal@wesleyan.edu
# A script that graphs the normalized residue hitting times for all simulations for each construct

# Input:  .csv files
# Output: .svg plot
# """

# import numpy as np
# import matplotlib.pyplot as plt
# import os

# # List of constructs and simulation counts
# constructs = ['P', 'PL', 'AP', 'APL']
# simulations = [100, 1000, 10000, 100000, 1000000]
# output_dir = ''  # Directory to save graphs as SVG

# # Loop through each construct
# for construct in constructs:
#     plt.figure(figsize=(10, 6))  # Create a new figure for each construct

#     # Loop through each simulation count for the current construct
#     for sim in simulations:
#         # Construct the file path for the current construct and simulation count
#         file_path = f'project_{sim}/analysis/{construct}/res_frequencies.csv'

#         # Load data using numpy
#         data = np.genfromtxt(file_path, delimiter=',', skip_header=1)
#         residues = data[:, 0]  # First column for Residue
#         frequencies = data[:, 1]  # Second column for Frequency

#         # Normalize frequencies by the number of simulations
#         normalized_frequencies = frequencies / sim

#         # Check if any normalized frequency is less than 1
#         for freq, i in normalized_frequencies:
#             if freq > 0:
#                 print("frequency for res {residue[i]} is {freq} (original {frequencies[i]}))")

#         # Plot each normalized dataset with a label
#         plt.plot(residues, normalized_frequencies, label=f'{sim} Simulations')

#     # Add titles, labels, and legend
#     plt.title(f'Normalized Hitting Times Across Simulations for {construct} Construct')
#     plt.xlabel('Residue')
#     plt.ylabel('Normalized Frequency (Proportion)')
#     plt.legend()  # Add a legend to distinguish between the different simulations

#     # Save the plot as an SVG file in the specified output directory
#     output_file = os.path.join(output_dir, f'{construct}_residue_comparison.svg')
#     plt.savefig(output_file, format='svg')
#     plt.close()  # Close the figure to prevent overlap when creating multiple plots

# print(f"Normalized graphs saved as SVG files in the directory: {output_dir}")
















# """
# Samvit Prem Singhal | spsinghal@wesleyan.edu
# A script that graphs the normalized residue hitting times for all simulations for each construct

# Input:  .csv files
# Output: .svg plot
# """

# import numpy as np
# import matplotlib.pyplot as plt
# import os

# # List of constructs and simulation counts
# constructs = ['P', 'PL', 'AP', 'APL']
# simulations = [100, 1000, 10000, 100000, 1000000]
# output_dir = ''  # Directory to save graphs as SVG

# # Loop through each construct
# for construct in constructs:
#     plt.figure(figsize=(10, 6))  # Create a new figure for each construct

#     # Loop through each simulation count for the current construct
#     for sim in simulations:
#         # Construct the file path for the current construct and simulation count
#         file_path = f'project_{sim}/analysis/{construct}/res_frequencies.csv'

#         # Load data using numpy
#         data = np.genfromtxt(file_path, delimiter=',', skip_header=1)
#         residues = data[:, 0]  # First column for Residue
#         frequencies = data[:, 1]  # Second column for Frequency

#         # Normalize frequencies by the number of simulations
#         normalized_frequencies = frequencies / sim

#         print(f"SIM: {sim}")
        
#         # Check if any normalized frequency is greater than 0
#         for i, freq in enumerate(normalized_frequencies):
#             if freq < 1:
#                 print(f"Frequency for residue {residues[i]} is {freq:.4f} (original {frequencies[i]})")

#         # Plot each normalized dataset with a label
#         plt.plot(residues, normalized_frequencies, label=f'{sim} Simulations')

#     # Add titles, labels, and legend
#     plt.title(f'Normalized Hitting Times Across Simulations for {construct} Construct')
#     plt.xlabel('Residue')
#     plt.ylabel('Normalized Frequency (Proportion)')
#     plt.legend()  # Add a legend to distinguish between the different simulations

#     # Save the plot as an SVG file in the specified output directory
#     output_file = os.path.join(output_dir, f'{construct}_residue_comparison.svg')
#     plt.savefig(output_file, format='svg')
#     plt.close()  # Close the figure to prevent overlap when creating multiple plots

# print(f"Normalized graphs saved as SVG files in the directory: {output_dir}")










# # """
# # Samvit Prem Singhal | spsinghal@wesleyan.edu
# # A script that graphs the normalized residue hitting times for all simulations for each construct

# # Input:  .csv files
# # Output: .svg plot
# # """

# # import numpy as np
# # import matplotlib.pyplot as plt
# # import os

# # # List of constructs and simulation counts
# # constructs = ['P', 'PL', 'AP', 'APL']
# # simulations = [100, 1000, 10000, 100000, 1000000]
# # output_dir = ''  # Directory to save graphs as SVG

# # # Loop through each construct
# # for construct in constructs:
# #     plt.figure(figsize=(10, 6))  # Create a new figure for each construct

# #     # Loop through each simulation count for the current construct
# #     for sim in simulations:
# #         # Construct the file path for the current construct and simulation count
# #         file_path = f'project_{sim}/analysis/{construct}/res_frequencies.csv'

# #         # Load data using numpy
# #         data = np.genfromtxt(file_path, delimiter=',', skip_header=1)
# #         residues = data[:, 0]  # First column for Residue
# #         frequencies = data[:, 1]  # Second column for Frequency

# #         # Normalize frequencies by the number of simulations
# #         normalized_frequencies = frequencies / sim

# #         # Check if any normalized frequency is less than 1
# #         for freq, i in normalized_frequencies:
# #             if freq > 0:
# #                 print("frequency for res {residue[i]} is {freq} (original {frequencies[i]}))")

# #         # Plot each normalized dataset with a label
# #         plt.plot(residues, normalized_frequencies, label=f'{sim} Simulations')

# #     # Add titles, labels, and legend
# #     plt.title(f'Normalized Hitting Times Across Simulations for {construct} Construct')
# #     plt.xlabel('Residue')
# #     plt.ylabel('Normalized Frequency (Proportion)')
# #     plt.legend()  # Add a legend to distinguish between the different simulations

# #     # Save the plot as an SVG file in the specified output directory
# #     output_file = os.path.join(output_dir, f'{construct}_residue_comparison.svg')
# #     plt.savefig(output_file, format='svg')
# #     plt.close()  # Close the figure to prevent overlap when creating multiple plots

# # print(f"Normalized graphs saved as SVG files in the directory: {output_dir}")