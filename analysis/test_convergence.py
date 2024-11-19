"""
Samvit Prem Singhal | spsinghal@wesleyan.edu
A script to check if the walks are converging. Basically checking if the values in the how_end.csv files are changing as the number of simulations decreases.

Input:  .csv files
Output: .svg plot saved in analysis/{construct} directory.
"""

import numpy as np
import matplotlib.pyplot as plt
import re
import os
import sys

# Ensure the script is called with the required arguments
if len(sys.argv) != 2:
    print("Usage: python check_convergence.py <construct>")
    sys.exit(1)

# Construct name from command-line argument
construct = sys.argv[1]

# Define file paths for each simulation count and category
files = {
    'AP': [
        'project_100/analysis/AP/how_end.csv',
        'project_1000/analysis/AP/how_end.csv',
        'project_10000/analysis/AP/how_end.csv',
        'project_100000/analysis/AP/how_end.csv',
        'project_1000000/analysis/AP/how_end.csv',
    ],
    'APL': [
        'project_100/analysis/APL/how_end.csv',
        'project_1000/analysis/APL/how_end.csv',
        'project_10000/analysis/APL/how_end.csv',
        'project_100000/analysis/APL/how_end.csv',
        'project_1000000/analysis/APL/how_end.csv',
    ],
    'P': [
        'project_100/analysis/P/how_end.csv',
        'project_1000/analysis/P/how_end.csv',
        'project_10000/analysis/P/how_end.csv',
        'project_100000/analysis/P/how_end.csv',
        'project_1000000/analysis/P/how_end.csv',
    ],
    'PL': [
        'project_100/analysis/PL/how_end.csv',
        'project_1000/analysis/PL/how_end.csv',
        'project_10000/analysis/PL/how_end.csv',
        'project_100000/analysis/PL/how_end.csv',
        'project_1000000/analysis/PL/how_end.csv',
    ],
}

# Check if the construct is valid
if construct not in files:
    print(f"Invalid construct: {construct}. Choose from: {', '.join(files.keys())}")
    sys.exit(1)

# Initialize a plot
plt.figure(figsize=(12, 8))

# Define a color map for different categories
colors = {
    'AP': 'orange',
    'APL': 'red',
    'P': 'blue',
    'PL': 'green',
}

# Define a marker for different metrics
markers = {
    'Captured paths': 'D',
    'Circular paths': 'o',
    'Dead ends': 's',
    'Max steps reached': '^',
}

# Get the file paths for the chosen construct
file_paths = files[construct]

# Loop through each file and load the data
x_labels = []  # To store the number of walks for the x-axis
# Initialize a dictionary to hold percentages for different metrics
percentages = {key: [] for key in markers.keys()}

for file_path in file_paths:
    try:
        data = np.genfromtxt(file_path, delimiter=',', skip_header=1, dtype=None, encoding='utf-8')

        # Extract categories and percentages
        categories = data[:, 0]
        percentages_data = data[:, 2]

        # Extract the number of walks from the file path using regex
        walks = int(re.search(r'project_(\d+)', file_path).group(1))
        x_labels.append(walks)  # Append the number of walks for the x-axis

        # Loop through each metric and extract its percentage
        for metric in percentages.keys():
            if metric in categories:
                index = np.where(categories == metric)[0][0]
                percentage_value = percentages_data[index]
                if percentage_value != 'N/A':
                    percentages[metric].append(float(percentage_value.replace('%', '')))
                else:
                    percentages[metric].append(np.nan)  # Append NaN for invalid values
            else:
                percentages[metric].append(np.nan)  # Append NaN if metric is missing

    except Exception as e:
        print(f"Error loading {file_path}: {e}")

# Plot the percentages for each metric
for metric, marker in markers.items():
    plt.plot(x_labels, percentages[metric], label=f"{construct} - {metric}", color=colors[construct], marker=marker)

# Customize the plot
plt.xscale('log')  # Set x-axis to logarithmic scale for better visualization
plt.xlabel('Number of Walks (Log Scale)')
plt.ylabel('Percentage (%)')
plt.title(f'Percentage of Various Metrics for {construct} Across Different Simulation Runs')
plt.grid(True)
plt.legend()
plt.tight_layout()

# Create output directory if it doesn't exist
output_dir = f"analysis/{construct}"
os.makedirs(output_dir, exist_ok=True)

# Save the plot as SVG
output_file = f"{output_dir}/convergence_plot_{construct}.svg"
plt.savefig(output_file)
plt.show()

print(f"Plot saved to {output_file}")














# """
# Samvit Prem Singhal | spsinghal@wesleyan.edu
# A script to check if the walks are converging. Basically checking if the values in the how_end.csv files are changing as the number of simulations decreases

# Input:  .csv files
# Output: .svg plot
# """

# import numpy as np
# import matplotlib.pyplot as plt
# import re

# # Define file paths for each simulation count and category
# files = {
#     'AP': [
#         'project_100/analysis/AP/how_end.csv',
#         'project_1000/analysis/AP/how_end.csv',
#         'project_10000/analysis/AP/how_end.csv',
#         'project_100000/analysis/AP/how_end.csv',
#         'project_1000000/analysis/AP/how_end.csv',
#     ],
#     'APL': [
#         'project_100/analysis/APL/how_end.csv',
#         'project_1000/analysis/APL/how_end.csv',
#         'project_10000/analysis/APL/how_end.csv',
#         'project_100000/analysis/APL/how_end.csv',
#         'project_1000000/analysis/APL/how_end.csv',
#     ],
#     'P': [
#         'project_100/analysis/P/how_end.csv',
#         'project_1000/analysis/P/how_end.csv',
#         'project_10000/analysis/P/how_end.csv',
#         'project_100000/analysis/P/how_end.csv',
#         'project_1000000/analysis/P/how_end.csv',
#     ],
#     'PL': [
#         'project_100/analysis/PL/how_end.csv',
#         'project_1000/analysis/PL/how_end.csv',
#         'project_10000/analysis/PL/how_end.csv',
#         'project_100000/analysis/PL/how_end.csv',
#         'project_1000000/analysis/PL/how_end.csv',
#     ],
# }

# # Initialize a plot
# plt.figure(figsize=(12, 8))

# # Define a color map for different categories
# colors = {
#     'AP': 'orange',
#     'APL': 'red',
#     'P': 'blue',
#     'PL': 'green',
# }

# # Define a marker for different metrics
# markers = {
#     'Captured paths': 'D',
#     'Circular paths': 'o',
#     'Dead ends': 's',
#     'Max steps reached': '^',
# }

# # Loop through each category and load the corresponding CSV files
# for category, file_paths in files.items():
#     x_labels = []  # To store the number of walks for the x-axis
#     # Initialize a dictionary to hold percentages for different metrics
#     percentages = {key: [] for key in markers.keys()}

#     for file_path in file_paths:
#         try:
#             data = np.genfromtxt(file_path, delimiter=',', skip_header=1, dtype=None, encoding='utf-8')

#             # Extract categories and percentages
#             categories = data[:, 0]
#             percentages_data = data[:, 2]

#             # Extract the number of walks from the file path using regex
#             walks = int(re.search(r'project_(\d+)', file_path).group(1))
#             x_labels.append(walks)  # Append the number of walks for the x-axis

#             # Loop through each metric and extract its percentage
#             for metric in percentages.keys():
#                 if metric in categories:
#                     index = np.where(categories == metric)[0][0]
#                     percentage_value = percentages_data[index]
#                     if percentage_value != 'N/A':
#                         percentages[metric].append(float(percentage_value.replace('%', '')))
#                     else:
#                         percentages[metric].append(np.nan)  # Append NaN for invalid values
#                 else:
#                     percentages[metric].append(np.nan)  # Append NaN if metric is missing

#         except Exception as e:
#             print(f"Error loading {file_path}: {e}")

#     # Plot the percentages for each metric
#     for metric, marker in markers.items():
#         plt.plot(x_labels, percentages[metric], label=f"{category} - {metric}", color=colors[category], marker=marker)

# # Customize the plot
# plt.xscale('log')  # Set x-axis to logarithmic scale for better visualization
# plt.xlabel('Number of Walks (Log Scale)')
# plt.ylabel('Percentage (%)')
# plt.title('Percentage of Various Metrics Across Different Simulation Runs')
# plt.grid(True)
# plt.legend()
# plt.tight_layout()

# # Save the plot as SVG
# plt.savefig('test_convergence.svg')
# plt.show()