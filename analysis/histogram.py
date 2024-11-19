import csv
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import sys

if len(sys.argv) != 2:
    print("Error!\nUsage: python3 histogram.py <construct>")
    sys.exit(1)

CONSTRUCT = sys.argv[1]

# Function to process a single CSV file
def process_csv_file(file_path):
    log_odd_scores = {
        -1: [],  # hit binding pocket residue
        -2: [],  # dead end
        -3: [],  # circular path
        -4: []   # max steps reached
    }
    
    with open(file_path, mode='r') as file:
        reader = csv.reader(file, delimiter=';')
        next(reader)  # skip header
        
        for row in reader:
            # Extract the 'Full Path' column for termination code (last atom in the path)
            full_path = row[2].split(',')  # Full Path column
            termination_code = int(full_path[-1])  # Last value in the 'Full Path' column
            
            # Extract the log probabilities (log) and sum them
            log_odds_str = row[3].split(',')  # 'Probabilities (Log)' column
            log_odds = list(map(float, log_odds_str))
            total_log_odd = sum(log_odds)
            
            # Assign the total log-odd score to the correct termination code category
            if termination_code in log_odd_scores:
                log_odd_scores[termination_code].append(total_log_odd)
    
    return log_odd_scores

# Function to aggregate log-odd scores from multiple CSV files
def aggregate_log_odd_scores(folder_path):
    aggregated_scores = {
        -1: [],  # hit binding pocket residue
        -2: [],  # dead end
        -3: [],  # circular path
        -4: []   # max steps reached
    }
    
    for file_name in os.listdir(folder_path):
        if file_name.endswith('.csv'):
            file_path = os.path.join(folder_path, file_name)
            scores = process_csv_file(file_path)
            
            for code in aggregated_scores:
                aggregated_scores[code].extend(scores[code])
    
    return aggregated_scores

# Function to plot histograms and save as SVG
def plot_log_odd_distributions(aggregated_scores, construct_name, output_folder):
    plt.figure(figsize=(10, 6))
    termination_labels = {
        -1: 'Hit Binding Pocket',
        -2: 'Dead End',
        -3: 'Circular Path',
        -4: 'Max Steps Reached'
    }

    # Updated soothing color palette
    colors = {
        -1: '#ff8a80',  # Soft red (flipped with dead end)
        -2: '#82b1ff',  # Soft blue (flipped with binding pocket)
        -3: '#b9f6ca',  # Soft green
        -4: '#b388ff'   # Soft purple
    }

    # Darker shades for line borders
    border_colors = {
        -1: '#b22222',  # Darker red (flipped with dead end)
        -2: '#0056b3',  # Darker blue (flipped with binding pocket)
        -3: '#00695c',  # Darker green
        -4: '#7b1fa2'   # Darker purple
    }

    # Set up logarithmic bins
    min_score = min([min(scores) for scores in aggregated_scores.values() if scores])
    max_score = max([max(scores) for scores in aggregated_scores.values() if scores])
    bins = np.logspace(np.log10(max(min_score, 1e-3)), np.log10(max_score), 30)

    # Calculate total scores for percentage calculation
    total_scores = sum(len(scores) for scores in aggregated_scores.values())
    
    # Sort aggregated scores by the length of scores (so higher frequencies are plotted first)
    sorted_scores = sorted(aggregated_scores.items(), key=lambda x: len(x[1]), reverse=True)
    
    for code, scores in sorted_scores:
        if scores:  # Plot only if there are scores for this termination code
            n = len(scores)
            percentage = (n / total_scores) * 100 if total_scores > 0 else 0
            plt.hist(scores, bins=bins, label=f"{termination_labels[code]} (n={n}, {percentage:.1f}%)", color=colors[code])
            
            # Calculate mean and median
            mean_score = np.mean(scores)
            median_score = np.median(scores)
            
            # Plot mean and median lines with border colors
            # Add a slight offset to the mean and median lines
            offset = 0.05  # Adjust as needed
            mean_line = plt.axvline(mean_score, color=border_colors[code], linestyle='dashed', linewidth=2,
                                    label=f"{termination_labels[code]} Mean")
            median_line = plt.axvline(median_score, color=border_colors[code], linestyle='solid', linewidth=2,
                                    label=f"{termination_labels[code]} Median")

            # Optional: Add markers at the end of the lines for clarity
            plt.plot(mean_score, plt.ylim()[1] * (1 - offset), marker='o', color=border_colors[code], markersize=5)
            plt.plot(median_score, plt.ylim()[1] * (1 - offset), marker='o', color=border_colors[code], markersize=5)

    plt.xscale('log')

    # Set specific ticks and labels for x-axis
    tick_values = [1e-3, 1e-2, 1e-1, 1, 10, 100, 1e3]
    tick_labels = [r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$']
    plt.xticks(tick_values, tick_labels, rotation=45)
    plt.xlim(1e-3, 1e2)  # Adjust based on the expected range of your data
    
    # Additional plot formatting
    plt.title(f"Log-Odd Score Distribution by Termination Code for {construct_name}")
    plt.xlabel("Log-Odd Score")
    plt.ylabel("Frequency")
    plt.legend()
    plt.grid(True)

    # Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Save the plot as an SVG file
    output_file = os.path.join(output_folder, f"{construct_name}_log_odd_distribution.svg")
    plt.savefig(output_file, format='svg')
    
    # Clear the plot to avoid overlapping in the next iteration
    plt.clf()


# Folder containing the CSV files for each construct
data_folder = f"../../data/{CONSTRUCT}/random_walk"
output_folder = "./"

# Aggregate log-odd scores for the current construct
aggregated_scores = aggregate_log_odd_scores(data_folder)

# Plot log-odd score distributions and save as SVG for the current construct
plot_log_odd_distributions(aggregated_scores, CONSTRUCT, output_folder)
print(f"Saved file in {output_folder}")