# import plotly.graph_objects as go
# import plotly.express as px
# import numpy as np
# import pandas as pd
# from utils.target import get_target_nodes
# from linalg.absorption_probs.absorption_probs import (
#     filter_by_abs_corr_threshold,
#     soften_threshold,
#     make_markov_transition,
#     absorption_dataframe,
# )

# # Function to compute node-wise mean/median/max absorption probabilities
# def compute_absorption_stats_nodewise(matrix_file: str, absorbing_nodes: set[int], thresholds: np.ndarray) -> pd.DataFrame:
#     stats = {"threshold": [], "mean": [], "median": [], "max": []}
#     raw_matrix = np.loadtxt(matrix_file)
#     raw_matrix = np.abs(raw_matrix)

#     for t in thresholds:
#         filtered = filter_by_abs_corr_threshold(raw_matrix, t)
#         filtered = soften_threshold(filtered, t)
#         markov_matrix = make_markov_transition(filtered, absorbing_nodes)
#         df = absorption_dataframe(markov_matrix, absorbing_nodes)

#         # Compute row-wise mean (over absorbing nodes)
#         row_means = df.mean(axis=1).values

#         stats["threshold"].append(t)
#         stats["mean"].append(np.mean(row_means))    # mean over transient nodes
#         stats["median"].append(np.median(row_means))  # median over transient nodes
#         stats["max"].append(np.max(row_means))      # max row mean

#     return pd.DataFrame(stats)

# # Config
# THRESHOLDS = np.arange(0.5, 0.95, 0.05)
# MATRIX_FILES = [
#     "/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/dccm/P.dat",
#     "/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/dccm/PL.dat",
#     "/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/dccm/AP.dat",
#     "/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/dccm/APL.dat"
# ]
# LABELS = ["P", "PL", "AP", "APL"]
# _, absorbing_nodes = get_target_nodes()

# # Plotly interactive plot
# colors = px.colors.qualitative.T10[:len(MATRIX_FILES)]
# line_styles = {"mean": "solid", "median": "dash", "max": "dot"}

# fig = go.Figure()

# for i, (matrix_file, label) in enumerate(zip(MATRIX_FILES, LABELS)):
#     stats_df = compute_absorption_stats_nodewise(matrix_file, absorbing_nodes, THRESHOLDS)
#     color = colors[i]

#     # Mean
#     fig.add_trace(go.Scatter(
#         x=stats_df["threshold"],
#         y=stats_df["mean"],
#         mode="lines+markers",
#         name=f"{label} - Mean",
#         line=dict(color=color, dash=line_styles["mean"])
#     ))

#     # Median
#     fig.add_trace(go.Scatter(
#         x=stats_df["threshold"],
#         y=stats_df["median"],
#         mode="lines+markers",
#         name=f"{label} - Median",
#         line=dict(color=color, dash=line_styles["median"])
#     ))

#     # Max
#     fig.add_trace(go.Scatter(
#         x=stats_df["threshold"],
#         y=stats_df["max"],
#         mode="lines+markers",
#         name=f"{label} - Max",
#         line=dict(color=color, dash=line_styles["max"])
#     ))

# fig.update_layout(
#     title="Node-wise Absorption Probability Statistics vs. Thresholds",
#     xaxis_title="Correlation Threshold",
#     yaxis_title="Absorption Probability",
#     legend_title="Protein & Statistic",
#     template="plotly_white"
# )

# fig.show()











"""
Samvit Prem Singhal | spsinghal@wesleyan.edu

Visualize absorption percentages for the four constructs at various thresholds.
Generates:
1. Line plots of summary stats (min, max, mean, median)
2. Heatmaps of node-wise absorption probabilities per protein
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from utils.target import get_target_nodes
from linalg.absorption_probs.absorption_probs import (
    filter_by_abs_corr_threshold,
    soften_threshold,
    make_markov_transition,
    absorption_dataframe,
)

# ============================================================
# CONFIGURATION
# ============================================================

MATRIX_FILES = [
    "/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/dccm/P.dat",
    "/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/dccm/PL.dat",
    "/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/dccm/AP.dat",
    "/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/dccm/APL.dat"
]
LABELS = ["P", "PL", "AP", "APL"]

THRESHOLDS = np.arange(0.5, 0.95, 0.05)
_, absorbing_nodes = get_target_nodes()

# ============================================================
# FUNCTIONS
# ============================================================

def compute_absorption_stats(matrix_file: str, absorbing_nodes: set[int], thresholds: np.ndarray) -> pd.DataFrame:
    """Compute min, max, mean, median absorption probabilities for multiple thresholds."""
    stats = {"threshold": [], "min": [], "max": [], "mean": [], "median": []}
    raw_matrix = np.loadtxt(matrix_file)
    raw_matrix = np.abs(raw_matrix)

    for t in thresholds:
        filtered = filter_by_abs_corr_threshold(raw_matrix, t)
        filtered = soften_threshold(filtered, t)
        markov_matrix = make_markov_transition(filtered, absorbing_nodes)

        df = absorption_dataframe(markov_matrix, absorbing_nodes)
        values = df.to_numpy().flatten()

        stats["threshold"].append(t)
        stats["min"].append(np.min(values))
        stats["max"].append(np.max(values))
        stats["mean"].append(np.mean(values))
        stats["median"].append(np.median(values))

    return pd.DataFrame(stats)

def compute_nodewise_absorption(matrix_file: str, absorbing_nodes: set[int], thresholds: np.ndarray) -> pd.DataFrame:
    """
    Returns a DataFrame: rows = thresholds, columns = node indices,
    values = mean absorption probability across absorbing nodes.
    """
    raw_matrix = np.loadtxt(matrix_file)
    raw_matrix = np.abs(raw_matrix)

    heatmap_data = []

    for t in thresholds:
        filtered = filter_by_abs_corr_threshold(raw_matrix, t)
        filtered = soften_threshold(filtered, t)
        markov_matrix = make_markov_transition(filtered, absorbing_nodes)

        df = absorption_dataframe(markov_matrix, absorbing_nodes)

        # Average across absorbing nodes for visualization
        mean_probs = df.mean(axis=1).values
        heatmap_data.append(mean_probs)

    # Convert to DataFrame: rows=thresholds, cols=nodes
    heatmap_df = pd.DataFrame(heatmap_data, index=[f"{t:.2f}" for t in thresholds])
    return heatmap_df

# ============================================================
# MAIN EXECUTION
# ============================================================





import plotly.graph_objects as go
import plotly.express as px

# Define colors for each protein
colors = px.colors.qualitative.T10[:len(MATRIX_FILES)]

fig = go.Figure()

line_styles = {
    "mean": "solid",
    "median": "dash",
    "max": "dot"
}

for i, (matrix_file, label) in enumerate(zip(MATRIX_FILES, LABELS)):
    stats_df = compute_absorption_stats(matrix_file, absorbing_nodes, THRESHOLDS)
    color = colors[i]

    # Plot mean
    fig.add_trace(go.Scatter(
        x=stats_df["threshold"],
        y=stats_df["mean"],
        mode="lines+markers",
        name=f"{label} - Mean",
        line=dict(color=color, dash=line_styles["mean"])
    ))

    # Plot median
    fig.add_trace(go.Scatter(
        x=stats_df["threshold"],
        y=stats_df["median"],
        mode="lines+markers",
        name=f"{label} - Median",
        line=dict(color=color, dash=line_styles["median"])
    ))

    # Plot max
    fig.add_trace(go.Scatter(
        x=stats_df["threshold"],
        y=stats_df["max"],
        mode="lines+markers",
        name=f"{label} - Max",
        line=dict(color=color, dash=line_styles["max"])
    ))

# Layout
fig.update_layout(
    title="Absorption Probability Statistics vs. Thresholds (Interactive)",
    xaxis_title="Correlation Threshold",
    yaxis_title="Absorption Probability",
    legend_title="Protein & Statistic",
    template="plotly_white"
)

out_dir = Path(__file__).resolve().parent / ".fds"
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / "absorption_stats.html"
fig.write_html(str(out_path))
print(f"Saved to {out_path}")






# plt.figure(figsize=(9, 7))

# # Define a consistent color palette
# colors = sns.color_palette("tab10", len(MATRIX_FILES))

# for i, (matrix_file, label) in enumerate(zip(MATRIX_FILES, LABELS)):
#     stats_df = compute_absorption_stats(matrix_file, absorbing_nodes, THRESHOLDS)
#     color = colors[i]

#     # Plot mean as solid line
#     plt.plot(
#         stats_df["threshold"], stats_df["mean"],
#         marker='o', linestyle='-', color=color, label=f"{label} - Mean"
#     )
    
#     # Plot median as dashed line with same color
#     plt.plot(
#         stats_df["threshold"], stats_df["median"],
#         marker='s', linestyle='--', color=color, label=f"{label} - Median"
#     )

# plt.xlabel("Correlation Threshold")
# plt.ylabel("Absorption Probability")
# plt.title("Absorption Probability Statistics vs. Thresholds")
# plt.legend(fontsize=8, ncol=2)
# plt.grid(True)
# plt.tight_layout()
# plt.show()

# # # ============================================================
# # # NODE-WISE HEATMAPS
# # # ============================================================
# 
# # for matrix_file, label in zip(MATRIX_FILES, LABELS):
# #     node_heatmap = compute_nodewise_absorption(matrix_file, absorbing_nodes, THRESHOLDS)
# 
# #     plt.figure(figsize=(12, 6))
# #     sns.heatmap(
# #         node_heatmap,
# #         cmap="magma",
# #         annot=False,
# #         cbar_kws={"label": "Absorption Probability"},
# #     )
# #     plt.title(f"Node-wise Absorption Probabilities for {label}")
# #     plt.xlabel("Node Index (Transient Nodes)")
# #     plt.ylabel("Correlation Threshold")
# #     plt.tight_layout()
# #     plt.show()




# """
# Samvit Prem Singhal | spsinghal@wesleyan.edu

# Script to visualize absorption percentages for the four constructs at various thresholds.
# """

# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns  # used for heatmap
# from utils.target import get_target_nodes
# from linalg.absorption_probs.absorption_probs import (
#     filter_by_abs_corr_threshold,
#     soften_threshold,
#     make_markov_transition,
#     absorption_dataframe,
# )

# # ============================================================
# # CONFIGURATION
# # ============================================================

# MATRIX_FILES = [
#     "/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/dccm/P.dat",
#     "/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/dccm/PL.dat",
#     "/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/dccm/AP.dat",
#     "/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/dccm/APL.dat"
# ]
# LABELS = ["P", "PL", "AP", "APL"]

# THRESHOLDS = np.arange(0.5, 0.95, 0.05)
# _, absorbing_nodes = get_target_nodes()

# # ============================================================
# # FUNCTIONS
# # ============================================================

# def compute_absorption_stats(matrix_file: str, absorbing_nodes: set[int], thresholds: np.ndarray) -> pd.DataFrame:
#     """Compute min, max, mean, median absorption probabilities for multiple thresholds."""
#     stats = {
#         "threshold": [],
#         "min": [],
#         "max": [],
#         "mean": [],
#         "median": [],
#     }

#     raw_matrix = np.loadtxt(matrix_file)
#     raw_matrix = np.abs(raw_matrix)

#     for t in thresholds:
#         filtered = filter_by_abs_corr_threshold(raw_matrix, t)
#         filtered = soften_threshold(filtered, t)
#         markov_matrix = make_markov_transition(filtered, absorbing_nodes)

#         df = absorption_dataframe(markov_matrix, absorbing_nodes)
#         values = df.to_numpy().flatten()

#         stats["threshold"].append(t)
#         stats["min"].append(np.min(values))
#         stats["max"].append(np.max(values))
#         stats["mean"].append(np.mean(values))
#         stats["median"].append(np.median(values))

#     return pd.DataFrame(stats)

# # ============================================================
# # MAIN EXECUTION
# # ============================================================

# plt.figure(figsize=(9, 7))

# # For heatmap data storage
# heatmap_data = pd.DataFrame(index=LABELS, columns=[f"{t:.2f}" for t in THRESHOLDS])

# for matrix_file, label in zip(MATRIX_FILES, LABELS):
#     stats_df = compute_absorption_stats(matrix_file, absorbing_nodes, THRESHOLDS)

#     # Overlay the 4 curves for each dataset with the same color hue
#     plt.plot(stats_df["threshold"], stats_df["mean"], marker='o', label=f"{label} - Mean")
#     plt.plot(stats_df["threshold"], stats_df["median"], marker='s', linestyle='--', label=f"{label} - Median")
#     plt.plot(stats_df["threshold"], stats_df["max"], linestyle=':', label=f"{label} - Max")
#     plt.plot(stats_df["threshold"], stats_df["min"], linestyle='-.', label=f"{label} - Min")

#     # Store mean probabilities for the heatmap
#     heatmap_data.loc[label] = stats_df["mean"].values

# plt.xlabel("Correlation Threshold")
# plt.ylabel("Absorption Probability")
# plt.title("Absorption Probability Statistics vs. Thresholds (Multiple Matrices)")
# plt.legend(fontsize=8, ncol=2)
# plt.grid(True)
# plt.tight_layout()
# plt.show()
