"""
Script to visualize and analyze a Markov chain absorption probability matrix.
Usage:
python3 -m linalg.absorption_probs.visualize /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/linalg/absorption_probs/P_trans_abs_probs.csv /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/linalg/absorption_probs

python3 -m linalg.absorption_probs.visualize /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/linalg/absorption_probs/PL_trans_abs_probs.csv /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/linalg/absorption_probs

python3 -m linalg.absorption_probs.visualize /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/linalg/absorption_probs/AP_trans_abs_probs.csv /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/linalg/absorption_probs

python3 -m linalg.absorption_probs.visualize /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/linalg/absorption_probs/APL_trans_abs_probs.csv /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/linalg/absorption_probs

"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go

if len(sys.argv) != 3:
    print("Usage: python3 analyze_absorption.py <absorption_probs_csv> <output_dir>")
    sys.exit(1)

csv_path = Path(sys.argv[1])
output_dir = Path(sys.argv[2])
output_dir.mkdir(parents=True, exist_ok=True)

# Extract prefix from filename (before first underscore)
prefix = csv_path.stem.split("_")[0]
log_file = output_dir / f"{prefix}_analysis_summary.txt"

if not csv_path.exists():
    print(f"File {csv_path} not found.")
    sys.exit(1)

with open(log_file, "w") as log:
    # Load CSV into DataFrame
    df = pd.read_csv(csv_path, index_col=0)
    
    log.write(f"Construct: {prefix}\n\n")
    log.write(f"Matrix shape: {df.shape} (transient nodes x absorbing nodes)\n\n")
    log.write("Example absorption profile (first row):\n")
    log.write(df.iloc[0].to_string() + "\n\n")

    # Sum over rows: how much influence each absorbing node has
    total_influence = df.sum(axis=0)
    total_influence /= total_influence.sum()
    log.write("Absorbing state influence (normalized):\n")
    log.write(total_influence.sort_values(ascending=False).to_string() + "\n\n")

    # Entropy of each transient node’s absorption profile
    entropy = -np.sum(df.values * np.log(df.values + 1e-12), axis=1)
    df['entropy'] = entropy

    log.write("Top 10 most certain transient nodes (lowest entropy):\n")
    log.write(df.sort_values(by='entropy').head(10).to_string() + "\n")

# Plotly heatmap of absorption probabilities (excluding entropy)
fig = px.imshow(
    df.drop(columns='entropy'),
    labels=dict(color="Absorption Probability"),
    x=df.columns[:-1],
    y=df.index,
    color_continuous_scale='Viridis',
    aspect="auto",
    title=f"{prefix} Absorption Probability Heatmap"
)
fig.update_layout(
    xaxis_title="Absorbing Nodes",
    yaxis_title="Transient Nodes",
    autosize=True,
    margin=dict(l=40, r=40, t=40, b=40),
    template="plotly_dark"
)
fig.write_html(str(output_dir / f"{prefix}_heatmap_absorption.html"))

# Bar chart of absorbing node influence
fig2 = go.Figure(data=[
    go.Bar(x=total_influence.index, y=total_influence.values)
])
fig2.update_layout(
    title=f"{prefix} Absorbing Node Influence (Normalized)",
    xaxis_title="Absorbing Nodes",
    yaxis_title="Total Absorption Probability",
    template="plotly_dark"
)
fig2.write_html(str(output_dir / f"{prefix}_bar_absorbing_influence.html"))

# Scatterplot of entropy across transient nodes
fig3 = px.scatter(
    df.reset_index(),
    x='index',
    y='entropy',
    title=f"{prefix} Entropy of Absorption Profiles",
    labels={"index": "Transient Node", "entropy": "Entropy"},
    template="plotly_dark"
)
fig3.write_html(str(output_dir / f"{prefix}_scatter_entropy.html"))
