import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from collections import Counter

def parse_walks(csv_file):
    """
    Parse the random walk CSV file and extract paths where termination_code == 1.
    Returns a list of paths (each path is a list of node IDs).
    """
    df = pd.read_csv(csv_file, sep=';')
    
    # Filter for termination_code == 1
    valid_walks = df[df['termination_code'] == 1]
    
    paths = []
    for path_str in valid_walks['path']:
        # Parse the path string - assuming comma or space separated integers
        if isinstance(path_str, str):
            # Try different separators
            if ',' in path_str:
                nodes = [int(x.strip()) for x in path_str.split(',') if x.strip()]
            else:
                nodes = [int(x.strip()) for x in path_str.split() if x.strip()]
            paths.append(nodes)
    
    return paths

def create_heatmap(csv_file, output_file='heatmap.html'):
    """
    Create an interactive 2D heatmap of node visit frequencies using Plotly.
    """
    # Parse walks and count node visits
    paths = parse_walks(csv_file)
    
    # Count how many times each node is visited
    all_nodes = []
    for path in paths:
        all_nodes.extend(path)
    
    node_counts = Counter(all_nodes)
    
    if not node_counts:
        print("No valid paths found!")
        return
    
    # Sort by node ID
    sorted_nodes = sorted(node_counts.items())
    node_ids = [n[0] for n in sorted_nodes]
    counts = [n[1] for n in sorted_nodes]
    
    # Create a matrix for the heatmap
    # Reshape into a grid (approximately square)
    n_nodes = len(node_ids)
    n_cols = int(np.ceil(np.sqrt(n_nodes)))
    n_rows = int(np.ceil(n_nodes / n_cols))
    
    # Pad with zeros if needed
    padded_counts = counts + [0] * (n_rows * n_cols - n_nodes)
    padded_node_ids = node_ids + [None] * (n_rows * n_cols - n_nodes)
    
    heatmap_matrix = np.array(padded_counts).reshape(n_rows, n_cols)
    node_id_matrix = np.array(padded_node_ids).reshape(n_rows, n_cols)
    
    # Create hover text with node IDs
    hover_text = []
    for i in range(n_rows):
        row_text = []
        for j in range(n_cols):
            if node_id_matrix[i, j] is not None:
                row_text.append(f'Node: {node_id_matrix[i, j]}<br>Visits: {heatmap_matrix[i, j]}')
            else:
                row_text.append('No data')
        hover_text.append(row_text)
    
    # Create the Plotly heatmap
    fig = go.Figure(data=go.Heatmap(
        z=heatmap_matrix,
        text=hover_text,
        hoverinfo='text',
        colorscale='Viridis',
        colorbar=dict(title='Visit Count'),
        showscale=True
    ))
    
    # Update layout
    fig.update_layout(
        title=dict(
            text=f'Node Visit Frequency Heatmap<br><sub>Total Walks: {len(paths)}, Unique Nodes: {len(node_counts)}</sub>',
            x=0.5,
            xanchor='center'
        ),
        xaxis=dict(title='Column Index'),
        yaxis=dict(title='Row Index'),
        width=1000,
        height=800,
        plot_bgcolor='white'
    )
    
    # Save to HTML file
    fig.write_html(output_file)
    print(f"Interactive heatmap saved to {output_file}")
    
    # Print statistics
    print(f"\nStatistics:")
    print(f"Total valid walks: {len(paths)}")
    print(f"Unique nodes visited: {len(node_counts)}")
    print(f"Most visited node: Node {max(node_counts, key=node_counts.get)} ({max(node_counts.values())} visits)")
    print(f"Least visited node: {min(node_counts.values())} visits")
    print(f"Average visits per node: {sum(node_counts.values()) / len(node_counts):.2f}")
    print(f"Median visits per node: {np.median(counts):.2f}")
    
    # Show top 10 most visited nodes
    print(f"\nTop 10 most visited nodes:")
    top_10 = sorted(node_counts.items(), key=lambda x: x[1], reverse=True)[:10]
    for node_id, count in top_10:
        print(f"  Node {node_id}: {count} visits")

if __name__ == "__main__":
    default_output = Path(__file__).resolve().parent / "figures" / "heatmap.html"
    parser = argparse.ArgumentParser(description='Create an interactive 2D heatmap of node visit frequencies')
    parser.add_argument('csv_file', help='Path to CSV file containing walks')
    parser.add_argument('--output', default=str(default_output), help='Output HTML file path')
    args = parser.parse_args()
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    create_heatmap(args.csv_file, args.output)