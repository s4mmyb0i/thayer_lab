"""
Visualize random walk paths as "highways" by combining walks that go through the same edges.

Usage:
    python3 analysis/visualize_highways.py <pdb_file> <csv_file> [max_paths] [min_edge]

Example:
    python3 analysis/visualize_highways.py /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/pdb/P.pdb "/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/random_walk/P_allos->lig_walks.csv" 1000 1
"""

import sys
from pathlib import Path
import pandas as pd
import plotly.graph_objects as go
from typing import Optional, Tuple
from collections import Counter
from utils import parse_path, create_atom_coord_map, get_residue_atom_sets


def visualize_highways_3d(
    pdb_file: str,
    csv_file: str,
    max_paths: Optional[int] = None,
    min_edge: int = 1,
    output_file: Optional[str] = None
):
    """
    Visualize random walk paths as "highways" by combining walks that go through the same edges.
    Edge thickness and opacity are proportional to how many paths traverse each edge.
    Only shows paths where termination_code == 1.
    Includes a slider to filter edges by minimum frequency (p).
    
    Args:
        pdb_file (str): Path to the PDB file
        csv_file (str): Path to the CSV file with random walk data
        max_paths (int, optional): Maximum number of paths to analyze, sorted by sum_neglog 
            (descending - highest values first). If None, analyzes all paths.
        min_edge (int): Minimum number of repetitions required for an edge. Only paths where
            all edges occur more than min_edge times will be shown. Default is 1.
        output_file (str, optional): Path to save the HTML figure. Defaults to analysis/figures/highways.html.
    """
    # Load data
    print("Loading CSV file...")
    df = pd.read_csv(csv_file, sep=';')
    
    # Filter to only show paths with termination_code == 1
    print("Filtering paths with termination_code == 1...")
    df = df[df['termination_code'] == 1].copy()
    print(f"Found {len(df)} paths with termination_code == 1")
    
    print("Extracting coordinates from PDB file...")
    coord_map = create_atom_coord_map(pdb_file)
    
    # Get residue to atom mapping for START_RES and END_RES
    print("Mapping residues to atoms for START_RES and END_RES...")
    start_atoms, end_atoms = get_residue_atom_sets(pdb_file)
    
    print(f"START_RES atoms: {len(start_atoms)}")
    print(f"END_RES atoms: {len(end_atoms)}")
    
    # Sort by sum_neglog (descending - highest values first)
    df_sorted = df.sort_values('sum_neglog', ascending=False).reset_index(drop=True)  # type: ignore
    
    # Parse paths
    print("Parsing paths...")
    df_sorted['parsed_path'] = df_sorted['path'].apply(parse_path)
    
    # Filter out empty paths
    df_sorted = df_sorted[df_sorted['parsed_path'].apply(len) > 0].copy()
    
    total_paths = len(df_sorted)
    if max_paths is None:
        max_paths = total_paths
    
    # Limit to max_paths
    df_sorted = df_sorted.head(max_paths).copy()
    
    print(f"Analyzing {len(df_sorted)} paths for highway visualization...")
    
    # First pass: Count edge frequencies (edges are consecutive node pairs)
    edge_frequencies: Counter[Tuple[int, int]] = Counter()
    
    for _, row in df_sorted.iterrows():
        path = row['parsed_path']
        # Count each edge (consecutive pair) in the path
        for i in range(len(path) - 1):
            node1 = int(path[i])
            node2 = int(path[i + 1])
            edge = (node1, node2)
            edge_frequencies[edge] += 1
    
    print(f"Found {len(edge_frequencies)} unique edges")
    if edge_frequencies:
        max_freq = max(edge_frequencies.values())
        min_freq = min(edge_frequencies.values())
        print(f"Edge frequency range: {min_freq} to {max_freq}")
    
    # Filter paths: only keep paths where all edges have frequency > min_edge
    print(f"Filtering paths to only include those where all edges occur more than {min_edge} times...")
    filtered_paths = []
    
    for _, row in df_sorted.iterrows():
        path = row['parsed_path']
        # Check if all edges in this path occur more than min_edge times
        path_valid = True
        for i in range(len(path) - 1):
            node1 = int(path[i])
            node2 = int(path[i + 1])
            edge = (node1, node2)
            if edge_frequencies[edge] <= min_edge:
                path_valid = False
                break
        
        if path_valid:
            filtered_paths.append(row)
    
    print(f"After filtering: {len(filtered_paths)} paths remain (out of {len(df_sorted)})")
    
    # Recalculate edge frequencies based only on filtered paths
    edge_frequencies_filtered: Counter[Tuple[int, int]] = Counter()
    
    for row in filtered_paths:
        path = row['parsed_path']
        for i in range(len(path) - 1):
            node1 = int(path[i])
            node2 = int(path[i + 1])
            edge = (node1, node2)
            edge_frequencies_filtered[edge] += 1
    
    print(f"Found {len(edge_frequencies_filtered)} unique edges in filtered paths")
    if edge_frequencies_filtered:
        max_freq = max(edge_frequencies_filtered.values())
        min_freq = min(edge_frequencies_filtered.values())
        print(f"Filtered edge frequency range: {min_freq} to {max_freq}")
    
    # Use filtered edge frequencies for visualization
    edge_frequencies = edge_frequencies_filtered
    
    # Create figure
    fig = go.Figure()
    
    # Store edge traces and their frequencies for slider control
    edge_trace_indices = []
    edge_freq_list = []  # Frequency for each edge trace
    trace_counter = 0
    
    if edge_frequencies:
        max_freq = max(edge_frequencies.values())
        min_freq = min(edge_frequencies.values())
        
        # Add each edge as a line segment
        for (node1, node2), freq in edge_frequencies.items():
            if node1 in coord_map and node2 in coord_map:
                x1, y1, z1 = coord_map[node1]
                x2, y2, z2 = coord_map[node2]
                
                # Normalize frequency to [0, 1] for color and width
                if max_freq > min_freq:
                    normalized_freq = (freq - min_freq) / (max_freq - min_freq)
                else:
                    normalized_freq = 1.0
                
                # Color: blue (low) to red (high)
                r = int(255 * normalized_freq)
                g = int(255 * normalized_freq * 0.3)
                b = int(255 * (1 - normalized_freq))
                color = f'rgb({r}, {g}, {b})'
                
                # Line width: 1 (thin) to 10 (thick) based on frequency
                line_width = 1 + normalized_freq * 9
                
                # Opacity: 0.3 (low) to 1.0 (high)
                opacity = 0.3 + normalized_freq * 0.7
                
                fig.add_trace(
                    go.Scatter3d(
                        x=[x1, x2],
                        y=[y1, y2],
                        z=[z1, z2],
                        mode='lines',
                        name=f'Edge ({node1}→{node2})',
                        line=dict(
                            color=color,
                            width=line_width
                        ),
                        opacity=opacity,
                        showlegend=False,
                        visible=True,  # All visible initially
                        hovertemplate=f'Edge: {node1} → {node2}<br>Frequency: {freq}<br>Traversed by {freq} path(s)<extra></extra>'
                    )
                )
                edge_trace_indices.append(trace_counter)
                edge_freq_list.append(freq)
                trace_counter += 1
    
    # Add START_RES atoms as special markers
    start_x, start_y, start_z = [], [], []
    for atom_id in start_atoms:
        if atom_id in coord_map:
            x, y, z = coord_map[atom_id]
            start_x.append(x)
            start_y.append(y)
            start_z.append(z)
    
    start_res_trace_idx = None
    end_res_trace_idx = None
    
    if start_x:
        fig.add_trace(
            go.Scatter3d(
                x=start_x,
                y=start_y,
                z=start_z,
                mode='markers',
                name='START_RES (Allosteric Binding Pocket)',
                marker=dict(size=10, color='green', symbol='diamond'),
                visible=True,
                hovertemplate='START_RES atom<extra></extra>'
            )
        )
        start_res_trace_idx = trace_counter
        trace_counter += 1
    
    # Add END_RES atoms as special markers
    end_x, end_y, end_z = [], [], []
    for atom_id in end_atoms:
        if atom_id in coord_map:
            x, y, z = coord_map[atom_id]
            end_x.append(x)
            end_y.append(y)
            end_z.append(z)
    
    if end_x:
        fig.add_trace(
            go.Scatter3d(
                x=end_x,
                y=end_y,
                z=end_z,
                mode='markers',
                name='END_RES (Ligand Binding Pocket)',
                marker=dict(size=10, color='red', symbol='square'),
                visible=True,
                hovertemplate='END_RES atom<extra></extra>'
            )
        )
        end_res_trace_idx = trace_counter
        trace_counter += 1
    
    # Get total number of traces
    num_traces = trace_counter
    
    # Create slider steps for minimum edge frequency (p)
    # p represents the minimum number of paths that must traverse an edge
    if edge_frequencies:
        max_freq = max(edge_frequencies.values())
        min_freq = min(edge_frequencies.values())
        
        # Create slider steps from min_freq to max_freq
        # p=0 means show all edges, p=max_freq means show only most frequent edges
        freq_steps = []
        for p in range(int(min_freq), int(max_freq) + 1):
            # Determine visibility: show edges with frequency > p
            visible_values = []
            for trace_idx in range(num_traces):
                if trace_idx in edge_trace_indices:
                    edge_idx = edge_trace_indices.index(trace_idx)
                    if edge_idx < len(edge_freq_list):
                        # Show edge if its frequency is greater than p
                        visible_values.append(edge_freq_list[edge_idx] > p)
                    else:
                        visible_values.append(False)
                else:
                    # Keep START_RES and END_RES traces visible
                    visible_values.append(True)
            
            step = dict(
                args=[{
                    "visible": visible_values
                }],
                label=f"{p}",
                method="restyle"
            )
            freq_steps.append(step)
        
        # Add slider
        sliders = [dict(
            active=0,  # Start with p=min_freq (show all edges)
            currentvalue={"prefix": "Min paths (p): "},
            pad={"t": 50},
            steps=freq_steps,
        )]
    else:
        sliders = []
    
    # Update layout
    fig.update_layout(
        title='3D Highway Visualization of Random Walk Paths<br><sub>Use slider to filter edges by minimum path frequency (p). Thickness and color indicate frequency.</sub>',
        scene=dict(
            xaxis_title='X (Angstrom)',
            yaxis_title='Y (Angstrom)',
            zaxis_title='Z (Angstrom)',
            aspectmode='data'
        ),
        sliders=sliders,
        height=800,
        width=1200
    )
    
    out_path = output_file or str(Path(__file__).resolve().parent / "figures" / "highways.html")
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    print(f"Saving highway visualization to {out_path}...")
    fig.write_html(out_path)


if __name__ == "__main__":
    import argparse
    default_output = str(Path(__file__).resolve().parent / "figures" / "highways.html")
    parser = argparse.ArgumentParser(description="Visualize random walk paths as highways in 3D")
    parser.add_argument("pdb_file", help="Path to PDB file")
    parser.add_argument("csv_file", help="Path to walk CSV file")
    parser.add_argument("max_paths", nargs="?", type=int, default=None, help="Max number of paths to analyze")
    parser.add_argument("min_edge", nargs="?", type=int, default=1, help="Minimum edge frequency to show")
    parser.add_argument("--output", "-o", dest="output_file", default=default_output, help="Output HTML file path")
    args = parser.parse_args()
    visualize_highways_3d(args.pdb_file, args.csv_file, args.max_paths, args.min_edge, args.output_file)

