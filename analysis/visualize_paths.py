"""
Visualize random walk paths in 3D space with interactive sliders.

Usage:
    python3 analysis/visualize_paths.py <pdb_file> <csv_file> [max_paths] [--output <path>]

Examples:
    python3 analysis/visualize_paths.py pdb/P.pdb random_walk/P_allos->lig_walks.csv 1000
    python3 analysis/visualize_paths.py pdb/P.pdb random_walk/P_allos->lig_walks.csv --output figures/P_paths.html
"""

import sys
from pathlib import Path
import pandas as pd
import plotly.graph_objects as go
from typing import Optional
from utils import parse_path, create_atom_coord_map, get_residue_atom_sets


def get_color_for_index(i: int, total: int) -> str:
    """
    Generate a color for path index i out of total paths.
    Creates a gradient from blue to red.
    """
    if total <= 1:
        return 'rgb(0, 0, 255)'
    ratio = i / max(total - 1, 1)
    r = int(255 * (1 - ratio))
    g = int(255 * ratio * 0.5)
    b = int(255 * ratio)
    return f'rgb({r}, {g}, {b})'


def visualize_paths_3d(
    pdb_file: str,
    csv_file: str,
    max_paths: Optional[int] = None,
    output_file: Optional[str] = None
):
    """
    Visualize random walk paths in 3D space with interactive sliders.
    Only shows paths where termination_code == 1.
    
    Args:
        pdb_file (str): Path to the PDB file
        csv_file (str): Path to the CSV file with random walk data
        max_paths (int, optional): Maximum number of paths to show, sorted by sum_neglog 
            (descending - highest values first). If None, shows all paths.
        output_file (str, optional): If set, save the figure to this HTML path instead of showing.
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
    
    # Get path length range for slider
    min_length = int(df_sorted['length'].min())
    max_length = int(df_sorted['length'].max())
    print(f"Path length range: {min_length} to {max_length}")
    
    total_paths = len(df_sorted)
    if max_paths is None:
        max_paths = total_paths
    
    print(f"Total paths: {total_paths}")
    print(f"Will show up to {max_paths} paths")
    
    # Create figure
    fig = go.Figure()
    
    n_paths_to_show = min(max_paths, total_paths)
    
    # Add all paths as traces (we'll control visibility with sliders)
    path_trace_indices = []
    trace_counter = 0
    for i in range(n_paths_to_show):
        row = df_sorted.iloc[i]
        path = row['parsed_path']
        
        # Get coordinates for this path
        x_coords = []
        y_coords = []
        z_coords = []
        
        for atom_id in path:
            if atom_id in coord_map:
                x, y, z = coord_map[atom_id]
                x_coords.append(x)
                y_coords.append(y)
                z_coords.append(z)
        
        if len(x_coords) > 0:
            # Initial color (will be updated by slider)
            color = get_color_for_index(i, n_paths_to_show)
            
            trace = go.Scatter3d(
                x=x_coords,
                y=y_coords,
                z=z_coords,
                mode='lines+markers',
                name=f'Path {i+1} (sum_neglog: {row["sum_neglog"]:.2f}, length: {row["length"]})',
                line=dict(color=color, width=3),
                marker=dict(size=4, color=color),
                visible=True,
                hovertemplate=f'Path {i+1}<br>sum_neglog: {row["sum_neglog"]:.2f}<br>length: {row["length"]}<extra></extra>'
            )
            fig.add_trace(trace)
            path_trace_indices.append(trace_counter)
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
    if start_x:
        fig.add_trace(
            go.Scatter3d(
                x=start_x,
                y=start_y,
                z=start_z,
                mode='markers',
                name='START_RES (Ligand Binding Pocket)',
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
    
    end_res_trace_idx = None
    if end_x:
        fig.add_trace(
            go.Scatter3d(
                x=end_x,
                y=end_y,
                z=end_z,
                mode='markers',
                name='END_RES (Allosteric Binding Pocket)',
                marker=dict(size=10, color='red', symbol='square'),
                visible=True,
                hovertemplate='END_RES atom<extra></extra>'
            )
        )
        end_res_trace_idx = trace_counter
        trace_counter += 1
    
    # Get total number of traces (paths + START_RES + END_RES)
    num_traces = trace_counter
    
    # Create slider steps for top K paths with proportional coloring
    # The slider value k means "show top k paths"
    top_k_steps = []
    for k in range(1, n_paths_to_show + 1):
        # Show top k paths (indices 0 to k-1)
        visibility_list = [True] * k + [False] * (n_paths_to_show - k)
        # Update colors proportionally for visible paths
        line_colors = []
        marker_colors = []
        for i in range(n_paths_to_show):
            if i < k:
                color = get_color_for_index(i, k)
            else:
                # Keep original color for hidden paths
                color = get_color_for_index(i, n_paths_to_show)
            line_colors.append(color)
            marker_colors.append(color)
        
        # Create args dict - properties as keys, lists as values (one per trace)
        # For restyle, we need to provide lists with one value per trace
        visible_values = []
        line_colors_values = []
        marker_colors_values = []
        
        for trace_idx in range(num_traces):
            if trace_idx in path_trace_indices:
                path_idx = path_trace_indices.index(trace_idx)
                if path_idx < len(visibility_list):
                    visible_values.append(visibility_list[path_idx])
                    line_colors_values.append(line_colors[path_idx])
                    marker_colors_values.append(marker_colors[path_idx])
                else:
                    visible_values.append(False)
                    # Use a default color for hidden paths
                    line_colors_values.append(get_color_for_index(0, n_paths_to_show))
                    marker_colors_values.append(get_color_for_index(0, n_paths_to_show))
            else:
                # Keep START_RES and END_RES traces visible
                # Use their existing colors (green for START, red for END)
                visible_values.append(True)
                if trace_idx == start_res_trace_idx:
                    line_colors_values.append('green')
                    marker_colors_values.append('green')
                elif trace_idx == end_res_trace_idx:
                    line_colors_values.append('red')
                    marker_colors_values.append('red')
                else:
                    line_colors_values.append('blue')  # fallback
                    marker_colors_values.append('blue')
        
        step = dict(
            args=[{
                "visible": visible_values,
                "line.color": line_colors_values,
                "marker.color": marker_colors_values
            }],
            label=f"{k}",
            method="restyle"
        )
        top_k_steps.append(step)
    
    # Create slider steps for path length filter
    length_steps = []
    for length_threshold in range(min_length, max_length + 1):
        # Filter paths by length <= threshold
        visibility_list = []
        for i in range(n_paths_to_show):
            row = df_sorted.iloc[i]
            visibility_list.append(row['length'] <= length_threshold)
        
        # Create args dict - properties as keys, lists as values (one per trace)
        visible_values = []
        
        for trace_idx in range(num_traces):
            if trace_idx in path_trace_indices:
                path_idx = path_trace_indices.index(trace_idx)
                if path_idx < len(visibility_list):
                    visible_values.append(visibility_list[path_idx])
                else:
                    visible_values.append(False)
            else:
                # Keep START_RES and END_RES traces visible
                visible_values.append(True)
        
        step = dict(
            args=[{
                "visible": visible_values
            }],
            label=f"{length_threshold}",
            method="restyle"
        )
        length_steps.append(step)
    
    # Add sliders
    sliders = [
        dict(
            active=n_paths_to_show - 1,  # Start with all paths visible
            currentvalue={"prefix": "Top K paths: "},
            pad={"t": 50},
            steps=top_k_steps,
        ),
        dict(
            active=len(length_steps) - 1,  # Start with all lengths
            currentvalue={"prefix": "Max path length: "},
            pad={"t": 100},
            steps=length_steps,
        )
    ]
    
    # Update layout
    fig.update_layout(
        title='3D Visualization of Random Walk Paths',
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
    
    output_path = output_file or str(Path(__file__).resolve().parent / "figures" / "paths.html")
    print(f"Saving visualization to {output_path}...")
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(output_path)


if __name__ == "__main__":
    import argparse
    default_output = str(Path(__file__).resolve().parent / "figures" / "paths.html")
    parser = argparse.ArgumentParser(description="Visualize random walk paths in 3D")
    parser.add_argument("pdb_file", help="Path to PDB file")
    parser.add_argument("csv_file", help="Path to walk CSV file")
    parser.add_argument("max_paths", nargs="?", type=int, default=None, help="Max number of paths to show (default: all)")
    parser.add_argument("--output", "-o", dest="output_file", default=default_output, help="Output HTML file path")
    args = parser.parse_args()

    visualize_paths_3d(args.pdb_file, args.csv_file, args.max_paths, args.output_file)

