import pandas as pd
import numpy as np
import plotly.graph_objects as go
from collections import Counter
import argparse
from utils import extract_xyz_dataframe

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


def create_heatmap(pdb_file, csv_file, output_html='heatmap_3d.html'):
    """
    Create an interactive 3D heatmap of atoms visited in random walks.
    """
    # Extract coordinates from PDB - keep all atoms, not grouped by residue
    coord_df = extract_xyz_dataframe(pdb_file)
    
    # Parse walks and count node visits
    paths = parse_walks(csv_file)
    
    # Count how many times each node (atom_id) is visited
    all_nodes = []
    for path in paths:
        all_nodes.extend(path)
    
    node_counts = Counter(all_nodes)
    
    # Create a mapping of atom_id to coordinates
    # Assuming the dataframe has an atom_id column (or index represents atom_id)
    # Adjust column name based on your actual DataFrame structure
    if 'atom_id' not in coord_df.columns:
        # If atom_id is the index, reset it
        coord_df = coord_df.reset_index()
        if 'index' in coord_df.columns:
            coord_df = coord_df.rename(columns={'index': 'atom_id'})
    
    # Merge counts with coordinates
    heatmap_data = []
    for atom_id, count in node_counts.items():
        coords = coord_df[coord_df['atom_id'] == atom_id]
        if not coords.empty:
            heatmap_data.append({
                'atom_id': atom_id,
                'count': count,
                'x': coords['x'].values[0],
                'y': coords['y'].values[0],
                'z': coords['z'].values[0],
                'res_id': coords['res_id'].values[0] if 'res_id' in coords.columns else 'N/A',
                'atom_name': coords['atom_name'].values[0] if 'atom_name' in coords.columns else 'N/A'
            })
    
    heatmap_df = pd.DataFrame(heatmap_data)
    
    if heatmap_df.empty:
        print("No matching atoms found between PDB and walk data!")
        return
    
    # Normalize counts for color scaling
    max_count = heatmap_df['count'].max()
    min_count = heatmap_df['count'].min()
    
    # Create traces for different visit count thresholds
    # Generate threshold values (max 50 steps for performance)
    count_range = max_count - min_count
    max_steps = min(50, int(count_range) + 1)
    
    if count_range <= max_steps:
        threshold_values = list(range(int(min_count), int(max_count) + 1))
    else:
        threshold_values = [int(min_count + i * count_range / (max_steps - 1)) 
                           for i in range(max_steps)]
        if threshold_values[-1] != int(max_count):
            threshold_values[-1] = int(max_count)
    
    # Create figure
    fig = go.Figure()
    
    # Create a trace for each threshold
    for i, threshold in enumerate(threshold_values):
        # Filter atoms by threshold
        filtered_df = heatmap_df[heatmap_df['count'] >= threshold]
        
        trace = go.Scatter3d(
            x=filtered_df['x'],
            y=filtered_df['y'],
            z=filtered_df['z'],
            mode='markers',
            marker=dict(
                size=6,
                color=filtered_df['count'],
                colorscale='Viridis',
                colorbar=dict(title="Visit Count"),
                cmin=min_count,
                cmax=max_count,
                showscale=True,
                line=dict(width=0.5, color='white')
            ),
            text=[f"Atom ID: {aid}<br>Atom: {aname}<br>Residue: {rid}<br>Visits: {cnt}<br>X: {x:.2f}<br>Y: {y:.2f}<br>Z: {z:.2f}" 
                  for aid, aname, rid, cnt, x, y, z in zip(
                      filtered_df['atom_id'], 
                      filtered_df['atom_name'],
                      filtered_df['res_id'],
                      filtered_df['count'],
                      filtered_df['x'],
                      filtered_df['y'],
                      filtered_df['z'])],
            hoverinfo='text',
            name='Atoms',
            visible=(i == 0)  # Only first trace visible initially
        )
        fig.add_trace(trace)
    
    # Create slider steps
    steps = []
    for i, threshold in enumerate(threshold_values):
        filtered_count = len(heatmap_df[heatmap_df['count'] >= threshold])
        step = dict(
            method="update",
            args=[{"visible": [j == i for j in range(len(threshold_values))]}],
            label=str(threshold),
            value=str(threshold)
        )
        steps.append(step)
    
    sliders = [dict(
        active=0,
        currentvalue={"prefix": "Min visit count: ", "visible": True},
        pad={"t": 50, "b": 10},
        steps=steps
    )]
    
    # Update layout
    fig.update_layout(
        title=f'3D Heatmap of Random Walk Atom Visits<br>Total Walks: {len(paths)}, Unique Atoms: {len(node_counts)}',
        scene=dict(
            xaxis_title='X Coordinate (Å)',
            yaxis_title='Y Coordinate (Å)',
            zaxis_title='Z Coordinate (Å)',
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=1.5)
            )
        ),
        width=1000,
        height=800,
        hovermode='closest',
        sliders=sliders
    )
    
    # Save to HTML file
    fig.write_html(output_html)
    print(f"Interactive 3D heatmap saved to {output_html}")
    
    # Print statistics
    print(f"\nStatistics:")
    print(f"Total valid walks: {len(paths)}")
    print(f"Unique atoms visited: {len(node_counts)}")
    print(f"Most visited atom: Atom ID {max(node_counts, key=node_counts.get)} ({max(node_counts.values())} visits)")
    print(f"Least visited atoms: {min(node_counts.values())} visits")
    print(f"Average visits per atom: {sum(node_counts.values()) / len(node_counts):.2f}")
    
    return fig


if __name__ == "__main__":
    from pathlib import Path
    default_output = Path(__file__).resolve().parent / "figures" / "heatmap_3d.html"
    parser = argparse.ArgumentParser(description='Create 3D heatmap of random walk atoms')
    parser.add_argument('pdb_file', help='Path to PDB file')
    parser.add_argument('csv_file', help='Path to CSV file containing walks')
    parser.add_argument('--output', default=str(default_output), help='Output HTML file path')
    args = parser.parse_args()
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    create_heatmap(args.pdb_file, args.csv_file, args.output)









# """

# Create a 3D heatmap visualization of nodes hit in random walks.

# Takes a PDB file (for xyz coordinates) and a CSV file (contains the walks),
# and displays a Plotly heatmap showing which nodes are visited most frequently.

# Usage:
#     python3 analysis/node_heatmap.py <pdb_file> <csv_file>

# Example:
#     python3 analysis/node_heatmap.py /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/pdb/P.pdb "/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/random_walk/P_allos->lig_walks.csv"
# """

# import sys
# import pandas as pd
# import numpy as np
# import plotly.graph_objects as go
# from collections import Counter
# from itertools import chain
# from typing import Dict, List, Tuple, Set
# from utils import parse_path, create_atom_coord_map, get_residue_atom_sets


# def count_node_visits(csv_file: str) -> Counter:
#     """
#     Count how many times each node (atom) is visited across all walks.
    
#     Args:
#         csv_file (str): Path to the CSV file with random walk data
        
#     Returns:
#         Counter: Dictionary-like object with atom_id as key and visit count as value
#     """
#     print("Loading CSV file...")
#     df = pd.read_csv(csv_file, sep=';')
    
#     print(f"Found {len(df)} total walks")
    
#     # Filter to only show paths with termination_code == 1 FIRST (before parsing)
#     print("Filtering paths with termination_code == 1...")
#     df = df[df['termination_code'] == 1]
#     print(f"Found {len(df)} paths with termination_code == 1")
    
#     # Filter out empty paths before parsing
#     df = df[df['path'].notna() & (df['path'] != '')]
    
#     # Parse only the paths we need
#     print("Parsing paths...")
#     df['parsed_path'] = df['path'].apply(parse_path)
    
#     # Filter out empty parsed paths
#     df = df[df['parsed_path'].apply(len) > 0]
    
#     # Count visits to each node using itertools.chain for efficiency
#     print("Counting node visits...")
#     node_counts = Counter(chain.from_iterable(df['parsed_path']))
    
#     print(f"Found {len(node_counts)} unique nodes visited")
#     print(f"Total node visits: {sum(node_counts.values())}")
    
#     return node_counts


# def create_heatmap(pdb_file: str, csv_file: str):
#     """
#     Create and display a 3D heatmap of nodes hit in walks.
    
#     Args:
#         pdb_file (str): Path to the PDB file
#         csv_file (str): Path to the CSV file with random walk data
#     """
#     # Count node visits
#     node_counts = count_node_visits(csv_file)
    
#     # Get coordinates for all nodes
#     print("Extracting coordinates from PDB file...")
#     coord_map = create_atom_coord_map(pdb_file)
    
#     # Get residue to atom mapping for START_RES and END_RES
#     print("Mapping residues to atoms for START_RES and END_RES...")
#     start_atoms, end_atoms = get_residue_atom_sets(pdb_file)
#     print(f"START_RES atoms: {len(start_atoms)}")
#     print(f"END_RES atoms: {len(end_atoms)}")
    
#     # Pre-allocate sets for faster lookup
#     start_atoms_set = set(start_atoms)
#     end_atoms_set = set(end_atoms)
    
#     # Separate nodes into categories using list comprehensions (faster than appending in loop)
#     start_data = {'x': [], 'y': [], 'z': [], 'counts': [], 'atom_ids': []}
#     end_data = {'x': [], 'y': [], 'z': [], 'counts': [], 'atom_ids': []}
#     other_data = {'x': [], 'y': [], 'z': [], 'counts': [], 'atom_ids': []}
    
#     # Filter nodes that have coordinates and categorize in one pass
#     for atom_id, count in node_counts.items():
#         if atom_id not in coord_map:
#             continue
#         x, y, z = coord_map[atom_id]
#         if atom_id in start_atoms_set:
#             start_data['x'].append(x)
#             start_data['y'].append(y)
#             start_data['z'].append(z)
#             start_data['counts'].append(count)
#             start_data['atom_ids'].append(atom_id)
#         elif atom_id in end_atoms_set:
#             end_data['x'].append(x)
#             end_data['y'].append(y)
#             end_data['z'].append(z)
#             end_data['counts'].append(count)
#             end_data['atom_ids'].append(atom_id)
#         else:
#             other_data['x'].append(x)
#             other_data['y'].append(y)
#             other_data['z'].append(z)
#             other_data['counts'].append(count)
#             other_data['atom_ids'].append(atom_id)
    
#     total_nodes = len(start_data['x']) + len(end_data['x']) + len(other_data['x'])
#     if total_nodes == 0:
#         print("Error: No nodes with coordinates found!")
#         return
    
#     # Get visit count range for consistent colorbar scaling
#     all_counts = start_data['counts'] + end_data['counts'] + other_data['counts']
#     min_count = min(all_counts) if all_counts else 1
#     max_count = max(all_counts) if all_counts else 1
    
#     print(f"\nTotal nodes: {total_nodes}")
#     print(f"START_RES nodes: {len(start_data['x'])}")
#     print(f"END_RES nodes: {len(end_data['x'])}")
#     print(f"Other nodes: {len(other_data['x'])}")
#     print(f"Visit count range: {min_count} to {max_count}")
    
#     # Define custom colorscales
#     # Green to white for START_RES
#     green_to_white = [[0, 'rgb(0, 128, 0)'], [1, 'rgb(255, 255, 255)']]
    
#     # Red to white for END_RES
#     red_to_white = [[0, 'rgb(128, 0, 0)'], [1, 'rgb(255, 255, 255)']]
    
#     # Light blue to dark blue for other nodes
#     light_blue_to_dark_blue = [[0, 'rgb(173, 216, 230)'], [1, 'rgb(0, 0, 139)']]
    
#     # Create figure
#     fig = go.Figure()
    
#     # First, add START_RES and END_RES traces (always visible, not affected by slider)
#     trace_counter = 0
#     start_trace_idx = None
#     end_trace_idx = None
    
#     # Add START_RES trace (green to white) - always visible
#     if start_data['x']:
#         marker_dict = {
#             'size': 10,
#             'color': start_data['counts'],
#             'colorscale': green_to_white,
#             'colorbar': dict(
#                 title=dict(text="Visit Count", side="right"),
#                 x=1.15
#             ),
#             'cmin': min_count,
#             'cmax': max_count,
#             'opacity': 0.9,
#             'line': dict(width=1, color='darkgreen')
#         }
        
#         trace = go.Scatter3d(
#             x=start_data['x'],
#             y=start_data['y'],
#             z=start_data['z'],
#             mode='markers',
#             name='START_RES (Allosteric Binding Pocket)',
#             marker=marker_dict,
#             text=[f'Atom {int(aid)}<br>Visits: {int(cnt)}<br>Type: START_RES' 
#                   for aid, cnt in zip(start_data['atom_ids'], start_data['counts'])],
#             hovertemplate='%{text}<extra></extra>',
#             visible=True,
#             showlegend=True
#         )
#         fig.add_trace(trace)
#         start_trace_idx = trace_counter
#         trace_counter += 1
    
#     # Add END_RES trace (red to white) - always visible
#     if end_data['x']:
#         marker_dict = {
#             'size': 10,
#             'color': end_data['counts'],
#             'colorscale': red_to_white,
#             'colorbar': dict(
#                 title=dict(text="Visit Count", side="right"),
#                 x=1.25
#             ),
#             'cmin': min_count,
#             'cmax': max_count,
#             'opacity': 0.9,
#             'line': dict(width=1, color='darkred')
#         }
        
#         trace = go.Scatter3d(
#             x=end_data['x'],
#             y=end_data['y'],
#             z=end_data['z'],
#             mode='markers',
#             name='END_RES (Ligand Binding Pocket)',
#             marker=marker_dict,
#             text=[f'Atom {int(aid)}<br>Visits: {int(cnt)}<br>Type: END_RES' 
#                   for aid, cnt in zip(end_data['atom_ids'], end_data['counts'])],
#             hovertemplate='%{text}<extra></extra>',
#             visible=True,
#             showlegend=True
#         )
#         fig.add_trace(trace)
#         end_trace_idx = trace_counter
#         trace_counter += 1
    
#     # Create filtered traces for different frequency thresholds (only for other nodes)
#     # Limit to max 50 traces for better performance (reduced from 100)
#     count_range = max_count - min_count
#     max_traces = min(50, int(max_count - min_count + 1))
    
#     # Create threshold values - use fewer traces but allow all integer values in slider
#     if count_range <= max_traces:
#         # Small range: create a trace for each integer value
#         threshold_values = list(range(int(min_count), int(max_count) + 1))
#     else:
#         # Large range: create evenly spaced thresholds
#         threshold_values = [int(min_count + i * count_range / (max_traces - 1)) 
#                            for i in range(max_traces)]
#         # Ensure max_count is included
#         if threshold_values[-1] != int(max_count):
#             threshold_values[-1] = int(max_count)
    
#     # Convert to numpy arrays for faster filtering (only if we have other nodes)
#     traces_per_threshold = []
#     threshold_to_trace_map = {}  # Map any threshold value to the nearest trace index
    
#     if other_data['x']:
#         other_x = np.array(other_data['x'])
#         other_y = np.array(other_data['y'])
#         other_z = np.array(other_data['z'])
#         other_counts = np.array(other_data['counts'])
#         other_atom_ids = np.array(other_data['atom_ids'])
#     else:
#         # No other nodes, create empty arrays
#         other_x = np.array([])
#         other_y = np.array([])
#         other_z = np.array([])
#         other_counts = np.array([])
#         other_atom_ids = np.array([])
    
#     for threshold in threshold_values:
#         # Filter using numpy for speed
#         if len(other_counts) > 0:
#             mask = other_counts >= threshold
#             other_filtered_x = other_x[mask]
#             other_filtered_y = other_y[mask]
#             other_filtered_z = other_z[mask]
#             other_filtered_counts = other_counts[mask]
#             other_filtered_atom_ids = other_atom_ids[mask]
#         else:
#             other_filtered_x = np.array([])
#             other_filtered_y = np.array([])
#             other_filtered_z = np.array([])
#             other_filtered_counts = np.array([])
#             other_filtered_atom_ids = np.array([])
        
#         threshold_traces = []
        
#         # Add filtered other nodes trace (light blue to dark blue)
#         if len(other_filtered_x) > 0 and len(other_data['x']) > 0:
#             marker_dict = {
#                 'size': 8,
#                 'color': other_filtered_counts.tolist(),
#                 'colorscale': light_blue_to_dark_blue,
#                 'cmin': min_count,
#                 'cmax': max_count,
#                 'opacity': 0.8,
#                 'line': dict(width=0.5, color='rgba(0,0,0,0.3)'),
#                 'colorbar': dict(
#                     title=dict(text="Visit Count", side="right"),
#                     x=1.35
#                 ),
#                 'showscale': True  # Show colorbar for all traces (only visible one will display)
#             }
            
#             # Create hover text efficiently using list comprehension
#             hover_text = [f'Atom {int(aid)}<br>Visits: {int(cnt)}<br>Type: Other' 
#                          for aid, cnt in zip(other_filtered_atom_ids, other_filtered_counts)]
            
#             trace = go.Scatter3d(
#                 x=other_filtered_x.tolist(),
#                 y=other_filtered_y.tolist(),
#                 z=other_filtered_z.tolist(),
#                 mode='markers',
#                 name='Other Nodes',
#                 marker=marker_dict,
#                 text=hover_text,
#                 hovertemplate='%{text}<extra></extra>',
#                 visible=(threshold == min_count),  # Only first threshold is visible initially
#                 showlegend=(threshold == min_count)  # Only show legend for first threshold
#             )
#             fig.add_trace(trace)
#             threshold_traces.append(trace_counter)
#             trace_counter += 1
        
#         traces_per_threshold.append(threshold_traces)
#         threshold_to_trace_map[threshold] = len(traces_per_threshold) - 1
    
#     # Create mapping for all integer values to nearest trace
#     all_thresholds = list(range(int(min_count), int(max_count) + 1))
#     for thresh in all_thresholds:
#         if thresh not in threshold_to_trace_map:
#             # Find nearest threshold value
#             nearest = min(threshold_values, key=lambda x: abs(x - thresh))
#             threshold_to_trace_map[thresh] = threshold_to_trace_map[nearest]
    
#     # Create slider steps for all integer values (smooth selection)
#     frequency_steps = []
#     total_traces = trace_counter
#     all_thresholds = list(range(int(min_count), int(max_count) + 1))
    
#     for threshold in all_thresholds:
#         visible_list = [False] * total_traces
        
#         # Always show START_RES and END_RES traces (not affected by slider)
#         if start_trace_idx is not None:
#             visible_list[start_trace_idx] = True
#         if end_trace_idx is not None:
#             visible_list[end_trace_idx] = True
        
#         # Show other node traces for this threshold (use mapping to find correct trace)
#         trace_idx_for_threshold = threshold_to_trace_map[threshold]
#         for trace_idx in traces_per_threshold[trace_idx_for_threshold]:
#             visible_list[trace_idx] = True
        
#         # Create step for this threshold
#         step = dict(
#             args=[{
#                 "visible": visible_list
#             }],
#             label=f"{threshold}",
#             method="restyle"
#         )
#         frequency_steps.append(step)
    
#     # Update layout with slider
#     fig.update_layout(
#         title='3D Heatmap of Nodes Hit in Random Walks',
#         scene=dict(
#             xaxis_title='X (Angstrom)',
#             yaxis_title='Y (Angstrom)',
#             zaxis_title='Z (Angstrom)',
#             aspectmode='data'
#         ),
#         height=800,
#         width=1200,
#         sliders=[dict(
#             active=0,  # Start with all nodes visible (min threshold)
#             currentvalue={"prefix": "Min visit count: "},
#             pad={"t": 50},
#             steps=frequency_steps,
#         )]
#     )
    
#     print(f"\nDisplaying heatmap with {total_nodes} nodes...")
#     fig.show()


# if __name__ == "__main__":
#     if len(sys.argv) != 3:
#         print("Usage: python node_heatmap.py <pdb_file> <csv_file>")
#         sys.exit(1)
    
#     pdb_file = sys.argv[1]
#     csv_file = sys.argv[2]
    
#     create_heatmap(pdb_file, csv_file)

