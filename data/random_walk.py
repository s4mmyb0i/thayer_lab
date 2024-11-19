import os
import sys
import networkx as nx
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

from res_atom_dict import residue_atoms_dict

if len(sys.argv) != 6:
    print("Error!\nUsage: python3 random_walk.py <construct> <steps> <number_of_walks> <threshold_val_lower> <project_number>")
    sys.exit(1)

# Constants from input arguments
CONSTRUCT = sys.argv[1]
STEPS = int(sys.argv[2])
NUM_WALKS = int(sys.argv[3])
THRESHOLD_VAL_LOWER = float(sys.argv[4])
PROJECT_NUMBER = str(sys.argv[5])
THRESHOLD_VAL_UPPER = 0.999

# # Define binding pocket residues for each construct
# binding_pocket_residues_dict = {
#     "P": [
#     158, 159, 160, 161, 162, 163, 164, 165,
#     166, 167, 168, 169, 170, 171, 172, 173,
#     174, 175, 176, 177, 178, 179, 180, 181,
#     182, 183, 184, 185,
#     186, 187, 188, 189,
#     190, 191, 192, 193, 194, 195, 196, 197, 198,
#     281, 282, 283, 284, 285,
#     517, 518, 519, 520, 521,
#     565, 566, 567, 568],
#     "PL": [
#     196, 197, 198, 199, 200, 201, 202, 203,
#     204, 205, 206, 207, 208, 209, 210, 211,
#     212, 213, 214, 215, 216, 217, 218, 219,
#     220, 221, 222, 223,
#     224, 225, 226, 227,
#     228, 229, 230, 231, 232, 233, 234, 235, 236,
#     319, 320, 321, 322, 323,
#     555, 556, 557, 558, 559,
#     603, 604, 605, 606],
#     "AP": [
#     196, 197, 198, 199, 200, 201, 202, 203,
#     204, 205, 206, 207, 208, 209, 210, 211,
#     212, 213, 214, 215, 216, 217, 218, 219,
#     220, 221, 222, 223,
#     224, 225, 226, 227,
#     228, 229, 230, 231, 232, 233, 234, 235, 236,
#     319, 320, 321, 322, 323,
#     555, 556, 557, 558, 559,
#     603, 604, 605, 606],
#     "AL": [
#     196, 197, 198, 199, 200, 201, 202, 203,
#     204, 205, 206, 207, 208, 209, 210, 211,
#     212, 213, 214, 215, 216, 217, 218, 219,
#     220, 221, 222, 223,
#     224, 225, 226, 227,
#     228, 229, 230, 231, 232, 233, 234, 235, 236,
#     319, 320, 321, 322, 323,
#     555, 556, 557, 558, 559,
#     603, 604, 605, 606],
#     "default": [26, 27, 28, 29, 30, 31, 43, 76, 83]
# }

# # Update the global BINDING_POCKET_RESIDUES based on the construct
# BINDING_POCKET_RESIDUES = binding_pocket_residues_dict.get(CONSTRUCT, binding_pocket_residues_dict["default"])

# Get dictionary of (residue : int, atoms : list[int])
pdb_file = ""
res_atoms = {}
if CONSTRUCT == "P":
    pdb_file = "../../../../from_leticia/P_1BFE_RIP/1bfe.pdb"
    res_atoms = residue_atoms_dict(pdb_file)
if CONSTRUCT == "PL":
    pdb_file = "../../../../from_leticia/PL_1BE9_RIP/1be9.pdb"
    res_atoms = residue_atoms_dict(pdb_file)
if CONSTRUCT == "AP":
    pdb_file = "../../../../from_leticia/AP_PDZ3CDC42_RIP/pdz3cdc42_3.pdb"
    res_atoms = residue_atoms_dict(pdb_file)
if CONSTRUCT == "APL":
    pdb_file = "../../../../from_leticia/APL_PDZ3_RIP/pdz3lc_3.pdb"
    res_atoms = residue_atoms_dict(pdb_file)


# Define the list of binding pocket residues
binding_pocket_residues = [326, 327, 328, 329, 330, 331, 343, 376, 383]

# Populate BINDING_POCKET_RESIDUES with atoms for each residue in the binding pocket
BINDING_POCKET_RESIDUES = [
    atom for residue in binding_pocket_residues if residue in res_atoms
    for atom in res_atoms[residue]
]

"""
Termination codes:
-1: hits binding pocket residue
-2: dead end
-3: circular path
-4: max steps
"""

# Directories for GraphML files and output
graphml_dir = f'../../graphml/{CONSTRUCT}/'  # Change this to your GraphML directory
output_dir = f'../../{PROJECT_NUMBER}/data/{CONSTRUCT}/random_walk'  # Change this to your desired output directory
os.makedirs(output_dir, exist_ok=True)

def load_graph_from_graphml(file_path: str):
    """ Load graphml file """
    return nx.read_graphml(file_path, node_type=int, edge_key_type=float)

def random_walk_weighted(graph, start_node, steps):
    """ Perform random walks based on edge weights """
    path = [start_node]
    current_node = start_node
    chosen_prob = []
    termination_code = None

    for step in range(steps):
        if current_node in BINDING_POCKET_RESIDUES:
            termination_code = -1
            break
        
        neighbors = list(graph.neighbors(current_node))
        if not neighbors:
            termination_code = -2
            break
        
        weights = np.array([graph[current_node][neighbor]['weight'] for neighbor in neighbors])
        possible_prob = ((weights - THRESHOLD_VAL_LOWER) / (THRESHOLD_VAL_UPPER - THRESHOLD_VAL_LOWER) + 0.001)
        possible_prob /= possible_prob.sum()

        next_node = np.random.choice(neighbors, p=possible_prob)
        if next_node in path:
            termination_code = -3
            break

        log_prob = np.log(possible_prob[neighbors.index(next_node)])
        chosen_prob.append(-log_prob if log_prob != 0 else 0.0)
        path.append(int(next_node))
        current_node = next_node
    
    sum_log_prob = np.sum(chosen_prob)
    if len(chosen_prob) == 0:
        chosen_prob.append("0.0")

    if termination_code is None:
        termination_code = -4
    
    path.append(termination_code)
    
    return path, chosen_prob, sum_log_prob

# def write_walks_to_file(all_paths, file_path):
#     """ Writes random walks to a .csv file """
#     try:
#         with open(file_path, 'w') as f:
#             f.write("Start Node; Full Path; Probabilities (Log); Sum of Log Probabilities\n")
#             for start_node, paths in all_paths.items():
#                 for full_path, probabilities, sum_log_prob in paths:
#                     full_path_str = ", ".join(map(str, full_path))
#                     chosen_probabilities_str = ", ".join(map(str, probabilities))
#                     f.write(f"{start_node}; {full_path_str}; {chosen_probabilities_str}; {sum_log_prob}\n")
#     except Exception as e:
#         print(f"Error writing to file {file_path}: {e}")

def write_walks_to_file_for_residue(residue, all_paths_for_residue, output_dir):
    """Writes all random walks for each atom in the residue to one file."""
    file_path = os.path.join(output_dir, f"{residue}.csv")  # Filename for each residue
    
    try:
        with open(file_path, 'w') as f:
            # Write headers to the file
            f.write("Residue; Atom; Full Path; Probabilities (Log); Sum of Log Probabilities\n")
            
            # Write all the random walk data for each atom in the residue
            for atom, paths in all_paths_for_residue:
                for full_path, probabilities, sum_log_prob in paths:
                    full_path_str = ", ".join(map(str, full_path))
                    chosen_probabilities_str = ", ".join(map(str, probabilities))
                    f.write(f"{residue}; {atom}; {full_path_str}; {chosen_probabilities_str}; {sum_log_prob}\n")
                    
        print(f"Random walks for residue {residue} written to {file_path}")
    except Exception as e:
        print(f"Error writing to file {file_path}: {e}")

# def perform_random_walks(G, start_node):
#     """ Perform random walks for a single start node """
#     paths = []
#     for _ in range(NUM_WALKS):
#         full_path, chosen_probabilities, sum_log_prob = random_walk_weighted(G, start_node, STEPS)
#         paths.append((full_path, chosen_probabilities, sum_log_prob))
#     return start_node, paths

def perform_random_walks(G, residue, atoms):
    """ Perform random walks for every atom in a given residue """
    all_paths_for_residue = []
    for atom in atoms:
        paths_for_atom = []
        for _ in range(NUM_WALKS):
            full_path, chosen_probabilities, sum_log_prob = random_walk_weighted(G, atom, STEPS)
            paths_for_atom.append((full_path, chosen_probabilities, sum_log_prob))
        all_paths_for_residue.append((atom, paths_for_atom))
    return residue, all_paths_for_residue

# def process_single_graph(graphml_file, output_dir):
#     """ Process a single GraphML file: load graph, perform random walks, and save results """
#     G = load_graph_from_graphml(graphml_file)
    
#     # Extract start node from filename (assuming it's numeric)
#     start_node = int(os.path.splitext(os.path.basename(graphml_file))[0])
    
#     # Perform random walks
#     start_node, paths = perform_random_walks(G, start_node)
    
#     # Store paths
#     all_paths = {start_node: paths}
    
#     # Write results to file
#     output_file = os.path.join(output_dir, f"{start_node}.csv")
#     write_walks_to_file(all_paths, output_file)
#     print(f"Finished processing {graphml_file}")
    
#     return output_file

def process_single_graph(graphml_file, output_dir):
    """ Process a single GraphML file: load graph, perform random walks, and save results """
    G = load_graph_from_graphml(graphml_file)
    
    # Extract residue number from filename (assuming it's numeric)
    residue = int(os.path.splitext(os.path.basename(graphml_file))[0]) + 300  # Adjust residue number as needed

    # Get the atoms for the current residue from res_atoms
    atoms_for_residue = res_atoms.get(residue, [])
    
    if not atoms_for_residue:
        print(f"Skipping residue {residue} because no atoms found.")
        return None

    # Perform random walks for all atoms in the residue
    residue, all_paths_for_residue = perform_random_walks(G, residue, atoms_for_residue)

    # Write all random walks for the residue to a single CSV file
    write_walks_to_file_for_residue(residue, all_paths_for_residue, output_dir)
    
    return residue  # Optionally return residue number for logging or further use

def process_graphmls_in_parallel(graphml_dir, output_dir):
    """ Process all GraphML files in the directory in parallel """
    graphml_files = [os.path.join(graphml_dir, file) for file in os.listdir(graphml_dir) if file.endswith('.graphml')]
    
    with ProcessPoolExecutor() as executor:
        future_to_file = {executor.submit(process_single_graph, graphml_file, output_dir): graphml_file for graphml_file in graphml_files}
        for future in as_completed(future_to_file):
            graphml_file = future_to_file[future]
            try:
                result = future.result()
            except Exception as exc:
                print(f"Error processing {graphml_file}: {exc}")

def main():
    print("1")
    # Process GraphML files in parallel
    process_graphmls_in_parallel(graphml_dir, output_dir)

if __name__ == "__main__":
    main()