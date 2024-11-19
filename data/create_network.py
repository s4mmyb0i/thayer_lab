"""
Samvit Prem Singhal | spsinghal@wesleyan.edu
Creates a weighted network where nodes are amino acids and edges are created when above a DCCM coefficient threshold. Edges are given a weight by the absolute value of the coefficient. These networks are also exported to a .graphml file.

The script takes in 2 system inputs: <construct> <number_of_residues> <threshold_val_lower>

Input:  dccm in .dat format
Output: graphml networks
"""

import os
import sys
import numpy as np
import networkx as nx
from concurrent.futures import ProcessPoolExecutor, as_completed

if len(sys.argv) != 4:
    print("Usage: python3 create_network.py <construct> <number_of_residues> <threshold_val_lower>")
    sys.exit(1)

# Constants from input arguments
CONSTRUCT = sys.argv[1]
NUM_RESIDUES = int(sys.argv[2])
THRESHOLD_VAL_LOWER = float(sys.argv[3])
THRESHOLD_VAL_UPPER = 0.999

DCCM_DIR = os.path.join('../../../dccms', CONSTRUCT)
EXPORT_DIR = os.path.join('../../graphml', CONSTRUCT)
os.makedirs(EXPORT_DIR, exist_ok=True)

def load_dccm(file_path: str):
    """Load the DCCM data from a given file."""
    return np.loadtxt(file_path)

def create_networkx(dccm_data):
    """Create a NetworkX DiGraph from the DCCM data."""
    G = nx.DiGraph()
    
    for i in range(1, NUM_RESIDUES):
        G.add_node(i)
    
    for i in range(1, NUM_RESIDUES):
        for j in range(1, NUM_RESIDUES):
            if i != j:
                weight = abs(dccm_data[i, j])
                if THRESHOLD_VAL_LOWER < weight < THRESHOLD_VAL_UPPER:
                    G.add_edge(i, j, weight=weight)
    
    return G

def export_to_graphml(G, file_path: str):
    """Export NetworkX DiGraph to GraphML format."""
    nx.write_graphml(G, file_path)
    print(f"Exported NetworkX graph to {file_path}")

def process_single_dccm(dccm_file: str):
    """Process a single DCCM file."""
    file_index = os.path.splitext(os.path.basename(dccm_file))[0]
    export_file = os.path.join(EXPORT_DIR, f'{file_index}.graphml')
    
    dccm_data = load_dccm(dccm_file)
    G = create_networkx(dccm_data)
    export_to_graphml(G, export_file)
    return f"Processed {dccm_file}"

def main():
    dccm_files = [os.path.join(DCCM_DIR, f) for f in os.listdir(DCCM_DIR) if f.endswith('.dat')]
    
    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(process_single_dccm, file): file for file in dccm_files}
        
        for future in as_completed(futures):
            result = future.result()
            print(result)

if __name__ == '__main__':
    main()