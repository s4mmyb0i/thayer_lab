import sys
import os
import networkx as nx
import scipy.sparse

if len(sys.argv) != 3:
    print("Error!\nUsage: python3 graphml_to_sparse_matrix.py <graphml_dir> <scipy_dir>")
    sys.exit(1)

GRAPHML_DIR = sys.argv[1]
SCIPY_DIR = sys.argv[2]

# Create output directory if it doesn't exist
os.makedirs(SCIPY_DIR, exist_ok=True)

def convert_to_sparse_matrix(G: nx.Graph):
    """Convert a NetworkX graph to a sparse adjacency matrix (CSR format)."""
    nodes = list(G.nodes())
    node_to_index = {node: i for i, node in enumerate(nodes)}
    
    # Convert to sparse CSR matrix
    sparse_adj_matrix = nx.to_scipy_sparse_array(G, nodelist=nodes, format="csr")

    return sparse_adj_matrix, node_to_index, nodes

# Iterate over constructs (subdirectories) in the directory
for construct in os.listdir(GRAPHML_DIR):
    construct_path = os.path.join(GRAPHML_DIR, construct)
    
    if not os.path.isdir(construct_path):
        continue  # Skip non-directory files

    for filename in os.listdir(construct_path):
        if filename.endswith(".graphml"):
            file_path = os.path.join(construct_path, filename)
            print(f"Processing: {file_path}")

            # Load graph
            G = nx.read_graphml(file_path)
            print(f"Nodes: {len(G.nodes())}, Edges: {len(G.edges())}")

            # Convert to sparse matrix
            sparse_matrix, node_to_index, nodes = convert_to_sparse_matrix(G)

            # Create subdirectory for this construct in SCIPY_DIR
            construct_output_dir = os.path.join(SCIPY_DIR, construct)
            os.makedirs(construct_output_dir, exist_ok=True)

            # Save sparse matrix
            output_file = os.path.join(construct_output_dir, filename.replace(".graphml", ".npz"))
            scipy.sparse.save_npz(output_file, sparse_matrix)

            print(f"Saved sparse matrix to {output_file} (shape: {sparse_matrix.shape})")
