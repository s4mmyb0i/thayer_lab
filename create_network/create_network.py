"""
Samvit Prem Singhal | spsinghal@wesleyan.edu
Unified script: Creates a weighted graph from a DCCM matrix, builds a sparse transition matrix, and sets specified nodes as absorbing for Markov analyses and random walks.

Takes 7 command-line arguments:
    <matrix> <pdb_file> <threshold_val_lower> <graphml_dir> <adj_matrix_dir> <trans_matrix_dir> <trans_abs_matrix_dir>

Outputs:
    - .graphml weighted directed graph
    - adjacency matrix (.npz)
    - transition matrix (.npz)
    - transition matrix with absorbing nodes (.npz)
"""

import os
import sys
import numpy as np
import networkx as nx
from pathlib import Path
from scipy.sparse import csr_matrix, diags, save_npz
from networkx.convert_matrix import to_scipy_sparse_array
from utils.network_utils import find_num_atoms
from utils.target import get_target_nodes

if len(sys.argv) != 8:
    print("Usage: python3 create_network.py <matrix> <pdb_file> <threshold_val_lower> <graphml_dir> <adj_matrix_dir> <trans_matrix_dir> <trans_abs_matrix_dir>")
    sys.exit(1)

matrix_file         = str(sys.argv[1])
pdb_file            = str(sys.argv[2])
threshold_val_lower = float(sys.argv[3])
threshold_val_upper = 0.999
graphml_dir         = str(sys.argv[4])
adj_matrix_dir      = str(sys.argv[5])
trans_matrix_dir    = str(sys.argv[6])
abs_dir             = str(sys.argv[7])

os.makedirs(graphml_dir, exist_ok=True)
os.makedirs(adj_matrix_dir, exist_ok=True)
os.makedirs(trans_matrix_dir, exist_ok=True)
os.makedirs(abs_dir, exist_ok=True)

_, absorption_nodes = get_target_nodes()

def load_matrix(matrix: str, pdb_file: str) -> np.ndarray:
    """
    Load a matrix from a .dat file and validate it.

    Args:
        matrix (str): Path to the matrix file.
        pdb_file (str): Path to the PDB file for validation and checks.

    Returns:
        np.ndarray: The loaded matrix data.
    """
    data = np.loadtxt(matrix)
    if data.shape[0] != data.shape[1]:
        raise ValueError("Matrix is not square.")
    if data.shape[0] != find_num_atoms(pdb_file):
        raise ValueError("Matrix size doesn't match number of atoms in PDB.")
    return data

def create_graph(matrix: np.ndarray, lower: float, upper: float) -> tuple[nx.DiGraph, csr_matrix]:
    """
    Create a weighted directed graph from a matrix and threshold.

    Args:
        matrix (np.ndarray): The input matrix.
        lower (float): Lower bound of the threshold range (inclusive).
        upper (float): Upper bound of the threshold range (inclusive).

    Returns:
        nx.DiGraph: The created graph.
        csr_matrix: The sparse adjacency matrix of the graph.
    """
    G = nx.DiGraph()
    # num_weakly_connected_components = nx.number_weakly_connected_components(G)

    # print("components: ", nx.number_strongly_connected_components(G))
    N = matrix.shape[0]
    for i in range(N):
        for j in range(N):
            if i != j:
                w = abs(matrix[i, j])
                if lower <= w <= upper:
                    G.add_edge(i, j, weight=w)
    G = nx.convert_node_labels_to_integers(G, ordering="sorted")
    A = to_scipy_sparse_array(G, weight="weight", format="csr")
    return G, A  # pyright: ignore[reportReturnType]

def to_transition_matrix(A: csr_matrix, self_loops: bool = True) -> csr_matrix:
    """
    Convert a weighted adjacency matrix into a transition matrix.

    This function takes a weighted adjacency matrix, A, and converts it into a
    transition matrix, P, by normalizing each row so that it sums to 1. If a
    row sums to 0, a self-loop is added with a weight of 1.0. The resulting
    transition matrix can be used to simulate random walks on the graph.

    Parameters
    ----------
    A : csr_matrix
        The weighted adjacency matrix to be converted.
    self_loops : bool, default=True
        If True, a self-loop is added to any row that sums to 0. If False,
        the row is left unchanged.

    Returns
    -------
    csr_matrix
        The transition matrix corresponding to the input matrix.
    """
    P = A.copy().astype(float)
    row_sum = np.asarray(P.sum(axis=1)).ravel()
    if self_loops:
        dangling = np.where(row_sum == 0)[0]
        P[dangling, dangling] = 1.0
        row_sum[dangling] = 1.0
    inv = np.zeros_like(row_sum)
    inv[row_sum > 0] = 1.0 / row_sum[row_sum > 0]
    return diags(inv) @ P

def set_absorbing(P: csr_matrix, absorbing: set[int]) -> csr_matrix:
    """
    Convert specified nodes in a transition matrix into absorbing states.

    Args:
        P (csr_matrix): The transition matrix to be modified.
        absorbing (set[int]): Indices of nodes to be set as absorbing.

    Returns:
        csr_matrix: The modified transition matrix with specified nodes set as absorbing.
    """

    P = P.tolil(copy=True) # type: ignore
    for node in absorbing:
        P.rows[node] = [node] # type: ignore
        P.data[node] = [1.0]
    return P.tocsr()

def set_pseudo_transitions(P: csr_matrix, absorbing: set[int], epsilon: float = 1e-4) -> csr_matrix:
    """
    Assign extremely low probabilities to nodes that cannot reach any absorbing state.
    This ensures that all rows remain stochastic, so Markovian analyses can work

    Args:
        P (csr_matrix): The transition matrix with absorbing nodes set.
        absorbing (set[int]): Indices of absorbing nodes.
        epsilon (float, optional): Very small probability assigned to unreachable nodes.

    Returns:
        csr_matrix: Modified transition matrix.
    """
    P_mod = P.tolil(copy=True)  # LIL is easier for row-wise manipulation
    G = nx.from_scipy_sparse_array(P_mod, create_using=nx.DiGraph)

    transient_nodes = set(range(P.shape[0])) - absorbing  # pyright: ignore[reportOptionalSubscript]
    for node in transient_nodes:
        if not any(desc in absorbing for desc in nx.descendants(G, node)):
            # Assign very low probability to all transitions except self-loop
            num_cols = P_mod.shape[1]  # pyright: ignore[reportOptionalSubscript]
            P_mod.rows[node] = [node]
            P_mod.data[node] = [1.0 - epsilon]
            
            # Optional: distribute epsilon equally among all other nodes
            other_nodes = list(range(P_mod.shape[1]))  # pyright: ignore[reportOptionalSubscript]
            other_nodes.remove(node)
            P_mod.rows[node].extend(other_nodes)
            P_mod.data[node].extend([epsilon / (num_cols - 1)] * (num_cols - 1))

    return P_mod.tocsr()  # pyright: ignore[reportReturnType]

def sanity_check(P: csr_matrix, absorbing: set[int]) -> bool:
    """
    Sanity checks for a transition matrix.

    Checks:
    - Rows sum to 1.0
    - Absorbing nodes have only self-loop of 1.0
    - Every transient state can reach at least one absorbing node

    Args:
        P (csr_matrix): The transition matrix to be checked.
        absorbing (set[int]): Indices of absorbing nodes.
    """
    is_safe = True
    row_sums = P.sum(axis=1).A1
    assert np.all(np.abs(row_sums - 1.0) < 1e-8), "Rows are not stochastic."
    for node in absorbing:
        if not (P[node, node] == 1.0 and P[node].count_nonzero() == 1): # type: ignore
            print(f"Warning: Node {node} not properly absorbing.")
            is_safe = False
    G = nx.from_scipy_sparse_array(P, create_using=nx.DiGraph)
    for t in set(range(P.shape[0])) - absorbing: # type: ignore
        if not any(r in absorbing for r in nx.descendants(G, t)):
            print(f"Transient node {t} can't reach absorbing states.")
            is_safe = False
    return is_safe

def export_all(base: str, G: nx.DiGraph, A: csr_matrix, P: csr_matrix, P_abs: csr_matrix):
    """
    Export a graph and its sparse matrices to the correct directories.

    Args:
        base (str): The base name of the files to be saved.
        G (nx.DiGraph): The graph to be saved.
        A (csr_matrix): The adjacency matrix to be saved.
        P (csr_matrix): The transition matrix to be saved.
        P_abs (csr_matrix): The absorbing transition matrix to be saved.

    """
    graphml_path = os.path.join(graphml_dir,       f"{base}.graphml")
    sparse_path  = os.path.join(adj_matrix_dir,    f"{base}_adj.npz")
    trans_path   = os.path.join(trans_matrix_dir,  f"{base}_trans.npz")
    absorb_path  = os.path.join(abs_dir,           f"{base}_trans_abs.npz")

    nx.write_graphml(G, graphml_path)
    save_npz(sparse_path, A)
    save_npz(trans_path, P)
    save_npz(absorb_path, P_abs)

    print(f"✅ Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    print(f"📁 Saved to: {graphml_path}, {sparse_path}, {trans_path}, {absorb_path}\n")

def main():
    base = Path(matrix_file).stem
    matrix = load_matrix(matrix_file, pdb_file)
    G, A = create_graph(matrix, threshold_val_lower, threshold_val_upper)
    P = to_transition_matrix(A, self_loops=True)
    P_abs = set_absorbing(P.copy(), absorption_nodes)
    P_abs = set_pseudo_transitions(P_abs, absorption_nodes)
    if not sanity_check(P_abs, absorption_nodes): print("This network is not Markov ready!")
    export_all(base, G, A, P, P_abs)

if __name__ == '__main__':
    main()
