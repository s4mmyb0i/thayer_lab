"""
This module provides functions for computing absorption probabilities for a Markov chain.

Usage:
python3 -m linalg.absorption_probs.absorption_probs \
  /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/matrices/trans_abs_matrix/P_trans_abs.npz \
  /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/linalg/absorption_probs

python3 -m linalg.absorption_probs.absorption_probs \
  /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/matrices/trans_abs_matrix/PL_trans_abs.npz \
  /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/linalg/absorption_probs

python3 -m linalg.absorption_probs.absorption_probs \
  /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/matrices/trans_abs_matrix/AP_trans_abs.npz \
  /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/linalg/absorption_probs

python3 -m linalg.absorption_probs.absorption_probs \
  /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/matrices/trans_abs_matrix/APL_trans_abs.npz \
  /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/linalg/absorption_probs

"""


from utils.target import get_target_nodes
import numpy as np
import pandas as pd

def filter_by_abs_corr_threshold(matrix: np.ndarray, threshold: float) -> np.ndarray:
    filtered = np.copy(matrix)
    filtered[np.abs(filtered) < threshold] = 0
    return filtered

def soften_threshold(matrix: np.ndarray, threshold: float, epsilon: float = 1e-6) -> np.ndarray:
    """Keep high correlations, but replace low ones with epsilon."""
    filtered = np.copy(matrix)
    low_mask = np.abs(filtered) < threshold
    filtered[low_mask] = epsilon
    return filtered

def create_np_matrix(file: str) -> np.ndarray:
    matrix = np.loadtxt(file)
    matrix = np.abs(matrix)
    threshold = 0.85
    filtered = filter_by_abs_corr_threshold(matrix, threshold)
    filtered = soften_threshold(filtered, threshold)
    return filtered

def make_markov_transition(matrix: np.ndarray, absorbing_nodes: set[int]) -> np.ndarray:

    def make_diagonals_zero(matrix: np.ndarray) -> np.ndarray:
        np.fill_diagonal(matrix, 0)
        return matrix

    def make_transition_matrix(matrix: np.ndarray) -> np.ndarray:
        return  matrix / matrix.sum(axis=1, keepdims=True)

    def set_absorbing(matrix: np.ndarray, absorbing_nodes: set[int]) -> np.ndarray:
        for absorbing_node in absorbing_nodes:
            matrix[absorbing_node, :] = 0
            matrix[absorbing_node][absorbing_node] = 1
        return matrix

    matrix = make_diagonals_zero(matrix)
    matrix = make_transition_matrix(matrix)
    matrix = set_absorbing(matrix, absorbing_nodes)
    return matrix

def calculate_absorption_prob(matrix: np.ndarray, absorbing_nodes: set[int]) -> np.ndarray:
    absorbing_nodes = absorbing_nodes
    all_nodes = np.arange(matrix.shape[0])
    transient_nodes = [i for i in all_nodes if i not in absorbing_nodes]
    
    Q = matrix[np.ix_(transient_nodes, transient_nodes)]  # transition to transition
    R = matrix[np.ix_(transient_nodes, list(absorbing_nodes))]  # transition to absorbing
    I = np.eye(len(Q))  # identity matrix of size Q
    F = np.linalg.inv(I - Q)  # fundamental matrix

    result = F @ R  # percentage of absorption probabilities
    return result

def absorption_dataframe(matrix, absorbing_nodes):
    all_nodes = np.arange(matrix.shape[0])
    transient_nodes = [i for i in all_nodes if i not in absorbing_nodes]
    probs = calculate_absorption_prob(matrix, absorbing_nodes)
    df = pd.DataFrame(
        probs,
        index=[f"Start {i}" for i in transient_nodes],  # pyright: ignore[reportArgumentType]
        columns=[f"Absorb {j}" for j in absorbing_nodes]  # pyright: ignore[reportArgumentType]
    )
    return df

matrix = create_np_matrix("/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/dccm/P.dat")
_, absorbing_nodes = get_target_nodes()

# matrix = np.array([
#     [1.0, 0.8, 0.0, 0.8, 0.0, 0.0],
#     [0.8, 1.0, 0.6, 0.0, 0.0, 0.0],
#     [0.0, 0.6, 1.0, 0.2, 0.0, 0.0],
#     [0.8, 0.0, 0.2, 1.0, 0.4, 0.0],
#     [0.0, 0.0, 0.0, 0.4, 1.0, 0.8],
#     [0.0, 0.0, 0.0, 0.0, 0.8, 1.0]
# ])

# absorbing_nodes = [1, 2, 5] # one indexed
# absorbing_nodes = set([0, 1, 3]) # zero indexed

matrix = make_markov_transition(matrix, absorbing_nodes)
df = absorption_dataframe(matrix, absorbing_nodes)

# Save DataFrame to CSV
output_path = "./linalg/absorption_probs/absorption_probabilities.csv"
df.to_csv(output_path, index=True)

print(f"\nAbsorption probability matrix saved to:\n{output_path}")


# Get index (row label) and column (absorbing node)
max_idx = df.stack().idxmax()
max_val = df.stack().max()

print(f"Largest value: {max_val:.6f}")
print(f"Located at: Row = {max_idx[0]}, Column = {max_idx[1]}")  # pyright: ignore[reportIndexIssue]


"""
The resulting answer:
[[0.21052632 0.78947368]
 [0.84210526 0.15789474]
 [0.84210526 0.15789474]
 [0.84210526 0.15789474]]
 
The column is the absorbing nodes
The row is the transient nodes
The ijth value is saying, starting from i, what is the probability we land in j,
where i is the transient node, and j is the absorbing node
"""










# import numpy as np
# import pandas as pd
# from scipy.sparse import load_npz, csr_matrix, eye, save_npz
# from scipy.sparse.linalg import inv
# import networkx as nx
# from utils.target import get_target_nodes
# from utils.network_utils import check_transition_matrix, get_reachable_transients
# import logging
# import sys
# from pathlib import Path

# if len(sys.argv) != 3:
#     logging.error("Error! Incorrect number of arguments.")
#     print("Error!\nUsage: python3 absorption_probs.py <matrix> <output_dir>")
#     sys.exit(1)

# matrix = str(sys.argv[1])
# output_dir = str(sys.argv[2])
# SCRIPT_DIR = Path(__file__).resolve().parent
# Path(output_dir).mkdir(parents=True, exist_ok=True)

# def setup_logger(matrix_file: str) -> None:
#     log_name = Path(matrix_file).stem + "_probs.log"
#     log_path = SCRIPT_DIR / log_name
#     logging.basicConfig(
#         filename=log_path,
#         level=logging.INFO,
#         format="%(asctime)s - %(levelname)s - %(message)s",
#         datefmt="%Y-%m-%d %H:%M:%S"
#     )

# def compute_absorption_probabilities(P: csr_matrix, absorbing: set[int]) -> tuple[csr_matrix, np.ndarray, np.ndarray]:
#     """
#     Compute the absorption probabilities of the given transition matrix.

#     The absorption probabilities are the probabilities of absorption into each
#     absorbing state, given that the random walk starts at a transient state.

#     Args:
#         P (csr_matrix): The transition matrix.
#         absorbing (set[int]): The indices of the absorbing nodes.

#     Returns:
#         tuple[csr_matrix, np.ndarray, np.ndarray]:
#             - absorption_probabilities (csr_matrix): The absorption probabilities.
#             - transient_nodes (np.ndarray): The indices of the transient nodes.
#             - absorbing_nodes (np.ndarray): The indices of the absorbing nodes.
#     """
#     absorbing_nodes = np.array(sorted(absorbing)) # absorbing nodes
#     # transient_nodes = get_reachable_transients(P, absorbing)
#     transient_nodes = np.setdiff1d(np.arange(P.shape[0]), list(absorbing))  # pyright: ignore[reportOptionalSubscript]
    
#     Q = P[transient_nodes][:, transient_nodes] # type: ignore
#     R = P[transient_nodes][:, absorbing_nodes] # type: ignore
#     I = eye(Q.shape[0], format='csr') # type: ignore
#     fundamental_matrix = inv(I - Q)
    
#     absorption_probabilities = fundamental_matrix @ R
#     return absorption_probabilities, transient_nodes, absorbing_nodes

# def main():
#     setup_logger(matrix)
#     logging.info("Starting simulations with the following parameters:")
#     logging.info(f"Matrix file: {matrix}")
#     logging.info(f"Output directory: {output_dir}")

#     P = load_npz(matrix)
#     _, absorbing_indices = get_target_nodes()
#     if not check_transition_matrix(P, absorbing=absorbing_indices):
#         print("❌ Matrix is not a valid transition matrix.")
#         return
    
#     absorption_probabilities, transient_nodes, absorbing_nodes = compute_absorption_probabilities(P, absorbing_indices)

#     df = pd.DataFrame(
#         absorption_probabilities.toarray(),
#         index=transient_nodes,
#         columns=absorbing_nodes
#     )

#     csv_out = Path(output_dir) / f"{Path(matrix).stem}_probs.csv"
#     npz_out = Path(output_dir) / f"{Path(matrix).stem}_probs.npz"
    
#     df.to_csv(csv_out)
#     save_npz(npz_out, absorption_probabilities)
    
#     logging.info(f"Saved absorption probabilities to {csv_out} and {npz_out}")
#     print(f"Absorption probabilities saved to {csv_out} and {npz_out}")

# if __name__ == "__main__":
#     main()