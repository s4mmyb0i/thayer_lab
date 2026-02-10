"""
Samvit Prem Singhal | spsinghal@wesleyan.edu
Simulated random walks starting at every node given a SciPy Sparse Matrix.

The script takes 5 command-line inputs:
    <matrix> <pdb_file> <max_steps> <n_walks> <random_walk_dir>

Input: Normalized Weight Transition Matrix in .npz format
Output: .csv file with simulated random walks

python3 -m simulate_walk.simulate_walk \
  transition_matrix/P_transition.npz \
  100 \
  100 \
  random_walk
"""

import os
import sys
from pathlib import Path
import logging
from tkinter import END
import numpy as np
import scipy.sparse as sp
from typing import Literal
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from utils.target import get_target_nodes

START_RES = {38, 42, 66, 102, 105} # ligand binding pocket
END_RES = {3, 4, 6, 9, 10, 11, 13, 14, 17, 19, 15, 18, 26, 78} # allosteric binding pocket

if len(sys.argv) != 5:
    logging.error("Error! Incorrect number of arguments.")
    print("Error!\nUsage: python3 simulate_walk.py <matrix> <max_steps> <n_walks> <random_walk_dir>")
    sys.exit(1)

matrix = str(sys.argv[1])
max_steps = int(sys.argv[2])
n_walks = int(sys.argv[3])
random_walk_dir = str(sys.argv[4])
SCRIPT_DIR = Path(__file__).resolve().parent
Path(random_walk_dir).mkdir(parents=True, exist_ok=True)

def setup_logger(matrix_file: str) -> None:
    """
    Set up a logger for the given matrix file.

    Args:
        matrix_file (str): Path to the matrix file.

    Returns:
        None
    """
    log_name = Path(matrix_file).stem + ".log"
    log_path = SCRIPT_DIR / log_name
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

termination_codes = {
    1: "hit target",
    2: "dead_end",
    3: "circular_walk",
    4: "max_steps"
}

def load_matrix(matrix_file: str) -> sp.csr_matrix:
    """Load the matrix data from a given file.
    
    Args:
        matrix_file (str): Path to the matrix file.
        
    Returns:
        sp.csr_matrix: The loaded matrix data.
    """
    return sp.load_npz(matrix_file)

def random_walk(P: sp.csr_matrix, start: int, max_steps: int, target_nodes: set[int], rng) -> tuple[list[int], Literal[1, 2, 3, 4], list[float]]:
    # print("# of target nodes:", len(target_nodes))
    # print("Min / Max target idx:", min(target_nodes), max(target_nodes))
    # assert max(target_nodes) < P.shape[0], "target_nodes has indices outside CSR matrix!"
    
    """
    Perform a random walk on a given transition matrix starting from a specified node.

    Args:
        P (sp.csr_matrix): A sparse matrix representing the transition probabilities between nodes.
        start (int): The starting node for the random walk.
        max_steps (int): The maximum number of steps to take in the walk.
        target_nodes (set[int]): A set of nodes that are considered targets for the walk.
        rng: A random number generator to use for sampling the next node.

    Returns:
        tuple: A tuple containing:
            - list[int]: The path taken during the walk as a list of node indices.
            - Literal[1, 2, 3, 4]: The termination code indicating the reason for walk termination.
                1: Reached target node
                2: Dead end (no outgoing edges)
                3: Circular walk (returned to a previously visited node)
                4: Max steps reached without hitting a target
            - list[float]: The list of log-transformed probabilities of each step in the walk.
    """

    path        = [start]
    current     = start
    probs_seen  = []
    term_code   = 0
    
    for _ in range(max_steps):
        start_ptr, end_ptr = P.indptr[current], P.indptr[current+1]
        neighs = P.indices[start_ptr:end_ptr]
        weights = P.data[start_ptr:end_ptr]
        
        next_node = rng.choice(neighs, p=weights)
        if len(neighs) == 0:
            term_code = 2; break
        
        if next_node in path:
            term_code = 3; break
            
        prob = weights[neighs.tolist().index(next_node)]
        probs_seen.append(-np.log(prob))
        path.append(next_node)
        current = next_node
        
        if current in target_nodes:
            term_code = 1; break
    
    else:
        term_code = 4
    
    return path, term_code, probs_seen

def worker(P_csr: sp.csr_matrix, atoms, max_steps: int, n_walks: int, end_nodes: set[int], atom2res: dict[int, int], rng_seed):
    """
    Worker function to perform multiple random walks starting from different nodes.

    Args:
        P_csr (sp.csr_matrix): The transition matrix.
        atoms (list[int]): A list of node indices to start the walks from.
        max_steps (int): The maximum number of steps to take in each walk.
        n_walks (int): The number of walks to perform starting from each node.
        target_nodes (set[int]): A set of target nodes to terminate the walks at.
        atom2res (dict[int, int]): A dictionary mapping node indices to residue numbers.
        rng_seed (int): The seed to use for the random number generator.

    Returns:
        list[tuple]: A list of tuples containing information about each walk.
            - int: The starting node index.
            - int: The residue number of the starting node.
            - list[int]: The path taken during the walk as a list of node indices.
            - Literal[1, 2, 3, 4]: The termination code indicating the reason for walk termination.
            - list[float]: The list of log-transformed probabilities of each step in the walk.
            - float: The sum of the log-transformed probabilities of each step in the walk.
    """
    rng = np.random.default_rng(rng_seed)
    out = []
    for atom in atoms:
        residue = atom2res.get(atom, -1)
        for _ in range(n_walks):
            path, code, probs = random_walk(P_csr, atom, max_steps, end_nodes, rng)
            out.append((
                atom,
                residue,
                path,
                code,
                len(path),
                probs,
                float(np.sum(probs))
            ))
        logging.info(f"Finished walks starting from atom {atom}.")
    return out

def export_csv(results, matrix_file, out_dir):
    """
    Export the results of random walks to a CSV file.

    Args:
        results (list[tuple]): A list of tuples containing information about each walk.
            Each tuple includes:
            - int: The starting node index (atom).
            - int: The residue number of the starting node.
            - list[int]: The path taken during the walk as a list of node indices.
            - Literal[1, 2, 3, 4]: The termination code indicating the reason for walk termination.
            - int: The length of the path.
            - list[float]: The list of negative log-transformed probabilities of each step in the walk.
            - float: The sum of the negative log-transformed probabilities of each step in the walk.
        matrix_file (str): The path to the matrix file, used to derive the output filename.
        out_dir (str): The directory where the CSV file will be saved.

    Returns:
        None
    """

    output_path = Path(out_dir) / (Path(matrix_file).stem + "_walks.csv")
    with open(output_path, 'w') as f:
        f.write("start;residue;path;termination_code;length;neglog_probs;sum_neglog\n")
        for start, res, path, code, length, probs, total in results:
            if length == 1: probs = [0.0]
            f.write(
                f"{start};{res};{','.join(map(str, path))};{code};{length};"
                f"{','.join(map(str, probs))};{total}\n"
            )
    logging.info(f"Written random walks to {output_path}")

def main():
    setup_logger(matrix)

    logging.info("Starting simulations with the following parameters:")
    logging.info(f"Matrix file: {matrix}")
    logging.info(f"Max steps: {max_steps}")
    logging.info(f"Number of walks: {n_walks}")
    logging.info(f"Random walk directory: {random_walk_dir}")

    P_csr = load_matrix(matrix)
    atom2res, start_nodes = get_target_nodes(START_RES)
    atom2res, end_nodes = get_target_nodes(END_RES)
    
    # res2atoms = get_res_atom_map(pdb_file)
    
    # atom2res = {atom: res for res, atoms in res2atoms.items() for atom in atoms}
    # assert 0 in atom2res, "Atom 0 not found in atom2res – index mismatch?"
    
    # target_nodes = {
    #     a for res in TARGET_RES
    #     for a in res2atoms.get(res, [])
    # }
    # all_atoms = np.arange(P_csr.shape[0])  # pyright: ignore[reportOptionalSubscript]
    all_atoms = np.array(sorted(start_nodes))
    n_workers = os.cpu_count() or 4
    chunks = np.array_split(all_atoms, n_workers)
    
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = [
            pool.submit(worker, P_csr, chunk, max_steps, n_walks, end_nodes, atom2res, 1234+i)
            for i, chunk in enumerate(chunks)
        ]
        results = [r for f in as_completed(futures) for r in f.result()]

    export_csv(results, matrix, random_walk_dir)
    logging.info("Simulation completed.")

if __name__ == "__main__":
    main()