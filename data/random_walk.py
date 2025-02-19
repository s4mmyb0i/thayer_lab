import os
import sys
import logging
import numpy as np
import scipy.sparse as sp
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

from res_atom_dict import residue_atoms_dict

def setup_logger(construct:str, walk_dir:str):
    log_file = os.path.join(walk_dir, f"{construct}_walks.log")  # Correctly join path
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

# Setup logger
construct = sys.argv[1].split('/')[-2] if len(sys.argv) > 1 else "unknown_construct"
walk_dir = sys.argv[3] if len(sys.argv) > 3 else "unknown_walk_dir"
setup_logger(construct, walk_dir)

# Ensure correct number of arguments
if len(sys.argv) != 7:
    logging.error("Error! Incorrect number of arguments.")
    print("Error!\nUsage: python3 random_walk.py <pdb_file> <network_dir> <walk_dir> <steps> <number_of_walks> <threshold_val_lower>")
    sys.exit(1)

# Log the start of the script
logging.info("Starting random walk script with the following parameters:")
logging.info(f"Walk directory: {sys.argv[3]}")
logging.info(f"Number of walks: {sys.argv[5]}")
logging.info(f"Steps: {sys.argv[4]}")


PDB_FILE = str(sys.argv[1])
NETWORK_DIR = str(sys.argv[2])
WALK_DIR = str(sys.argv[3])
STEPS = int(sys.argv[4])
NUM_WALKS = int(sys.argv[5])
THRESHOLD_VAL_UPPER = 0.999
THRESHOLD_VAL_LOWER = float(sys.argv[6])
BINDING_POCKET_RESIDUES = [326, 327, 328, 329, 330, 331, 343, 376, 383]

RESIDUE_ATOMS = residue_atoms_dict(PDB_FILE)  # returns a dict (residue : atoms) of type (int : int list) for every residue

BINDING_POCKET_ATOMS = [
    atom for residue in BINDING_POCKET_RESIDUES if residue in RESIDUE_ATOMS
    for atom in RESIDUE_ATOMS[residue]
]

TERMINATION_CODES = {
    1: "hit_target",
    2: "dead_end",
    3: "circular_path",
    4: "max_steps"
}

def load_network(network_file:str) -> sp.csr_matrix:
    matrix = sp.load_npz(network_file)
    return matrix

def get_probs(matrix:sp.csr_matrix, current_node:int, neighborhood:list[int]) -> list[float]:
    """
    Returns all possible probabilities in log space from current_node to neighborhood nodes
    """
    # Extract weights from the matrix
    weights = matrix[current_node, neighborhood].toarray().flatten()
    
    # Apply threshold scaling and use logarithm directly
    log_weights = np.log(weights - THRESHOLD_VAL_LOWER + 0.001)  # Ensure positive values before taking the log
    
    # Normalize the log probabilities (equivalent to dividing in normal space)
    log_prob_sum = np.logaddexp.reduce(log_weights)  # sum in log space
    log_probabilities = log_weights - log_prob_sum  # subtract log(sum) for normalization
    
    # print(f"got log-probabilities {log_probabilities}")
    return log_probabilities


def random_walk(matrix:sp.csr_matrix, start_node:int) -> tuple[int, list[int], int, list[float], float]:
    """
    Termination Codes:
        1: hit_target
        2: dead_end
        3: circular_path
        4: max_steps
    Returns:
    start_node, path, termination_code, probabilities, sum_probability
    """
    path = [start_node]  # int list initialized with starting node
    current_node = start_node
    probabilities = []  # float list of each probability
    termination_code = 0  # initialized to None value
    
    for i in range(STEPS):
        if current_node in BINDING_POCKET_ATOMS:
            termination_code = 1  # Hit target
            break
            
        # Ensure current_node is within the matrix bounds
        if current_node >= matrix.shape[0]:
            termination_code = 2  # Invalid node
            break

        # Handle boundary case when current_node is the last node in the matrix
        if current_node < matrix.shape[0]-1:
            neighborhood = matrix.indices[matrix.indptr[current_node]:matrix.indptr[current_node+1]]
        else:
            neighborhood = matrix.indices[matrix.indptr[current_node]:]

        if len(neighborhood) == 0:
            termination_code = 2
            break
        
        possible_prob = get_probs(matrix, current_node, neighborhood)
        possible_prob = np.exp(possible_prob)  # convert back to normal space for sampling
        next_node = np.random.choice(neighborhood, p=possible_prob/possible_prob.sum())  # chose next node based on probability
        if next_node in path:
            termination_code = 3
            break
        
        prob = possible_prob[neighborhood.tolist().index(next_node)]
        
        # Checking for edge case where probability is zero
        if prob > 0: probabilities.append(float(-np.log(prob)))
        else: probabilities.append(0.0)  # Or handle zero probability differently
        
        probabilities.append(float(-np.log(prob) if prob != 0 else 0.0))
        path.append(int(next_node))
        current_node = next_node
        
        if i == STEPS-1:
            termination_code = 4
            break
        
        # logging.debug(f"Step {i}: Current node: {current_node}, Path so far: {path}")
        # logging.debug(f"Neighborhood: {neighborhood}, Probabilities: {possible_prob}")
        
    sum_probability = np.sum(probabilities)
    # logging.info(f"Random walk from {start_node} to {path[-1]} completed.")
    return start_node, path, termination_code, probabilities, sum_probability

def process_walks_for_atom(matrix:sp.csr_matrix, atom:int, residue:int):
    path_counter = defaultdict(int)
    path_probabilities = {}  # New dictionary to store probabilities per unique path
    data_to_write = []

    for walk_counter in range(1, NUM_WALKS+1):
        start_node, path, termination_code, probabilities, sum_probability = random_walk(matrix, atom)
        path_tuple = tuple(path)
        path_counter[path_tuple] += 1

        # Store the probabilities for this specific path
        path_probabilities[path_tuple] = (probabilities, sum_probability)

    for path_tuple, count in path_counter.items():
        path_str = ",".join(map(str, path_tuple))
        probabilities, sum_probability = path_probabilities[path_tuple]  # Retrieve correct probabilities
        probabilities_str = ",".join(map(str, probabilities))
        path_length = len(path_tuple)
        data_to_write.append((atom, residue, path_str, termination_code, path_length, count, probabilities_str, sum_probability))

    # logging.info(f"Completed {NUM_WALKS} walks for atom {atom} (Residue {residue}).")
    return data_to_write, path_counter

def simulate_walks():
    # print("Hello")
    # logging.info(f"Starting random walks for {construct} in directory {WALK_DIR}...")
    
    os.makedirs(WALK_DIR, exist_ok=True)

    for network_file in os.listdir(NETWORK_DIR):
        permutation = os.path.splitext(network_file)[0]
        
        construct = permutation.split('_')[0]
        
        walk_residue = int(permutation)+305  # Residue number from which the walk will start. 306 is the first residue number in the P/PL/AP/APL constructs
        if walk_residue not in RESIDUE_ATOMS:
            logging.info(f"Residue {walk_residue} for starting walk not in pdb file. Skipping.")
            continue
        
        output_file = os.path.join(WALK_DIR, f"{permutation}.csv")
        network_file_path = os.path.join(NETWORK_DIR, network_file)  # Full path to the network file
        # logging.info(f"Output file is {output_file} in {WALK_DIR}")
        # print(network_file)
        # Load the network once
        matrix = load_network(network_file_path)

        with open(output_file, 'w') as f:
            f.write("atom;residue;path;termination_code;length;counter;probabilities;sum_probability\n")

            # Prepare all (network, atom, residue) tuples before parallel execution
            tasks = [(matrix, atom, walk_residue) for atom in RESIDUE_ATOMS[walk_residue]] # residue, atoms in RESIDUE_ATOMS.items() for atom in atoms]

            # Use ProcessPoolExecutor with optimized `map()`
            with ProcessPoolExecutor() as executor:
                for data_to_write, _ in executor.map(process_walks_for_atom, *zip(*tasks), chunksize=10):
                    # Use list comprehension for batch writing
                    f.writelines([
                        f"{atom};{walk_residue};{path_str};{termination_code};{length};{counter};{probabilities_str};{sum_probability}\n"
                        for atom, walk_residue, path_str, termination_code, length, counter, probabilities_str, sum_probability in data_to_write
                    ])

        logging.info(f"Finished processing network file {network_file_path}. Output written to {output_file}")

    logging.info("Simulation complete.")

if __name__ == "__main__":
    # print("check")
    simulate_walks()
