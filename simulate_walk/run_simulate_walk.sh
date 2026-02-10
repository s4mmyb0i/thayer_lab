#!/bin/bash

# A script to run create_network.py for all matrix/PDB file pairs
# Usage:
# $ simulate_walk/./run_simulate_walk.sh /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/transition_matrix /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/pdb 100 100 /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/random_walk

# $ simulate_walk/./run_simulate_walk.sh /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/matrices/trans_abs_matrix 100 100 /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/random_walk

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <matrix_dir> <max_steps> <n_walks> <random_walk_dir>"
    exit 1
fi

MATRIX_DIR="$1"
MAX_STEPS="$2"
N_WALKS="$3"
RANDOM_WALK_DIR="$4"

LOG_FILE="./simulate_walk/simulation.log"
exec >"$LOG_FILE" 2>&1  # Redirect all stdout and stderr to the log file

echo "==== Batch Network Creation Started ===="
echo "Matrix directory: $MATRIX_DIR"
echo "Max steps: $MAX_STEPS"
echo "Number of walks: $N_WALKS"
echo "Random walk output directory: $RANDOM_WALK_DIR"
echo "Logging to: $LOG_FILE"
echo ""

mkdir -p "$RANDOM_WALK_DIR"

if [ ! -d "$MATRIX_DIR" ]; then
    echo "Matrix directory $MATRIX_DIR does not exist."
    exit 1
fi

shopt -s nullglob
matrix_files=("$MATRIX_DIR"/*.npz)
if [ ${#matrix_files[@]} -eq 0 ]; then
    echo "No .dat files found in $MATRIX_DIR"
    exit 1
fi

for matrix_file in "${matrix_files[@]}"; do
    echo "Processing $matrix_file"
    base_name=$(basename "$matrix_file" .npz)

    python3 -m simulate_walk.simulate_walk "$matrix_file" "$MAX_STEPS" "$N_WALKS" "$RANDOM_WALK_DIR"
    echo ""
done

echo "==== Batch Random Simulated Walks Completed ===="