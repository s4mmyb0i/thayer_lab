#!/bin/bash

# A script to generate random walks for all four constructs

# Define the PDB files to use
PDB_FILES=(
    "/home/samvit/Documents/2024/project/from_leticia/P_1BFE/1bfe.pdb"
    "/home/samvit/Documents/2024/project/from_leticia/PL_1BE9/1be9.pdb"
    "/home/samvit/Documents/2024/project/from_leticia/AP_PDZ3CDC42/pdz3cdc42.pdb"
    "/home/samvit/Documents/2024/project/from_leticia/APL_PDZ-3_2/pdz3lc_3.pdb"
)

# Define the network directories (strictly these 4 directories)
NETWORK_DIRS=(
    "/home/samvit/Documents/2024/project/all_atom/project_0.8/scipy_matrix/P"
    "/home/samvit/Documents/2024/project/all_atom/project_0.8/scipy_matrix/PL"
    "/home/samvit/Documents/2024/project/all_atom/project_0.8/scipy_matrix/AP"
    "/home/samvit/Documents/2024/project/all_atom/project_0.8/scipy_matrix/APL"
)

# Define the output directories for each construct
OUTPUT_BASE_DIR="/home/samvit/Documents/2024/project/all_atom/project_0.8"
PROJECT_DIRS=(
    # "$OUTPUT_BASE_DIR/project_100"
    "$OUTPUT_BASE_DIR/project_1000"
    "$OUTPUT_BASE_DIR/project_10000"
    "$OUTPUT_BASE_DIR/project_100000"
    "$OUTPUT_BASE_DIR/project_1000000"
)

# Define other variables
THRESHOLD_VAL_LOWER=0.8

# Loop through each PDB file and network directory
for idx in ${!PDB_FILES[@]}; do
    PDB_FILE="${PDB_FILES[$idx]}"
    NETWORK_DIR="${NETWORK_DIRS[$idx]}"

    # Extract construct name from the PDB file path (P, PL, AP, APL)
    construct=$(basename $(dirname $PDB_FILE))

    # Loop through each project directory and corresponding construct subdirectory
    for project_dir in "${PROJECT_DIRS[@]}"; do
        WALK_DIR="${project_dir}/data/${construct}/random_walk"
        
        # Ensure that the walk directory exists
        mkdir -p "$WALK_DIR"

        # Extract the project number (e.g., 1000 from project_1000)
        project_number=$(basename "$project_dir" | sed 's/project_//')

        # Set the number of walks based on the project number
        NUM_WALKS=$project_number

        # Run the Python script for each PDB file and corresponding network directory
        python3 /home/samvit/Documents/2024/project/all_atom/project_0.8/scripts/data/random_walk.py \
            "$PDB_FILE" \
            "$NETWORK_DIR" \
            "$WALK_DIR" \
            100 "$NUM_WALKS" $THRESHOLD_VAL_LOWER
    done
done
