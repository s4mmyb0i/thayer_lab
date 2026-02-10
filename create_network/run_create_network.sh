#!/bin/bash

# A script to run create_network.py for all matrix/PDB file pairs
# Usage:
# $ create_network/./run_create_network.sh /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/dccm /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/pdb 0.8 /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/matrices/graphml /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/matrices/adj_matrix /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/matrices/trans_matrix /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/matrices/trans_abs_matrix

if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <matrix_dir> <pdb_dir> <threshold_val_lower> <graphml_dir> <adj_matrix_dir> <trans_matrix_dir> <trans_absorbing_matrix_dir>"
    exit 1
fi

MATRIX_DIR="$1"
PDB_DIR="$2"
THRESHOLD_VAL_LOWER="$3"
GRAPHML_DIR="$4"
ADJ_MATRIX_DIR="$5"
TRANS_MATRIX_DIR="$6"
ABS_TRANS_MATRIX_DIR="$7"

LOG_FILE="./create_network/network_creation.log"
exec >"$LOG_FILE" 2>&1  # Redirect stdout and stderr to log

echo "==== Batch Network Creation Started ===="
echo "Matrix directory: $MATRIX_DIR"
echo "PDB directory: $PDB_DIR"
echo "Threshold lower bound: $THRESHOLD_VAL_LOWER"
echo "GraphML output directory: $GRAPHML_DIR"
echo "Adjacency matrix output directory: $ADJ_MATRIX_DIR"
echo "Transition matrix output directory: $TRANS_MATRIX_DIR"
echo "Absorbing transition matrix output directory: $ABS_TRANS_MATRIX_DIR"
echo "Logging to: $LOG_FILE"
echo ""

mkdir -p "$GRAPHML_DIR" "$ADJ_MATRIX_DIR" "$TRANS_MATRIX_DIR" "$ABS_TRANS_MATRIX_DIR"

if [ ! -d "$MATRIX_DIR" ]; then
    echo "Matrix directory $MATRIX_DIR does not exist."
    exit 1
fi

if [ ! -d "$PDB_DIR" ]; then
    echo "PDB directory $PDB_DIR does not exist."
    exit 1
fi

shopt -s nullglob
matrix_files=("$MATRIX_DIR"/*.dat)
if [ ${#matrix_files[@]} -eq 0 ]; then
    echo "No .dat files found in $MATRIX_DIR"
    exit 1
fi

for matrix_file in "${matrix_files[@]}"; do
    echo "Processing $matrix_file"
    base_name=$(basename "$matrix_file" .dat)
    pdb_file="$PDB_DIR/$base_name.pdb"

    if [ ! -f "$pdb_file" ]; then
        echo "⚠️ PDB file $pdb_file does not exist. Skipping."
        continue
    fi

    echo "Using PDB file: $pdb_file"
    python3 -m create_network.create_network "$matrix_file" "$pdb_file" "$THRESHOLD_VAL_LOWER" "$GRAPHML_DIR" "$ADJ_MATRIX_DIR" "$TRANS_MATRIX_DIR" "$ABS_TRANS_MATRIX_DIR"
    echo ""
done

echo "==== Batch Network Creation Completed ===="