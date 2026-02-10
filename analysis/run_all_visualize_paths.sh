#!/usr/bin/env bash
# Run visualize_paths.py for all 8 walk CSVs and save HTML in analysis/figures/.
# Run from project root: ./analysis/run_all_visualize_paths.sh
# Optional: set MAX_PATHS to limit paths per figure (default 1000 for smaller files).

set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
FIGURES_DIR="$SCRIPT_DIR/figures"
mkdir -p "$FIGURES_DIR"
cd "$PROJECT_ROOT"

# Limit paths per visualization (set empty for all paths; 1000 keeps HTML manageable)
MAX_PATHS="${MAX_PATHS:-1000}"

PDB_DIR="pdb"
WALK_DIR="random_walk"
SCRIPT="analysis/visualize_paths.py"

for csv in \
  "$WALK_DIR/P_allos->lig_walks.csv" \
  "$WALK_DIR/PL_allos->lig_walks.csv" \
  "$WALK_DIR/AP_allos->lig_walks.csv" \
  "$WALK_DIR/APL_allos->lig_walks.csv" \
  "$WALK_DIR/P_lig->allos_walks.csv" \
  "$WALK_DIR/PL_lig->allos_walks.csv" \
  "$WALK_DIR/AP_lig->allos_walks.csv" \
  "$WALK_DIR/APL_lig->allos_walks.csv"; do
  base=$(basename "$csv" _walks.csv)
  pdb_name="${base%%_*}"
  pdb_file="$PDB_DIR/${pdb_name}.pdb"
  out="$FIGURES_DIR/${base}_paths.html"
  echo "Running: $pdb_file + $csv -> $out"
  if [ -n "$MAX_PATHS" ]; then
    python3 "$SCRIPT" "$pdb_file" "$csv" "$MAX_PATHS" --output "$out"
  else
    python3 "$SCRIPT" "$pdb_file" "$csv" --output "$out"
  fi
done

echo "Done. Figures saved in $FIGURES_DIR"
