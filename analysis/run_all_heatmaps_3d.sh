#!/usr/bin/env bash
# Run node_heatmap_3d.py for all 8 walk CSVs and save HTML in analysis/figures/.
# Run from project root: ./analysis/run_all_heatmaps_3d.sh

set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
FIGURES_DIR="$SCRIPT_DIR/figures"
mkdir -p "$FIGURES_DIR"
cd "$PROJECT_ROOT"

PDB_DIR="pdb"
WALK_DIR="random_walk"
HEATMAP_SCRIPT="analysis/node_heatmap_3d.py"

# Each CSV: base is e.g. P_allos->lig or PL_lig->allos; PDB is first segment (P, PL, AP, APL)
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
  # PDB name is first component before underscore (P, PL, AP, APL)
  pdb_name="${base%%_*}"
  pdb_file="$PDB_DIR/${pdb_name}.pdb"
  out="$FIGURES_DIR/${base}_heatmap_3d.html"
  echo "Running: $pdb_file + $csv -> $out"
  python3 "$HEATMAP_SCRIPT" "$pdb_file" "$csv" --output "$out"
done

echo "Done. Figures saved in $FIGURES_DIR"
