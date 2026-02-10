#!/usr/bin/env bash
# Run node_heatmap.py for all 8 walk CSVs and save HTML in analysis/figures/.
# Run from project root: ./analysis/run_all_heatmaps.sh

set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
FIGURES_DIR="$SCRIPT_DIR/figures"
mkdir -p "$FIGURES_DIR"
cd "$PROJECT_ROOT"

WALK_DIR="random_walk"
HEATMAP_SCRIPT="analysis/node_heatmap.py"

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
  out="$FIGURES_DIR/${base}_heatmap.html"
  echo "Running: $csv -> $out"
  python3 "$HEATMAP_SCRIPT" "$csv" --output "$out"
done

echo "Done. Figures saved in $FIGURES_DIR"
