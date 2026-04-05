#!/usr/bin/env bash
# Jaccard maze pilot: 3 seeds × 8 topologies × 500 generations
# Parameters match prior pilots: 8 islands, pop=50, migInterval=10, migrants=5
set -euo pipefail

export PATH="$HOME/.ghcup/bin:$PATH"

OUTDIR="results/jaccard_maze_pilot"
mkdir -p "$OUTDIR"

TOPOS=(disconnected ring star complete hypercube barbell watts-strogatz random-regular)
SEEDS=(42 137 2025)
ISLANDS=8
POP=50
MIG_INT=10
MIG_N=5
GENS=500

TOTAL=$((${#TOPOS[@]} * ${#SEEDS[@]}))
COUNT=0

echo "=== Jaccard Maze Pilot ==="
echo "Topologies: ${TOPOS[*]}"
echo "Seeds: ${SEEDS[*]}"
echo "Total runs: $TOTAL"
echo ""

for topo in "${TOPOS[@]}"; do
  for seed in "${SEEDS[@]}"; do
    COUNT=$((COUNT + 1))
    OUTFILE="$OUTDIR/${topo}_seed${seed}.csv"
    echo "[$COUNT/$TOTAL] $topo seed=$seed -> $OUTFILE"
    cabal run topology-sim -- maze "$topo" $ISLANDS $POP $MIG_INT $MIG_N $GENS $seed \
      > "$OUTFILE" 2>/dev/null
  done
done

echo ""
echo "=== All $TOTAL runs complete ==="
echo "Results in $OUTDIR/"
ls -la "$OUTDIR/"
