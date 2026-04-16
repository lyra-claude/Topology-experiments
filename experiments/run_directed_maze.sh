#!/usr/bin/env bash
# Directed cycle experiment on Maze domain (8x8 grid for speed)
# 8 directed topologies x 1 domain (maze) x 30 seeds = 240 runs
# Tests whether directed cycle count (kappa) predicts diversity on maze navigation.

set -euo pipefail

OUTDIR="results/directed-maze"

# 8 directed topologies ordered by simple cycle count (kappa):
# dag-layer (0), dag-wide (0), lowcyc-1 (3), bidir-ring (10),
# two-cliques (14), mesh-cyclic (20), dense-triangles (29), ring-skip2 (47)
TOPOLOGIES="dag-layer dag-wide lowcyc-1 bidir-ring two-cliques mesh-cyclic dense-triangles ring-skip2"

# Same 30 seeds as original directed cycle experiment
SEEDS="42 137 2718 314 1618 7 99 256 512 1024 2048 4096 8192 16384 32768 11 23 37 53 71 89 97 113 131 151 173 191 199 211 223"

ISLANDS=8
POP=50
MIG_INTERVAL=10
MIGRANTS=5
GENS=500

total=240  # 8 topologies x 30 seeds
count=0

mkdir -p "${OUTDIR}"

for topo in $TOPOLOGIES; do
  for seed in $SEEDS; do
    count=$((count + 1))
    outfile="${OUTDIR}/${topo}_seed${seed}.csv"
    if [ -f "$outfile" ]; then
      echo "[$count/$total] SKIP (exists): maze $topo seed=$seed"
      continue
    fi
    echo "[$count/$total] Running: maze $topo seed=$seed -> $outfile"
    cabal run topology-sim -- --domain maze --grid 8 "$topo" $ISLANDS $POP $MIG_INTERVAL $MIGRANTS $GENS "$seed" \
      > "$outfile" 2>/dev/null
  done
done

echo ""
echo "All $total runs complete. Results in $OUTDIR/"
