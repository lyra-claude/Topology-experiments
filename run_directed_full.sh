#!/usr/bin/env bash
# Full directed cycle experiment
# 30 runs x 8 directed topologies x 500 generations
# 8 islands, pop 50 per island, migration every 10 gens, 5 migrants
# Same parameters as the undirected OneMax pilot for comparability.

set -euo pipefail

OUTDIR="results/directed_full"
mkdir -p "$OUTDIR"

# 8 directed topologies ordered by simple cycle count:
# dag-layer (0), dag-wide (0), lowcyc-1 (3), bidir-ring (10),
# two-cliques (14), mesh-cyclic (20), dense-triangles (29), ring-skip2 (47)
TOPOLOGIES="dag-layer dag-wide lowcyc-1 bidir-ring two-cliques mesh-cyclic dense-triangles ring-skip2"

# 30 seeds for statistical power
SEEDS="42 137 2718 314 1618 7 99 256 512 1024 2048 4096 8192 16384 32768 11 23 37 53 71 89 97 113 131 151 173 191 199 211 223"

ISLANDS=8
POP=50
MIG_INTERVAL=10
MIGRANTS=5
GENS=500

total=240  # 8 topologies x 30 seeds
count=0

for topo in $TOPOLOGIES; do
  for seed in $SEEDS; do
    count=$((count + 1))
    outfile="${OUTDIR}/${topo}_seed${seed}.csv"
    if [ -f "$outfile" ]; then
      echo "[$count/$total] SKIP (exists): $topo seed=$seed"
      continue
    fi
    echo "[$count/$total] Running: $topo seed=$seed -> $outfile"
    cabal run topology-sim -- --domain onemax "$topo" $ISLANDS $POP $MIG_INTERVAL $MIGRANTS $GENS "$seed" \
      > "$outfile" 2>/dev/null
  done
done

echo ""
echo "All $total runs complete. Results in $OUTDIR/"
