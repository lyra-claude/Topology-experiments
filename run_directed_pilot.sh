#!/usr/bin/env bash
# Pilot study: Directed cycle experiment
# 5 runs x 3 topologies (DAG-Layer, Bidir-Ring, Ring-Skip2) x 500 generations
# Tests 0, 10, and 47 cycle families to verify signal before full run.
# 8 islands, pop 50 per island, migration every 10 gens, 5 migrants

set -euo pipefail

OUTDIR="results/directed_pilot"
mkdir -p "$OUTDIR"

TOPOLOGIES="dag-layer bidir-ring ring-skip2"
SEEDS="42 137 2718 314 1618"

ISLANDS=8
POP=50
MIG_INTERVAL=10
MIGRANTS=5
GENS=500

total=15  # 3 topologies x 5 seeds
count=0

for topo in $TOPOLOGIES; do
  for seed in $SEEDS; do
    count=$((count + 1))
    outfile="${OUTDIR}/${topo}_seed${seed}.csv"
    echo "[$count/$total] Running: $topo seed=$seed -> $outfile"
    cabal run topology-sim -- --domain onemax "$topo" $ISLANDS $POP $MIG_INTERVAL $MIGRANTS $GENS "$seed" \
      > "$outfile" 2>/dev/null
  done
done

echo ""
echo "All $total runs complete. Results in $OUTDIR/"
