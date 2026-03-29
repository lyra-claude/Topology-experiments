#!/usr/bin/env bash
# Pilot study: OneMax domain
# 3 runs × 8 topologies × 500 generations
# 8 islands, pop 50 per island, migration every 10 gens, 5 migrants

set -euo pipefail

OUTDIR="results/pilot_onemax"
mkdir -p "$OUTDIR"

TOPOLOGIES="disconnected ring star complete hypercube barbell watts-strogatz random-regular"
SEEDS="42 137 2718"

ISLANDS=8
POP=50
MIG_INTERVAL=10
MIGRANTS=5
GENS=500

total=24  # 8 topologies × 3 seeds
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
