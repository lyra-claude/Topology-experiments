#!/usr/bin/env bash
# Transient study: OneMax domain, 10 runs per topology
# 10 runs × 8 topologies × 100 generations
# 8 islands, pop 50 per island, migration every 10 gens, 5 migrants

set -euo pipefail

BIN=$(cabal list-bin topology-sim 2>/dev/null)
OUTDIR="results/transient_onemax_10runs"
mkdir -p "$OUTDIR"

TOPOLOGIES="disconnected ring star complete hypercube barbell watts-strogatz random-regular"
SEEDS="42 137 2718 314 1618 271 577 997 1729 4242"

ISLANDS=8
POP=50
MIG_INTERVAL=10
MIGRANTS=5
GENS=100

total=80  # 8 topologies × 10 seeds
count=0
running=0
MAX_PARALLEL=8

for topo in $TOPOLOGIES; do
  for seed in $SEEDS; do
    count=$((count + 1))
    outfile="${OUTDIR}/${topo}_seed${seed}.csv"
    echo "[$count/$total] Running: $topo seed=$seed -> $outfile"
    "$BIN" --domain onemax "$topo" $ISLANDS $POP $MIG_INTERVAL $MIGRANTS $GENS "$seed" \
      > "$outfile" 2>/dev/null &
    running=$((running + 1))
    if [ $running -ge $MAX_PARALLEL ]; then
      wait
      running=0
    fi
  done
done

wait
echo ""
echo "All $total runs complete. Results in $OUTDIR/"
