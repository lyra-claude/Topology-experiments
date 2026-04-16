#!/usr/bin/env bash
# Interference experiment on Maze domain
# 3 topologies x 1 domain (maze) x 30 seeds = 90 runs

set -euo pipefail

OUTDIR="results/interference-maze"

TOPOLOGIES="interference-adjacent interference-separated interference-nested"

SEEDS="42 137 2718 314 1618 7 99 256 512 1024 2048 4096 8192 16384 32768 11 23 37 53 71 89 97 113 131 151 173 191 199 211 223"

ISLANDS=8
POP=50
MIG_INTERVAL=10
MIGRANTS=5
GENS=500

total=90
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
