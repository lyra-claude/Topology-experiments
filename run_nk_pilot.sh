#!/usr/bin/env bash
# NK landscape pilot experiment
# 3 domains (nk0, nk4, nk6) x 3 directed topologies x 10 seeds = 90 runs
# Tests whether ruggedness (K) modulates the cycle-count effect.
#
# Topologies (constant density: n=8, m=16):
#   dag-layer     (0 directed simple cycles)
#   bidir-ring    (10 directed simple cycles)
#   ring-skip2    (47 directed simple cycles)

set -euo pipefail

OUTDIR="results/nk_pilot"
mkdir -p "$OUTDIR"

DOMAINS="nk0 nk4 nk6"
TOPOLOGIES="dag-layer bidir-ring ring-skip2"
SEEDS="42 137 2718 314 1618 7 99 256 512 1024"

ISLANDS=8
POP=50
MIG_INTERVAL=10
MIGRANTS=5
GENS=500

total=90  # 3 domains x 3 topologies x 10 seeds
count=0

for domain in $DOMAINS; do
  for topo in $TOPOLOGIES; do
    for seed in $SEEDS; do
      count=$((count + 1))
      outfile="${OUTDIR}/${domain}_${topo}_${seed}.csv"
      if [ -f "$outfile" ]; then
        echo "[$count/$total] SKIP (exists): $domain $topo seed=$seed"
        continue
      fi
      echo "[$count/$total] Running: $domain $topo seed=$seed -> $outfile"
      cabal run topology-sim -- --domain "$domain" "$topo" $ISLANDS $POP $MIG_INTERVAL $MIGRANTS $GENS "$seed" \
        > "$outfile" 2>/dev/null
    done
  done
done

echo ""
echo "All $total runs complete. Results in $OUTDIR/"
