#!/usr/bin/env bash
# Multi-domain directed cycle experiment
# 8 directed topologies x 4 NK domains x 30 seeds = 960 runs
# Tests whether directed cycle count (kappa) predicts diversity
# across landscapes of varying ruggedness.

set -euo pipefail

OUTDIR="results/directed_multidomain"

# 8 directed topologies ordered by simple cycle count (kappa):
# dag-layer (0), dag-wide (0), lowcyc-1 (3), bidir-ring (10),
# two-cliques (14), mesh-cyclic (20), dense-triangles (29), ring-skip2 (47)
TOPOLOGIES="dag-layer dag-wide lowcyc-1 bidir-ring two-cliques mesh-cyclic dense-triangles ring-skip2"

# 4 NK domains: smooth to very rugged
DOMAINS="nk0 nk2 nk4 nk6"

# Same 30 seeds as original directed cycle experiment
SEEDS="42 137 2718 314 1618 7 99 256 512 1024 2048 4096 8192 16384 32768 11 23 37 53 71 89 97 113 131 151 173 191 199 211 223"

ISLANDS=8
POP=50
MIG_INTERVAL=10
MIGRANTS=5
GENS=500

total=960  # 8 topologies x 4 domains x 30 seeds
count=0

for domain in $DOMAINS; do
  mkdir -p "${OUTDIR}/${domain}"
  for topo in $TOPOLOGIES; do
    for seed in $SEEDS; do
      count=$((count + 1))
      outfile="${OUTDIR}/${domain}/${topo}_seed${seed}.csv"
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
