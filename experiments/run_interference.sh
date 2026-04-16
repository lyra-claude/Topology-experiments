#!/usr/bin/env bash
# Interference experiment: does cycle ARRANGEMENT matter at constant beta_1?
#
# 3 directed topologies, all with n=8, |E|=9, beta_1=2, kappa=2:
#   1. interference-adjacent  — two cycles sharing a vertex (figure-eight)
#   2. interference-separated — two cycles connected by a path (no shared vertices)
#   3. interference-nested    — one cycle inside another (concentric)
#
# NK4 domain (where topology effects are strongest), 30 seeds, 500 generations.
# Motivated by Clio's confirmation that holonomy group is NON-ABELIAN (cactus J_n).
# If arrangement matters, eta-squared between graph types should be significant.

set -euo pipefail

export PATH="/home/lyra/.ghcup/bin:$PATH"

OUTDIR="results/interference"
DOMAIN="nk4"

TOPOLOGIES="interference-adjacent interference-separated interference-nested"

# Same 30 seeds as directed cycle experiments for comparability
SEEDS="42 137 2718 314 1618 7 99 256 512 1024 2048 4096 8192 16384 32768 11 23 37 53 71 89 97 113 131 151 173 191 199 211 223"

ISLANDS=8
POP=50
MIG_INTERVAL=10
MIGRANTS=5
GENS=500

total=90  # 3 topologies x 30 seeds
count=0

mkdir -p "$OUTDIR"

for topo in $TOPOLOGIES; do
  for seed in $SEEDS; do
    count=$((count + 1))
    outfile="${OUTDIR}/${topo}_seed${seed}.csv"
    if [ -f "$outfile" ]; then
      echo "[$count/$total] SKIP (exists): $topo seed=$seed"
      continue
    fi
    echo "[$count/$total] Running: $topo seed=$seed -> $outfile"
    cabal run topology-sim -- --domain "$DOMAIN" "$topo" $ISLANDS $POP $MIG_INTERVAL $MIGRANTS $GENS "$seed" \
      > "$outfile" 2>/dev/null
  done
done

echo ""
echo "All $total runs complete. Results in $OUTDIR/"
