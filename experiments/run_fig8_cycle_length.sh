#!/usr/bin/env bash
# Figure-eight cycle-length experiment
# Tests whether mean cycle length predicts evolutionary dynamics at constant
# node count (n=12), edge count (|E|=13), and cycle rank (beta_1=2).
#
# Three directed figure-eight topologies from Claudius (2026-04-25):
#   fig8-m1: C3+C5 at v0, pendant chain (mean cycle length 4.0)
#   fig8-m2: C5+C6 at v0, short pendant  (mean cycle length 5.5)
#   fig8-m3: C7+C6 at v0, no pendant     (mean cycle length 6.5)
#
# Expected: shorter cycles -> faster mixing -> larger lambda_2
#           M1 > M2 > M3 in spectral gap
#
# Domains:
#   onemax — negative control (smooth, topology should be irrelevant)
#   nk4    — rugged landscape (topology effects strongest here)
#   maze   — deceptive/path-dependent (real-world proxy)
#
# 3 topologies x 3 domains x 30 seeds = 270 runs

set -euo pipefail

export PATH="/home/lyra/.ghcup/bin:$PATH"

OUTDIR="results/fig8_cycle_length"

TOPOLOGIES="fig8-m1 fig8-m2 fig8-m3"
DOMAINS="onemax nk4 maze"

# Same 30 seeds as all other experiments for comparability
SEEDS="42 137 2718 314 1618 7 99 256 512 1024 2048 4096 8192 16384 32768 11 23 37 53 71 89 97 113 131 151 173 191 199 211 223"

ISLANDS=12
POP=50
MIG_INTERVAL=10
MIGRANTS=5
GENS=500

total=270  # 3 topologies x 3 domains x 30 seeds
count=0

for domain in $DOMAINS; do
  mkdir -p "${OUTDIR}/${domain}"

  # Domain-specific flags
  DOMAIN_FLAGS="--domain $domain"
  if [ "$domain" = "maze" ]; then
    DOMAIN_FLAGS="--domain maze --grid 15"
  fi

  for topo in $TOPOLOGIES; do
    for seed in $SEEDS; do
      count=$((count + 1))
      outfile="${OUTDIR}/${domain}/${topo}_seed${seed}.csv"
      if [ -f "$outfile" ]; then
        echo "[$count/$total] SKIP (exists): $domain $topo seed=$seed"
        continue
      fi
      echo "[$count/$total] Running: $domain $topo seed=$seed -> $outfile"
      cabal run topology-sim -- $DOMAIN_FLAGS "$topo" $ISLANDS $POP $MIG_INTERVAL $MIGRANTS $GENS "$seed" \
        > "$outfile" 2>/dev/null
    done
  done
done

echo ""
echo "All $total runs complete. Results in $OUTDIR/"
echo ""
echo "Conditions (all n=12, |E|=13, beta_1=2):"
echo "  fig8-m1: C3+C5, mean cycle length 4.0, pendant chain (5 nodes)"
echo "  fig8-m2: C5+C6, mean cycle length 5.5, pendant chain (2 nodes)"
echo "  fig8-m3: C7+C6, mean cycle length 6.5, no pendant"
echo ""
echo "Domains: onemax (neg. control), nk4 (rugged), maze (deceptive)"
