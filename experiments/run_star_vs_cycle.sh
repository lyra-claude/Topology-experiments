#!/usr/bin/env bash
# Star_n vs Cycle_{n-1} experiment: edge-controlled topology comparison
#
# KEY INSIGHT: star(n) has n-1 edges, cycle(n-1) has n-1 edges.
# Same edge count, but star has β₁=0 and cycle has β₁=1.
# This isolates the topology effect (cycle rank) from edge count.
#
# Conditions:
#   star_8 (8 nodes, 7 edges, β₁=0) vs cycle_7 (7 nodes, 7 edges, β₁=1)
#   star_12 (12 nodes, 11 edges, β₁=0) vs cycle_11 (11 nodes, 11 edges, β₁=1)
#
# NK landscapes: K=0 (smooth), K=2 (moderate), K=4 (rugged)
# 30 seeds each
# Total: 2 sizes × 2 topologies × 3 K values × 30 seeds = 360 runs

set -euo pipefail

export PATH="/home/lyra/.ghcup/bin:$PATH"

OUTDIR="results/star_vs_cycle"
mkdir -p "$OUTDIR"

# 3 NK domains
DOMAINS="nk0 nk2 nk4"

# 30 seeds (same as other experiments for comparability)
SEEDS="42 137 2718 314 1618 7 99 256 512 1024 2048 4096 8192 16384 32768 11 23 37 53 71 89 97 113 131 151 173 191 199 211 223"

POP=50
MIG_INTERVAL=10
MIGRANTS=5
GENS=500

# Two size conditions: (star_n, cycle_{n-1})
# Size 8: star 8 nodes vs cycle 7 nodes (both 7 edges)
# Size 12: star 12 nodes vs cycle 11 nodes (both 11 edges)

total=360  # 2 sizes x 2 topos x 3 domains x 30 seeds
count=0

for domain in $DOMAINS; do
  # --- Size condition: 8/7 (7 edges) ---
  for seed in $SEEDS; do
    count=$((count + 1))
    outfile="${OUTDIR}/star8_${domain}_seed${seed}.csv"
    if [ -f "$outfile" ]; then
      echo "[$count/$total] SKIP (exists): star8 $domain seed=$seed"
      continue
    fi
    echo "[$count/$total] Running: star8 $domain seed=$seed -> $outfile"
    cabal run topology-sim -- --domain "$domain" star 8 $POP $MIG_INTERVAL $MIGRANTS $GENS "$seed" \
      > "$outfile" 2>/dev/null
  done

  for seed in $SEEDS; do
    count=$((count + 1))
    outfile="${OUTDIR}/cycle7_${domain}_seed${seed}.csv"
    if [ -f "$outfile" ]; then
      echo "[$count/$total] SKIP (exists): cycle7 $domain seed=$seed"
      continue
    fi
    echo "[$count/$total] Running: cycle7 $domain seed=$seed -> $outfile"
    cabal run topology-sim -- --domain "$domain" ring 7 $POP $MIG_INTERVAL $MIGRANTS $GENS "$seed" \
      > "$outfile" 2>/dev/null
  done

  # --- Size condition: 12/11 (11 edges) ---
  for seed in $SEEDS; do
    count=$((count + 1))
    outfile="${OUTDIR}/star12_${domain}_seed${seed}.csv"
    if [ -f "$outfile" ]; then
      echo "[$count/$total] SKIP (exists): star12 $domain seed=$seed"
      continue
    fi
    echo "[$count/$total] Running: star12 $domain seed=$seed -> $outfile"
    cabal run topology-sim -- --domain "$domain" star 12 $POP $MIG_INTERVAL $MIGRANTS $GENS "$seed" \
      > "$outfile" 2>/dev/null
  done

  for seed in $SEEDS; do
    count=$((count + 1))
    outfile="${OUTDIR}/cycle11_${domain}_seed${seed}.csv"
    if [ -f "$outfile" ]; then
      echo "[$count/$total] SKIP (exists): cycle11 $domain seed=$seed"
      continue
    fi
    echo "[$count/$total] Running: cycle11 $domain seed=$seed -> $outfile"
    cabal run topology-sim -- --domain "$domain" ring 11 $POP $MIG_INTERVAL $MIGRANTS $GENS "$seed" \
      > "$outfile" 2>/dev/null
  done
done

echo ""
echo "All $total runs complete. Results in $OUTDIR/"
echo ""
echo "Conditions:"
echo "  star8:   8 nodes, 7 edges, β₁=0"
echo "  cycle7:  7 nodes, 7 edges, β₁=1"
echo "  star12: 12 nodes, 11 edges, β₁=0"
echo "  cycle11: 11 nodes, 11 edges, β₁=1"
