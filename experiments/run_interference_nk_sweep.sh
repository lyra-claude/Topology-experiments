#!/usr/bin/env bash
# Interference experiment across NK roughness levels.
#
# Tests whether the 10:1 ratio (primary beta_1 effect vs arrangement effect)
# holds across NK roughness levels. Clio's cactus group prediction:
# arrangement effect should GROW with landscape ruggedness.
#
# 3 directed topologies, all with n=8, |E|=9, beta_1=2, kappa=2:
#   1. interference-adjacent  — two cycles sharing a vertex (figure-eight)
#   2. interference-separated — two cycles connected by a path (no shared vertices)
#   3. interference-nested    — one cycle inside another (concentric)
#
# NK0, NK2, NK4, NK6: 3 topologies x 30 seeds x 4 NK levels = 360 runs
# (NK4 already done in results/interference/ — skip or redo as needed)

set -euo pipefail

export PATH="/usr/local/.ghcup/bin:$PATH"

TOPOLOGIES="interference-adjacent interference-separated interference-nested"

# Same 30 seeds as all other experiments
SEEDS="42 137 2718 314 1618 7 99 256 512 1024 2048 4096 8192 16384 32768 11 23 37 53 71 89 97 113 131 151 173 191 199 211 223"

ISLANDS=8
POP=50
MIG_INTERVAL=10
MIGRANTS=5
GENS=500

# Which NK levels to run (skip nk4 since it exists in results/interference/)
NK_LEVELS="${1:-nk0 nk2 nk6}"

for nk in $NK_LEVELS; do
  OUTDIR="results/interference-${nk}"
  mkdir -p "$OUTDIR"

  total=90  # 3 topologies x 30 seeds
  count=0

  echo "=== Starting $nk interference runs ==="

  for topo in $TOPOLOGIES; do
    for seed in $SEEDS; do
      count=$((count + 1))
      outfile="${OUTDIR}/${topo}_seed${seed}.csv"
      if [ -f "$outfile" ]; then
        echo "[$count/$total] SKIP (exists): $nk $topo seed=$seed"
        continue
      fi
      echo "[$count/$total] Running: $nk $topo seed=$seed -> $outfile"
      cabal run topology-sim -- --domain "$nk" "$topo" $ISLANDS $POP $MIG_INTERVAL $MIGRANTS $GENS "$seed" \
        > "$outfile" 2>/dev/null
    done
  done

  echo "=== $nk complete: $total runs in $OUTDIR/ ==="
  echo ""
done

echo "All interference NK sweep runs complete."
