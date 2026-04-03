#!/usr/bin/env bash
# Foster census cubic symmetric graph sweep
# 30 runs x 13 topologies x 500 generations
# Each topology uses its natural vertex count as island count.
# All graphs are 3-regular (cubic symmetric).

set -euo pipefail

export PATH="$HOME/.ghcup/bin:$PATH"

OUTDIR="results/foster_sweep"
mkdir -p "$OUTDIR"

# 13 Foster census cubic symmetric graphs, ordered by vertex count.
# Format: "name:islands"
TOPOLOGIES=(
  "k4:4"
  "k33:6"
  "cube:8"
  "petersen:10"
  "heawood:14"
  "mobius-kantor:16"
  "pappus:18"
  "dodecahedron:20"
  "desargues:20"
  "nauru:24"
  "f26a:26"
  "coxeter:28"
  "tutte-coxeter:30"
)

# 30 seeds for statistical power (same as directed experiment)
SEEDS="42 137 2718 314 1618 7 99 256 512 1024 2048 4096 8192 16384 32768 11 23 37 53 71 89 97 113 131 151 173 191 199 211 223"

POP=100
MIG_INTERVAL=10
MIGRANTS=1
GENS=500

total=$((13 * 30))  # 13 topologies x 30 seeds = 390
count=0

for entry in "${TOPOLOGIES[@]}"; do
  topo="${entry%%:*}"
  islands="${entry##*:}"
  for seed in $SEEDS; do
    count=$((count + 1))
    outfile="${OUTDIR}/${topo}_seed${seed}.csv"
    if [ -f "$outfile" ]; then
      echo "[$count/$total] SKIP (exists): $topo seed=$seed"
      continue
    fi
    echo "[$count/$total] Running: $topo (${islands} islands) seed=$seed -> $outfile"
    cabal run topology-sim -- --domain onemax "$topo" "$islands" $POP $MIG_INTERVAL $MIGRANTS $GENS "$seed" \
      > "$outfile" 2>/dev/null
  done
done

echo ""
echo "All $total runs complete. Results in $OUTDIR/"
