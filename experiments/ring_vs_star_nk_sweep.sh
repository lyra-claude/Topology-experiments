#!/usr/bin/env bash
# Ring vs Star NK Sweep
# Runs ring and star topologies across NK landscapes (K=0,2,4,6)
# with 20 seeds each for statistical power.
#
# Parameters: 8 islands, pop 50, migration every 10 gens,
#             5 migrants per event, 500 generations.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BINARY="$(source /home/lyra/.ghcup/env && cd "$PROJECT_DIR" && cabal list-bin topology-sim 2>/dev/null)"
OUTDIR="$PROJECT_DIR/results/ring_vs_star_nk"

mkdir -p "$OUTDIR"

TOPOLOGIES="ring star"
DOMAINS="nk0 nk2 nk4 nk6"
SEEDS=$(seq 1 20)
N_ISLANDS=8
POP_SIZE=50
MIG_INTERVAL=10
N_MIGRANTS=5
TOTAL_GENS=500

TOTAL_RUNS=$(echo "$TOPOLOGIES" | wc -w)
TOTAL_RUNS=$((TOTAL_RUNS * 4 * 20))
RUN=0

for domain in $DOMAINS; do
    for topo in $TOPOLOGIES; do
        for seed in $SEEDS; do
            RUN=$((RUN + 1))
            OUTFILE="$OUTDIR/${domain}_${topo}_seed${seed}.csv"
            if [ -f "$OUTFILE" ] && [ -s "$OUTFILE" ]; then
                echo "[$RUN/$TOTAL_RUNS] SKIP (exists): $domain $topo seed=$seed"
                continue
            fi
            echo "[$RUN/$TOTAL_RUNS] Running: $domain $topo seed=$seed"
            $BINARY --domain "$domain" "$topo" $N_ISLANDS $POP_SIZE $MIG_INTERVAL $N_MIGRANTS $TOTAL_GENS "$seed" > "$OUTFILE" 2>/dev/null
        done
    done
done

echo ""
echo "All runs complete. Results in $OUTDIR/"
echo "Total CSV files: $(ls "$OUTDIR"/*.csv 2>/dev/null | wc -l)"
