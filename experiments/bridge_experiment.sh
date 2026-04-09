#!/usr/bin/env bash
# Bridge Experiment: Disentangling beta_1 from lambda_2
#
# Two iso-spectral families (constant lambda_2, varying beta_1):
#   Family 1 (lambda_2=0.5858): ring, ring-chord1, ring-chord2, ring-chord3
#   Family 2 (lambda_2=1.0):    star, star-leaf1, star-leaf2, star-leaf3
#
# 8 topologies x 2 NK domains (K=0, K=4) x 20 seeds = 320 runs.
# Parameters: 8 islands, pop 50, migration every 10 gens, 5 migrants, 500 gens.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BINARY="$(source /home/lyra/.ghcup/env && cd "$PROJECT_DIR" && cabal list-bin topology-sim 2>/dev/null)"
OUTDIR="$PROJECT_DIR/results/bridge_experiment"

mkdir -p "$OUTDIR"

TOPOLOGIES="ring ring-chord1 ring-chord2 ring-chord3 star star-leaf1 star-leaf2 star-leaf3"
DOMAINS="nk0 nk4"
SEEDS=$(seq 1 20)
N_ISLANDS=8
POP_SIZE=50
MIG_INTERVAL=10
N_MIGRANTS=5
TOTAL_GENS=500

TOTAL_RUNS=$((8 * 2 * 20))
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
