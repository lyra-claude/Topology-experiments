#!/usr/bin/env bash
# Cycle Length Experiment: Does mean cycle length predict diversity independently of beta_1?
#
# Three two-cycle-bridge graphs, all beta_1=2, varying mean cycle length:
#   G1: cycle-bridge-3-3  (mean_cycle=3.0,  n=6,  lambda_2=0.4384)
#   G2: cycle-bridge-3-9  (mean_cycle=6.0,  n=12, lambda_2=0.1907)
#   G3: cycle-bridge-9-9  (mean_cycle=9.0,  n=18, lambda_2=0.0822)
#
# 3 topologies x 2 NK domains (K=0, K=4) x 10 seeds = 60 runs.
# Parameters: pop 50 per island, migration every 10 gens, 5 migrants, 500 gens.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BINARY="$(source /home/lyra/.ghcup/env && cd "$PROJECT_DIR" && cabal list-bin topology-sim 2>/dev/null)"
OUTDIR="$PROJECT_DIR/results/cycle_length"

mkdir -p "$OUTDIR"

# Topology name -> island count (each graph has a different number of nodes)
declare -A TOPO_ISLANDS
TOPO_ISLANDS["cycle-bridge-3-3"]=6
TOPO_ISLANDS["cycle-bridge-3-9"]=12
TOPO_ISLANDS["cycle-bridge-9-9"]=18

TOPOLOGIES="cycle-bridge-3-3 cycle-bridge-3-9 cycle-bridge-9-9"
DOMAINS="nk0 nk4"
SEEDS=$(seq 1 10)
POP_SIZE=50
MIG_INTERVAL=10
N_MIGRANTS=5
TOTAL_GENS=500

TOTAL_RUNS=$((3 * 2 * 10))
RUN=0

for domain in $DOMAINS; do
    for topo in $TOPOLOGIES; do
        N_ISLANDS=${TOPO_ISLANDS[$topo]}
        for seed in $SEEDS; do
            RUN=$((RUN + 1))
            OUTFILE="$OUTDIR/${domain}_${topo}_seed${seed}.csv"
            if [ -f "$OUTFILE" ] && [ -s "$OUTFILE" ]; then
                echo "[$RUN/$TOTAL_RUNS] SKIP (exists): $domain $topo seed=$seed"
                continue
            fi
            echo "[$RUN/$TOTAL_RUNS] Running: $domain $topo (n=$N_ISLANDS) seed=$seed"
            $BINARY --domain "$domain" "$topo" $N_ISLANDS $POP_SIZE $MIG_INTERVAL $N_MIGRANTS $TOTAL_GENS "$seed" > "$OUTFILE" 2>/dev/null
        done
    done
done

echo ""
echo "All runs complete. Results in $OUTDIR/"
echo "Total CSV files: $(ls "$OUTDIR"/*.csv 2>/dev/null | wc -l)"
