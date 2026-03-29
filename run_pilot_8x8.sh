#!/usr/bin/env bash
# Pilot study: 3 runs x 8 topologies x 500 generations (8x8 maze)
# BFS-only fitness, 8 islands, pop=50, migInterval=10, migrants=5

set -euo pipefail

RESULTS_DIR="/home/lyra/projects/Topology-experiments/results/pilot_batch1_8x8"
mkdir -p "$RESULTS_DIR"

# Find the binary
BINARY=$(cabal list-bin topology-sim 2>/dev/null)
echo "Binary: $BINARY"
echo "Results directory: $RESULTS_DIR"

# Topologies
TOPOS=(disconnected ring star complete hypercube barbell watts-strogatz random-regular)

# Seeds for 3 independent runs
SEEDS=(42 137 2718)

# GA parameters
ISLANDS=8
POP=50
MIG_INTERVAL=10
MIGRANTS=5
GENS=500
GRID=8

echo "Starting 8x8 pilot study: ${#TOPOS[@]} topologies x ${#SEEDS[@]} runs x $GENS generations"
echo "Parameters: grid=${GRID}x${GRID} islands=$ISLANDS pop=$POP migInterval=$MIG_INTERVAL migrants=$MIGRANTS"
echo ""

OVERALL_START=$(date +%s)

# Run all topology/seed combinations in parallel (up to 24 at a time)
PIDS=()
for TOPO in "${TOPOS[@]}"; do
    for SEED in "${SEEDS[@]}"; do
        OUTFILE="$RESULTS_DIR/${TOPO}_seed${SEED}.csv"
        TIMEFILE="$RESULTS_DIR/${TOPO}_seed${SEED}.time"

        (
            START=$(date +%s%N)
            "$BINARY" --grid "$GRID" "$TOPO" "$ISLANDS" "$POP" "$MIG_INTERVAL" "$MIGRANTS" "$GENS" "$SEED" \
                > "$OUTFILE" 2> "$RESULTS_DIR/${TOPO}_seed${SEED}.log"
            END=$(date +%s%N)
            ELAPSED_MS=$(( (END - START) / 1000000 ))
            echo "$ELAPSED_MS" > "$TIMEFILE"
            echo "  Done: $TOPO seed=$SEED (${ELAPSED_MS}ms)"
        ) &
        PIDS+=($!)
    done
done

echo "Launched ${#PIDS[@]} parallel jobs. Waiting..."
echo ""

# Wait for all jobs
FAILED=0
for PID in "${PIDS[@]}"; do
    if ! wait "$PID"; then
        FAILED=$((FAILED + 1))
    fi
done

OVERALL_END=$(date +%s)
OVERALL_ELAPSED=$((OVERALL_END - OVERALL_START))

echo ""
echo "All jobs complete. Failed: $FAILED / ${#PIDS[@]}"
echo "Total wall-clock time: ${OVERALL_ELAPSED}s"
echo ""

# Summary: list all output files with line counts
echo "Output files:"
for TOPO in "${TOPOS[@]}"; do
    for SEED in "${SEEDS[@]}"; do
        OUTFILE="$RESULTS_DIR/${TOPO}_seed${SEED}.csv"
        TIMEFILE="$RESULTS_DIR/${TOPO}_seed${SEED}.time"
        LINES=$(wc -l < "$OUTFILE" 2>/dev/null || echo "0")
        TIME_MS=$(cat "$TIMEFILE" 2>/dev/null || echo "?")
        echo "  $TOPO seed=$SEED: $LINES lines, ${TIME_MS}ms"
    done
done
