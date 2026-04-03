#!/usr/bin/env python3
"""
Analyze directed cycle experiment pilot: 3 topologies x 5 seeds.
Compares DAG-Layer (0 cycles), Bidir-Ring (10 cycles), Ring-Skip2 (47 cycles).
"""

import os
import csv
from collections import defaultdict

RESULTS_DIR = "results/directed_pilot"

CYCLE_COUNTS = {
    "dag-layer": 0,
    "bidir-ring": 10,
    "ring-skip2": 47,
}

def load_csv(path):
    """Load CSV, return list of dicts. Skip non-CSV lines (e.g. 'Up to date')."""
    with open(path) as f:
        lines = [l for l in f if l.strip() and l.strip().startswith(("g", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))]
    reader = csv.DictReader(lines)
    return list(reader)

def main():
    # Collect results by topology
    results = defaultdict(list)  # topo -> list of run data

    for fname in sorted(os.listdir(RESULTS_DIR)):
        if not fname.endswith(".csv"):
            continue
        topo = fname.rsplit("_seed", 1)[0]
        rows = load_csv(os.path.join(RESULTS_DIR, fname))
        results[topo].append(rows)

    print("DIRECTED CYCLE EXPERIMENT — PILOT ANALYSIS")
    print("=" * 70)
    print(f"{'Topology':<20} {'Cycles':>7} {'Seeds':>6} {'Final Mean':>11} {'Final Best':>11} {'Final Div':>10} {'Gen50 Div':>10}")
    print("-" * 70)

    for topo in sorted(results.keys(), key=lambda t: CYCLE_COUNTS.get(t, 999)):
        runs = results[topo]
        cycles = CYCLE_COUNTS.get(topo, "?")

        final_means = []
        final_bests = []
        final_divs = []
        gen50_divs = []

        for run in runs:
            if run:
                last = run[-1]
                final_means.append(float(last["meanFitness"]))
                final_bests.append(float(last["bestFitness"]))
                final_divs.append(float(last["diversity"]))

                # Find gen ~50
                for row in run:
                    if int(row["generation"]) == 50:
                        gen50_divs.append(float(row["diversity"]))
                        break

        n = len(runs)
        avg_mean = sum(final_means) / n if n else 0
        avg_best = sum(final_bests) / n if n else 0
        avg_div = sum(final_divs) / n if n else 0
        avg_g50_div = sum(gen50_divs) / len(gen50_divs) if gen50_divs else 0

        print(f"{topo:<20} {cycles:>7} {n:>6} {avg_mean:>11.6f} {avg_best:>11.6f} {avg_div:>10.6f} {avg_g50_div:>10.6f}")

    print()

    # Detailed per-gen comparison at key checkpoints
    print("\nPER-GENERATION COMPARISON (averaged over seeds)")
    print(f"{'Gen':>5}", end="")
    for topo in sorted(results.keys(), key=lambda t: CYCLE_COUNTS.get(t, 999)):
        cycles = CYCLE_COUNTS.get(topo, "?")
        print(f"  {topo}({cycles}cy) mean/div", end="")
    print()
    print("-" * 120)

    checkpoints = [0, 10, 20, 30, 50, 100, 200, 300, 500]
    for gen in checkpoints:
        print(f"{gen:>5}", end="")
        for topo in sorted(results.keys(), key=lambda t: CYCLE_COUNTS.get(t, 999)):
            runs = results[topo]
            means = []
            divs = []
            for run in runs:
                for row in run:
                    if int(row["generation"]) == gen:
                        means.append(float(row["meanFitness"]))
                        divs.append(float(row["diversity"]))
                        break
            if means:
                avg_m = sum(means) / len(means)
                avg_d = sum(divs) / len(divs)
                print(f"  {avg_m:.4f}/{avg_d:.4f}         ", end="")
            else:
                print(f"  {'N/A':>20}         ", end="")
        print()

    print()
    print("KEY QUESTION: Does cycle count predict diversity retention at constant density?")
    print("If DAG has different diversity trajectory than Ring-Skip2, the experiment works.")

if __name__ == "__main__":
    main()
