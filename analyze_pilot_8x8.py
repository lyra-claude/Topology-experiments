#!/usr/bin/env python3
"""Analyze pilot_batch1_8x8 results: diversity ordering vs lambda_2."""

import csv
import os
import statistics

RESULTS_DIR = "/home/lyra/projects/Topology-experiments/results/pilot_batch1_8x8"

TOPOS = [
    "disconnected", "ring", "star", "complete",
    "hypercube", "barbell", "watts-strogatz", "random-regular"
]
SEEDS = [42, 137, 2718]

# Known lambda_2 values for 8-island topologies (algebraic connectivity)
# disconnected: 0 (no edges)
# ring: 2 - 2*cos(2*pi/8) = 2 - 2*cos(pi/4) = 2 - sqrt(2) ~ 0.586
# star: 1.0
# complete: 8.0 (lambda_2 = n for complete graph)
# hypercube (k=3): 2.0
# barbell: small (bottleneck) ~ 0.235 (from Cheeger bound; two K4 + bridge)
# watts-strogatz: depends on rewiring, typically between ring and small-world
# random-regular (3-regular): typically ~0.5-1.5

LAMBDA2 = {
    "disconnected": 0.0,
    "ring": 0.586,
    "star": 1.0,
    "complete": 8.0,
    "hypercube": 2.0,
    "barbell": 0.235,
    "watts-strogatz": 1.2,  # approximate
    "random-regular": 0.8,  # approximate
}


def load_csv(filepath):
    """Load a CSV and return list of dicts."""
    rows = []
    with open(filepath) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                "generation": int(row["generation"]),
                "meanFitness": float(row["meanFitness"]),
                "bestFitness": float(row["bestFitness"]),
                "diversity": float(row["diversity"]),
            })
    return rows


def get_diversity_at_gen(rows, gen):
    """Get diversity at a specific generation."""
    for row in rows:
        if row["generation"] == gen:
            return row["diversity"]
    return None


def main():
    print("=" * 80)
    print("8x8 Maze Pilot Study Results (3 runs x 8 topologies x 500 gens)")
    print("=" * 80)
    print()

    # Collect final diversity for each topology across seeds
    topo_data = {}
    for topo in TOPOS:
        diversities_0 = []
        diversities_100 = []
        diversities_250 = []
        diversities_500 = []
        min_diversities = []
        fitness_500 = []
        for seed in SEEDS:
            filepath = os.path.join(RESULTS_DIR, f"{topo}_seed{seed}.csv")
            rows = load_csv(filepath)
            d0 = get_diversity_at_gen(rows, 0)
            d100 = get_diversity_at_gen(rows, 100)
            d250 = get_diversity_at_gen(rows, 250)
            d500 = get_diversity_at_gen(rows, 500)
            min_d = min(r["diversity"] for r in rows)
            f500 = get_diversity_at_gen(rows, 500)  # using diversity field
            # Get actual fitness at 500
            for r in rows:
                if r["generation"] == 500:
                    f500_real = r["bestFitness"]
            diversities_0.append(d0)
            diversities_100.append(d100)
            diversities_250.append(d250)
            diversities_500.append(d500)
            min_diversities.append(min_d)
            fitness_500.append(f500_real)

        topo_data[topo] = {
            "div_0_mean": statistics.mean(diversities_0),
            "div_100_mean": statistics.mean(diversities_100),
            "div_250_mean": statistics.mean(diversities_250),
            "div_500_mean": statistics.mean(diversities_500),
            "div_500_std": statistics.stdev(diversities_500) if len(diversities_500) > 1 else 0,
            "div_min_mean": statistics.mean(min_diversities),
            "fitness_500_mean": statistics.mean(fitness_500),
            "lambda2": LAMBDA2.get(topo, "?"),
            "div_decay": statistics.mean(diversities_0) - statistics.mean(diversities_500),
        }

    # Print diversity at gen 500, sorted by diversity (ascending = most decay)
    print("Diversity at gen 500 (mean +/- std across 3 runs), sorted by diversity:")
    print("-" * 80)
    print(f"{'Topology':<20} {'lambda_2':>8} {'div@0':>8} {'div@100':>8} {'div@250':>8} {'div@500':>8} {'std':>8} {'decay':>8}")
    print("-" * 80)
    sorted_topos = sorted(TOPOS, key=lambda t: topo_data[t]["div_500_mean"])
    for topo in sorted_topos:
        d = topo_data[topo]
        print(f"{topo:<20} {d['lambda2']:>8.3f} {d['div_0_mean']:>8.4f} {d['div_100_mean']:>8.4f} {d['div_250_mean']:>8.4f} {d['div_500_mean']:>8.4f} {d['div_500_std']:>8.4f} {d['div_decay']:>8.4f}")
    print()

    # Print diversity decay ordering vs lambda_2 ordering
    print("Diversity decay ordering (most decay first):")
    sorted_by_decay = sorted(TOPOS, key=lambda t: topo_data[t]["div_decay"], reverse=True)
    for i, topo in enumerate(sorted_by_decay, 1):
        d = topo_data[topo]
        print(f"  {i}. {topo:<20} decay={d['div_decay']:.4f}  lambda_2={d['lambda2']:.3f}")
    print()

    print("Lambda_2 ordering (highest first -- should lose diversity fastest):")
    sorted_by_lambda = sorted(TOPOS, key=lambda t: topo_data[t]["lambda2"], reverse=True)
    for i, topo in enumerate(sorted_by_lambda, 1):
        d = topo_data[topo]
        print(f"  {i}. {topo:<20} lambda_2={d['lambda2']:.3f}  decay={d['div_decay']:.4f}")
    print()

    # Compare 8x8 vs 15x15 (load 15x15 if available)
    pilot15_dir = "/home/lyra/projects/Topology-experiments/results/pilot_batch1"
    if os.path.isdir(pilot15_dir):
        print("=" * 80)
        print("Comparison: 8x8 vs 15x15 (mean diversity @ gen 500)")
        print("-" * 80)
        print(f"{'Topology':<20} {'8x8 div@500':>12} {'15x15 div@500':>14} {'8x8 decay':>10} {'15x15 decay':>12}")
        print("-" * 80)
        for topo in TOPOS:
            div15_list = []
            d0_15_list = []
            for seed in SEEDS:
                fp15 = os.path.join(pilot15_dir, f"{topo}_seed{seed}.csv")
                if os.path.exists(fp15):
                    rows15 = load_csv(fp15)
                    d15_500 = get_diversity_at_gen(rows15, 500)
                    d15_0 = get_diversity_at_gen(rows15, 0)
                    if d15_500 is not None:
                        div15_list.append(d15_500)
                    if d15_0 is not None:
                        d0_15_list.append(d15_0)
            d8 = topo_data[topo]
            if div15_list:
                d15_mean = statistics.mean(div15_list)
                decay15 = statistics.mean(d0_15_list) - d15_mean if d0_15_list else 0
                print(f"{topo:<20} {d8['div_500_mean']:>12.4f} {d15_mean:>14.4f} {d8['div_decay']:>10.4f} {decay15:>12.4f}")
            else:
                print(f"{topo:<20} {d8['div_500_mean']:>12.4f} {'N/A':>14} {d8['div_decay']:>10.4f} {'N/A':>12}")
        print()

    # Minimum diversity reached (proxy for how much the topology can compress)
    print("=" * 80)
    print("Minimum diversity reached during 500 gens (mean across 3 runs):")
    print("-" * 80)
    sorted_by_min = sorted(TOPOS, key=lambda t: topo_data[t]["div_min_mean"])
    for topo in sorted_by_min:
        d = topo_data[topo]
        print(f"  {topo:<20} min_div={d['div_min_mean']:.4f}  lambda_2={d['lambda2']:.3f}")
    print()

    # Best fitness at gen 500
    print("Best fitness at gen 500 (mean across 3 runs):")
    print("-" * 80)
    sorted_by_fit = sorted(TOPOS, key=lambda t: topo_data[t]["fitness_500_mean"], reverse=True)
    for topo in sorted_by_fit:
        d = topo_data[topo]
        print(f"  {topo:<20} bestFit@500={d['fitness_500_mean']:.4f}  lambda_2={d['lambda2']:.3f}")


if __name__ == "__main__":
    main()
