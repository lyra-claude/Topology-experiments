#!/usr/bin/env python3
"""Analyze NK pilot experiment results.

Key prediction: eta-squared (effect size of topology) increases with K.
Higher K = more rugged landscape = topology matters more.

Computes one-way ANOVA (topology as factor) for each K value at multiple
generations, reports eta-squared.
"""

import os
import csv
import numpy as np
from collections import defaultdict
from scipy import stats

RESULTS_DIR = "/home/lyra/projects/Topology-experiments/results/nk_pilot"

DOMAINS = ["nk0", "nk4", "nk6"]
TOPOLOGIES = ["dag-layer", "bidir-ring", "ring-skip2"]
SEEDS = [42, 137, 2718, 314, 1618, 7, 99, 256, 512, 1024]

# Generations to analyze
ANALYSIS_GENS = [50, 100, 200, 300, 500]


def load_run(domain, topo, seed):
    """Load a single run CSV. Returns dict: gen -> {meanFitness, bestFitness, diversity}."""
    path = os.path.join(RESULTS_DIR, f"{domain}_{topo}_{seed}.csv")
    data = {}
    with open(path) as f:
        # Skip any non-CSV lines (e.g., "Up to date" from cabal)
        lines = f.readlines()
        csv_lines = [l for l in lines if not l.startswith("Up to date")]
        reader = csv.DictReader(csv_lines)
        for row in reader:
            gen = int(row["generation"])
            data[gen] = {
                "meanFitness": float(row["meanFitness"]),
                "bestFitness": float(row["bestFitness"]),
                "diversity": float(row["diversity"]),
            }
    return data


def compute_eta_squared(groups):
    """Compute eta-squared from a list of arrays (one per group)."""
    # One-way ANOVA
    f_stat, p_val = stats.f_oneway(*groups)
    # eta^2 = SS_between / SS_total
    all_vals = np.concatenate(groups)
    grand_mean = np.mean(all_vals)
    ss_total = np.sum((all_vals - grand_mean) ** 2)
    ss_between = sum(len(g) * (np.mean(g) - grand_mean) ** 2 for g in groups)
    eta_sq = ss_between / ss_total if ss_total > 0 else 0.0
    return f_stat, p_val, eta_sq


def main():
    # Load all data
    all_data = {}  # (domain, topo, seed) -> gen_data
    for domain in DOMAINS:
        for topo in TOPOLOGIES:
            for seed in SEEDS:
                all_data[(domain, topo, seed)] = load_run(domain, topo, seed)

    print("=" * 80)
    print("NK PILOT EXPERIMENT — TOPOLOGY EFFECT ANALYSIS")
    print("=" * 80)
    print(f"\nDomains: {DOMAINS}")
    print(f"Topologies: {TOPOLOGIES} (0, 10, 47 directed cycles)")
    print(f"Seeds per condition: {len(SEEDS)}")
    print(f"Total runs: {len(all_data)}")

    # For each metric, for each domain, for each generation: ANOVA across topologies
    for metric_name, metric_key in [("Best Fitness", "bestFitness"),
                                     ("Mean Fitness", "meanFitness"),
                                     ("Diversity", "diversity")]:
        print(f"\n{'=' * 80}")
        print(f"METRIC: {metric_name}")
        print(f"{'=' * 80}")

        # Table header
        print(f"\n{'Gen':>5} | ", end="")
        for domain in DOMAINS:
            print(f"  {domain:>6} eta²  {'p':>8}  F", end="  | ")
        print()
        print("-" * 80)

        for gen in ANALYSIS_GENS:
            print(f"{gen:>5} | ", end="")
            for domain in DOMAINS:
                groups = []
                for topo in TOPOLOGIES:
                    vals = []
                    for seed in SEEDS:
                        d = all_data[(domain, topo, seed)]
                        if gen in d:
                            vals.append(d[gen][metric_key])
                    groups.append(np.array(vals))

                f_stat, p_val, eta_sq = compute_eta_squared(groups)
                sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
                print(f"  {eta_sq:>6.4f}  {p_val:>8.2e} {f_stat:>6.2f}{sig:<3}", end="  | ")
            print()

    # Summary: eta-squared at gen 500 for best fitness
    print(f"\n{'=' * 80}")
    print("SUMMARY: eta² for Best Fitness at gen 500 (key prediction test)")
    print(f"{'=' * 80}")
    print(f"\nPrediction: eta² should INCREASE with K (nk0 < nk4 < nk6)")
    print()

    eta_values = {}
    for domain in DOMAINS:
        groups = []
        for topo in TOPOLOGIES:
            vals = [all_data[(domain, topo, seed)][500]["bestFitness"] for seed in SEEDS]
            groups.append(np.array(vals))
        f_stat, p_val, eta_sq = compute_eta_squared(groups)
        eta_values[domain] = eta_sq
        print(f"  {domain}: eta² = {eta_sq:.4f}, F = {f_stat:.3f}, p = {p_val:.2e}")

    # Test the prediction
    print()
    if eta_values["nk0"] < eta_values["nk4"] < eta_values["nk6"]:
        print("  PREDICTION CONFIRMED: eta² increases monotonically with K")
    elif eta_values["nk0"] < eta_values["nk6"]:
        print("  PREDICTION PARTIALLY CONFIRMED: nk6 > nk0, but not strictly monotonic")
    else:
        print("  PREDICTION NOT CONFIRMED: eta² does NOT increase with K")

    # Per-topology means at gen 500
    print(f"\n{'=' * 80}")
    print("TOPOLOGY MEANS at gen 500")
    print(f"{'=' * 80}")
    for metric_name, metric_key in [("Best Fitness", "bestFitness"),
                                     ("Diversity", "diversity")]:
        print(f"\n{metric_name}:")
        print(f"  {'Topology':<15}", end="")
        for domain in DOMAINS:
            print(f"  {domain:>10}", end="")
        print()
        for topo in TOPOLOGIES:
            print(f"  {topo:<15}", end="")
            for domain in DOMAINS:
                vals = [all_data[(domain, topo, seed)][500][metric_key] for seed in SEEDS]
                print(f"  {np.mean(vals):>10.6f}", end="")
            print(f"  (std: {np.std(vals):.4f})")

    # Kruskal-Wallis as non-parametric check
    print(f"\n{'=' * 80}")
    print("KRUSKAL-WALLIS (non-parametric) at gen 500, Best Fitness")
    print(f"{'=' * 80}")
    for domain in DOMAINS:
        groups = []
        for topo in TOPOLOGIES:
            vals = [all_data[(domain, topo, seed)][500]["bestFitness"] for seed in SEEDS]
            groups.append(vals)
        h_stat, p_val = stats.kruskal(*groups)
        print(f"  {domain}: H = {h_stat:.3f}, p = {p_val:.2e}")

    # Temporal dynamics: eta² over all generations for best fitness
    print(f"\n{'=' * 80}")
    print("TEMPORAL DYNAMICS: eta² (best fitness) over generations")
    print(f"{'=' * 80}")
    all_gens = sorted(all_data[(DOMAINS[0], TOPOLOGIES[0], SEEDS[0])].keys())
    print(f"\n{'Gen':>5}", end="")
    for domain in DOMAINS:
        print(f"  {domain:>8}", end="")
    print()
    print("-" * 40)
    for gen in all_gens[::5]:  # every 50 gens
        print(f"{gen:>5}", end="")
        for domain in DOMAINS:
            groups = []
            for topo in TOPOLOGIES:
                vals = [all_data[(domain, topo, seed)][gen]["bestFitness"] for seed in SEEDS]
                groups.append(np.array(vals))
            _, _, eta_sq = compute_eta_squared(groups)
            print(f"  {eta_sq:>8.4f}", end="")
        print()


if __name__ == "__main__":
    main()
