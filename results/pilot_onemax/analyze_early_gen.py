#!/usr/bin/env python3
"""
Analyze early-generation dynamics of OneMax pilot study.
Look for topology effects in the transient window (gen 0-50).
"""

import csv
import os
from collections import defaultdict
from scipy import stats
import math

DATA_DIR = "/home/lyra/projects/Topology-experiments/results/pilot_onemax/"

TOPOLOGIES = [
    "disconnected", "barbell", "random-regular", "ring",
    "star", "watts-strogatz", "hypercube", "complete"
]

SEEDS = [42, 137, 2718]

# Lambda_2 values (algebraic connectivity)
LAMBDA2 = {
    "disconnected": 0.0,
    "barbell": 0.07,
    "random-regular": 0.27,
    "ring": 0.59,
    "star": 1.0,
    "watts-strogatz": 1.5,
    "hypercube": 2.0,
    "complete": 8.0,
}

def read_csv(filepath):
    """Read a CSV file into a list of dicts with numeric values."""
    rows = []
    with open(filepath, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                "generation": int(row["generation"]),
                "meanFitness": float(row["meanFitness"]),
                "bestFitness": float(row["bestFitness"]),
                "diversity": float(row["diversity"]),
            })
    return rows

def get_gen_value(rows, gen, field):
    """Get a specific field value at a specific generation."""
    for row in rows:
        if row["generation"] == gen:
            return row[field]
    return None

def first_gen_reaching(rows, field, threshold):
    """Find the first generation where field >= threshold."""
    for row in rows:
        if row[field] >= threshold:
            return row["generation"]
    return None  # Never reached

# ============================================================
# Step 1: Load all data
# ============================================================
data = {}  # data[(topology, seed)] = list of rows
for topo in TOPOLOGIES:
    for seed in SEEDS:
        fname = f"{topo}_seed{seed}.csv"
        fpath = os.path.join(DATA_DIR, fname)
        if os.path.exists(fpath):
            data[(topo, seed)] = read_csv(fpath)
        else:
            print(f"WARNING: Missing file {fname}")

# ============================================================
# Step 2: Compute per-topology metrics (averaged across seeds)
# ============================================================
print("=" * 80)
print("EARLY-GENERATION ANALYSIS: OneMax Pilot Study")
print("=" * 80)

# 2a: Generation to reach best fitness 1.0
print("\n--- Generation to First Best Fitness = 1.0 ---")
print(f"{'Topology':<20} {'Seed42':>8} {'Seed137':>8} {'Seed2718':>8} {'Mean':>8} {'Lambda2':>8}")
gen_best_1 = {}
for topo in TOPOLOGIES:
    vals = []
    row_strs = []
    for seed in SEEDS:
        g = first_gen_reaching(data[(topo, seed)], "bestFitness", 1.0)
        vals.append(g)
        row_strs.append(f"{g:>8}")
    mean_val = sum(vals) / len(vals)
    gen_best_1[topo] = mean_val
    print(f"{topo:<20} {'  '.join([str(v) for v in vals]):>26} {mean_val:>8.1f} {LAMBDA2[topo]:>8.2f}")

# 2b: Generation to reach mean fitness 0.9
print("\n--- Generation to First Mean Fitness >= 0.9 ---")
print(f"{'Topology':<20} {'Seed42':>8} {'Seed137':>8} {'Seed2718':>8} {'Mean':>8} {'Lambda2':>8}")
gen_mean_09 = {}
for topo in TOPOLOGIES:
    vals = []
    for seed in SEEDS:
        g = first_gen_reaching(data[(topo, seed)], "meanFitness", 0.9)
        vals.append(g)
    mean_val = sum(vals) / len(vals)
    gen_mean_09[topo] = mean_val
    print(f"{topo:<20} {vals[0]:>8} {vals[1]:>8} {vals[2]:>8} {mean_val:>8.1f} {LAMBDA2[topo]:>8.2f}")

# 2c: Diversity at gen 10, 20, 30, 50
div_at = {gen: {} for gen in [10, 20, 30, 50]}
print("\n--- Diversity at Key Generations (averaged over 3 runs) ---")
print(f"{'Topology':<20} {'Gen10':>10} {'Gen20':>10} {'Gen30':>10} {'Gen50':>10} {'Lambda2':>8}")
for topo in TOPOLOGIES:
    for gen in [10, 20, 30, 50]:
        vals = [get_gen_value(data[(topo, seed)], gen, "diversity") for seed in SEEDS]
        div_at[gen][topo] = sum(vals) / len(vals)
    print(f"{topo:<20} {div_at[10][topo]:>10.4f} {div_at[20][topo]:>10.4f} {div_at[30][topo]:>10.4f} {div_at[50][topo]:>10.4f} {LAMBDA2[topo]:>8.2f}")

# 2d: Mean fitness at gen 10, 20
mf_at = {gen: {} for gen in [10, 20]}
print("\n--- Mean Fitness at Key Generations (averaged over 3 runs) ---")
print(f"{'Topology':<20} {'Gen10':>10} {'Gen20':>10} {'Lambda2':>8}")
for topo in TOPOLOGIES:
    for gen in [10, 20]:
        vals = [get_gen_value(data[(topo, seed)], gen, "meanFitness") for seed in SEEDS]
        mf_at[gen][topo] = sum(vals) / len(vals)
    print(f"{topo:<20} {mf_at[10][topo]:>10.4f} {mf_at[20][topo]:>10.4f} {LAMBDA2[topo]:>8.2f}")

# ============================================================
# Step 3: Spearman correlations with lambda_2
# ============================================================
print("\n" + "=" * 80)
print("SPEARMAN RANK CORRELATIONS WITH LAMBDA_2")
print("=" * 80)

lambda2_vals = [LAMBDA2[t] for t in TOPOLOGIES]

def compute_spearman(metric_dict, metric_name, expected_sign):
    """Compute Spearman correlation and print results."""
    metric_vals = [metric_dict[t] for t in TOPOLOGIES]
    rho, pval = stats.spearmanr(lambda2_vals, metric_vals)
    sign = "POSITIVE" if rho > 0 else "NEGATIVE"
    match = "YES" if (expected_sign == "negative" and rho < 0) or (expected_sign == "positive" and rho > 0) else "NO"
    sig = "***" if pval < 0.01 else "**" if pval < 0.05 else "*" if pval < 0.1 else ""
    print(f"\n{metric_name}:")
    print(f"  rho = {rho:.4f}, p = {pval:.4f} {sig}")
    print(f"  Direction: {sign} (expected: {expected_sign.upper()}, match: {match})")

    # Also show the rank ordering
    ranked = sorted(TOPOLOGIES, key=lambda t: metric_dict[t])
    print(f"  Rank order (low to high): {' < '.join(ranked)}")
    return rho, pval

print("\nH1: Higher lambda_2 → faster convergence (NEGATIVE correlation with gen to reach best=1.0)")
rho1, p1 = compute_spearman(gen_best_1, "Gen to best fitness = 1.0", "negative")

print("\nH2: Higher lambda_2 → faster convergence (NEGATIVE correlation with gen to reach mean>=0.9)")
rho2, p2 = compute_spearman(gen_mean_09, "Gen to mean fitness >= 0.9", "negative")

print("\nH3: Higher lambda_2 → more information sharing → LOWER diversity at gen 10")
rho3, p3 = compute_spearman(div_at[10], "Diversity at gen 10", "negative")

print("\nH4: Higher lambda_2 → LOWER diversity at gen 20")
rho4, p4 = compute_spearman(div_at[20], "Diversity at gen 20", "negative")

print("\nH5: Higher lambda_2 → LOWER diversity at gen 30")
rho5, p5 = compute_spearman(div_at[30], "Diversity at gen 30", "negative")

print("\nH6: Higher lambda_2 → HIGHER mean fitness at gen 10")
rho6, p6 = compute_spearman(mf_at[10], "Mean fitness at gen 10", "positive")

print("\nH7: Higher lambda_2 → HIGHER mean fitness at gen 20")
rho7, p7 = compute_spearman(mf_at[20], "Mean fitness at gen 20", "positive")

# ============================================================
# Step 4: Additional analysis — per-seed consistency
# ============================================================
print("\n" + "=" * 80)
print("PER-SEED SPEARMAN CORRELATIONS (gen to best=1.0 vs lambda_2)")
print("=" * 80)
for seed in SEEDS:
    vals = [first_gen_reaching(data[(t, seed)], "bestFitness", 1.0) for t in TOPOLOGIES]
    rho, pval = stats.spearmanr(lambda2_vals, vals)
    print(f"  Seed {seed}: rho = {rho:.4f}, p = {pval:.4f}")

# ============================================================
# Step 5: Check if gen 0 and gen 10 are truly identical across topologies
# ============================================================
print("\n" + "=" * 80)
print("IDENTITY CHECK: Are gen 0 and gen 10 identical across topologies?")
print("=" * 80)
for seed in SEEDS:
    print(f"\n  Seed {seed}:")
    for gen in [0, 10]:
        mean_fits = set()
        best_fits = set()
        divs = set()
        for topo in TOPOLOGIES:
            rows = data[(topo, seed)]
            row = [r for r in rows if r["generation"] == gen][0]
            mean_fits.add(f"{row['meanFitness']:.10f}")
            best_fits.add(f"{row['bestFitness']:.10f}")
            divs.add(f"{row['diversity']:.10f}")
        all_same_mean = len(mean_fits) == 1
        all_same_best = len(best_fits) == 1
        all_same_div = len(divs) == 1
        print(f"    Gen {gen}: meanFit same={all_same_mean}, bestFit same={all_same_best}, div same={all_same_div}")
        if not all_same_mean:
            print(f"      meanFitness values: {mean_fits}")

# ============================================================
# Step 6: Detailed gen 20 comparison
# ============================================================
print("\n" + "=" * 80)
print("DETAILED GEN 20 COMPARISON (the critical divergence point)")
print("=" * 80)
print(f"{'Topology':<20} {'meanFit_42':>12} {'meanFit_137':>12} {'meanFit_2718':>12} {'bestFit_42':>10} {'bestFit_137':>10} {'bestFit_2718':>10}")
for topo in TOPOLOGIES:
    vals = []
    for seed in SEEDS:
        row = [r for r in data[(topo, seed)] if r["generation"] == 20][0]
        vals.append(row)
    print(f"{topo:<20} {vals[0]['meanFitness']:>12.6f} {vals[1]['meanFitness']:>12.6f} {vals[2]['meanFitness']:>12.6f} {vals[0]['bestFitness']:>10.2f} {vals[1]['bestFitness']:>10.2f} {vals[2]['bestFitness']:>10.2f}")

# ============================================================
# Step 7: Summary statistics
# ============================================================
print("\n" + "=" * 80)
print("SUMMARY TABLE: All Metrics by Topology")
print("=" * 80)
print(f"{'Topology':<18} {'λ₂':>6} {'GenBest1':>9} {'GenMean9':>9} {'Div@10':>8} {'Div@20':>8} {'Div@30':>8} {'Div@50':>8} {'MF@10':>8} {'MF@20':>8}")
for topo in TOPOLOGIES:
    print(f"{topo:<18} {LAMBDA2[topo]:>6.2f} {gen_best_1[topo]:>9.1f} {gen_mean_09[topo]:>9.1f} {div_at[10][topo]:>8.4f} {div_at[20][topo]:>8.4f} {div_at[30][topo]:>8.4f} {div_at[50][topo]:>8.4f} {mf_at[10][topo]:>8.4f} {mf_at[20][topo]:>8.4f}")

print("\n\nDONE.")
