#!/usr/bin/env python3
"""Comprehensive analysis of Jaccard maze pilot experiment.

8 topologies x 3 seeds x 500 generations.
Metrics: meanFitness, bestFitness, diversity.
"""

import os
import glob
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from itertools import combinations

# ── Config ──────────────────────────────────────────────────────────────
DATA_DIR = "/home/lyra/projects/Topology-experiments/results/jaccard_maze_pilot"
FIG_DIR = os.path.join(DATA_DIR, "figures")
METRICS = ["meanFitness", "bestFitness", "diversity"]
KEY_GENS = [30, 50, 100, 200, 500]
SEEDS = [42, 137, 2025]

# Topology display order and colors
TOPO_ORDER = ["ring", "star", "complete", "hypercube", "barbell",
              "random-regular", "watts-strogatz", "disconnected"]
TOPO_COLORS = {
    "ring": "#1f77b4", "star": "#ff7f0e", "complete": "#2ca02c",
    "hypercube": "#d62728", "barbell": "#9467bd", "random-regular": "#8c564b",
    "watts-strogatz": "#e377c2", "disconnected": "#7f7f7f"
}
TOPO_SHORT = {
    "ring": "Ring", "star": "Star", "complete": "Complete",
    "hypercube": "Hypercube", "barbell": "Barbell",
    "random-regular": "Rand-Reg", "watts-strogatz": "Watts-Strog",
    "disconnected": "Disconn"
}

# ── Load data ───────────────────────────────────────────────────────────
def load_all():
    """Load all CSVs into a single DataFrame with topology and seed columns."""
    frames = []
    for path in sorted(glob.glob(os.path.join(DATA_DIR, "*.csv"))):
        fname = os.path.basename(path).replace(".csv", "")
        parts = fname.rsplit("_seed", 1)
        topo = parts[0]
        seed = int(parts[1])
        df = pd.read_csv(path, comment='U')
        # Skip any non-numeric header lines
        if 'generation' not in df.columns:
            df = pd.read_csv(path, skiprows=1)
        df["topology"] = topo
        df["seed"] = seed
        frames.append(df)
    return pd.concat(frames, ignore_index=True)

data = load_all()
topos = sorted(data["topology"].unique())
gens = sorted(data["generation"].unique())
print(f"Loaded {len(data)} rows: {len(topos)} topologies, {len(SEEDS)} seeds, {len(gens)} generations")
print(f"Topologies: {topos}")
print(f"Generation range: {min(gens)} - {max(gens)}")
print()

# ── Eta-squared over time ──────────────────────────────────────────────
def compute_eta_sq(data, metric, gen):
    """One-way ANOVA eta-squared for topology effect at a given generation."""
    subset = data[data["generation"] == gen]
    groups = [g[metric].values for _, g in subset.groupby("topology")]
    # Filter out groups with < 2 observations
    groups = [g for g in groups if len(g) >= 1]
    if len(groups) < 2:
        return np.nan, np.nan

    # Compute SS_between and SS_total
    all_vals = np.concatenate(groups)
    grand_mean = np.mean(all_vals)
    ss_total = np.sum((all_vals - grand_mean) ** 2)
    ss_between = sum(len(g) * (np.mean(g) - grand_mean) ** 2 for g in groups)

    if ss_total == 0:
        return 0.0, 1.0

    eta_sq = ss_between / ss_total

    # Also get ANOVA p-value
    try:
        f_stat, p_val = stats.f_oneway(*groups)
    except:
        p_val = np.nan

    return eta_sq, p_val

print("=" * 70)
print("ETA-SQUARED (η²) OVER TIME — Topology Effect Size")
print("=" * 70)

eta_results = {}
for metric in METRICS:
    eta_results[metric] = {"gens": [], "eta_sq": [], "p_val": []}
    for gen in gens:
        eta, p = compute_eta_sq(data, metric, gen)
        eta_results[metric]["gens"].append(gen)
        eta_results[metric]["eta_sq"].append(eta)
        eta_results[metric]["p_val"].append(p)

# Print at key generations
for gen in KEY_GENS:
    if gen in gens:
        print(f"\nGeneration {gen}:")
        for metric in METRICS:
            idx = gens.index(gen)
            eta = eta_results[metric]["eta_sq"][idx]
            p = eta_results[metric]["p_val"][idx]
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            print(f"  {metric:15s}: η² = {eta:.4f}, p = {p:.4e} {sig}")

# Peak eta-squared
print(f"\n{'=' * 70}")
print("PEAK η² VALUES")
print("=" * 70)
for metric in METRICS:
    arr = np.array(eta_results[metric]["eta_sq"])
    peak_idx = np.argmax(arr)
    peak_gen = eta_results[metric]["gens"][peak_idx]
    peak_eta = arr[peak_idx]
    peak_p = eta_results[metric]["p_val"][peak_idx]
    print(f"  {metric:15s}: peak η² = {peak_eta:.4f} at gen {peak_gen} (p = {peak_p:.4e})")

# ── Statistical tests at key generations ────────────────────────────────
print(f"\n{'=' * 70}")
print("STATISTICAL TESTS AT KEY GENERATIONS")
print("=" * 70)

significant_results = []

for gen in KEY_GENS:
    if gen not in gens:
        continue
    print(f"\n--- Generation {gen} ---")
    subset = data[data["generation"] == gen]

    for metric in METRICS:
        groups = [g[metric].values for name, g in subset.groupby("topology")]
        group_names = [name for name, g in subset.groupby("topology")]

        # ANOVA
        try:
            f_stat, p_anova = stats.f_oneway(*groups)
        except:
            f_stat, p_anova = np.nan, np.nan

        # Kruskal-Wallis
        try:
            h_stat, p_kw = stats.kruskal(*groups)
        except:
            h_stat, p_kw = np.nan, np.nan

        print(f"\n  {metric}:")
        print(f"    ANOVA:          F = {f_stat:.4f}, p = {p_anova:.4e}")
        print(f"    Kruskal-Wallis: H = {h_stat:.4f}, p = {p_kw:.4e}")

        if p_anova < 0.05:
            significant_results.append((gen, metric, p_anova))
            # Post-hoc pairwise t-tests (Bonferroni)
            n_comparisons = len(group_names) * (len(group_names) - 1) // 2
            print(f"    Post-hoc pairwise t-tests (Bonferroni, {n_comparisons} comparisons):")
            sig_pairs = []
            for (i, name_i), (j, name_j) in combinations(enumerate(group_names), 2):
                t_stat, p_pair = stats.ttest_ind(groups[i], groups[j])
                p_adj = min(p_pair * n_comparisons, 1.0)
                if p_adj < 0.05:
                    sig_pairs.append((name_i, name_j, p_adj))
                    print(f"      {TOPO_SHORT.get(name_i, name_i):10s} vs {TOPO_SHORT.get(name_j, name_j):10s}: p_adj = {p_adj:.4e}")
            if not sig_pairs:
                print(f"      No pairwise comparisons significant after Bonferroni correction")

# ── Topology rankings ──────────────────────────────────────────────────
print(f"\n{'=' * 70}")
print("TOPOLOGY RANKINGS")
print("=" * 70)

for metric in METRICS:
    print(f"\n  {metric} (mean over all gens and seeds):")
    ranking = data.groupby("topology")[metric].mean().sort_values(ascending=False)
    for i, (topo, val) in enumerate(ranking.items()):
        print(f"    {i+1}. {TOPO_SHORT.get(topo, topo):12s}: {val:.6f}")

# Rankings at final generation
print(f"\n  Rankings at generation {max(gens)}:")
final = data[data["generation"] == max(gens)]
for metric in METRICS:
    print(f"\n  {metric}:")
    ranking = final.groupby("topology")[metric].mean().sort_values(ascending=False)
    for i, (topo, val) in enumerate(ranking.items()):
        std = final[final["topology"] == topo][metric].std()
        print(f"    {i+1}. {TOPO_SHORT.get(topo, topo):12s}: {val:.6f} (±{std:.6f})")

# ── Comparison with OneMax ──────────────────────────────────────────────
print(f"\n{'=' * 70}")
print("COMPARISON WITH ONEMAX")
print("=" * 70)
print("\nOneMax reference (from directed cycle experiment):")
print("  η² at gen 30: fitness ~0.88, diversity ~0.76 (transient, peaks early)")
print()
for metric in METRICS:
    if 30 in gens:
        idx30 = gens.index(30)
        eta30 = eta_results[metric]["eta_sq"][idx30]
    else:
        eta30 = "N/A"
    arr = np.array(eta_results[metric]["eta_sq"])
    peak = np.max(arr)
    print(f"  Maze {metric:15s}: η² at gen 30 = {eta30:.4f}, peak η² = {peak:.4f}")

# ── Figure A: Fitness trajectories ──────────────────────────────────────
print(f"\nGenerating figures...")

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.suptitle("Jaccard Maze Pilot — Trajectories by Topology", fontsize=14, fontweight='bold')

for ax, metric in zip(axes, METRICS):
    for topo in TOPO_ORDER:
        topo_data = data[data["topology"] == topo]
        means = topo_data.groupby("generation")[metric].mean()
        stds = topo_data.groupby("generation")[metric].std()

        ax.plot(means.index, means.values, label=TOPO_SHORT[topo],
                color=TOPO_COLORS[topo], linewidth=1.5)
        ax.fill_between(means.index, means.values - stds.values,
                       means.values + stds.values,
                       color=TOPO_COLORS[topo], alpha=0.15)

    ax.set_xlabel("Generation", fontsize=11)
    ax.set_ylabel(metric, fontsize=11)
    ax.set_title(metric, fontsize=12)
    ax.grid(True, alpha=0.3)

# Single legend outside
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=8, fontsize=9,
           bbox_to_anchor=(0.5, -0.02))
plt.tight_layout(rect=[0, 0.05, 1, 0.95])

pdf_path_a = os.path.join(FIG_DIR, "fig_maze_trajectories.pdf")
plt.savefig(pdf_path_a, bbox_inches='tight', dpi=300)
plt.savefig(pdf_path_a.replace('.pdf', '.png'), bbox_inches='tight', dpi=150)
plt.close()
print(f"  Saved: {pdf_path_a}")

# ── Figure B: η² over time ─────────────────────────────────────────────
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
fig.suptitle("Jaccard Maze Pilot — Topology Effect Size (η²) Over Generations",
             fontsize=14, fontweight='bold')

metric_styles = {
    "meanFitness": ("-", "#2ca02c", "Mean Fitness"),
    "bestFitness": ("--", "#1f77b4", "Best Fitness"),
    "diversity": ("-.", "#d62728", "Diversity")
}

for metric in METRICS:
    style, color, label = metric_styles[metric]
    ax.plot(eta_results[metric]["gens"], eta_results[metric]["eta_sq"],
            linestyle=style, color=color, linewidth=2, label=label)

# Reference lines
ax.axhline(y=0.14, color='gray', linestyle=':', alpha=0.5, label="Large effect (0.14)")
ax.axhline(y=0.06, color='gray', linestyle=':', alpha=0.3, label="Medium effect (0.06)")

# Mark significant generations
for metric in METRICS:
    style, color, label = metric_styles[metric]
    for i, (gen, eta, p) in enumerate(zip(
            eta_results[metric]["gens"],
            eta_results[metric]["eta_sq"],
            eta_results[metric]["p_val"])):
        if p < 0.05:
            ax.scatter([gen], [eta], color=color, s=20, zorder=5, marker='o')

ax.set_xlabel("Generation", fontsize=12)
ax.set_ylabel("η² (ANOVA effect size)", fontsize=12)
ax.legend(fontsize=10, loc='upper right')
ax.grid(True, alpha=0.3)
ax.set_ylim(bottom=0)

plt.tight_layout()
pdf_path_b = os.path.join(FIG_DIR, "fig_maze_eta_squared.pdf")
plt.savefig(pdf_path_b, bbox_inches='tight', dpi=300)
plt.savefig(pdf_path_b.replace('.pdf', '.png'), bbox_inches='tight', dpi=150)
plt.close()
print(f"  Saved: {pdf_path_b}")

# ── Summary statistics for report ───────────────────────────────────────
print(f"\n{'=' * 70}")
print("SUMMARY FOR REPORT")
print("=" * 70)

any_sig = len(significant_results) > 0
print(f"\nSignificant ANOVA results (p < 0.05): {len(significant_results)}")
for gen, metric, p in significant_results:
    print(f"  gen {gen}, {metric}: p = {p:.4e}")

# Overall spread at final generation
final = data[data["generation"] == max(gens)]
for metric in METRICS:
    vals = final.groupby("topology")[metric].mean()
    spread = vals.max() - vals.min()
    mean_val = vals.mean()
    pct_spread = 100 * spread / mean_val if mean_val > 0 else 0
    print(f"\n  {metric} spread at gen {max(gens)}: {spread:.6f} ({pct_spread:.2f}% of mean)")
    print(f"    Best:  {vals.idxmax()} = {vals.max():.6f}")
    print(f"    Worst: {vals.idxmin()} = {vals.min():.6f}")

print("\nDone.")
