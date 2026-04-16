#!/usr/bin/env python3
"""
Analyze interference experiment: does cycle ARRANGEMENT matter at constant beta_1?

3 directed topologies, all with n=8, |E|=9, beta_1=2, kappa=2:
  1. interference-adjacent  — two cycles sharing a vertex (figure-eight)
  2. interference-separated — two cycles connected by a path
  3. interference-nested    — one cycle inside another (concentric)

NK4 domain, 30 seeds each, 500 generations.

Key question: If the holonomy group is non-abelian (cactus J_n),
different cycle arrangements should produce different diversity dynamics,
even at identical beta_1. This would show up as significant eta-squared
between graph types.
"""

import os
import csv
import math
from collections import defaultdict

RESULTS_DIR = "results/interference"

TOPOLOGIES = [
    "interference-adjacent",
    "interference-separated",
    "interference-nested",
]

TOPO_LABELS = {
    "interference-adjacent": "Adjacent (figure-eight)",
    "interference-separated": "Separated (path-connected)",
    "interference-nested": "Nested (concentric)",
}

TOPO_COLORS = {
    "interference-adjacent": "#e74c3c",     # red
    "interference-separated": "#2ecc71",    # green
    "interference-nested": "#3498db",       # blue
}

# Report at every 10 generations (matching migration interval)
CHECKPOINTS = list(range(0, 510, 10))


def load_csv(path):
    """Load CSV, skip non-CSV header lines."""
    with open(path) as f:
        lines = [l for l in f if l.strip() and
                 (l.strip()[0].isdigit() or l.strip().startswith("generation"))]
    reader = csv.DictReader(lines)
    return list(reader)


def mean(xs):
    return sum(xs) / len(xs) if xs else 0.0


def std(xs):
    if len(xs) < 2:
        return 0.0
    m = mean(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))


def stderr(xs):
    return std(xs) / math.sqrt(len(xs)) if xs else 0.0


def eta_squared(groups):
    """One-way ANOVA eta-squared: fraction of variance explained by group."""
    all_vals = []
    for g in groups:
        all_vals.extend(g)
    if not all_vals:
        return 0.0
    grand_mean = mean(all_vals)
    ss_total = sum((x - grand_mean) ** 2 for x in all_vals)
    if ss_total == 0:
        return 0.0
    ss_between = sum(len(g) * (mean(g) - grand_mean) ** 2 for g in groups if g)
    return ss_between / ss_total


def f_statistic(groups):
    """One-way ANOVA F-statistic."""
    k = len(groups)
    n_total = sum(len(g) for g in groups)
    if k < 2 or n_total <= k:
        return 0.0

    grand_mean = mean([x for g in groups for x in g])
    ss_between = sum(len(g) * (mean(g) - grand_mean) ** 2 for g in groups if g)
    ss_within = sum(sum((x - mean(g)) ** 2 for x in g) for g in groups if g)

    df_between = k - 1
    df_within = n_total - k

    if ss_within == 0 or df_within == 0:
        return float('inf') if ss_between > 0 else 0.0

    ms_between = ss_between / df_between
    ms_within = ss_within / df_within
    return ms_between / ms_within


def cohens_d(xs, ys):
    """Cohen's d effect size between two groups."""
    if len(xs) < 2 or len(ys) < 2:
        return 0.0
    mx, my = mean(xs), mean(ys)
    sx, sy = std(xs), std(ys)
    pooled = math.sqrt(((len(xs) - 1) * sx**2 + (len(ys) - 1) * sy**2) /
                       (len(xs) + len(ys) - 2))
    if pooled == 0:
        return 0.0
    return (mx - my) / pooled


def main():
    # Load all data: topo -> gen -> [diversity values across seeds]
    data = defaultdict(lambda: defaultdict(list))
    fitness_data = defaultdict(lambda: defaultdict(list))
    file_counts = defaultdict(int)

    for topo in TOPOLOGIES:
        for fname in os.listdir(RESULTS_DIR):
            if not fname.startswith(topo + "_seed"):
                continue
            path = os.path.join(RESULTS_DIR, fname)
            rows = load_csv(path)
            file_counts[topo] += 1
            for row in rows:
                gen = int(row["generation"])
                div_val = float(row["diversity"])
                fit_val = float(row["bestFitness"])
                data[topo][gen].append(div_val)
                fitness_data[topo][gen].append(fit_val)

    print("=" * 70)
    print("INTERFERENCE EXPERIMENT — Cycle Arrangement at Constant beta_1")
    print("=" * 70)
    print()
    print("Graph properties (all identical):")
    print("  n=8, |E|=9, beta_1=2, kappa=2 (directed simple cycles)")
    print()
    for topo in TOPOLOGIES:
        print(f"  {TOPO_LABELS[topo]}: {file_counts[topo]} runs")
    print()

    # --- Table 1: Diversity at key checkpoints ---
    key_gens = [0, 30, 50, 100, 200, 300, 500]
    print("-" * 70)
    print("Table 1: Mean diversity (SE) at key generations")
    print("-" * 70)
    header = f"{'Gen':>5}"
    for topo in TOPOLOGIES:
        header += f"  {TOPO_LABELS[topo]:>22}"
    header += f"  {'eta^2':>8}  {'F':>8}"
    print(header)
    print("-" * 70)

    for gen in key_gens:
        row = f"{gen:>5}"
        groups = []
        for topo in TOPOLOGIES:
            vals = data[topo].get(gen, [])
            groups.append(vals)
            if vals:
                row += f"  {mean(vals):>17.4f} ({stderr(vals):.4f})"
            else:
                row += f"  {'N/A':>22}"
        eta2 = eta_squared(groups)
        f_val = f_statistic(groups)
        row += f"  {eta2:>8.4f}  {f_val:>8.2f}"
        print(row)
    print()

    # --- Table 2: Best fitness at key checkpoints ---
    print("-" * 70)
    print("Table 2: Mean best fitness (SE) at key generations")
    print("-" * 70)
    header = f"{'Gen':>5}"
    for topo in TOPOLOGIES:
        header += f"  {TOPO_LABELS[topo]:>22}"
    header += f"  {'eta^2':>8}  {'F':>8}"
    print(header)
    print("-" * 70)

    for gen in key_gens:
        row = f"{gen:>5}"
        groups = []
        for topo in TOPOLOGIES:
            vals = fitness_data[topo].get(gen, [])
            groups.append(vals)
            if vals:
                row += f"  {mean(vals):>17.4f} ({stderr(vals):.4f})"
            else:
                row += f"  {'N/A':>22}"
        eta2 = eta_squared(groups)
        f_val = f_statistic(groups)
        row += f"  {eta2:>8.4f}  {f_val:>8.2f}"
        print(row)
    print()

    # --- Eta-squared over time (every 10 gens) ---
    print("-" * 70)
    print("Eta-squared (diversity) over time — every 50 generations")
    print("-" * 70)
    eta2_diversity = {}
    eta2_fitness = {}
    for gen in CHECKPOINTS:
        groups_div = [data[topo].get(gen, []) for topo in TOPOLOGIES]
        groups_fit = [fitness_data[topo].get(gen, []) for topo in TOPOLOGIES]
        if all(len(g) > 0 for g in groups_div):
            eta2_diversity[gen] = eta_squared(groups_div)
            eta2_fitness[gen] = eta_squared(groups_fit)

    for gen in sorted(eta2_diversity.keys()):
        if gen % 50 == 0:
            print(f"  Gen {gen:>4}: eta^2(diversity) = {eta2_diversity[gen]:.4f}, "
                  f"eta^2(fitness) = {eta2_fitness[gen]:.4f}")
    print()

    # --- Pairwise Cohen's d at peak effect ---
    peak_gen = max(eta2_diversity, key=eta2_diversity.get)
    print(f"Peak eta^2(diversity) = {eta2_diversity[peak_gen]:.4f} at generation {peak_gen}")
    print()
    print(f"Pairwise Cohen's d at generation {peak_gen}:")
    for i, t1 in enumerate(TOPOLOGIES):
        for t2 in TOPOLOGIES[i+1:]:
            d = cohens_d(data[t1][peak_gen], data[t2][peak_gen])
            print(f"  {TOPO_LABELS[t1]} vs {TOPO_LABELS[t2]}: d = {d:.3f}")
    print()

    # --- Max eta-squared ---
    max_eta2 = max(eta2_diversity.values())
    max_gen = max(eta2_diversity, key=eta2_diversity.get)
    print(f"Maximum eta^2(diversity) = {max_eta2:.4f} at generation {max_gen}")

    min_eta2 = min(v for g, v in eta2_diversity.items() if g > 0)
    min_gen = min((g for g, v in eta2_diversity.items() if g > 0 and v == min_eta2))
    print(f"Minimum eta^2(diversity) = {min_eta2:.4f} at generation {min_gen}")
    print()

    # --- Interpretation ---
    print("=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    if max_eta2 > 0.14:
        print("LARGE EFFECT: Cycle arrangement matters significantly at constant beta_1.")
        print("This is consistent with non-abelian holonomy (cactus group J_n).")
        print("Different cycle arrangements create different diversity dynamics")
        print("even when beta_1, node count, edge count, and cycle count are identical.")
    elif max_eta2 > 0.06:
        print("MEDIUM EFFECT: Cycle arrangement has a moderate influence at constant beta_1.")
        print("Suggestive of non-abelian holonomy effects, but not overwhelming.")
    elif max_eta2 > 0.01:
        print("SMALL EFFECT: Cycle arrangement has a small but detectable influence.")
        print("Weak evidence for holonomy interference.")
    else:
        print("NEGLIGIBLE EFFECT: Cycle arrangement does NOT matter at constant beta_1.")
        print("This would be consistent with abelian holonomy (arrangement is irrelevant).")
    print()

    # --- Generate plot ---
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np

        fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

        # Top panel: diversity trajectories
        ax1 = axes[0]
        for topo in TOPOLOGIES:
            gens = sorted(data[topo].keys())
            means = [mean(data[topo][g]) for g in gens]
            sems = [stderr(data[topo][g]) for g in gens]
            color = TOPO_COLORS[topo]
            ax1.plot(gens, means, color=color, linewidth=2,
                     label=TOPO_LABELS[topo])
            ax1.fill_between(gens,
                             [m - s for m, s in zip(means, sems)],
                             [m + s for m, s in zip(means, sems)],
                             color=color, alpha=0.15)

        ax1.set_ylabel("Diversity (mean pairwise distance)", fontsize=12)
        ax1.set_title("Interference Experiment: Cycle Arrangement at Constant $\\beta_1 = 2$\n"
                       "(n=8, |E|=9, $\\kappa$=2, NK4 domain, 30 seeds each)",
                       fontsize=13)
        ax1.legend(fontsize=10, loc='upper right')
        ax1.grid(True, alpha=0.3)

        # Bottom panel: eta-squared over time
        ax2 = axes[1]
        gens_eta = sorted(eta2_diversity.keys())
        eta_vals = [eta2_diversity[g] for g in gens_eta]
        ax2.plot(gens_eta, eta_vals, color='#8e44ad', linewidth=2,
                 label='$\\eta^2$ (diversity)')

        eta_fit_vals = [eta2_fitness[g] for g in gens_eta]
        ax2.plot(gens_eta, eta_fit_vals, color='#e67e22', linewidth=2,
                 linestyle='--', label='$\\eta^2$ (best fitness)')

        # Reference lines
        ax2.axhline(y=0.14, color='gray', linestyle=':', alpha=0.5, label='Large effect (0.14)')
        ax2.axhline(y=0.06, color='gray', linestyle='--', alpha=0.3, label='Medium effect (0.06)')

        ax2.set_xlabel("Generation", fontsize=12)
        ax2.set_ylabel("$\\eta^2$ (effect size)", fontsize=12)
        ax2.legend(fontsize=10, loc='upper right')
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim(bottom=0)

        plt.tight_layout()
        plot_path = os.path.join(RESULTS_DIR, "interference_experiment.png")
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        print(f"Plot saved to: {plot_path}")

        # Also save a pairwise comparison plot
        fig2, axes2 = plt.subplots(1, 3, figsize=(15, 5))
        pairs = [
            ("interference-adjacent", "interference-separated"),
            ("interference-adjacent", "interference-nested"),
            ("interference-separated", "interference-nested"),
        ]
        for ax, (t1, t2) in zip(axes2, pairs):
            gens = sorted(set(data[t1].keys()) & set(data[t2].keys()))
            d_vals = [cohens_d(data[t1][g], data[t2][g]) for g in gens]
            ax.plot(gens, d_vals, linewidth=1.5)
            ax.axhline(y=0, color='gray', linestyle='-', alpha=0.3)
            ax.axhline(y=0.8, color='red', linestyle=':', alpha=0.3, label='Large (0.8)')
            ax.axhline(y=-0.8, color='red', linestyle=':', alpha=0.3)
            ax.set_title(f"{TOPO_LABELS[t1].split()[0]} vs {TOPO_LABELS[t2].split()[0]}",
                         fontsize=10)
            ax.set_xlabel("Generation")
            ax.set_ylabel("Cohen's d")
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=8)

        fig2.suptitle("Pairwise Cohen's d (diversity) over time", fontsize=13)
        plt.tight_layout()
        pairwise_path = os.path.join(RESULTS_DIR, "interference_pairwise.png")
        fig2.savefig(pairwise_path, dpi=150, bbox_inches='tight')
        print(f"Pairwise plot saved to: {pairwise_path}")

    except ImportError:
        print("matplotlib not available; skipping plot generation.")


if __name__ == "__main__":
    main()
