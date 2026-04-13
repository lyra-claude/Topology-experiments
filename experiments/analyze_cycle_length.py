#!/usr/bin/env python3
"""
Cycle Length Experiment Analysis: Does mean cycle length predict diversity
independently of beta_1?

Three two-cycle-bridge graphs, all beta_1=2, varying mean cycle length:
  G1: cycle-bridge-3-3  (mean_cycle=3.0,  n=6,  lambda_2=0.4384)
  G2: cycle-bridge-3-9  (mean_cycle=6.0,  n=12, lambda_2=0.1907)
  G3: cycle-bridge-9-9  (mean_cycle=9.0,  n=18, lambda_2=0.0822)

Analyses:
  1. Diversity at gen 10, 30, 100, 500 for each topology
  2. Mann-Whitney U test: G1 vs G3 at gen 30
  3. Cohen's d effect size
  4. Regression: diversity ~ mean_cycle_length + lambda_2
  5. One-way ANOVA (eta^2) across all three topologies
  6. Figures: diversity trajectory (mean +/- SD)

CAVEAT: These graphs have different node counts (6, 12, 18) and different
lambda_2 values. Any diversity effect could be driven by:
  (a) mean cycle length (the hypothesis),
  (b) lambda_2 / algebraic connectivity, or
  (c) island count / population size.
The regression in Analysis 4 attempts to tease apart (a) and (b).
"""

import os
import re
import numpy as np
import pandas as pd
from scipy import stats as sp_stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
RESULTS_DIR = Path(__file__).parent.parent / "results" / "cycle_length"
TIMEPOINTS = [10, 30, 100, 500]

# Topology metadata: (mean_cycle_length, beta_1, lambda_2, n_islands)
TOPO_META = {
    "cycle-bridge-3-3": (3.0,  2, 0.4384, 6),
    "cycle-bridge-3-9": (6.0,  2, 0.1907, 12),
    "cycle-bridge-9-9": (9.0,  2, 0.0822, 18),
}

TOPO_SHORT = {
    "cycle-bridge-3-3": "G1 (3+3)",
    "cycle-bridge-3-9": "G2 (3+9)",
    "cycle-bridge-9-9": "G3 (9+9)",
}

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

def load_all_results():
    """Load all CSV files into a single DataFrame with metadata."""
    rows = []
    pattern = re.compile(r"(nk\d+)_(.+)_seed(\d+)\.csv")

    for f in sorted(RESULTS_DIR.glob("*.csv")):
        m = pattern.match(f.name)
        if not m:
            continue
        domain = m.group(1)
        topo = m.group(2)
        seed = int(m.group(3))

        if topo not in TOPO_META:
            continue

        mean_cycle, beta1, lambda2, n_islands = TOPO_META[topo]

        df = pd.read_csv(f)
        df["domain"] = domain
        df["topology"] = topo
        df["mean_cycle"] = mean_cycle
        df["beta1"] = beta1
        df["lambda2"] = lambda2
        df["n_islands"] = n_islands
        df["seed"] = seed
        rows.append(df)

    if not rows:
        raise FileNotFoundError(f"No CSV files found in {RESULTS_DIR}")

    return pd.concat(rows, ignore_index=True)


# ---------------------------------------------------------------------------
# Statistical helpers
# ---------------------------------------------------------------------------

def eta_squared_groups(*groups):
    """Compute eta-squared for k groups (one-way ANOVA effect size)."""
    all_vals = np.concatenate(groups)
    grand_mean = np.mean(all_vals)
    ss_between = sum(len(g) * (np.mean(g) - grand_mean)**2 for g in groups)
    ss_total = np.sum((all_vals - grand_mean)**2)
    if ss_total == 0:
        return 0.0
    return ss_between / ss_total


def cohens_d(g1, g2):
    """Cohen's d for two groups (pooled SD)."""
    n1, n2 = len(g1), len(g2)
    m1, m2 = np.mean(g1), np.mean(g2)
    s1, s2 = np.std(g1, ddof=1), np.std(g2, ddof=1)
    pooled_sd = np.sqrt(((n1-1)*s1**2 + (n2-1)*s2**2) / (n1+n2-2))
    if pooled_sd == 0:
        return 0.0
    return (m1 - m2) / pooled_sd


# ---------------------------------------------------------------------------
# Analysis 1: Diversity at timepoints
# ---------------------------------------------------------------------------

def diversity_summary(data):
    """Diversity at key timepoints for each topology."""
    print("=" * 80)
    print("ANALYSIS 1: Diversity at key timepoints")
    print("=" * 80)

    for domain in sorted(data["domain"].unique()):
        print(f"\n  Domain: {domain}")
        print(f"  {'Gen':>5} | ", end="")
        for topo in TOPO_META:
            print(f"  {TOPO_SHORT[topo]:>12}", end="")
        print(f"  | {'eta^2':>8} | {'F':>8} | {'p':>10}")
        print(f"  {'─'*5}-+-{'─'*42}-+-{'─'*8}-+-{'─'*8}-+-{'─'*10}")

        for gen in TIMEPOINTS:
            snapshot = data[(data["domain"] == domain) & (data["generation"] == gen)]
            print(f"  {gen:>5} | ", end="")

            groups = []
            for topo in TOPO_META:
                vals = snapshot[snapshot["topology"] == topo]["diversity"].values
                groups.append(vals)
                m = np.mean(vals) if len(vals) > 0 else float('nan')
                s = np.std(vals) if len(vals) > 0 else float('nan')
                print(f"  {m:.4f}±{s:.3f}", end="")

            if len(groups) >= 2 and all(len(g) > 0 for g in groups):
                eta2 = eta_squared_groups(*groups)
                f_stat, p_val = sp_stats.f_oneway(*groups)
                sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
                print(f"  | {eta2:>8.4f} | {f_stat:>8.2f} | {p_val:>9.2e} {sig}")
            else:
                print(f"  | {'N/A':>8} | {'N/A':>8} | {'N/A':>10}")


# ---------------------------------------------------------------------------
# Analysis 2: Mann-Whitney U test (G1 vs G3 at gen 30)
# ---------------------------------------------------------------------------

def mannwhitney_test(data):
    """Mann-Whitney U test: G1 vs G3 at gen 30."""
    print("\n\n" + "=" * 80)
    print("ANALYSIS 2: Mann-Whitney U (G1 vs G3) at key generations")
    print("=" * 80)

    for domain in sorted(data["domain"].unique()):
        print(f"\n  Domain: {domain}")
        print(f"  {'Gen':>5} | {'G1 mean':>10} | {'G3 mean':>10} | {'U':>10} | {'p':>10} | {'Cohen d':>10}")
        print(f"  {'─'*5}-+-{'─'*10}-+-{'─'*10}-+-{'─'*10}-+-{'─'*10}-+-{'─'*10}")

        for gen in TIMEPOINTS:
            snapshot = data[(data["domain"] == domain) & (data["generation"] == gen)]
            g1_vals = snapshot[snapshot["topology"] == "cycle-bridge-3-3"]["diversity"].values
            g3_vals = snapshot[snapshot["topology"] == "cycle-bridge-9-9"]["diversity"].values

            if len(g1_vals) == 0 or len(g3_vals) == 0:
                continue

            u_stat, p_val = sp_stats.mannwhitneyu(g1_vals, g3_vals, alternative='two-sided')
            d = cohens_d(g1_vals, g3_vals)
            sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"

            print(f"  {gen:>5} | {np.mean(g1_vals):>10.4f} | {np.mean(g3_vals):>10.4f} | {u_stat:>10.1f} | {p_val:>9.2e} {sig} | {d:>10.3f}")


# ---------------------------------------------------------------------------
# Analysis 3: Regression
# ---------------------------------------------------------------------------

def regression_analysis(data):
    """Regression: diversity ~ mean_cycle_length + lambda_2."""
    print("\n\n" + "=" * 80)
    print("ANALYSIS 3: Regression (diversity ~ mean_cycle + lambda_2)")
    print("=" * 80)
    print("  CAVEAT: Only 3 topology levels — treat as exploratory, not confirmatory.")

    for domain in sorted(data["domain"].unique()):
        print(f"\n  Domain: {domain}")

        for gen in TIMEPOINTS:
            snapshot = data[(data["domain"] == domain) & (data["generation"] == gen)]
            if snapshot.empty:
                continue

            y = snapshot["diversity"].values
            X_cycle = snapshot["mean_cycle"].values
            X_lambda = snapshot["lambda2"].values

            # Simple correlations
            r_cycle, p_cycle = sp_stats.pearsonr(X_cycle, y)
            r_lambda, p_lambda = sp_stats.pearsonr(X_lambda, y)

            print(f"\n  Gen {gen}:")
            print(f"    Pearson r(diversity, mean_cycle): {r_cycle:+.4f} (p={p_cycle:.2e})")
            print(f"    Pearson r(diversity, lambda_2):   {r_lambda:+.4f} (p={p_lambda:.2e})")

            # Multiple regression via OLS (manual — avoid statsmodels dependency)
            X = np.column_stack([np.ones_like(X_cycle), X_cycle, X_lambda])
            try:
                beta = np.linalg.lstsq(X, y, rcond=None)[0]
                y_hat = X @ beta
                ss_res = np.sum((y - y_hat)**2)
                ss_tot = np.sum((y - np.mean(y))**2)
                r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
                print(f"    OLS: diversity = {beta[0]:.4f} + {beta[1]:.4f}*mean_cycle + {beta[2]:.4f}*lambda_2  (R²={r_squared:.4f})")
            except np.linalg.LinAlgError:
                print("    OLS: singular matrix, cannot compute")


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def make_figures(data):
    """Publication-quality figures."""
    figdir = RESULTS_DIR / "figures"
    figdir.mkdir(exist_ok=True)

    colors = {
        "cycle-bridge-3-3": "#2196F3",  # blue
        "cycle-bridge-3-9": "#9C27B0",  # purple
        "cycle-bridge-9-9": "#FF5722",  # red-orange
    }

    # Figure 1: Diversity trajectories (mean +/- SD) per domain
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    for ax, domain in zip(axes, ["nk0", "nk4"]):
        dom_data = data[data["domain"] == domain]

        for topo in TOPO_META:
            topo_data = dom_data[dom_data["topology"] == topo]
            if topo_data.empty:
                continue

            grouped = topo_data.groupby("generation")["diversity"]
            mean_div = grouped.mean()
            std_div = grouped.std()

            color = colors[topo]
            label = TOPO_SHORT[topo]
            mcl, _, lam2, n = TOPO_META[topo]
            label_full = f"{label}, MCL={mcl:.0f}, λ₂={lam2:.2f}, n={n}"

            ax.plot(mean_div.index, mean_div.values, color=color,
                    linewidth=2, label=label_full)
            ax.fill_between(mean_div.index,
                           (mean_div - std_div).values,
                           (mean_div + std_div).values,
                           color=color, alpha=0.15)

        ax.set_xlabel("Generation", fontsize=11)
        ax.set_title(f"Domain: {domain.upper()}", fontsize=12, fontweight='bold')
        ax.legend(fontsize=8, loc='best')
        ax.grid(True, alpha=0.3)

    axes[0].set_ylabel("Diversity", fontsize=11)
    fig.suptitle("Cycle Length Experiment: Diversity Trajectories (all β₁=2)",
                 fontsize=13, fontweight='bold')
    plt.tight_layout()

    for ext in ["png", "pdf"]:
        fig.savefig(figdir / f"cycle_length_diversity_trajectories.{ext}",
                   dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {figdir / 'cycle_length_diversity_trajectories.png'}")

    # Figure 2: Diversity at gen 30 boxplot
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    for ax, domain in zip(axes, ["nk0", "nk4"]):
        snapshot = data[(data["domain"] == domain) & (data["generation"] == 30)]
        box_data = []
        labels = []
        for topo in TOPO_META:
            vals = snapshot[snapshot["topology"] == topo]["diversity"].values
            box_data.append(vals)
            labels.append(TOPO_SHORT[topo])

        bp = ax.boxplot(box_data, labels=labels, patch_artist=True)
        for patch, topo in zip(bp['boxes'], TOPO_META):
            patch.set_facecolor(colors[topo])
            patch.set_alpha(0.6)

        ax.set_ylabel("Diversity at Gen 30", fontsize=11)
        ax.set_title(f"Domain: {domain.upper()}", fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')

    fig.suptitle("Cycle Length Experiment: Diversity at Generation 30 (all β₁=2)",
                 fontsize=13, fontweight='bold')
    plt.tight_layout()

    for ext in ["png", "pdf"]:
        fig.savefig(figdir / f"cycle_length_boxplot_gen30.{ext}",
                   dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {figdir / 'cycle_length_boxplot_gen30.png'}")


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def print_summary(data):
    """Print the key verdict."""
    print("\n\n" + "=" * 80)
    print("VERDICT")
    print("=" * 80)

    print("\n  Design confound: mean_cycle correlates with lambda_2 and n_islands.")
    print("  All three graphs have beta_1=2. The question: does cycle length add")
    print("  predictive power beyond what lambda_2 already explains?")

    for domain in sorted(data["domain"].unique()):
        print(f"\n  Domain: {domain}")
        for gen in [30, 100, 500]:
            snapshot = data[(data["domain"] == domain) & (data["generation"] == gen)]
            groups = [snapshot[snapshot["topology"] == t]["diversity"].values for t in TOPO_META]
            if all(len(g) > 0 for g in groups):
                eta2 = eta_squared_groups(*groups)
                d = cohens_d(groups[0], groups[2])
                print(f"    Gen {gen}: eta^2={eta2:.4f}, Cohen's d(G1 vs G3)={d:.3f}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    data = load_all_results()
    n_runs = len(data.groupby(["domain", "topology", "seed"]))
    print(f"Loaded {len(data)} rows from {n_runs} runs\n")

    diversity_summary(data)
    mannwhitney_test(data)
    regression_analysis(data)
    print_summary(data)

    print("\n\nGenerating figures...")
    make_figures(data)

    print(f"\nDone. Results in {RESULTS_DIR}/")
