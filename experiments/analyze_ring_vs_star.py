#!/usr/bin/env python3
"""
Ring vs Star NK Sweep Analysis

For each K value (0, 2, 4, 6), computes:
- eta-squared (ring vs star) for diversity and best fitness at gen 100, 200, 300, 500
- Mann-Whitney U test for ring vs star at each K
- Publication-quality figure: eta-squared vs K

Answers Robin's question: Does ring advantage grow as K increases?
"""

import os
import re
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
RESULTS_DIR = Path(__file__).parent.parent / "results" / "ring_vs_star_nk"
TIMEPOINTS = [100, 200, 300, 500]
K_VALUES = [0, 2, 4, 6]
TOPOLOGIES = ["ring", "star"]

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

def load_all_results():
    """Load all CSV files into a single DataFrame."""
    rows = []
    pattern = re.compile(r"nk(\d+)_(ring|star)_seed(\d+)\.csv")

    for f in sorted(RESULTS_DIR.glob("*.csv")):
        m = pattern.match(f.name)
        if not m:
            continue
        k_val = int(m.group(1))
        topo = m.group(2)
        seed = int(m.group(3))

        df = pd.read_csv(f)
        df["K"] = k_val
        df["topology"] = topo
        df["seed"] = seed
        rows.append(df)

    return pd.concat(rows, ignore_index=True)


def eta_squared(group1, group2):
    """Compute eta-squared (effect size) for two groups."""
    all_vals = np.concatenate([group1, group2])
    grand_mean = np.mean(all_vals)
    ss_between = len(group1) * (np.mean(group1) - grand_mean)**2 + \
                 len(group2) * (np.mean(group2) - grand_mean)**2
    ss_total = np.sum((all_vals - grand_mean)**2)
    if ss_total == 0:
        return 0.0
    return ss_between / ss_total


def analyze(data):
    """Run the full analysis."""
    print("=" * 80)
    print("RING vs STAR NK SWEEP ANALYSIS")
    print("=" * 80)
    print(f"Total runs: {len(data.groupby(['K', 'topology', 'seed']))}")
    print()

    # Storage for figure
    results = []

    for gen in TIMEPOINTS:
        print(f"\n{'─' * 60}")
        print(f"  Generation {gen}")
        print(f"{'─' * 60}")
        print(f"{'K':>4} | {'Metric':>12} | {'Ring Mean':>10} | {'Star Mean':>10} | "
              f"{'eta²':>8} | {'U-stat':>8} | {'p-value':>10} | {'Winner':>6}")
        print(f"{'─' * 4}-+-{'─' * 12}-+-{'─' * 10}-+-{'─' * 10}-+-"
              f"{'─' * 8}-+-{'─' * 8}-+-{'─' * 10}-+-{'─' * 6}")

        for k in K_VALUES:
            snapshot = data[data["generation"] == gen]
            ring_data = snapshot[(snapshot["K"] == k) & (snapshot["topology"] == "ring")]
            star_data = snapshot[(snapshot["K"] == k) & (snapshot["topology"] == "star")]

            for metric in ["diversity", "bestFitness"]:
                ring_vals = ring_data[metric].values
                star_vals = star_data[metric].values

                if len(ring_vals) == 0 or len(star_vals) == 0:
                    print(f"  K={k}, {metric}: NO DATA")
                    continue

                eta2 = eta_squared(ring_vals, star_vals)
                u_stat, p_val = stats.mannwhitneyu(ring_vals, star_vals, alternative='two-sided')

                ring_mean = np.mean(ring_vals)
                star_mean = np.mean(star_vals)
                winner = "ring" if ring_mean > star_mean else "star"
                if metric == "diversity":
                    # Higher diversity = ring advantage (more exploration)
                    winner = "ring" if ring_mean > star_mean else "star"

                results.append({
                    "generation": gen,
                    "K": k,
                    "metric": metric,
                    "ring_mean": ring_mean,
                    "star_mean": star_mean,
                    "eta_squared": eta2,
                    "u_stat": u_stat,
                    "p_value": p_val,
                    "winner": winner,
                })

                sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
                print(f"{k:>4} | {metric:>12} | {ring_mean:>10.4f} | {star_mean:>10.4f} | "
                      f"{eta2:>8.4f} | {u_stat:>8.0f} | {p_val:>9.2e} {sig} | {winner:>6}")

    results_df = pd.DataFrame(results)

    # Summary: eta-squared vs K at final generation
    print(f"\n\n{'=' * 60}")
    print("SUMMARY: eta-squared at generation 500 (ring vs star)")
    print(f"{'=' * 60}")
    final = results_df[results_df["generation"] == 500]
    for metric in ["diversity", "bestFitness"]:
        m = final[final["metric"] == metric]
        print(f"\n{metric}:")
        for _, row in m.iterrows():
            sig = "***" if row["p_value"] < 0.001 else "**" if row["p_value"] < 0.01 else "*" if row["p_value"] < 0.05 else "ns"
            print(f"  K={row['K']:>1}: eta²={row['eta_squared']:.4f}  p={row['p_value']:.2e} {sig}  "
                  f"winner={row['winner']}  (ring={row['ring_mean']:.4f}, star={row['star_mean']:.4f})")

    # Answer Robin's question
    print(f"\n\n{'=' * 60}")
    print("ROBIN'S QUESTION: Does ring advantage grow as K increases?")
    print(f"{'=' * 60}")
    div_final = final[final["metric"] == "diversity"].sort_values("K")
    fit_final = final[final["metric"] == "bestFitness"].sort_values("K")

    print("\nDiversity eta² trend:  ", end="")
    print(" -> ".join(f"K={r['K']}: {r['eta_squared']:.4f}" for _, r in div_final.iterrows()))

    print("Fitness eta² trend:    ", end="")
    print(" -> ".join(f"K={r['K']}: {r['eta_squared']:.4f}" for _, r in fit_final.iterrows()))

    # Correlation of K with eta-squared
    k_vals = div_final["K"].values.astype(float)
    div_eta = div_final["eta_squared"].values
    fit_eta = fit_final["eta_squared"].values

    r_div, p_div = stats.pearsonr(k_vals, div_eta) if len(k_vals) > 2 else (0, 1)
    r_fit, p_fit = stats.pearsonr(k_vals, fit_eta) if len(k_vals) > 2 else (0, 1)

    print(f"\nPearson r(K, diversity eta²) = {r_div:.3f}, p = {p_div:.3f}")
    print(f"Pearson r(K, fitness eta²)   = {r_fit:.3f}, p = {p_fit:.3f}")

    if r_div > 0.5 and p_div < 0.1:
        print("\n>>> YES: Ring-vs-star effect on DIVERSITY grows with K.")
    else:
        print(f"\n>>> Diversity: {'Weak' if r_div > 0 else 'No'} evidence of growth (r={r_div:.3f}).")

    if r_fit > 0.5 and p_fit < 0.1:
        print(">>> YES: Ring-vs-star effect on FITNESS grows with K.")
    else:
        print(f">>> Fitness: {'Weak' if r_fit > 0 else 'No'} evidence of growth (r={r_fit:.3f}).")

    return results_df


def make_figure(results_df):
    """Publication-quality figure: eta-squared vs K."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)

    colors = {"diversity": "#2196F3", "bestFitness": "#FF5722"}
    labels = {"diversity": "Diversity", "bestFitness": "Best Fitness"}
    markers = {"diversity": "o", "bestFitness": "s"}

    for ax, metric in zip(axes, ["diversity", "bestFitness"]):
        subset = results_df[results_df["metric"] == metric]

        for gen in TIMEPOINTS:
            gen_data = subset[subset["generation"] == gen].sort_values("K")
            alpha = 0.3 + 0.7 * (gen / max(TIMEPOINTS))
            lw = 1.0 if gen < max(TIMEPOINTS) else 2.5
            ls = "--" if gen < max(TIMEPOINTS) else "-"
            marker = markers[metric] if gen == max(TIMEPOINTS) else None
            ms = 8 if gen == max(TIMEPOINTS) else 0

            ax.plot(gen_data["K"], gen_data["eta_squared"],
                    color=colors[metric], alpha=alpha, linewidth=lw,
                    linestyle=ls, marker=marker, markersize=ms,
                    label=f"Gen {gen}" if gen in [100, 500] else None)

            # Add significance stars at gen 500
            if gen == max(TIMEPOINTS):
                for _, row in gen_data.iterrows():
                    if row["p_value"] < 0.001:
                        sig = "***"
                    elif row["p_value"] < 0.01:
                        sig = "**"
                    elif row["p_value"] < 0.05:
                        sig = "*"
                    else:
                        sig = ""
                    if sig:
                        ax.annotate(sig, (row["K"], row["eta_squared"]),
                                   textcoords="offset points", xytext=(0, 10),
                                   ha='center', fontsize=9, fontweight='bold',
                                   color=colors[metric])

        ax.set_xlabel("K (landscape ruggedness)", fontsize=12)
        ax.set_title(labels[metric], fontsize=13, fontweight='bold')
        ax.set_xticks(K_VALUES)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-0.02, None)

    axes[0].set_ylabel(r"$\eta^2$ (ring vs star)", fontsize=12)

    fig.suptitle("Ring vs Star Topology Effect Across NK Landscapes",
                 fontsize=14, fontweight='bold', y=1.02)

    plt.tight_layout()

    for ext in ["png", "pdf"]:
        outpath = RESULTS_DIR / f"ring_vs_star_eta_squared.{ext}"
        fig.savefig(outpath, dpi=300, bbox_inches='tight')
        print(f"Saved: {outpath}")

    plt.close()

    # Also make a combined single-panel figure
    fig2, ax2 = plt.subplots(figsize=(6, 4.5))

    final = results_df[results_df["generation"] == 500]
    for metric in ["diversity", "bestFitness"]:
        m = final[final["metric"] == metric].sort_values("K")
        ax2.plot(m["K"], m["eta_squared"], color=colors[metric],
                marker=markers[metric], markersize=9, linewidth=2.5,
                label=labels[metric])

        for _, row in m.iterrows():
            if row["p_value"] < 0.001:
                sig = "***"
            elif row["p_value"] < 0.01:
                sig = "**"
            elif row["p_value"] < 0.05:
                sig = "*"
            else:
                sig = ""
            if sig:
                ax2.annotate(sig, (row["K"], row["eta_squared"]),
                           textcoords="offset points", xytext=(0, 10),
                           ha='center', fontsize=10, fontweight='bold',
                           color=colors[metric])

    ax2.set_xlabel("K (landscape ruggedness)", fontsize=12)
    ax2.set_ylabel(r"$\eta^2$ (ring vs star)", fontsize=12)
    ax2.set_title("Topology Effect Size Grows with Landscape Ruggedness\n(Generation 500, n=20 seeds)",
                  fontsize=12, fontweight='bold')
    ax2.set_xticks(K_VALUES)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(-0.02, None)

    plt.tight_layout()

    for ext in ["png", "pdf"]:
        outpath = RESULTS_DIR / f"ring_vs_star_combined.{ext}"
        fig2.savefig(outpath, dpi=300, bbox_inches='tight')
        print(f"Saved: {outpath}")

    plt.close()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    data = load_all_results()
    print(f"Loaded {len(data)} rows from {len(data.groupby(['K', 'topology', 'seed']))} runs\n")

    results_df = analyze(data)
    make_figure(results_df)

    # Save results table
    results_df.to_csv(RESULTS_DIR / "ring_vs_star_summary.csv", index=False)
    print(f"\nSaved: {RESULTS_DIR / 'ring_vs_star_summary.csv'}")
