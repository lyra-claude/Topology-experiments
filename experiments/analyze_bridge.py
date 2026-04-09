#!/usr/bin/env python3
"""
Bridge Experiment Analysis: Disentangling beta_1 from lambda_2

Two iso-spectral families with constant lambda_2 and varying beta_1:
  Family 1 (lambda_2=0.5858): ring(b1=1), ring-chord1(b1=2), ring-chord2(b1=3), ring-chord3(b1=4)
  Family 2 (lambda_2=1.0):    star(b1=0), star-leaf1(b1=1), star-leaf2(b1=2), star-leaf3(b1=3)

Analyses:
  1. Within-family one-way ANOVA on diversity at gen 50, 100, 200, 500 → eta^2
  2. Cross-family comparison at matched beta_1 levels
  3. Two-way ANOVA: family x beta_1
  4. Publication figures
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
RESULTS_DIR = Path(__file__).parent.parent / "results" / "bridge_experiment"
TIMEPOINTS = [50, 100, 200, 500]

# Topology metadata: (name, family, beta_1, lambda_2)
TOPO_META = {
    "ring":         ("cycle", 1, 0.5858),
    "ring-chord1":  ("cycle", 2, 0.5858),
    "ring-chord2":  ("cycle", 3, 0.5858),
    "ring-chord3":  ("cycle", 4, 0.5858),
    "star":         ("star",  0, 1.0),
    "star-leaf1":   ("star",  1, 1.0),
    "star-leaf2":   ("star",  2, 1.0),
    "star-leaf3":   ("star",  3, 1.0),
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

        family, beta1, lambda2 = TOPO_META[topo]

        df = pd.read_csv(f)
        df["domain"] = domain
        df["topology"] = topo
        df["family"] = family
        df["beta1"] = beta1
        df["lambda2"] = lambda2
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


def eta_squared_two(g1, g2):
    """Convenience: eta-squared for exactly two groups."""
    return eta_squared_groups(g1, g2)


# ---------------------------------------------------------------------------
# Analysis 1: Within-family ANOVA
# ---------------------------------------------------------------------------

def within_family_analysis(data):
    """One-way ANOVA on diversity within each iso-spectral family."""
    print("=" * 80)
    print("ANALYSIS 1: Within-family ANOVA (constant lambda_2, varying beta_1)")
    print("=" * 80)

    results = []

    for domain in sorted(data["domain"].unique()):
        print(f"\n  Domain: {domain}")
        for family_name, lambda2 in [("cycle", 0.5858), ("star", 1.0)]:
            print(f"\n  Family: {family_name} (lambda_2 = {lambda2})")
            print(f"  {'Gen':>5} | {'Metric':>12} | {'eta^2':>8} | {'F-stat':>8} | {'p-value':>10} | {'beta_1 means':>40}")
            print(f"  {'─'*5}-+-{'─'*12}-+-{'─'*8}-+-{'─'*8}-+-{'─'*10}-+-{'─'*40}")

            for gen in TIMEPOINTS:
                snapshot = data[(data["domain"] == domain) &
                               (data["family"] == family_name) &
                               (data["generation"] == gen)]

                for metric in ["diversity", "bestFitness"]:
                    groups = []
                    beta1_levels = sorted(snapshot["beta1"].unique())
                    means = []
                    for b1 in beta1_levels:
                        vals = snapshot[snapshot["beta1"] == b1][metric].values
                        groups.append(vals)
                        means.append(f"b{b1}={np.mean(vals):.4f}")

                    if len(groups) < 2 or any(len(g) == 0 for g in groups):
                        continue

                    eta2 = eta_squared_groups(*groups)
                    f_stat, p_val = sp_stats.f_oneway(*groups)

                    sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
                    means_str = ", ".join(means)
                    print(f"  {gen:>5} | {metric:>12} | {eta2:>8.4f} | {f_stat:>8.2f} | {p_val:>9.2e} {sig} | {means_str}")

                    results.append({
                        "analysis": "within_family",
                        "domain": domain,
                        "family": family_name,
                        "lambda2": lambda2,
                        "generation": gen,
                        "metric": metric,
                        "eta_squared": eta2,
                        "f_stat": f_stat,
                        "p_value": p_val,
                    })

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Analysis 2: Cross-family comparison at matched beta_1
# ---------------------------------------------------------------------------

def cross_family_analysis(data):
    """Compare families at matched beta_1 levels."""
    print("\n\n" + "=" * 80)
    print("ANALYSIS 2: Cross-family comparison (matched beta_1, different lambda_2)")
    print("=" * 80)

    # Matched levels: beta_1 = 1, 2, 3
    matched_pairs = {
        1: ("ring", "star-leaf1"),
        2: ("ring-chord1", "star-leaf2"),
        3: ("ring-chord2", "star-leaf3"),
    }

    results = []

    for domain in sorted(data["domain"].unique()):
        print(f"\n  Domain: {domain}")
        print(f"  {'Gen':>5} | {'beta_1':>6} | {'Metric':>12} | {'Cycle Mean':>10} | {'Star Mean':>10} | {'eta^2':>8} | {'p-value':>10}")
        print(f"  {'─'*5}-+-{'─'*6}-+-{'─'*12}-+-{'─'*10}-+-{'─'*10}-+-{'─'*8}-+-{'─'*10}")

        for gen in TIMEPOINTS:
            for b1, (cycle_topo, star_topo) in sorted(matched_pairs.items()):
                snapshot = data[(data["domain"] == domain) & (data["generation"] == gen)]

                for metric in ["diversity", "bestFitness"]:
                    cycle_vals = snapshot[snapshot["topology"] == cycle_topo][metric].values
                    star_vals = snapshot[snapshot["topology"] == star_topo][metric].values

                    if len(cycle_vals) == 0 or len(star_vals) == 0:
                        continue

                    eta2 = eta_squared_two(cycle_vals, star_vals)
                    u_stat, p_val = sp_stats.mannwhitneyu(cycle_vals, star_vals, alternative='two-sided')

                    sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
                    print(f"  {gen:>5} | {b1:>6} | {metric:>12} | {np.mean(cycle_vals):>10.4f} | {np.mean(star_vals):>10.4f} | {eta2:>8.4f} | {p_val:>9.2e} {sig}")

                    results.append({
                        "analysis": "cross_family",
                        "domain": domain,
                        "beta1": b1,
                        "generation": gen,
                        "metric": metric,
                        "cycle_mean": np.mean(cycle_vals),
                        "star_mean": np.mean(star_vals),
                        "eta_squared": eta2,
                        "p_value": p_val,
                    })

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Analysis 3: Two-way ANOVA (family x beta_1)
# ---------------------------------------------------------------------------

def twoway_anova(data):
    """Two-way ANOVA: family x beta_1 on diversity/fitness."""
    print("\n\n" + "=" * 80)
    print("ANALYSIS 3: Two-way ANOVA (family x beta_1)")
    print("=" * 80)
    print("  (Using matched beta_1 levels 1, 2, 3 that appear in both families)")

    results = []

    # Only use beta_1 levels present in both families: 1, 2, 3
    matched = data[data["beta1"].isin([1, 2, 3])].copy()

    for domain in sorted(matched["domain"].unique()):
        print(f"\n  Domain: {domain}")

        for gen in TIMEPOINTS:
            snapshot = matched[(matched["domain"] == domain) & (matched["generation"] == gen)]

            for metric in ["diversity", "bestFitness"]:
                # Main effect of beta_1 (pooling families)
                b1_groups = [snapshot[snapshot["beta1"] == b][metric].values
                             for b in [1, 2, 3]]
                eta2_beta1 = eta_squared_groups(*b1_groups)
                _, p_beta1 = sp_stats.f_oneway(*b1_groups)

                # Main effect of family (pooling beta_1)
                fam_groups = [snapshot[snapshot["family"] == f][metric].values
                              for f in ["cycle", "star"]]
                eta2_family = eta_squared_groups(*fam_groups)
                _, p_family = sp_stats.f_oneway(*fam_groups)

                sig_b = "***" if p_beta1 < 0.001 else "**" if p_beta1 < 0.01 else "*" if p_beta1 < 0.05 else "ns"
                sig_f = "***" if p_family < 0.001 else "**" if p_family < 0.01 else "*" if p_family < 0.05 else "ns"

                print(f"  Gen {gen:>3}, {metric:>12}: "
                      f"beta_1 eta^2={eta2_beta1:.4f} (p={p_beta1:.2e} {sig_b}), "
                      f"family eta^2={eta2_family:.4f} (p={p_family:.2e} {sig_f})")

                results.append({
                    "analysis": "twoway",
                    "domain": domain,
                    "generation": gen,
                    "metric": metric,
                    "eta2_beta1": eta2_beta1,
                    "p_beta1": p_beta1,
                    "eta2_family": eta2_family,
                    "p_family": p_family,
                })

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def print_summary(within_df, cross_df, twoway_df):
    """Print the key verdict."""
    print("\n\n" + "=" * 80)
    print("VERDICT")
    print("=" * 80)

    # Within-family eta^2 at gen 500 for diversity
    within_500 = within_df[(within_df["generation"] == 500) & (within_df["metric"] == "diversity")]
    if not within_500.empty:
        avg_within = within_500["eta_squared"].mean()
        print(f"\n  Within-family eta^2 (diversity, gen 500): {avg_within:.4f}")

    # Cross-family eta^2 at gen 500 for diversity
    cross_500 = cross_df[(cross_df["generation"] == 500) & (cross_df["metric"] == "diversity")]
    if not cross_500.empty:
        avg_cross = cross_500["eta_squared"].mean()
        print(f"  Cross-family eta^2 (diversity, gen 500):  {avg_cross:.4f}")

    # Two-way at gen 500
    twoway_500 = twoway_df[(twoway_df["generation"] == 500) & (twoway_df["metric"] == "diversity")]
    if not twoway_500.empty:
        for _, row in twoway_500.iterrows():
            print(f"  Two-way ({row['domain']}): beta_1 eta^2={row['eta2_beta1']:.4f}, family eta^2={row['eta2_family']:.4f}")

    if not within_500.empty and not cross_500.empty:
        if avg_within > 3 * avg_cross:
            print("\n  >>> beta_1 DOMINATES lambda_2. Case closed.")
        elif avg_within > avg_cross:
            print("\n  >>> beta_1 effect stronger than lambda_2, but not overwhelming.")
        else:
            print("\n  >>> lambda_2 shows comparable or stronger effect — needs more investigation.")

    # Compare with ring-vs-star baseline
    print("\n  Reference: Ring vs Star (same beta_1=1, different lambda_2) gave eta^2 ~ 0.08-0.13.")
    if not within_500.empty:
        if avg_within > 0.13:
            print(f"  Within-family beta_1 effect ({avg_within:.4f}) > ring-vs-star lambda_2 effect (0.13).")
        else:
            print(f"  Within-family beta_1 effect ({avg_within:.4f}) <= ring-vs-star lambda_2 effect (0.13).")


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def make_figures(data):
    """Publication-quality figures."""
    figdir = RESULTS_DIR / "figures"
    figdir.mkdir(exist_ok=True)

    colors_family = {"cycle": "#2196F3", "star": "#FF5722"}
    markers_family = {"cycle": "o", "star": "s"}

    # Figure 1: Diversity trajectories by topology, one panel per domain
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    for ax, domain in zip(axes, ["nk0", "nk4"]):
        dom_data = data[data["domain"] == domain]

        for topo, meta in TOPO_META.items():
            family, beta1, _ = meta
            topo_data = dom_data[dom_data["topology"] == topo]
            if topo_data.empty:
                continue

            mean_div = topo_data.groupby("generation")["diversity"].mean()
            color = colors_family[family]
            alpha = 0.4 + 0.2 * beta1  # darker = higher beta_1
            lw = 1.0 + 0.3 * beta1
            label = f"{topo} (b1={beta1})"

            ax.plot(mean_div.index, mean_div.values, color=color,
                    alpha=min(alpha, 1.0), linewidth=lw, label=label)

        ax.set_xlabel("Generation", fontsize=11)
        ax.set_title(f"Domain: {domain.upper()}", fontsize=12, fontweight='bold')
        ax.legend(fontsize=7, ncol=2, loc='best')
        ax.grid(True, alpha=0.3)

    axes[0].set_ylabel("Diversity", fontsize=11)
    fig.suptitle("Bridge Experiment: Diversity Trajectories by Topology",
                 fontsize=13, fontweight='bold')
    plt.tight_layout()

    for ext in ["png", "pdf"]:
        fig.savefig(figdir / f"bridge_diversity_trajectories.{ext}", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {figdir / 'bridge_diversity_trajectories.png'}")

    # Figure 2: eta^2 vs beta_1 within each family at gen 500
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

    for ax, domain in zip(axes, ["nk0", "nk4"]):
        for family_name, lambda2 in [("cycle", 0.5858), ("star", 1.0)]:
            snapshot = data[(data["domain"] == domain) &
                           (data["family"] == family_name) &
                           (data["generation"] == 500)]

            beta1_levels = sorted(snapshot["beta1"].unique())
            div_means = []
            div_sems = []
            for b1 in beta1_levels:
                vals = snapshot[snapshot["beta1"] == b1]["diversity"].values
                div_means.append(np.mean(vals))
                div_sems.append(np.std(vals) / np.sqrt(len(vals)))

            color = colors_family[family_name]
            marker = markers_family[family_name]
            ax.errorbar(beta1_levels, div_means, yerr=div_sems,
                       color=color, marker=marker, markersize=8, linewidth=2,
                       capsize=4, label=f"{family_name} (λ₂={lambda2})")

        ax.set_xlabel("β₁ (cycle rank)", fontsize=11)
        ax.set_title(f"Domain: {domain.upper()}", fontsize=12, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

    axes[0].set_ylabel("Diversity at Generation 500", fontsize=11)
    fig.suptitle("Bridge Experiment: Diversity vs β₁ (Constant λ₂)",
                 fontsize=13, fontweight='bold')
    plt.tight_layout()

    for ext in ["png", "pdf"]:
        fig.savefig(figdir / f"bridge_diversity_vs_beta1.{ext}", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {figdir / 'bridge_diversity_vs_beta1.png'}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    data = load_all_results()
    n_runs = len(data.groupby(["domain", "topology", "seed"]))
    print(f"Loaded {len(data)} rows from {n_runs} runs\n")

    within_df = within_family_analysis(data)
    cross_df = cross_family_analysis(data)
    twoway_df = twoway_anova(data)
    print_summary(within_df, cross_df, twoway_df)

    print("\n\nGenerating figures...")
    make_figures(data)

    # Save results
    within_df.to_csv(RESULTS_DIR / "within_family_results.csv", index=False)
    cross_df.to_csv(RESULTS_DIR / "cross_family_results.csv", index=False)
    twoway_df.to_csv(RESULTS_DIR / "twoway_results.csv", index=False)
    print(f"\nSaved CSV results to {RESULTS_DIR}/")
