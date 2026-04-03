#!/usr/bin/env python3
"""
Publication-quality figures for the directed cycle experiment.
Produces 4 figures saved as PNG (300 dpi) and PDF.
"""

import os
import csv
import math
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from scipy import stats

RESULTS_DIR = "results/directed_full"
FIG_DIR = "results/directed_full/figures"

CYCLE_COUNTS = {
    "dag-layer": 0,
    "dag-wide": 0,
    "lowcyc-1": 3,
    "bidir-ring": 10,
    "two-cliques": 14,
    "mesh-cyclic": 20,
    "dense-triangles": 29,
    "ring-skip2": 47,
}


def load_csv(path):
    with open(path) as f:
        lines = [l for l in f if l.strip() and
                 (l.strip()[0].isdigit() or l.strip().startswith("generation"))]
    reader = csv.DictReader(lines)
    return list(reader)


def mean(xs):
    return sum(xs) / len(xs) if xs else 0.0


def stderr(xs):
    if len(xs) < 2:
        return 0.0
    m = mean(xs)
    s = math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))
    return s / math.sqrt(len(xs))


def load_all_data():
    results = defaultdict(list)
    for fname in sorted(os.listdir(RESULTS_DIR)):
        if not fname.endswith(".csv"):
            continue
        topo = fname.rsplit("_seed", 1)[0]
        rows = load_csv(os.path.join(RESULTS_DIR, fname))
        if rows:
            results[topo].append(rows)
    return results


def get_all_generations(results):
    """Get sorted list of all generation numbers from data."""
    gens = set()
    for topo in results:
        for run in results[topo]:
            for row in run:
                gens.add(int(row["generation"]))
    return sorted(gens)


def get_metric_trajectories(results, metric="diversity"):
    """Return dict: topo -> (gens_array, means_array, se_array)."""
    topos_sorted = sorted(results.keys(), key=lambda t: CYCLE_COUNTS.get(t, 999))
    all_gens = get_all_generations(results)

    trajectories = {}
    for topo in topos_sorted:
        gens_out = []
        means_out = []
        ses_out = []
        for gen in all_gens:
            vals = []
            for run in results[topo]:
                for row in run:
                    if int(row["generation"]) == gen:
                        vals.append(float(row[metric]))
                        break
            if vals:
                gens_out.append(gen)
                means_out.append(mean(vals))
                ses_out.append(stderr(vals))
        trajectories[topo] = (np.array(gens_out), np.array(means_out), np.array(ses_out))
    return trajectories


def setup_style():
    """Set up clean academic plotting style."""
    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 11,
        'axes.linewidth': 0.8,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'figure.dpi': 150,
    })


def fig1_diversity_trajectories(results):
    """Line plot of diversity over generations, colored by cycle count."""
    trajectories = get_metric_trajectories(results, "diversity")
    topos_sorted = sorted(results.keys(), key=lambda t: CYCLE_COUNTS.get(t, 999))

    norm = Normalize(vmin=0, vmax=47)
    cmap = plt.cm.viridis

    fig, ax = plt.subplots(figsize=(8, 5))

    for topo in topos_sorted:
        gens, means, ses = trajectories[topo]
        cc = CYCLE_COUNTS[topo]
        color = cmap(norm(cc))
        ax.plot(gens, means, color=color, linewidth=1.5, label=f"{topo} ({cc})")
        ax.fill_between(gens, means - ses, means + ses, color=color, alpha=0.15)

    ax.set_xlabel("Generation")
    ax.set_ylabel("Mean Diversity")
    ax.set_title("Diversity Trajectories by Directed Cycle Count")
    ax.legend(fontsize=8, loc='upper right', frameon=False, ncol=2)

    # Colorbar
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02, aspect=30)
    cbar.set_label("Directed Cycle Count")

    fig.tight_layout()
    return fig


def fig2_correlation_over_time(results):
    """Pearson r over time with significance shading."""
    topos_sorted = sorted(results.keys(), key=lambda t: CYCLE_COUNTS.get(t, 999))
    all_gens = get_all_generations(results)

    rs = []
    ps = []
    valid_gens = []

    for gen in all_gens:
        cycle_counts_list = []
        mean_divs_list = []
        for topo in topos_sorted:
            vals = []
            for run in results[topo]:
                for row in run:
                    if int(row["generation"]) == gen:
                        vals.append(float(row["diversity"]))
                        break
            if vals:
                cycle_counts_list.append(CYCLE_COUNTS[topo])
                mean_divs_list.append(mean(vals))

        if len(cycle_counts_list) >= 3:
            r, p = stats.pearsonr(cycle_counts_list, mean_divs_list)
            rs.append(r)
            ps.append(p)
            valid_gens.append(gen)

    valid_gens = np.array(valid_gens)
    rs = np.array(rs)
    ps = np.array(ps)

    fig, ax = plt.subplots(figsize=(8, 4.5))

    # Shade significant regions
    sig_mask = ps < 0.05
    if np.any(sig_mask):
        for i in range(len(valid_gens)):
            if sig_mask[i]:
                # Shade a band around this generation
                if i == 0:
                    left = valid_gens[i]
                else:
                    left = (valid_gens[i-1] + valid_gens[i]) / 2
                if i == len(valid_gens) - 1:
                    right = valid_gens[i]
                else:
                    right = (valid_gens[i] + valid_gens[i+1]) / 2
                ax.axvspan(left, right, color='lightcoral', alpha=0.2)

    ax.plot(valid_gens, rs, 'k-', linewidth=1.5)
    ax.plot(valid_gens[sig_mask], rs[sig_mask], 'ro', markersize=5, label='p < 0.05')
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)

    ax.set_xlabel("Generation")
    ax.set_ylabel("Pearson r (cycle count vs diversity)")
    ax.set_title("Correlation Between Directed Cycle Count and Diversity Over Time")
    ax.legend(frameon=False)

    fig.tight_layout()
    return fig


def fig3_scatter_gen30(results):
    """Scatter: cycle count vs diversity at gen 30 with regression."""
    topos_sorted = sorted(results.keys(), key=lambda t: CYCLE_COUNTS.get(t, 999))

    xs = []
    ys = []
    labels = []

    for topo in topos_sorted:
        vals = []
        for run in results[topo]:
            for row in run:
                if int(row["generation"]) == 30:
                    vals.append(float(row["diversity"]))
                    break
        if vals:
            xs.append(CYCLE_COUNTS[topo])
            ys.append(mean(vals))
            labels.append(topo)

    xs = np.array(xs)
    ys = np.array(ys)

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.scatter(xs, ys, s=60, c='steelblue', edgecolors='black', linewidths=0.8, zorder=3)

    # Regression line
    slope, intercept, r, p, se = stats.linregress(xs, ys)
    x_line = np.linspace(xs.min() - 2, xs.max() + 2, 100)
    y_line = slope * x_line + intercept
    ax.plot(x_line, y_line, 'r--', linewidth=1.2,
            label=f"r = {r:.3f}, p = {p:.3f}")

    # Label points
    for xi, yi, lab in zip(xs, ys, labels):
        ax.annotate(lab, (xi, yi), textcoords="offset points",
                    xytext=(5, 6), fontsize=8, color='0.3')

    ax.set_xlabel("Directed Simple Cycle Count")
    ax.set_ylabel("Mean Diversity at Generation 30")
    ax.set_title("Cycle Count vs Diversity (n=8, m=16, constant density)")
    ax.legend(frameon=False, fontsize=10)

    fig.tight_layout()
    return fig


def fig4_effect_size(results):
    """Eta-squared for diversity over time with reference lines."""
    topos_sorted = sorted(results.keys(), key=lambda t: CYCLE_COUNTS.get(t, 999))
    all_gens = get_all_generations(results)

    eta2_div = []
    eta2_fit = []
    valid_gens = []

    for gen in all_gens:
        for metric, store in [("diversity", eta2_div), ("meanFitness", eta2_fit)]:
            all_vals = []
            group_means = []
            for topo in topos_sorted:
                vals = []
                for run in results[topo]:
                    for row in run:
                        if int(row["generation"]) == gen:
                            vals.append(float(row[metric]))
                            break
                all_vals.extend(vals)
                if vals:
                    group_means.append((mean(vals), len(vals)))

            if all_vals and len(group_means) > 1:
                grand_mean = mean(all_vals)
                ss_between = sum(n * (m - grand_mean)**2 for m, n in group_means)
                ss_total = sum((x - grand_mean)**2 for x in all_vals)
                eta2 = ss_between / ss_total if ss_total > 0 else 0
            else:
                eta2 = 0
            store.append(eta2)

        valid_gens.append(gen)

    valid_gens = np.array(valid_gens)
    eta2_div = np.array(eta2_div)
    eta2_fit = np.array(eta2_fit)

    fig, ax = plt.subplots(figsize=(8, 4.5))

    ax.plot(valid_gens, eta2_div, 'b-', linewidth=1.5, label=r'$\eta^2$ (diversity)')
    ax.plot(valid_gens, eta2_fit, 'g-', linewidth=1.5, label=r'$\eta^2$ (fitness)')

    # Reference lines
    ax.axhline(y=0.01, color='gray', linestyle=':', linewidth=0.8, alpha=0.7)
    ax.axhline(y=0.06, color='gray', linestyle=':', linewidth=0.8, alpha=0.7)
    ax.axhline(y=0.14, color='gray', linestyle=':', linewidth=0.8, alpha=0.7)

    ax.text(valid_gens[-1] + 5, 0.01, 'small', fontsize=8, color='gray', va='center')
    ax.text(valid_gens[-1] + 5, 0.06, 'medium', fontsize=8, color='gray', va='center')
    ax.text(valid_gens[-1] + 5, 0.14, 'large', fontsize=8, color='gray', va='center')

    ax.set_xlabel("Generation")
    ax.set_ylabel(r"Effect Size ($\eta^2$)")
    ax.set_title(r"Effect Size of Topology on Diversity and Fitness ($\eta^2$)")
    ax.legend(frameon=False)
    ax.set_ylim(bottom=-0.005)

    fig.tight_layout()
    return fig


def main():
    setup_style()
    os.makedirs(FIG_DIR, exist_ok=True)

    print("Loading data...")
    results = load_all_data()
    topos = sorted(results.keys(), key=lambda t: CYCLE_COUNTS.get(t, 999))
    print(f"  Loaded {sum(len(v) for v in results.values())} runs across {len(topos)} topologies")

    figures = [
        ("fig1_diversity_trajectories", fig1_diversity_trajectories),
        ("fig2_correlation_over_time", fig2_correlation_over_time),
        ("fig3_scatter_gen30", fig3_scatter_gen30),
        ("fig4_effect_size", fig4_effect_size),
    ]

    for name, func in figures:
        print(f"  Generating {name}...")
        fig = func(results)
        fig.savefig(os.path.join(FIG_DIR, f"{name}.png"), dpi=300, bbox_inches='tight')
        fig.savefig(os.path.join(FIG_DIR, f"{name}.pdf"), bbox_inches='tight')
        plt.close(fig)
        print(f"    Saved {name}.png and {name}.pdf")

    print(f"\nAll figures saved to {FIG_DIR}/")


if __name__ == "__main__":
    main()
