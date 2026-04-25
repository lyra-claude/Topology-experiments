#!/usr/bin/env python3
"""
Analyze figure-eight cycle-length experiment.

Three directed figure-eight topologies (n=12, |E|=13, beta_1=2)
varying mean cycle length:
  M1: C3+C5, mean cycle length 4.0
  M2: C5+C6, mean cycle length 5.5
  M3: C7+C6, mean cycle length 6.5

Hypothesis: shorter cycles -> faster mixing -> larger lambda_2 ->
higher diversity maintenance under rugged landscapes.

Tested across 3 domains:
  onemax — negative control (smooth, no topology effect expected)
  nk4    — rugged (strong topology effects)
  maze   — deceptive/path-dependent
"""

import os
import csv
import math
from collections import defaultdict

RESULTS_DIR = "results/fig8_cycle_length"

TOPOLOGIES = ["fig8-m1", "fig8-m2", "fig8-m3"]
DOMAINS = ["onemax", "nk4", "maze"]

TOPO_LABELS = {
    "fig8-m1": "M1 (C3+C5, mcl=4.0)",
    "fig8-m2": "M2 (C5+C6, mcl=5.5)",
    "fig8-m3": "M3 (C7+C6, mcl=6.5)",
}

DOMAIN_LABELS = {
    "onemax": "OneMax (smooth)",
    "nk4": "NK K=4 (rugged)",
    "maze": "Maze 15x15 (deceptive)",
}

TOPO_COLORS = {
    "fig8-m1": "#e74c3c",   # red   — shortest cycles
    "fig8-m2": "#f39c12",   # amber — medium cycles
    "fig8-m3": "#3498db",   # blue  — longest cycles
}

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


def spearman_rank(xs, ys):
    """Spearman rank correlation (ties handled by average rank)."""
    if len(xs) != len(ys) or len(xs) < 3:
        return 0.0

    def rank(vals):
        sorted_idx = sorted(range(len(vals)), key=lambda i: vals[i])
        ranks = [0.0] * len(vals)
        i = 0
        while i < len(sorted_idx):
            j = i
            while j < len(sorted_idx) and vals[sorted_idx[j]] == vals[sorted_idx[i]]:
                j += 1
            avg_rank = (i + j - 1) / 2.0 + 1
            for k in range(i, j):
                ranks[sorted_idx[k]] = avg_rank
            i = j
        return ranks

    rx, ry = rank(xs), rank(ys)
    n = len(xs)
    d2 = sum((a - b) ** 2 for a, b in zip(rx, ry))
    return 1 - 6 * d2 / (n * (n * n - 1))


def analyze_domain(domain, data, fitness_data, file_counts):
    """Analyze one domain and print results."""
    print(f"\n{'=' * 70}")
    print(f"  {DOMAIN_LABELS[domain]}")
    print(f"{'=' * 70}")
    print()

    for topo in TOPOLOGIES:
        print(f"  {TOPO_LABELS[topo]}: {file_counts.get(topo, 0)} runs")
    print()

    if all(file_counts.get(t, 0) == 0 for t in TOPOLOGIES):
        print("  NO DATA — skipping analysis")
        return None

    # Table: Diversity at key checkpoints
    key_gens = [0, 50, 100, 200, 300, 500]
    print("-" * 70)
    print("Diversity: mean (SE) at key generations")
    print("-" * 70)
    header = f"{'Gen':>5}"
    for topo in TOPOLOGIES:
        header += f"  {TOPO_LABELS[topo]:>22}"
    header += f"  {'eta^2':>8}  {'F':>8}"
    print(header)
    print("-" * 70)

    eta2_diversity = {}
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
        eta2_diversity[gen] = eta2
        row += f"  {eta2:>8.4f}  {f_val:>8.2f}"
        print(row)
    print()

    # Compute eta-squared at all checkpoints
    all_eta2 = {}
    for gen in CHECKPOINTS:
        groups = [data[topo].get(gen, []) for topo in TOPOLOGIES]
        if all(len(g) > 0 for g in groups):
            all_eta2[gen] = eta_squared(groups)

    if all_eta2:
        max_eta2 = max(all_eta2.values())
        max_gen = max(all_eta2, key=all_eta2.get)
        print(f"Peak eta^2(diversity) = {max_eta2:.4f} at generation {max_gen}")

        # Pairwise Cohen's d at peak
        print(f"\nPairwise Cohen's d (diversity) at generation {max_gen}:")
        for i, t1 in enumerate(TOPOLOGIES):
            for t2 in TOPOLOGIES[i+1:]:
                d = cohens_d(data[t1].get(max_gen, []), data[t2].get(max_gen, []))
                print(f"  {TOPO_LABELS[t1]} vs {TOPO_LABELS[t2]}: d = {d:.3f}")

        # Monotonicity check: is diversity ordered M1 > M2 > M3?
        print(f"\nOrdering check at generation {max_gen} (expect M1 > M2 > M3 diversity):")
        means_at_peak = {t: mean(data[t].get(max_gen, [])) for t in TOPOLOGIES}
        for t in TOPOLOGIES:
            print(f"  {TOPO_LABELS[t]}: {means_at_peak[t]:.4f}")
        if means_at_peak["fig8-m1"] > means_at_peak["fig8-m2"] > means_at_peak["fig8-m3"]:
            print("  --> CONFIRMED: M1 > M2 > M3 (shorter cycles = more diversity)")
        else:
            print("  --> NOT confirmed. Ordering deviates from spectral gap prediction.")

    return all_eta2


def main():
    print("=" * 70)
    print("FIGURE-EIGHT CYCLE-LENGTH EXPERIMENT")
    print("n=12, |E|=13, beta_1=2 — varying mean cycle length")
    print("=" * 70)

    all_domain_eta2 = {}

    for domain in DOMAINS:
        domain_dir = os.path.join(RESULTS_DIR, domain)
        if not os.path.isdir(domain_dir):
            print(f"\n  {DOMAIN_LABELS[domain]}: directory not found, skipping")
            continue

        data = defaultdict(lambda: defaultdict(list))
        fitness_data = defaultdict(lambda: defaultdict(list))
        file_counts = defaultdict(int)

        for topo in TOPOLOGIES:
            for fname in os.listdir(domain_dir):
                if not fname.startswith(topo + "_seed"):
                    continue
                path = os.path.join(domain_dir, fname)
                rows = load_csv(path)
                file_counts[topo] += 1
                for row in rows:
                    gen = int(row["generation"])
                    div_val = float(row["diversity"])
                    fit_val = float(row["bestFitness"])
                    data[topo][gen].append(div_val)
                    fitness_data[topo][gen].append(fit_val)

        eta2 = analyze_domain(domain, data, fitness_data, file_counts)
        if eta2:
            all_domain_eta2[domain] = eta2

    # Cross-domain comparison
    if len(all_domain_eta2) > 1:
        print(f"\n{'=' * 70}")
        print("CROSS-DOMAIN COMPARISON")
        print(f"{'=' * 70}")
        print()
        print("Peak eta^2(diversity) by domain:")
        for domain in DOMAINS:
            if domain in all_domain_eta2:
                eta2 = all_domain_eta2[domain]
                peak = max(eta2.values())
                peak_gen = max(eta2, key=eta2.get)
                print(f"  {DOMAIN_LABELS[domain]:>30}: {peak:.4f} (gen {peak_gen})")

        print()
        print("Prediction: onemax ~ 0 (no topology effect),")
        print("            nk4 >> 0 (strong topology effect),")
        print("            maze >= nk4 (deceptive landscapes amplify topology).")

    # Generate plot if matplotlib available
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        n_domains = len(all_domain_eta2)
        if n_domains == 0:
            return

        fig, axes = plt.subplots(n_domains, 2, figsize=(14, 5 * n_domains),
                                 squeeze=False)

        for idx, domain in enumerate(d for d in DOMAINS if d in all_domain_eta2):
            domain_dir = os.path.join(RESULTS_DIR, domain)
            data = defaultdict(lambda: defaultdict(list))
            for topo in TOPOLOGIES:
                for fname in os.listdir(domain_dir):
                    if not fname.startswith(topo + "_seed"):
                        continue
                    path = os.path.join(domain_dir, fname)
                    rows = load_csv(path)
                    for row in rows:
                        gen = int(row["generation"])
                        data[topo][gen].append(float(row["diversity"]))

            # Left: diversity trajectories
            ax1 = axes[idx][0]
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
            ax1.set_ylabel("Diversity", fontsize=11)
            ax1.set_title(f"{DOMAIN_LABELS[domain]} — Diversity", fontsize=12)
            ax1.legend(fontsize=9)
            ax1.grid(True, alpha=0.3)

            # Right: eta-squared over time
            ax2 = axes[idx][1]
            eta2 = all_domain_eta2[domain]
            gens_eta = sorted(eta2.keys())
            ax2.plot(gens_eta, [eta2[g] for g in gens_eta],
                     color='#8e44ad', linewidth=2)
            ax2.axhline(y=0.14, color='gray', linestyle=':', alpha=0.5,
                         label='Large (0.14)')
            ax2.axhline(y=0.06, color='gray', linestyle='--', alpha=0.3,
                         label='Medium (0.06)')
            ax2.set_ylabel("$\\eta^2$", fontsize=11)
            ax2.set_title(f"{DOMAIN_LABELS[domain]} — Effect size", fontsize=12)
            ax2.legend(fontsize=9)
            ax2.grid(True, alpha=0.3)
            ax2.set_ylim(bottom=0)

            if idx == n_domains - 1:
                ax1.set_xlabel("Generation", fontsize=11)
                ax2.set_xlabel("Generation", fontsize=11)

        fig.suptitle("Figure-Eight Cycle Length Experiment\n"
                     "n=12, |E|=13, $\\beta_1$=2 — varying mean cycle length",
                     fontsize=14, y=1.02)
        plt.tight_layout()
        plot_path = os.path.join(RESULTS_DIR, "fig8_cycle_length.png")
        fig.savefig(plot_path, dpi=150, bbox_inches='tight')
        print(f"\nPlot saved to: {plot_path}")

    except ImportError:
        print("\nmatplotlib not available; skipping plot generation.")


if __name__ == "__main__":
    main()
