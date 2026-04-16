#!/usr/bin/env python3
"""
Analyze directed cycle experiment on Maze domain.
8 topologies x 30 seeds = 240 runs.

Key questions:
1. Does directed cycle count (kappa) predict diversity on maze navigation?
2. How does the maze domain compare to NK landscapes?
3. Does temporal inversion occur (positive correlation early, negative late)?
"""

import os
import csv
import math
from collections import defaultdict

RESULTS_DIR = "results/directed-maze"

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

CHECKPOINTS = [30, 100, 200, 500]


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


def pearson_r(xs, ys):
    """Pearson correlation coefficient."""
    n = len(xs)
    if n < 3:
        return 0.0
    mx, my = mean(xs), mean(ys)
    sx = math.sqrt(sum((x - mx)**2 for x in xs))
    sy = math.sqrt(sum((y - my)**2 for y in ys))
    if sx == 0 or sy == 0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / (sx * sy)


def spearman_rho(xs, ys):
    """Spearman rank correlation."""
    def rank(vals):
        indexed = sorted(enumerate(vals), key=lambda p: p[1])
        ranks = [0.0] * len(vals)
        i = 0
        while i < len(indexed):
            j = i
            while j < len(indexed) and indexed[j][1] == indexed[i][1]:
                j += 1
            avg_rank = (i + j - 1) / 2.0 + 1
            for k in range(i, j):
                ranks[indexed[k][0]] = avg_rank
            i = j
        return ranks
    return pearson_r(rank(xs), rank(ys))


def eta_squared(groups):
    """One-way ANOVA eta-squared: fraction of variance explained by group."""
    all_vals = []
    for g in groups:
        all_vals.extend(g)
    if not all_vals:
        return 0.0
    grand_mean = mean(all_vals)
    ss_total = sum((x - grand_mean)**2 for x in all_vals)
    if ss_total == 0:
        return 0.0
    ss_between = sum(len(g) * (mean(g) - grand_mean)**2 for g in groups if g)
    return ss_between / ss_total


def load_data():
    """Load all maze results. Returns {topology: {seed: [rows]}}."""
    if not os.path.isdir(RESULTS_DIR):
        print(f"ERROR: Results directory not found: {RESULTS_DIR}")
        return {}
    data = defaultdict(dict)
    for fname in os.listdir(RESULTS_DIR):
        if not fname.endswith(".csv"):
            continue
        parts = fname.replace(".csv", "").rsplit("_seed", 1)
        if len(parts) != 2:
            continue
        topo, seed = parts[0], int(parts[1])
        rows = load_csv(os.path.join(RESULTS_DIR, fname))
        data[topo][seed] = rows
    return data


def get_metric_at_gen(rows, gen, metric="diversity"):
    """Get metric value at or near a specific generation."""
    for row in rows:
        if int(row["generation"]) == gen:
            return float(row[metric])
    closest = min(rows, key=lambda r: abs(int(r["generation"]) - gen))
    return float(closest[metric])


def analyze(data):
    """Full analysis of maze domain results."""
    print("Directed Cycle Experiment — Maze Domain (8x8)")
    print("=" * 70)

    topos = sorted(data.keys(), key=lambda t: CYCLE_COUNTS.get(t, 0))
    n_topos = len(topos)
    total_runs = sum(len(seeds) for seeds in data.values())
    print(f"Topologies: {n_topos}, Total runs: {total_runs}")

    eta2_summary = {}

    for gen in CHECKPOINTS:
        print(f"\n--- Generation {gen} ---")

        topo_diversities = {}
        topo_fitnesses = {}
        groups_div = []
        groups_fit = []

        for topo in topos:
            divs = []
            fits = []
            for seed, rows in data[topo].items():
                try:
                    divs.append(get_metric_at_gen(rows, gen, "diversity"))
                    fits.append(get_metric_at_gen(rows, gen, "meanFitness"))
                except (KeyError, ValueError, IndexError):
                    continue
            topo_diversities[topo] = divs
            topo_fitnesses[topo] = fits
            groups_div.append(divs)
            groups_fit.append(fits)

        # Print summary table
        print(f"{'Topology':<18} {'kappa':>6} {'Diversity':>12} {'SE':>8} {'MeanFit':>12} {'SE':>8} {'N':>4}")
        print(f"{'-'*18} {'-'*6} {'-'*12} {'-'*8} {'-'*12} {'-'*8} {'-'*4}")
        for topo in topos:
            kappa = CYCLE_COUNTS.get(topo, -1)
            divs = topo_diversities[topo]
            fits = topo_fitnesses[topo]
            print(f"{topo:<18} {kappa:>6} {mean(divs):>12.6f} {stderr(divs):>8.6f} "
                  f"{mean(fits):>12.6f} {stderr(fits):>8.6f} {len(divs):>4}")

        # Eta-squared
        eta2_div = eta_squared(groups_div)
        eta2_fit = eta_squared(groups_fit)
        print(f"\neta^2 (topology -> diversity): {eta2_div:.4f}")
        print(f"eta^2 (topology -> fitness):   {eta2_fit:.4f}")
        eta2_summary[gen] = (eta2_div, eta2_fit)

        # Correlation with kappa
        all_kappas = []
        all_divs = []
        kappas = []
        mean_divs = []
        for topo in topos:
            if topo_diversities[topo]:
                k = CYCLE_COUNTS.get(topo, 0)
                kappas.append(k)
                mean_divs.append(mean(topo_diversities[topo]))
                for d in topo_diversities[topo]:
                    all_kappas.append(k)
                    all_divs.append(d)

        if len(kappas) >= 3:
            r_div = pearson_r(all_kappas, all_divs)
            rho_div = spearman_rho(all_kappas, all_divs)
            r_mean = pearson_r(kappas, mean_divs)
            print(f"Pearson r(kappa, diversity):  {r_div:+.4f}  (per-run, N={len(all_kappas)})")
            print(f"Spearman rho(kappa, div):     {rho_div:+.4f}")
            print(f"Pearson r(kappa, mean_div):   {r_mean:+.4f}  (topology means, N={len(kappas)})")

    # Summary comparison table
    print(f"\n{'='*70}")
    print("SUMMARY: eta^2 across generations")
    print(f"{'Gen':>6} {'eta2_div':>10} {'eta2_fit':>10}")
    print(f"{'-'*6} {'-'*10} {'-'*10}")
    for gen in CHECKPOINTS:
        ed, ef = eta2_summary[gen]
        print(f"{gen:>6} {ed:>10.4f} {ef:>10.4f}")

    # Temporal inversion check
    print(f"\nTEMPORAL INVERSION CHECK")
    print("(Does kappa-diversity correlation flip sign over time?)")
    signs = []
    for gen in CHECKPOINTS:
        all_kappas = []
        all_divs = []
        for topo in topos:
            k = CYCLE_COUNTS.get(topo, 0)
            for seed, rows in data[topo].items():
                try:
                    d = get_metric_at_gen(rows, gen, "diversity")
                    all_kappas.append(k)
                    all_divs.append(d)
                except:
                    continue
        r = pearson_r(all_kappas, all_divs)
        sign = "+" if r > 0 else "-"
        signs.append(sign)
        print(f"  gen {gen:>3}: r = {r:+.4f} ({sign})")
    sign_str = " -> ".join(f"gen{g}:{s}" for g, s in zip(CHECKPOINTS, signs))
    print(f"  Pattern: {sign_str}")

    print(f"\n{'='*70}")
    print("DONE")


def generate_plots(data):
    """Generate plots for maze results."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available — skipping plots")
        return

    topos = sorted(data.keys(), key=lambda t: CYCLE_COUNTS.get(t, 0))
    plot_dir = "results/directed-maze/plots"
    os.makedirs(plot_dir, exist_ok=True)

    # --- Plot 1: Diversity over time by topology ---
    fig, ax = plt.subplots(figsize=(12, 7))
    colors = plt.cm.viridis([i / (len(topos) - 1) for i in range(len(topos))])

    for idx, topo in enumerate(topos):
        # Collect all generations
        all_gens = set()
        for seed, rows in data[topo].items():
            for row in rows:
                all_gens.add(int(row["generation"]))
        all_gens = sorted(all_gens)

        mean_divs = []
        se_divs = []
        for gen in all_gens:
            vals = []
            for seed, rows in data[topo].items():
                try:
                    vals.append(get_metric_at_gen(rows, gen, "diversity"))
                except:
                    continue
            mean_divs.append(mean(vals))
            se_divs.append(stderr(vals))

        kappa = CYCLE_COUNTS.get(topo, 0)
        ax.plot(all_gens, mean_divs, color=colors[idx],
                label=f"{topo} (kappa={kappa})", linewidth=1.5)

    ax.set_xlabel("Generation", fontsize=12)
    ax.set_ylabel("Mean Diversity", fontsize=12)
    ax.set_title("Directed Topology Effect on Diversity — Maze Domain (8x8)", fontsize=14)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/diversity_over_time.png", dpi=150)
    plt.close()
    print(f"Saved: {plot_dir}/diversity_over_time.png")

    # --- Plot 2: eta-squared over time ---
    fig, ax = plt.subplots(figsize=(10, 6))

    all_gens = set()
    for topo in topos:
        for seed, rows in data[topo].items():
            for row in rows:
                all_gens.add(int(row["generation"]))
    all_gens = sorted(all_gens)

    eta2_vals = []
    for gen in all_gens:
        groups = []
        for topo in topos:
            vals = []
            for seed, rows in data[topo].items():
                try:
                    vals.append(get_metric_at_gen(rows, gen, "diversity"))
                except:
                    continue
            groups.append(vals)
        eta2_vals.append(eta_squared(groups))

    ax.plot(all_gens, eta2_vals, 'b-', linewidth=2)
    ax.set_xlabel("Generation", fontsize=12)
    ax.set_ylabel("eta-squared (topology -> diversity)", fontsize=12)
    ax.set_title("Effect Size Over Time — Maze Domain (8x8)", fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/eta_squared_over_time.png", dpi=150)
    plt.close()
    print(f"Saved: {plot_dir}/eta_squared_over_time.png")

    # --- Plot 3: Kappa vs diversity at checkpoints ---
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    for ax_idx, gen in enumerate(CHECKPOINTS):
        ax = axes[ax_idx // 2][ax_idx % 2]
        for topo in topos:
            kappa = CYCLE_COUNTS.get(topo, 0)
            divs = []
            for seed, rows in data[topo].items():
                try:
                    divs.append(get_metric_at_gen(rows, gen, "diversity"))
                except:
                    continue
            # Jitter kappa slightly for visibility
            jitter = [kappa + (i - len(divs)/2) * 0.3 for i in range(len(divs))]
            ax.scatter(jitter, divs, alpha=0.4, s=15)
            ax.scatter([kappa], [mean(divs)], color='red', s=80, zorder=5,
                      edgecolors='black', linewidth=1)

        groups = []
        for topo in topos:
            vals = []
            for seed, rows in data[topo].items():
                try:
                    vals.append(get_metric_at_gen(rows, gen, "diversity"))
                except:
                    continue
            groups.append(vals)
        eta2 = eta_squared(groups)

        ax.set_title(f"Gen {gen} (eta^2 = {eta2:.4f})", fontsize=11)
        ax.set_xlabel("Directed cycle count (kappa)")
        ax.set_ylabel("Diversity")
        ax.grid(True, alpha=0.3)

    plt.suptitle("Kappa vs Diversity at Checkpoints — Maze Domain (8x8)", fontsize=14)
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/kappa_vs_diversity.png", dpi=150)
    plt.close()
    print(f"Saved: {plot_dir}/kappa_vs_diversity.png")

    # --- Plot 4: Fitness over time ---
    fig, ax = plt.subplots(figsize=(12, 7))
    for idx, topo in enumerate(topos):
        all_gens_t = set()
        for seed, rows in data[topo].items():
            for row in rows:
                all_gens_t.add(int(row["generation"]))
        all_gens_t = sorted(all_gens_t)

        mean_fits = []
        for gen in all_gens_t:
            vals = []
            for seed, rows in data[topo].items():
                try:
                    vals.append(get_metric_at_gen(rows, gen, "meanFitness"))
                except:
                    continue
            mean_fits.append(mean(vals))

        kappa = CYCLE_COUNTS.get(topo, 0)
        ax.plot(all_gens_t, mean_fits, color=colors[idx],
                label=f"{topo} (kappa={kappa})", linewidth=1.5)

    ax.set_xlabel("Generation", fontsize=12)
    ax.set_ylabel("Mean Fitness", fontsize=12)
    ax.set_title("Directed Topology Effect on Fitness — Maze Domain (8x8)", fontsize=14)
    ax.legend(fontsize=9, loc='lower right')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/fitness_over_time.png", dpi=150)
    plt.close()
    print(f"Saved: {plot_dir}/fitness_over_time.png")


if __name__ == "__main__":
    data = load_data()
    if data:
        analyze(data)
        print()
        generate_plots(data)
    else:
        print("No data found!")
