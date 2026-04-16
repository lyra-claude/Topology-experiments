#!/usr/bin/env python3
"""
Analyze star_n vs cycle_{n-1} experiment.

Edge-controlled comparison: star(n) and cycle(n-1) both have n-1 edges.
Star has β₁=0 (tree), cycle has β₁=1 (one independent loop).

Key questions:
1. Does β₁=1 (cycle) produce higher diversity than β₁=0 (star) at constant edge count?
2. Does the effect depend on landscape ruggedness (NK K)?
3. Does size (8/7 vs 12/11) modulate the effect?
4. Temporal dynamics: when does the topology effect emerge?
"""

import os
import csv
import math
from collections import defaultdict

RESULTS_DIR = "results/star_vs_cycle"

DOMAINS = ["nk0", "nk2", "nk4"]

# Topology conditions with their properties
TOPOS = {
    # name: (nodes, edges, β₁)
    "star8":   (8, 7, 0),
    "cycle7":  (7, 7, 1),
    "star12":  (12, 11, 0),
    "cycle11": (11, 11, 1),
}

# Group by size condition
SIZE_GROUPS = {
    "small": ["star8", "cycle7"],     # 7 edges
    "large": ["star12", "cycle11"],   # 11 edges
}

CHECKPOINTS = [10, 30, 50, 100, 200, 300, 500]


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


def cohens_d(xs, ys):
    """Cohen's d effect size (pooled SD)."""
    nx, ny = len(xs), len(ys)
    if nx < 2 or ny < 2:
        return 0.0
    mx, my = mean(xs), mean(ys)
    sx, sy = std(xs), std(ys)
    pooled = math.sqrt(((nx - 1) * sx**2 + (ny - 1) * sy**2) / (nx + ny - 2))
    if pooled == 0:
        return 0.0
    return (mx - my) / pooled


def eta_squared_oneway(groups):
    """One-way eta-squared: SS_between / SS_total."""
    all_vals = [v for g in groups for v in g]
    if not all_vals:
        return 0.0
    grand_mean = mean(all_vals)
    ss_total = sum((v - grand_mean)**2 for v in all_vals)
    if ss_total == 0:
        return 0.0
    ss_between = sum(len(g) * (mean(g) - grand_mean)**2 for g in groups if g)
    return ss_between / ss_total


def mann_whitney_u(xs, ys):
    """Mann-Whitney U test (approximate z for n >= 20)."""
    nx, ny = len(xs), len(ys)
    if nx == 0 or ny == 0:
        return 0.0, 1.0
    combined = [(v, 'x') for v in xs] + [(v, 'y') for v in ys]
    combined.sort(key=lambda p: p[0])
    # Assign ranks (handle ties)
    ranks = [0.0] * len(combined)
    i = 0
    while i < len(combined):
        j = i
        while j < len(combined) and combined[j][0] == combined[i][0]:
            j += 1
        avg_rank = (i + j + 1) / 2.0  # 1-based
        for k in range(i, j):
            ranks[k] = avg_rank
        i = j
    r1 = sum(ranks[i] for i in range(len(combined)) if combined[i][1] == 'x')
    u1 = r1 - nx * (nx + 1) / 2
    mu = nx * ny / 2
    sigma = math.sqrt(nx * ny * (nx + ny + 1) / 12)
    if sigma == 0:
        return u1, 1.0
    z = (u1 - mu) / sigma
    # Two-tailed p approximation
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return u1, p


def main():
    if not os.path.exists(RESULTS_DIR):
        print(f"ERROR: Results directory {RESULTS_DIR} not found.")
        print("Run experiments/run_star_vs_cycle.sh first.")
        return

    seeds = "42 137 2718 314 1618 7 99 256 512 1024 2048 4096 8192 16384 32768 11 23 37 53 71 89 97 113 131 151 173 191 199 211 223".split()

    # Load all data: data[topo][domain][seed] = {gen: {metric: value}}
    data = defaultdict(lambda: defaultdict(dict))
    missing = 0
    loaded = 0

    for topo in TOPOS:
        for domain in DOMAINS:
            for seed in seeds:
                fname = f"{topo}_{domain}_seed{seed}.csv"
                path = os.path.join(RESULTS_DIR, fname)
                if not os.path.exists(path):
                    missing += 1
                    continue
                rows = load_csv(path)
                gen_data = {}
                for row in rows:
                    gen = int(row["generation"])
                    gen_data[gen] = {
                        "meanFitness": float(row["meanFitness"]),
                        "bestFitness": float(row["bestFitness"]),
                        "diversity": float(row["diversity"]),
                    }
                data[topo][domain][seed] = gen_data
                loaded += 1

    print(f"Loaded {loaded} runs, {missing} missing")
    print(f"Expected: {len(TOPOS) * len(DOMAINS) * len(seeds)} = {len(TOPOS)}×{len(DOMAINS)}×{len(seeds)}")
    print()

    if loaded == 0:
        print("No data to analyze.")
        return

    # =========================================================================
    # 1. Summary table: mean diversity at each checkpoint
    # =========================================================================
    print("=" * 80)
    print("DIVERSITY SUMMARY (mean ± SE across 30 seeds)")
    print("=" * 80)

    for domain in DOMAINS:
        print(f"\n--- {domain.upper()} ---")
        header = f"{'Topology':<12} {'β₁':>3} {'nodes':>5} {'edges':>5}"
        for gen in CHECKPOINTS:
            header += f" {'gen'+str(gen):>10}"
        print(header)
        print("-" * len(header))

        for topo in TOPOS:
            nodes, edges, beta1 = TOPOS[topo]
            row = f"{topo:<12} {beta1:>3} {nodes:>5} {edges:>5}"
            for gen in CHECKPOINTS:
                vals = []
                for seed in seeds:
                    if seed in data[topo][domain] and gen in data[topo][domain][seed]:
                        vals.append(data[topo][domain][seed][gen]["diversity"])
                if vals:
                    row += f" {mean(vals):>6.4f}±{stderr(vals):.4f}"
                else:
                    row += f" {'N/A':>10}"
            print(row)

    # =========================================================================
    # 2. Effect sizes: star vs cycle at each checkpoint
    # =========================================================================
    print("\n" + "=" * 80)
    print("EFFECT SIZE: star vs cycle (Cohen's d, η², Mann-Whitney p)")
    print("Positive d means cycle has HIGHER diversity than star")
    print("=" * 80)

    for size_label, (star_name, cycle_name) in [("small (7 edges)", ("star8", "cycle7")),
                                                  ("large (11 edges)", ("star12", "cycle11"))]:
        print(f"\n--- {size_label} ---")
        header = f"{'Domain':<6} {'Gen':>5} {'d(div)':>8} {'η²(div)':>8} {'p(div)':>8} {'d(fit)':>8} {'η²(fit)':>8}"
        print(header)
        print("-" * len(header))

        for domain in DOMAINS:
            for gen in CHECKPOINTS:
                star_div = []
                cycle_div = []
                star_fit = []
                cycle_fit = []
                for seed in seeds:
                    if seed in data[star_name][domain] and gen in data[star_name][domain][seed]:
                        star_div.append(data[star_name][domain][seed][gen]["diversity"])
                        star_fit.append(data[star_name][domain][seed][gen]["meanFitness"])
                    if seed in data[cycle_name][domain] and gen in data[cycle_name][domain][seed]:
                        cycle_div.append(data[cycle_name][domain][seed][gen]["diversity"])
                        cycle_fit.append(data[cycle_name][domain][seed][gen]["meanFitness"])

                if not star_div or not cycle_div:
                    continue

                d_div = cohens_d(cycle_div, star_div)
                eta2_div = eta_squared_oneway([star_div, cycle_div])
                _, p_div = mann_whitney_u(star_div, cycle_div)

                d_fit = cohens_d(cycle_fit, star_fit)
                eta2_fit = eta_squared_oneway([star_fit, cycle_fit])

                sig = "***" if p_div < 0.001 else "**" if p_div < 0.01 else "*" if p_div < 0.05 else ""
                print(f"{domain:<6} {gen:>5} {d_div:>8.3f} {eta2_div:>8.3f} {p_div:>7.4f}{sig} {d_fit:>8.3f} {eta2_fit:>8.3f}")

    # =========================================================================
    # 3. Temporal dynamics: η² over time for diversity
    # =========================================================================
    print("\n" + "=" * 80)
    print("TEMPORAL DYNAMICS: η²(diversity) over generations")
    print("=" * 80)

    all_gens = sorted(set(
        gen for topo in data for domain in data[topo]
        for seed in data[topo][domain] for gen in data[topo][domain][seed]
    ))
    sample_gens = [g for g in all_gens if g % 10 == 0 or g in [1, 5]]

    for size_label, (star_name, cycle_name) in [("small (7 edges)", ("star8", "cycle7")),
                                                  ("large (11 edges)", ("star12", "cycle11"))]:
        print(f"\n--- {size_label} ---")
        header = f"{'Gen':>5}"
        for domain in DOMAINS:
            header += f" {'η²_'+domain:>10}"
        print(header)
        print("-" * len(header))

        for gen in sample_gens:
            row = f"{gen:>5}"
            for domain in DOMAINS:
                star_vals = [data[star_name][domain][seed][gen]["diversity"]
                            for seed in seeds
                            if seed in data[star_name][domain] and gen in data[star_name][domain][seed]]
                cycle_vals = [data[cycle_name][domain][seed][gen]["diversity"]
                             for seed in seeds
                             if seed in data[cycle_name][domain] and gen in data[cycle_name][domain][seed]]
                if star_vals and cycle_vals:
                    eta2 = eta_squared_oneway([star_vals, cycle_vals])
                    row += f" {eta2:>10.4f}"
                else:
                    row += f" {'N/A':>10}"
            print(row)

    # =========================================================================
    # 4. Node count confound check
    # =========================================================================
    print("\n" + "=" * 80)
    print("CONFOUND CHECK: star has more nodes than cycle")
    print("star8 has 8 nodes, cycle7 has 7 nodes (1 extra node)")
    print("star12 has 12 nodes, cycle11 has 11 nodes (1 extra node)")
    print("Both pairs have identical edge count — this is the key control")
    print("=" * 80)
    print()
    print("If star (more nodes, β₁=0) shows LOWER diversity than cycle (fewer nodes, β₁=1),")
    print("then the effect cannot be explained by node count — it must be topology (β₁).")
    print("Star having more total population (more nodes × same pop/island) should give it")
    print("a diversity ADVANTAGE from population size alone. If cycle still wins, β₁ is causal.")
    print()

    for domain in DOMAINS:
        for size_label, (star_name, cycle_name) in [("7-edge", ("star8", "cycle7")),
                                                      ("11-edge", ("star12", "cycle11"))]:
            star_vals = [data[star_name][domain][seed][500]["diversity"]
                        for seed in seeds
                        if seed in data[star_name][domain] and 500 in data[star_name][domain][seed]]
            cycle_vals = [data[cycle_name][domain][seed][500]["diversity"]
                         for seed in seeds
                         if seed in data[cycle_name][domain] and 500 in data[cycle_name][domain][seed]]
            if star_vals and cycle_vals:
                diff = mean(cycle_vals) - mean(star_vals)
                direction = "CYCLE > STAR (β₁ wins)" if diff > 0 else "STAR > CYCLE (β₁ loses)"
                _, p = mann_whitney_u(star_vals, cycle_vals)
                sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."
                print(f"  {domain} {size_label}: {direction} (Δ={diff:+.4f}, p={p:.4f} {sig})")

    # =========================================================================
    # 5. Size comparison: does effect scale with n?
    # =========================================================================
    print("\n" + "=" * 80)
    print("SIZE SCALING: does the β₁ effect grow with network size?")
    print("=" * 80)

    for domain in DOMAINS:
        print(f"\n{domain}: Cohen's d at gen 500:")
        for size_label, (star_name, cycle_name) in [("  7-edge (8/7)", ("star8", "cycle7")),
                                                      ("  11-edge (12/11)", ("star12", "cycle11"))]:
            star_vals = [data[star_name][domain][seed][500]["diversity"]
                        for seed in seeds
                        if seed in data[star_name][domain] and 500 in data[star_name][domain][seed]]
            cycle_vals = [data[cycle_name][domain][seed][500]["diversity"]
                         for seed in seeds
                         if seed in data[cycle_name][domain] and 500 in data[cycle_name][domain][seed]]
            if star_vals and cycle_vals:
                d = cohens_d(cycle_vals, star_vals)
                print(f"{size_label}: d={d:.3f}")


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/..")
    main()
