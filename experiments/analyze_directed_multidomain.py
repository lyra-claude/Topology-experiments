#!/usr/bin/env python3
"""
Analyze multi-domain directed cycle experiment.
8 topologies x 4 NK domains x 30 seeds = 960 runs.

Key questions:
1. Does directed cycle count (kappa) predict diversity across NK landscapes?
2. How does landscape ruggedness (K) modulate the topology effect?
3. Does the temporal inversion from OneMax replicate on NK?
   (OneMax: more cycles = LESS diversity at steady state)
"""

import os
import csv
import math
from collections import defaultdict

RESULTS_DIR = "results/directed_multidomain"

DOMAINS = ["nk0", "nk2", "nk4", "nk6"]

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


def load_domain_data(domain):
    """Load all results for a domain. Returns {topology: {seed: [rows]}}."""
    domain_dir = os.path.join(RESULTS_DIR, domain)
    if not os.path.isdir(domain_dir):
        return {}
    data = defaultdict(dict)
    for fname in os.listdir(domain_dir):
        if not fname.endswith(".csv"):
            continue
        # Parse: topology_seed123.csv
        parts = fname.replace(".csv", "").rsplit("_seed", 1)
        if len(parts) != 2:
            continue
        topo, seed = parts[0], int(parts[1])
        rows = load_csv(os.path.join(domain_dir, fname))
        data[topo][seed] = rows
    return data


def get_metric_at_gen(rows, gen, metric="diversity"):
    """Get metric value at or near a specific generation."""
    for row in rows:
        if int(row["generation"]) == gen:
            return float(row[metric])
    # Fallback: get closest
    closest = min(rows, key=lambda r: abs(int(r["generation"]) - gen))
    return float(closest[metric])


def analyze_domain(domain, data):
    """Analyze a single domain's results."""
    print(f"\n{'='*70}")
    print(f"  DOMAIN: {domain}")
    print(f"{'='*70}")

    topos = sorted(data.keys(), key=lambda t: CYCLE_COUNTS.get(t, 0))
    n_topos = len(topos)
    total_runs = sum(len(seeds) for seeds in data.values())
    print(f"  Topologies: {n_topos}, Total runs: {total_runs}")

    for gen in CHECKPOINTS:
        print(f"\n  --- Generation {gen} ---")

        # Collect per-topology mean diversity
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
        print(f"  {'Topology':<18} {'kappa':>6} {'Diversity':>12} {'SE':>8} {'MeanFit':>12} {'SE':>8} {'N':>4}")
        print(f"  {'-'*18} {'-'*6} {'-'*12} {'-'*8} {'-'*12} {'-'*8} {'-'*4}")
        for topo in topos:
            kappa = CYCLE_COUNTS.get(topo, -1)
            divs = topo_diversities[topo]
            fits = topo_fitnesses[topo]
            print(f"  {topo:<18} {kappa:>6} {mean(divs):>12.6f} {stderr(divs):>8.6f} "
                  f"{mean(fits):>12.6f} {stderr(fits):>8.6f} {len(divs):>4}")

        # Eta-squared
        eta2_div = eta_squared(groups_div)
        eta2_fit = eta_squared(groups_fit)
        print(f"\n  eta^2 (topology -> diversity): {eta2_div:.4f}")
        print(f"  eta^2 (topology -> fitness):   {eta2_fit:.4f}")

        # Correlation with kappa (using topology means)
        kappas = []
        mean_divs = []
        mean_fits = []
        for topo in topos:
            if topo_diversities[topo]:
                kappas.append(CYCLE_COUNTS.get(topo, 0))
                mean_divs.append(mean(topo_diversities[topo]))
                mean_fits.append(mean(topo_fitnesses[topo]))

        if len(kappas) >= 3:
            # Per-run correlation (all 240 datapoints)
            all_kappas = []
            all_divs = []
            for topo in topos:
                k = CYCLE_COUNTS.get(topo, 0)
                for d in topo_diversities[topo]:
                    all_kappas.append(k)
                    all_divs.append(d)

            r_div = pearson_r(all_kappas, all_divs)
            rho_div = spearman_rho(all_kappas, all_divs)
            r_mean = pearson_r(kappas, mean_divs)
            print(f"  Pearson r(kappa, diversity):  {r_div:+.4f}  (per-run, N={len(all_kappas)})")
            print(f"  Spearman rho(kappa, div):     {rho_div:+.4f}")
            print(f"  Pearson r(kappa, mean_div):   {r_mean:+.4f}  (topology means, N={len(kappas)})")

    return topos, topo_diversities


def cross_domain_comparison(all_results):
    """Compare effect sizes across domains."""
    print(f"\n{'='*70}")
    print(f"  CROSS-DOMAIN COMPARISON")
    print(f"{'='*70}")

    print(f"\n  {'Domain':<8}", end="")
    for gen in CHECKPOINTS:
        print(f"  {'eta2_g'+str(gen):>10}  {'r_g'+str(gen):>10}", end="")
    print()
    print(f"  {'-'*8}", end="")
    for gen in CHECKPOINTS:
        print(f"  {'-'*10}  {'-'*10}", end="")
    print()

    for domain in DOMAINS:
        data = all_results.get(domain, {})
        if not data:
            continue
        topos = sorted(data.keys(), key=lambda t: CYCLE_COUNTS.get(t, 0))
        print(f"  {domain:<8}", end="")

        for gen in CHECKPOINTS:
            groups = []
            all_kappas = []
            all_divs = []
            for topo in topos:
                divs = []
                for seed, rows in data[topo].items():
                    try:
                        d = get_metric_at_gen(rows, gen, "diversity")
                        divs.append(d)
                        all_kappas.append(CYCLE_COUNTS.get(topo, 0))
                        all_divs.append(d)
                    except:
                        continue
                groups.append(divs)

            eta2 = eta_squared(groups) if groups else 0
            r = pearson_r(all_kappas, all_divs) if len(all_kappas) >= 3 else 0
            print(f"  {eta2:>10.4f}  {r:>+10.4f}", end="")
        print()

    # Temporal inversion check
    print(f"\n  TEMPORAL INVERSION CHECK")
    print(f"  (OneMax showed: kappa correlates POSITIVELY with diversity early,")
    print(f"   NEGATIVELY at steady state. Does this replicate on NK?)")
    print()
    for domain in DOMAINS:
        data = all_results.get(domain, {})
        if not data:
            continue
        topos = sorted(data.keys(), key=lambda t: CYCLE_COUNTS.get(t, 0))
        signs = []
        for gen in CHECKPOINTS:
            all_kappas = []
            all_divs = []
            for topo in topos:
                for seed, rows in data[topo].items():
                    try:
                        d = get_metric_at_gen(rows, gen, "diversity")
                        all_kappas.append(CYCLE_COUNTS.get(topo, 0))
                        all_divs.append(d)
                    except:
                        continue
            r = pearson_r(all_kappas, all_divs) if len(all_kappas) >= 3 else 0
            signs.append("+" if r > 0 else "-")
        sign_str = " -> ".join(f"gen{g}:{s}" for g, s in zip(CHECKPOINTS, signs))
        print(f"  {domain}: {sign_str}")


def main():
    print("Multi-Domain Directed Cycle Experiment Analysis")
    print("=" * 70)

    all_results = {}
    for domain in DOMAINS:
        data = load_domain_data(domain)
        if data:
            all_results[domain] = data
            analyze_domain(domain, data)
        else:
            print(f"\n  WARNING: No data found for {domain}")

    if len(all_results) >= 2:
        cross_domain_comparison(all_results)

    print(f"\n{'='*70}")
    print("  DONE")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
