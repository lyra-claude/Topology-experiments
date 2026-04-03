#!/usr/bin/env python3
"""
Analyze full directed cycle experiment: 8 topologies x 30 seeds.
Tests whether simple directed cycle count predicts GA performance
at constant density (n=8, m=16 directed edges).
"""

import os
import csv
import math
from collections import defaultdict

RESULTS_DIR = "results/directed_full"

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


def main():
    # Load all results
    results = defaultdict(list)  # topo -> list of run data
    for fname in sorted(os.listdir(RESULTS_DIR)):
        if not fname.endswith(".csv"):
            continue
        topo = fname.rsplit("_seed", 1)[0]
        rows = load_csv(os.path.join(RESULTS_DIR, fname))
        if rows:
            results[topo].append(rows)

    topos_sorted = sorted(results.keys(), key=lambda t: CYCLE_COUNTS.get(t, 999))

    print("DIRECTED CYCLE EXPERIMENT — FULL ANALYSIS (8 topologies x 30 seeds)")
    print("=" * 90)
    print()

    # === Summary table ===
    print("SUMMARY: FINAL GENERATION (gen 500)")
    print(f"{'Topology':<20} {'Cycles':>7} {'N':>4} {'MeanFit':>9} {'±SE':>7} {'BestFit':>9} {'Diversity':>10} {'±SE':>7}")
    print("-" * 80)

    for topo in topos_sorted:
        runs = results[topo]
        cycles = CYCLE_COUNTS.get(topo, "?")
        final_means = [float(r[-1]["meanFitness"]) for r in runs]
        final_bests = [float(r[-1]["bestFitness"]) for r in runs]
        final_divs = [float(r[-1]["diversity"]) for r in runs]
        n = len(runs)
        print(f"{topo:<20} {cycles:>7} {n:>4} {mean(final_means):>9.6f} {stderr(final_means):>7.4f} "
              f"{mean(final_bests):>9.6f} {mean(final_divs):>10.6f} {stderr(final_divs):>7.4f}")

    print()

    # === Transient analysis (where topology signal is strongest) ===
    print("TRANSIENT ANALYSIS: DIVERSITY BY GENERATION")
    print(f"{'Gen':>5}", end="")
    for topo in topos_sorted:
        cycles = CYCLE_COUNTS.get(topo, "?")
        print(f"  {topo}({cycles})", end="")
    print()
    print("-" * 140)

    checkpoints = [0, 10, 20, 30, 40, 50, 100, 200, 500]
    gen_div_data = {}  # gen -> {topo -> [diversities]}

    for gen in checkpoints:
        print(f"{gen:>5}", end="")
        gen_div_data[gen] = {}
        for topo in topos_sorted:
            runs = results[topo]
            divs = []
            for run in runs:
                for row in run:
                    if int(row["generation"]) == gen:
                        divs.append(float(row["diversity"]))
                        break
            gen_div_data[gen][topo] = divs
            if divs:
                print(f"  {mean(divs):.4f}±{stderr(divs):.4f}     ", end="")
            else:
                print(f"  {'N/A':>14}     ", end="")
        print()

    print()

    # === Correlation: cycle count vs diversity at each generation ===
    print("CORRELATION: Cycle Count vs Mean Diversity (across topologies)")
    print(f"{'Gen':>5} {'r':>8} {'r²':>8} {'Direction':>12}")
    print("-" * 40)

    for gen in checkpoints:
        cycle_counts_list = []
        mean_divs_list = []
        for topo in topos_sorted:
            divs = gen_div_data[gen].get(topo, [])
            if divs:
                cycle_counts_list.append(CYCLE_COUNTS.get(topo, 0))
                mean_divs_list.append(mean(divs))

        if len(cycle_counts_list) >= 3:
            r = pearson_r(cycle_counts_list, mean_divs_list)
            direction = "more cycles → higher div" if r > 0 else "more cycles → lower div"
            print(f"{gen:>5} {r:>8.4f} {r*r:>8.4f} {direction:>12}")

    print()

    # === Correlation: cycle count vs mean fitness at each generation ===
    print("CORRELATION: Cycle Count vs Mean Fitness (across topologies)")
    print(f"{'Gen':>5} {'r':>8} {'r²':>8}")
    print("-" * 25)

    for gen in checkpoints:
        cycle_counts_list = []
        mean_fits_list = []
        for topo in topos_sorted:
            runs = results[topo]
            fits = []
            for run in runs:
                for row in run:
                    if int(row["generation"]) == gen:
                        fits.append(float(row["meanFitness"]))
                        break
            if fits:
                cycle_counts_list.append(CYCLE_COUNTS.get(topo, 0))
                mean_fits_list.append(mean(fits))

        if len(cycle_counts_list) >= 3:
            r = pearson_r(cycle_counts_list, mean_fits_list)
            print(f"{gen:>5} {r:>8.4f} {r*r:>8.4f}")

    print()

    # === ANOVA-like: between-group vs within-group variance ===
    print("EFFECT SIZE: Eta-squared (between-topology variance / total variance)")
    print(f"{'Gen':>5} {'eta²(div)':>10} {'eta²(fit)':>10}")
    print("-" * 30)

    for gen in checkpoints:
        # Diversity
        all_divs = []
        group_means_d = []
        for topo in topos_sorted:
            divs = gen_div_data[gen].get(topo, [])
            all_divs.extend(divs)
            if divs:
                group_means_d.append((mean(divs), len(divs)))

        if all_divs and len(group_means_d) > 1:
            grand_mean_d = mean(all_divs)
            ss_between_d = sum(n * (m - grand_mean_d)**2 for m, n in group_means_d)
            ss_total_d = sum((x - grand_mean_d)**2 for x in all_divs)
            eta2_d = ss_between_d / ss_total_d if ss_total_d > 0 else 0
        else:
            eta2_d = 0

        # Fitness
        all_fits = []
        group_means_f = []
        for topo in topos_sorted:
            runs = results[topo]
            fits = []
            for run in runs:
                for row in run:
                    if int(row["generation"]) == gen:
                        fits.append(float(row["meanFitness"]))
                        break
            all_fits.extend(fits)
            if fits:
                group_means_f.append((mean(fits), len(fits)))

        if all_fits and len(group_means_f) > 1:
            grand_mean_f = mean(all_fits)
            ss_between_f = sum(n * (m - grand_mean_f)**2 for m, n in group_means_f)
            ss_total_f = sum((x - grand_mean_f)**2 for x in all_fits)
            eta2_f = ss_between_f / ss_total_f if ss_total_f > 0 else 0
        else:
            eta2_f = 0

        print(f"{gen:>5} {eta2_d:>10.4f} {eta2_f:>10.4f}")

    print()

    # === Key finding ===
    print("=" * 90)
    print("KEY FINDINGS:")
    print()
    print("The independent variable is simple directed cycle count (0, 0, 3, 10, 14, 20, 29, 47)")
    print("at CONSTANT density (n=8, m=16 for all topologies).")
    print()
    print("If eta² for diversity is high in the transient (gen 10-50), then cycle COUNT")
    print("(not just density/cycle rank) affects GA dynamics. This would be the first")
    print("controlled evidence separating cycles from density.")


if __name__ == "__main__":
    main()
