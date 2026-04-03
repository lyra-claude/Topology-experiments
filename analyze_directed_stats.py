#!/usr/bin/env python3
"""
Statistical analysis of directed cycle experiment with p-values.
Extends analyze_directed_full.py with:
  - Pearson correlation p-values (scipy.stats.pearsonr)
  - One-way ANOVA (scipy.stats.f_oneway)
  - Kruskal-Wallis test (scipy.stats.kruskal)
"""

import os
import csv
import math
from collections import defaultdict
from scipy import stats

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


def stderr(xs):
    if len(xs) < 2:
        return 0.0
    m = mean(xs)
    s = math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))
    return s / math.sqrt(len(xs))


def main():
    # Load all results
    results = defaultdict(list)
    for fname in sorted(os.listdir(RESULTS_DIR)):
        if not fname.endswith(".csv"):
            continue
        topo = fname.rsplit("_seed", 1)[0]
        rows = load_csv(os.path.join(RESULTS_DIR, fname))
        if rows:
            results[topo].append(rows)

    topos_sorted = sorted(results.keys(), key=lambda t: CYCLE_COUNTS.get(t, 999))

    print("DIRECTED CYCLE EXPERIMENT — STATISTICAL ANALYSIS")
    print("=" * 110)
    print(f"Topologies: {len(topos_sorted)}, Seeds per topology: {len(results[topos_sorted[0]])}")
    print()

    checkpoints = [10, 20, 30, 50, 100, 200, 500]

    # Precompute: for each generation, extract per-topology diversity and fitness
    def get_values_at_gen(gen, metric="diversity"):
        """Return dict: topo -> [values across seeds]."""
        data = {}
        for topo in topos_sorted:
            vals = []
            for run in results[topo]:
                for row in run:
                    if int(row["generation"]) == gen:
                        vals.append(float(row[metric]))
                        break
            data[topo] = vals
        return data

    # ============================================================
    # Table 1: Diversity statistics with p-values
    # ============================================================
    print("TABLE 1: DIVERSITY — Correlation and Group Tests by Generation")
    print("-" * 110)
    print(f"{'Gen':>5}  {'r_pearson':>10}  {'p_pearson':>10}  {'r²':>8}  "
          f"{'F_anova':>10}  {'p_anova':>10}  {'H_kruskal':>10}  {'p_kruskal':>10}  {'sig':>5}")
    print("-" * 110)

    for gen in checkpoints:
        div_data = get_values_at_gen(gen, "diversity")

        # Pearson: cycle count vs mean diversity (8 data points, one per topology)
        cycle_counts_list = []
        mean_divs_list = []
        for topo in topos_sorted:
            if div_data[topo]:
                cycle_counts_list.append(CYCLE_COUNTS[topo])
                mean_divs_list.append(mean(div_data[topo]))

        if len(cycle_counts_list) >= 3:
            r, p_pearson = stats.pearsonr(cycle_counts_list, mean_divs_list)
            r_sq = r * r
        else:
            r, p_pearson, r_sq = 0, 1, 0

        # ANOVA: groups are topologies, values are individual seed diversities
        groups = [div_data[topo] for topo in topos_sorted if div_data[topo]]
        if len(groups) >= 2 and all(len(g) >= 2 for g in groups):
            F_anova, p_anova = stats.f_oneway(*groups)
        else:
            F_anova, p_anova = 0, 1

        # Kruskal-Wallis
        if len(groups) >= 2 and all(len(g) >= 2 for g in groups):
            H_kruskal, p_kruskal = stats.kruskal(*groups)
        else:
            H_kruskal, p_kruskal = 0, 1

        sig = ""
        if p_pearson < 0.001:
            sig = "***"
        elif p_pearson < 0.01:
            sig = "**"
        elif p_pearson < 0.05:
            sig = "*"

        print(f"{gen:>5}  {r:>10.4f}  {p_pearson:>10.4e}  {r_sq:>8.4f}  "
              f"{F_anova:>10.2f}  {p_anova:>10.4e}  {H_kruskal:>10.2f}  {p_kruskal:>10.4e}  {sig:>5}")

    print()

    # ============================================================
    # Table 2: Fitness statistics with p-values
    # ============================================================
    print("TABLE 2: FITNESS — Correlation and Group Tests by Generation")
    print("-" * 110)
    print(f"{'Gen':>5}  {'r_pearson':>10}  {'p_pearson':>10}  {'r²':>8}  "
          f"{'F_anova':>10}  {'p_anova':>10}  {'H_kruskal':>10}  {'p_kruskal':>10}  {'sig':>5}")
    print("-" * 110)

    for gen in checkpoints:
        fit_data = get_values_at_gen(gen, "meanFitness")

        cycle_counts_list = []
        mean_fits_list = []
        for topo in topos_sorted:
            if fit_data[topo]:
                cycle_counts_list.append(CYCLE_COUNTS[topo])
                mean_fits_list.append(mean(fit_data[topo]))

        if len(cycle_counts_list) >= 3:
            r, p_pearson = stats.pearsonr(cycle_counts_list, mean_fits_list)
            r_sq = r * r
        else:
            r, p_pearson, r_sq = 0, 1, 0

        groups = [fit_data[topo] for topo in topos_sorted if fit_data[topo]]
        if len(groups) >= 2 and all(len(g) >= 2 for g in groups):
            F_anova, p_anova = stats.f_oneway(*groups)
        else:
            F_anova, p_anova = 0, 1

        if len(groups) >= 2 and all(len(g) >= 2 for g in groups):
            H_kruskal, p_kruskal = stats.kruskal(*groups)
        else:
            H_kruskal, p_kruskal = 0, 1

        sig = ""
        if p_pearson < 0.001:
            sig = "***"
        elif p_pearson < 0.01:
            sig = "**"
        elif p_pearson < 0.05:
            sig = "*"

        print(f"{gen:>5}  {r:>10.4f}  {p_pearson:>10.4e}  {r_sq:>8.4f}  "
              f"{F_anova:>10.2f}  {p_anova:>10.4e}  {H_kruskal:>10.2f}  {p_kruskal:>10.4e}  {sig:>5}")

    print()

    # ============================================================
    # Table 3: Effect sizes (eta-squared) with ANOVA p-values
    # ============================================================
    print("TABLE 3: EFFECT SIZES — Eta-squared with ANOVA p-values")
    print("-" * 70)
    print(f"{'Gen':>5}  {'eta²(div)':>10}  {'p_anova(div)':>13}  {'eta²(fit)':>10}  {'p_anova(fit)':>13}")
    print("-" * 70)

    for gen in checkpoints:
        for metric, label in [("diversity", "div"), ("meanFitness", "fit")]:
            data = get_values_at_gen(gen, metric)
            all_vals = []
            group_means = []
            groups = []
            for topo in topos_sorted:
                vals = data[topo]
                all_vals.extend(vals)
                groups.append(vals)
                if vals:
                    group_means.append((mean(vals), len(vals)))

            if all_vals and len(group_means) > 1:
                grand_mean = mean(all_vals)
                ss_between = sum(n * (m - grand_mean)**2 for m, n in group_means)
                ss_total = sum((x - grand_mean)**2 for x in all_vals)
                eta2 = ss_between / ss_total if ss_total > 0 else 0
            else:
                eta2 = 0

            valid_groups = [g for g in groups if len(g) >= 2]
            if len(valid_groups) >= 2:
                _, p_anova = stats.f_oneway(*valid_groups)
            else:
                p_anova = 1.0

            if label == "div":
                eta2_d, p_d = eta2, p_anova
            else:
                eta2_f, p_f = eta2, p_anova

        print(f"{gen:>5}  {eta2_d:>10.4f}  {p_d:>13.4e}  {eta2_f:>10.4f}  {p_f:>13.4e}")

    print()

    # ============================================================
    # Post-hoc: pairwise comparisons at gen 30 (peak signal)
    # ============================================================
    print("POST-HOC: Pairwise Mann-Whitney U at gen 30 (diversity)")
    print("(Bonferroni-corrected alpha = 0.05/28 = 0.0018)")
    print("-" * 60)

    div30 = get_values_at_gen(30, "diversity")
    n_comparisons = len(topos_sorted) * (len(topos_sorted) - 1) // 2
    alpha_corrected = 0.05 / n_comparisons

    sig_pairs = []
    for i, t1 in enumerate(topos_sorted):
        for j, t2 in enumerate(topos_sorted):
            if j <= i:
                continue
            if div30[t1] and div30[t2]:
                U, p = stats.mannwhitneyu(div30[t1], div30[t2], alternative='two-sided')
                if p < alpha_corrected:
                    sig_pairs.append((t1, t2, U, p))

    if sig_pairs:
        for t1, t2, U, p in sig_pairs:
            c1, c2 = CYCLE_COUNTS[t1], CYCLE_COUNTS[t2]
            print(f"  {t1}({c1}) vs {t2}({c2}): U={U:.0f}, p={p:.4e} ***")
    else:
        print("  No pairwise comparisons significant after Bonferroni correction.")

    print()
    print("=" * 110)
    print("NOTES:")
    print("  - Pearson r computed on 8 topology-level means (N=8, df=6)")
    print("  - ANOVA/Kruskal-Wallis computed on all 240 individual seed observations")
    print("  - *p<0.05  **p<0.01  ***p<0.001 (uncorrected, for Pearson)")
    print("  - Eta-squared: 0.01=small, 0.06=medium, 0.14=large (Cohen)")


if __name__ == "__main__":
    main()
