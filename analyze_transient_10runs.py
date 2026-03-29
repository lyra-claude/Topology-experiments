#!/usr/bin/env python3
"""
Comprehensive transient analysis for OneMax 10-run study.
80 simulations: 10 seeds × 8 topologies × 100 generations.
"""

import csv
import os
import sys
import math
from collections import defaultdict

# ─── Configuration ───────────────────────────────────────────────────────────

RESULTS_DIR = "results/transient_onemax_10runs"
TOPOLOGIES = ["disconnected", "ring", "star", "complete", "hypercube", "barbell", "watts-strogatz", "random-regular"]
SEEDS = [42, 137, 2718, 314, 1618, 271, 577, 997, 1729, 4242]
CHECKPOINTS = [0, 10, 20, 30, 40, 50]  # Early transient window
ALL_GENS = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

# Lambda_2 values (algebraic connectivity)
LAMBDA2 = {
    "disconnected": 0.0,
    "barbell": 0.07,
    "random-regular": 0.27,
    "ring": 0.59,
    "star": 1.0,
    "watts-strogatz": 1.5,
    "hypercube": 2.0,
    "complete": 8.0,
}

# ─── Data Loading ────────────────────────────────────────────────────────────

def load_csv(filepath):
    """Load a single CSV file -> dict of gen -> {meanFitness, bestFitness, diversity}"""
    data = {}
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            gen = int(row['generation'])
            data[gen] = {
                'meanFitness': float(row['meanFitness']),
                'bestFitness': float(row['bestFitness']),
                'diversity': float(row['diversity']),
            }
    return data

def load_all():
    """Load all results -> dict of topo -> list of run data."""
    all_data = defaultdict(list)
    for topo in TOPOLOGIES:
        for seed in SEEDS:
            filepath = os.path.join(RESULTS_DIR, f"{topo}_seed{seed}.csv")
            if not os.path.exists(filepath):
                print(f"WARNING: Missing {filepath}", file=sys.stderr)
                continue
            all_data[topo].append(load_csv(filepath))
    return all_data

# ─── Statistics Helpers ──────────────────────────────────────────────────────

def mean(xs):
    return sum(xs) / len(xs)

def std(xs):
    m = mean(xs)
    return math.sqrt(sum((x - m)**2 for x in xs) / (len(xs) - 1)) if len(xs) > 1 else 0.0

def se(xs):
    return std(xs) / math.sqrt(len(xs)) if len(xs) > 0 else 0.0

def ci95(xs):
    """95% CI half-width using t-distribution approximation for n=10."""
    # t_{0.025, 9} = 2.262
    t_val = 2.262
    return t_val * se(xs)

def spearman_rank(xs, ys):
    """Compute Spearman rank correlation and approximate p-value."""
    n = len(xs)
    if n < 3:
        return 0.0, 1.0

    # Rank the values
    def rank_data(vals):
        indexed = sorted(enumerate(vals), key=lambda x: x[1])
        ranks = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j < n - 1 and indexed[j+1][1] == indexed[j][1]:
                j += 1
            avg_rank = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                ranks[indexed[k][0]] = avg_rank
            i = j + 1
        return ranks

    rx = rank_data(xs)
    ry = rank_data(ys)

    # Pearson correlation on ranks
    mx = mean(rx)
    my = mean(ry)
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    den_x = math.sqrt(sum((rx[i] - mx)**2 for i in range(n)))
    den_y = math.sqrt(sum((ry[i] - my)**2 for i in range(n)))

    if den_x == 0 or den_y == 0:
        return 0.0, 1.0

    rho = num / (den_x * den_y)

    # Approximate p-value using t-distribution
    if abs(rho) >= 1.0:
        return rho, 0.0
    t_stat = rho * math.sqrt((n - 2) / (1 - rho**2))
    # Approximate two-tailed p-value using normal distribution for n >= 8
    p_val = 2 * (1 - normal_cdf(abs(t_stat)))

    return rho, p_val

def normal_cdf(x):
    """Approximate standard normal CDF."""
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))

def t_test_independent(xs, ys):
    """Two-sample independent t-test (Welch's t-test)."""
    nx, ny = len(xs), len(ys)
    mx, my = mean(xs), mean(ys)
    vx = sum((x - mx)**2 for x in xs) / (nx - 1) if nx > 1 else 0
    vy = sum((y - my)**2 for y in ys) / (ny - 1) if ny > 1 else 0

    se_diff = math.sqrt(vx/nx + vy/ny)
    if se_diff == 0:
        return 0.0, 1.0

    t_stat = (mx - my) / se_diff

    # Welch-Satterthwaite degrees of freedom
    num = (vx/nx + vy/ny)**2
    den = (vx/nx)**2/(nx-1) + (vy/ny)**2/(ny-1)
    df = num / den if den > 0 else nx + ny - 2

    # Approximate p-value using normal for df >= 10
    p_val = 2 * (1 - normal_cdf(abs(t_stat)))

    return t_stat, p_val

def kruskal_wallis(groups):
    """Kruskal-Wallis H test."""
    # Combine all values with group labels
    all_vals = []
    for g_idx, g in enumerate(groups):
        for v in g:
            all_vals.append((v, g_idx))

    N = len(all_vals)
    # Rank all values
    sorted_vals = sorted(all_vals, key=lambda x: x[0])
    ranks = [0.0] * N
    i = 0
    while i < N:
        j = i
        while j < N - 1 and sorted_vals[j+1][0] == sorted_vals[j][0]:
            j += 1
        avg_rank = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[k] = avg_rank
        i = j + 1

    # Compute group rank sums
    group_ranks = defaultdict(list)
    for idx, (val, g_idx) in enumerate(sorted_vals):
        group_ranks[g_idx].append(ranks[idx])

    # H statistic
    H = 0.0
    for g_idx in range(len(groups)):
        ni = len(group_ranks[g_idx])
        if ni == 0:
            continue
        Ri = sum(group_ranks[g_idx])
        H += Ri**2 / ni

    H = (12 / (N * (N + 1))) * H - 3 * (N + 1)

    k = len(groups)
    df = k - 1

    # Approximate p-value using chi-squared (normal approximation for large df)
    # For df=7, use Wilson-Hilferty transformation
    z = (H / df)**(1/3) - (1 - 2/(9*df))
    z /= math.sqrt(2/(9*df))
    p_val = 1 - normal_cdf(z) if z > 0 else 1.0

    return H, p_val, df

# ─── Analysis Functions ──────────────────────────────────────────────────────

def compute_gen_to_optimal(runs):
    """For each run, find first generation where bestFitness >= 1.0."""
    results = []
    for run in runs:
        found = None
        for gen in sorted(run.keys()):
            if run[gen]['bestFitness'] >= 1.0:
                found = gen
                break
        results.append(found if found is not None else float('inf'))
    return results

def get_values_at_gen(runs, gen, metric):
    """Get metric values at a specific generation across all runs."""
    return [run[gen][metric] for run in runs if gen in run]

def print_section(title):
    print(f"\n{'='*80}")
    print(f"  {title}")
    print(f"{'='*80}\n")

# ─── Main Analysis ───────────────────────────────────────────────────────────

def main():
    data = load_all()

    # Verify we have all data
    for topo in TOPOLOGIES:
        n = len(data[topo])
        if n != 10:
            print(f"WARNING: {topo} has {n} runs instead of 10", file=sys.stderr)

    print("=" * 80)
    print("  ONEMAX TRANSIENT STUDY: 10-RUN ANALYSIS")
    print("  8 topologies × 10 seeds × 100 generations")
    print("=" * 80)

    # ─── Table 1: Mean Fitness at Each Checkpoint ────────────────────────────
    print_section("TABLE 1: Mean Fitness at Each Checkpoint (mean ± SE, n=10)")

    header = f"{'Topology':<18}" + "".join(f"{'Gen '+str(g):>14}" for g in CHECKPOINTS)
    print(header)
    print("-" * len(header))

    fitness_at_gen = {}  # (topo, gen) -> list of values
    for topo in TOPOLOGIES:
        row = f"{topo:<18}"
        for gen in CHECKPOINTS:
            vals = get_values_at_gen(data[topo], gen, 'meanFitness')
            fitness_at_gen[(topo, gen)] = vals
            m, s = mean(vals), se(vals)
            row += f"  {m:.4f}±{s:.4f}"
        print(row)

    # ─── Table 2: Diversity at Each Checkpoint ───────────────────────────────
    print_section("TABLE 2: Diversity at Each Checkpoint (mean ± SE, n=10)")

    header = f"{'Topology':<18}" + "".join(f"{'Gen '+str(g):>14}" for g in CHECKPOINTS)
    print(header)
    print("-" * len(header))

    diversity_at_gen = {}
    for topo in TOPOLOGIES:
        row = f"{topo:<18}"
        for gen in CHECKPOINTS:
            vals = get_values_at_gen(data[topo], gen, 'diversity')
            diversity_at_gen[(topo, gen)] = vals
            m, s = mean(vals), se(vals)
            row += f"  {m:.4f}±{s:.4f}"
        print(row)

    # ─── Table 3: Generation to Optimal ──────────────────────────────────────
    print_section("TABLE 3: Generation to Reach bestFitness = 1.0 (mean ± SE, n=10)")

    gen_to_opt = {}
    for topo in TOPOLOGIES:
        vals = compute_gen_to_optimal(data[topo])
        gen_to_opt[topo] = vals
        finite = [v for v in vals if v != float('inf')]
        if len(finite) == len(vals):
            print(f"  {topo:<18}  {mean(vals):6.1f} ± {se(vals):5.1f}  (all 10 runs converged)")
        else:
            print(f"  {topo:<18}  {mean(finite):6.1f} ± {se(finite):5.1f}  ({len(finite)}/10 runs converged by gen 100)")

    # ─── Table 4: Full Convergence Curves ────────────────────────────────────
    print_section("TABLE 4: Full Convergence Curves - Mean Fitness (mean across 10 runs)")

    header = f"{'Topology':<18}" + "".join(f"{'Gen '+str(g):>10}" for g in ALL_GENS)
    print(header)
    print("-" * len(header))

    for topo in TOPOLOGIES:
        row = f"{topo:<18}"
        for gen in ALL_GENS:
            vals = get_values_at_gen(data[topo], gen, 'meanFitness')
            row += f"  {mean(vals):.4f}  "
        print(row)

    print_section("TABLE 4b: Full Convergence Curves - Diversity (mean across 10 runs)")

    header = f"{'Topology':<18}" + "".join(f"{'Gen '+str(g):>10}" for g in ALL_GENS)
    print(header)
    print("-" * len(header))

    for topo in TOPOLOGIES:
        row = f"{topo:<18}"
        for gen in ALL_GENS:
            vals = get_values_at_gen(data[topo], gen, 'diversity')
            row += f"  {mean(vals):.4f}  "
        print(row)

    # ─── Spearman Correlations ───────────────────────────────────────────────
    print_section("SPEARMAN RANK CORRELATIONS: λ₂ vs Metrics")

    # Order topologies by lambda_2 for correlation
    topo_order = sorted(TOPOLOGIES, key=lambda t: LAMBDA2[t])
    lambda2_vals = [LAMBDA2[t] for t in topo_order]

    print("--- λ₂ vs DIVERSITY at each checkpoint ---")
    print(f"{'Generation':<12}  {'ρ':>8}  {'p-value':>10}  {'Significant?':>14}")
    print("-" * 50)

    for gen in CHECKPOINTS:
        div_means = [mean(diversity_at_gen[(t, gen)]) for t in topo_order]
        rho, p = spearman_rank(lambda2_vals, div_means)
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."
        print(f"  Gen {gen:<6}  {rho:>8.4f}  {p:>10.6f}  {sig:>14}")

    print("\n--- λ₂ vs MEAN FITNESS at each checkpoint ---")
    print(f"{'Generation':<12}  {'ρ':>8}  {'p-value':>10}  {'Significant?':>14}")
    print("-" * 50)

    for gen in CHECKPOINTS:
        fit_means = [mean(fitness_at_gen[(t, gen)]) for t in topo_order]
        rho, p = spearman_rank(lambda2_vals, fit_means)
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."
        print(f"  Gen {gen:<6}  {rho:>8.4f}  {p:>10.6f}  {sig:>14}")

    print("\n--- λ₂ vs GENERATION TO OPTIMAL ---")
    gen_opt_means = []
    for t in topo_order:
        vals = gen_to_opt[t]
        finite = [v for v in vals if v != float('inf')]
        gen_opt_means.append(mean(finite) if finite else 100)

    rho, p = spearman_rank(lambda2_vals, gen_opt_means)
    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."
    print(f"  ρ = {rho:.4f}, p = {p:.6f}  {sig}")
    print(f"  (Higher λ₂ = {'faster' if rho < 0 else 'slower'} convergence)")

    # ─── Pairwise t-tests: Complete vs Disconnected ──────────────────────────
    print_section("PAIRWISE t-TESTS: Complete vs Disconnected")

    print("--- Diversity ---")
    print(f"{'Generation':<12}  {'Complete':>10}  {'Disconn':>10}  {'t-stat':>8}  {'p-value':>10}  {'Sig?':>6}")
    print("-" * 65)

    for gen in [10, 20, 30, 40, 50]:
        comp_vals = diversity_at_gen[("complete", gen)]
        disc_vals = diversity_at_gen[("disconnected", gen)]
        t, p = t_test_independent(comp_vals, disc_vals)
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."
        print(f"  Gen {gen:<6}  {mean(comp_vals):>10.4f}  {mean(disc_vals):>10.4f}  {t:>8.3f}  {p:>10.6f}  {sig:>6}")

    print("\n--- Mean Fitness ---")
    print(f"{'Generation':<12}  {'Complete':>10}  {'Disconn':>10}  {'t-stat':>8}  {'p-value':>10}  {'Sig?':>6}")
    print("-" * 65)

    for gen in [10, 20, 30, 40, 50]:
        comp_vals = fitness_at_gen[("complete", gen)]
        disc_vals = fitness_at_gen[("disconnected", gen)]
        t, p = t_test_independent(comp_vals, disc_vals)
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."
        print(f"  Gen {gen:<6}  {mean(comp_vals):>10.4f}  {mean(disc_vals):>10.4f}  {t:>8.3f}  {p:>10.6f}  {sig:>6}")

    # ─── Kruskal-Wallis ──────────────────────────────────────────────────────
    print_section("KRUSKAL-WALLIS: Does Topology Affect Diversity?")

    for gen in [10, 20, 30]:
        groups = [diversity_at_gen[(t, gen)] for t in TOPOLOGIES]
        H, p, df = kruskal_wallis(groups)
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."
        print(f"  Gen {gen}: H({df}) = {H:.3f}, p = {p:.6f}  {sig}")

    # ─── Anomaly Investigation ───────────────────────────────────────────────
    print_section("ANOMALY INVESTIGATION")

    # Star anomaly
    print("--- STAR: Does it behave anomalously? ---")
    print("  Star's λ₂ = 1.0 (moderate connectivity)")
    print("  Expected: moderate fitness, moderate diversity")
    print()

    for gen in [10, 20, 30]:
        star_fit = fitness_at_gen[("star", gen)]
        star_div = diversity_at_gen[("star", gen)]

        # Rank among all topologies
        topo_fits = [(t, mean(fitness_at_gen[(t, gen)])) for t in TOPOLOGIES]
        topo_divs = [(t, mean(diversity_at_gen[(t, gen)])) for t in TOPOLOGIES]

        fit_rank = sorted(topo_fits, key=lambda x: -x[1])
        div_rank = sorted(topo_divs, key=lambda x: -x[1])

        star_fit_rank = [i for i, (t, _) in enumerate(fit_rank) if t == "star"][0] + 1
        star_div_rank = [i for i, (t, _) in enumerate(div_rank) if t == "star"][0] + 1

        print(f"  Gen {gen}: fitness={mean(star_fit):.4f} (rank {star_fit_rank}/8), "
              f"diversity={mean(star_div):.4f} (rank {star_div_rank}/8)")

    # Barbell anomaly
    print("\n--- BARBELL: Does it outperform its λ₂ prediction? ---")
    print("  Barbell's λ₂ = 0.07 (very low — near disconnected)")
    print("  If λ₂ is predictive, Barbell should be near disconnected in fitness/diversity.")
    print()

    for gen in [10, 20, 30]:
        barbell_fit = fitness_at_gen[("barbell", gen)]
        disc_fit = fitness_at_gen[("disconnected", gen)]

        topo_fits = [(t, mean(fitness_at_gen[(t, gen)])) for t in TOPOLOGIES]
        fit_rank = sorted(topo_fits, key=lambda x: -x[1])
        barbell_rank = [i for i, (t, _) in enumerate(fit_rank) if t == "barbell"][0] + 1
        disc_rank = [i for i, (t, _) in enumerate(fit_rank) if t == "disconnected"][0] + 1

        t_stat, p = t_test_independent(barbell_fit, disc_fit)
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."

        print(f"  Gen {gen}: barbell={mean(barbell_fit):.4f} (rank {barbell_rank}), "
              f"disconn={mean(disc_fit):.4f} (rank {disc_rank}), "
              f"t={t_stat:.3f}, p={p:.6f} {sig}")

    # ─── Topology Rankings at Gen 20 ─────────────────────────────────────────
    print_section("TOPOLOGY RANKINGS AT GEN 20 (with 95% CI)")

    print("--- By Mean Fitness (descending) ---")
    rankings = []
    for topo in TOPOLOGIES:
        vals = fitness_at_gen[(topo, 20)]
        m = mean(vals)
        c = ci95(vals)
        rankings.append((topo, m, c, LAMBDA2[topo]))

    rankings.sort(key=lambda x: -x[1])
    print(f"  {'Rank':<6}  {'Topology':<18}  {'Mean Fitness':>14}  {'95% CI':>16}  {'λ₂':>6}")
    print("  " + "-" * 65)
    for i, (topo, m, c, l2) in enumerate(rankings):
        print(f"  {i+1:<6}  {topo:<18}  {m:>14.4f}  [{m-c:.4f}, {m+c:.4f}]  {l2:>6.2f}")

    print("\n--- By Diversity at Gen 20 (descending) ---")
    rankings_div = []
    for topo in TOPOLOGIES:
        vals = diversity_at_gen[(topo, 20)]
        m = mean(vals)
        c = ci95(vals)
        rankings_div.append((topo, m, c, LAMBDA2[topo]))

    rankings_div.sort(key=lambda x: -x[1])
    print(f"  {'Rank':<6}  {'Topology':<18}  {'Mean Diversity':>14}  {'95% CI':>16}  {'λ₂':>6}")
    print("  " + "-" * 65)
    for i, (topo, m, c, l2) in enumerate(rankings_div):
        print(f"  {i+1:<6}  {topo:<18}  {m:>14.4f}  [{m-c:.4f}, {m+c:.4f}]  {l2:>6.2f}")

    # ─── Lambda_2 vs Rank Check ──────────────────────────────────────────────
    print_section("λ₂ vs OBSERVED RANKING COMPARISON")

    print("  Expected ranking by λ₂ (most → least connected):")
    l2_rank = sorted(TOPOLOGIES, key=lambda t: -LAMBDA2[t])
    for i, t in enumerate(l2_rank):
        print(f"    {i+1}. {t} (λ₂ = {LAMBDA2[t]})")

    print("\n  Observed ranking by mean fitness at gen 20:")
    obs_rank = [x[0] for x in rankings]
    for i, t in enumerate(obs_rank):
        l2_pos = l2_rank.index(t) + 1
        delta = l2_pos - (i + 1)
        arrow = f"({'=' if delta == 0 else '+' + str(delta) if delta > 0 else str(delta)})"
        print(f"    {i+1}. {t} (λ₂ = {LAMBDA2[t]}) {arrow}")

    # ─── Summary Statistics ──────────────────────────────────────────────────
    print_section("SUMMARY: KEY FINDINGS")

    # Overall Spearman at gen 20
    div_means_20 = [mean(diversity_at_gen[(t, 20)]) for t in topo_order]
    rho_div_20, p_div_20 = spearman_rank(lambda2_vals, div_means_20)

    fit_means_20 = [mean(fitness_at_gen[(t, 20)]) for t in topo_order]
    rho_fit_20, p_fit_20 = spearman_rank(lambda2_vals, fit_means_20)

    print(f"  1. λ₂ vs Diversity at gen 20: ρ = {rho_div_20:.4f}, p = {p_div_20:.6f}")
    print(f"  2. λ₂ vs Fitness at gen 20: ρ = {rho_fit_20:.4f}, p = {p_fit_20:.6f}")

    # Complete vs Disconnected at gen 20
    t20, p20 = t_test_independent(
        diversity_at_gen[("complete", 20)],
        diversity_at_gen[("disconnected", 20)]
    )
    print(f"  3. Complete vs Disconnected diversity at gen 20: t = {t20:.3f}, p = {p20:.6f}")

    # Kruskal-Wallis at gen 20
    groups_20 = [diversity_at_gen[(t, 20)] for t in TOPOLOGIES]
    H20, pH20, df20 = kruskal_wallis(groups_20)
    print(f"  4. Kruskal-Wallis on diversity at gen 20: H({df20}) = {H20:.3f}, p = {pH20:.6f}")

    print()

if __name__ == "__main__":
    main()
