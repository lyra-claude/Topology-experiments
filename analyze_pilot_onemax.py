#!/usr/bin/env python3
"""Analyze OneMax pilot study results: 3 runs x 8 topologies x 500 generations.

Compares diversity ordering at gen 500 with lambda_2 and computes Spearman rho.
"""

import csv
import os
import statistics

RESULTS_DIR = "results/pilot_onemax"

TOPOLOGIES = [
    "disconnected",
    "ring",
    "star",
    "complete",
    "hypercube",
    "barbell",
    "watts-strogatz",
    "random-regular",
]

SEEDS = [42, 137, 2718]

# Lambda_2 values (algebraic connectivity of each topology with 8 islands)
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


def read_csv(filepath):
    """Read a CSV file and return list of dicts with numeric values.
    Handles cabal 'Up to date' prefix by scanning for the header line."""
    rows = []
    with open(filepath) as f:
        lines = f.readlines()
    # Find the header line (skip any cabal output)
    header_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith("generation,"):
            header_idx = i
            break
    if header_idx is None:
        raise ValueError(f"No CSV header found in {filepath}")
    # Parse from header onward
    import io
    csv_text = "".join(lines[header_idx:])
    reader = csv.DictReader(io.StringIO(csv_text))
    for row in reader:
        rows.append({
            "generation": int(row["generation"]),
            "meanFitness": float(row["meanFitness"]),
            "bestFitness": float(row["bestFitness"]),
            "diversity": float(row["diversity"]),
        })
    return rows


def analyze():
    """Main analysis."""
    # Collect all data
    all_data = {}  # topo -> seed -> [rows]

    for topo in TOPOLOGIES:
        all_data[topo] = {}
        for seed in SEEDS:
            csv_path = os.path.join(RESULTS_DIR, f"{topo}_seed{seed}.csv")
            all_data[topo][seed] = read_csv(csv_path)

    # ===================================================================
    # SECTION 1: Summary table at key checkpoints
    # ===================================================================
    checkpoints = [0, 50, 100, 250, 500]

    print("=" * 100)
    print("ONEMAX PILOT STUDY RESULTS: 3 runs x 8 topologies x 500 generations")
    print("OneMax (N=100) | 8 islands | pop=50 | migInterval=10 | migrants=5")
    print("=" * 100)
    print()

    for gen_target in checkpoints:
        print(f"--- Generation {gen_target} ---")
        print(f"{'Topology':<18} {'lambda2':>8} {'MeanFit':>10} {'BestFit':>10} {'Diversity':>12} {'StdDev(Div)':>12}")
        print("-" * 72)

        for topo in TOPOLOGIES:
            mean_fits = []
            best_fits = []
            divs = []
            for seed in SEEDS:
                rows = all_data[topo][seed]
                row = None
                for r in rows:
                    if r["generation"] == gen_target:
                        row = r
                        break
                if row is None:
                    row = min(rows, key=lambda r: abs(r["generation"] - gen_target))
                mean_fits.append(row["meanFitness"])
                best_fits.append(row["bestFitness"])
                divs.append(row["diversity"])

            avg_mean_fit = statistics.mean(mean_fits)
            avg_best_fit = statistics.mean(best_fits)
            avg_div = statistics.mean(divs)
            std_div = statistics.stdev(divs) if len(divs) > 1 else 0.0

            print(f"{topo:<18} {LAMBDA2[topo]:>8.3f} {avg_mean_fit:>10.4f} {avg_best_fit:>10.4f} {avg_div:>12.6f} {std_div:>12.6f}")
        print()

    # ===================================================================
    # SECTION 2: Final generation (500) detailed analysis
    # ===================================================================
    print()
    print("=" * 100)
    print("FINAL GENERATION (500) -- DETAILED ANALYSIS")
    print("=" * 100)
    print()

    final_results = []
    for topo in TOPOLOGIES:
        mean_fits = []
        best_fits = []
        divs = []
        for seed in SEEDS:
            rows = all_data[topo][seed]
            final = rows[-1]
            mean_fits.append(final["meanFitness"])
            best_fits.append(final["bestFitness"])
            divs.append(final["diversity"])

        final_results.append({
            "topo": topo,
            "lambda2": LAMBDA2[topo],
            "mean_fit": statistics.mean(mean_fits),
            "best_fit_avg": statistics.mean(best_fits),
            "best_fit_max": max(best_fits),
            "diversity_avg": statistics.mean(divs),
            "diversity_std": statistics.stdev(divs) if len(divs) > 1 else 0.0,
            "mean_fits": mean_fits,
            "best_fits": best_fits,
            "divs": divs,
        })

    # Sort by diversity (descending)
    final_sorted = sorted(final_results, key=lambda x: x["diversity_avg"], reverse=True)

    print("Sorted by DIVERSITY (descending):")
    print(f"{'Rank':<5} {'Topology':<18} {'lambda2':>8} {'AvgDiv':>12} {'StdDiv':>10} {'AvgBestFit':>12}")
    print("-" * 67)
    for rank, r in enumerate(final_sorted, 1):
        print(f"{rank:<5} {r['topo']:<18} {r['lambda2']:>8.3f} {r['diversity_avg']:>12.6f} {r['diversity_std']:>10.6f} {r['best_fit_avg']:>12.4f}")

    # Sort by lambda2 (ascending) -- predicted diversity order
    predicted = sorted(final_results, key=lambda x: x["lambda2"])

    print()
    print("Sorted by LAMBDA_2 (ascending) -- PREDICTED diversity order:")
    print(f"{'Rank':<5} {'Topology':<18} {'lambda2':>8} {'AvgDiv':>12} {'DivRank':>8}")
    print("-" * 55)
    div_ranking = {r["topo"]: rank for rank, r in enumerate(final_sorted, 1)}
    for rank, r in enumerate(predicted, 1):
        print(f"{rank:<5} {r['topo']:<18} {r['lambda2']:>8.3f} {r['diversity_avg']:>12.6f} {div_ranking[r['topo']]:>8}")

    # ===================================================================
    # SECTION 3: Diversity trajectory (sampled generations)
    # ===================================================================
    print()
    print("=" * 100)
    print("DIVERSITY TRAJECTORIES (mean across 3 runs)")
    print("=" * 100)
    print()

    sample_gens = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
    header = f"{'Gen':>5}"
    for topo in TOPOLOGIES:
        header += f"  {topo[:12]:>12}"
    print(header)
    print("-" * (5 + 14 * len(TOPOLOGIES)))

    for gen_target in sample_gens:
        line = f"{gen_target:>5}"
        for topo in TOPOLOGIES:
            divs = []
            for seed in SEEDS:
                rows = all_data[topo][seed]
                row = None
                for r in rows:
                    if r["generation"] == gen_target:
                        row = r
                        break
                if row is None:
                    row = min(rows, key=lambda r: abs(r["generation"] - gen_target))
                divs.append(row["diversity"])
            avg_div = statistics.mean(divs)
            line += f"  {avg_div:>12.6f}"
        print(line)

    # ===================================================================
    # SECTION 4: Fitness trajectory
    # ===================================================================
    print()
    print("=" * 100)
    print("BEST FITNESS TRAJECTORIES (mean across 3 runs)")
    print("=" * 100)
    print()

    header = f"{'Gen':>5}"
    for topo in TOPOLOGIES:
        header += f"  {topo[:12]:>12}"
    print(header)
    print("-" * (5 + 14 * len(TOPOLOGIES)))

    for gen_target in sample_gens:
        line = f"{gen_target:>5}"
        for topo in TOPOLOGIES:
            fits = []
            for seed in SEEDS:
                rows = all_data[topo][seed]
                row = None
                for r in rows:
                    if r["generation"] == gen_target:
                        row = r
                        break
                if row is None:
                    row = min(rows, key=lambda r: abs(r["generation"] - gen_target))
                fits.append(row["bestFitness"])
            avg_fit = statistics.mean(fits)
            line += f"  {avg_fit:>12.4f}"
        print(line)

    # ===================================================================
    # SECTION 5: Diversity drop analysis
    # ===================================================================
    print()
    print("=" * 100)
    print("DIVERSITY DROP FROM INITIAL (gen 0) TO FINAL (gen 500)")
    print("=" * 100)
    print()

    drops = []
    print(f"{'Topology':<18} {'lambda2':>8} {'InitDiv':>12} {'FinalDiv':>12} {'Drop':>10} {'Drop%':>8}")
    print("-" * 70)
    for topo in TOPOLOGIES:
        init_divs = []
        final_divs = []
        for seed in SEEDS:
            rows = all_data[topo][seed]
            init_divs.append(rows[0]["diversity"])
            final_divs.append(rows[-1]["diversity"])
        init_avg = statistics.mean(init_divs)
        final_avg = statistics.mean(final_divs)
        drop = init_avg - final_avg
        drop_pct = (drop / init_avg) * 100 if init_avg > 0 else 0
        drops.append((topo, LAMBDA2[topo], init_avg, final_avg, drop, drop_pct))
        print(f"{topo:<18} {LAMBDA2[topo]:>8.3f} {init_avg:>12.6f} {final_avg:>12.6f} {drop:>10.6f} {drop_pct:>7.2f}%")

    # Sort by drop (descending)
    drops_sorted = sorted(drops, key=lambda x: x[4], reverse=True)
    print()
    print("Sorted by diversity drop (largest first):")
    for i, (topo, l2, init_d, final_d, drop, drop_pct) in enumerate(drops_sorted, 1):
        print(f"  {i}. {topo:<18} drop={drop:.6f} ({drop_pct:.2f}%)  lambda2={l2:.3f}")

    # ===================================================================
    # SECTION 6: Lambda_2 vs Diversity Correlation
    # ===================================================================
    print()
    print("=" * 100)
    print("LAMBDA_2 vs DIVERSITY CORRELATION ANALYSIS")
    print("=" * 100)
    print()

    # Compute Spearman rank correlation
    lambda2_order = sorted(final_results, key=lambda x: x["lambda2"])
    diversity_order = sorted(final_results, key=lambda x: x["diversity_avg"], reverse=True)

    lambda2_ranks = {r["topo"]: i for i, r in enumerate(lambda2_order)}
    diversity_ranks = {r["topo"]: i for i, r in enumerate(diversity_order)}

    n = len(TOPOLOGIES)
    d_squared_sum = sum(
        (lambda2_ranks[t] - diversity_ranks[t]) ** 2
        for t in TOPOLOGIES
    )
    spearman = 1 - (6 * d_squared_sum) / (n * (n ** 2 - 1))

    print(f"Spearman rank correlation (lambda2 ascending vs diversity descending): {spearman:.4f}")
    print()
    print("Interpretation:")
    if spearman > 0.7:
        print("  STRONG positive correlation: lower lambda2 -> higher diversity")
        print("  This SUPPORTS the theoretical prediction.")
    elif spearman > 0.3:
        print("  MODERATE positive correlation: some evidence for the prediction.")
    elif spearman > -0.3:
        print("  WEAK or no correlation: inconclusive.")
    else:
        print("  NEGATIVE correlation: evidence AGAINST the prediction.")

    print()
    print("Rank comparison (lower lambda2 should predict higher diversity):")
    print(f"{'Topology':<18} {'L2_rank':>8} {'Div_rank':>9} {'Diff':>6}")
    print("-" * 43)
    for topo in sorted(TOPOLOGIES, key=lambda t: lambda2_ranks[t]):
        lr = lambda2_ranks[topo]
        dr = diversity_ranks[topo]
        diff = lr - dr
        print(f"{topo:<18} {lr:>8} {dr:>9} {diff:>6}")

    # ===================================================================
    # SECTION 7: Per-run details
    # ===================================================================
    print()
    print("=" * 100)
    print("PER-RUN DETAILS (Generation 500)")
    print("=" * 100)
    print()

    print(f"{'Topology':<18} {'Seed':>6} {'MeanFit':>10} {'BestFit':>10} {'Diversity':>12}")
    print("-" * 58)
    for topo in TOPOLOGIES:
        for seed in SEEDS:
            rows = all_data[topo][seed]
            final = rows[-1]
            print(f"{topo:<18} {seed:>6} {final['meanFitness']:>10.4f} {final['bestFitness']:>10.4f} {final['diversity']:>12.6f}")
        print()

    # ===================================================================
    # SECTION 8: Convergence speed (generation to reach best fitness = 1.0)
    # ===================================================================
    print()
    print("=" * 100)
    print("CONVERGENCE SPEED (first generation where bestFitness = 1.0)")
    print("=" * 100)
    print()

    print(f"{'Topology':<18} {'lambda2':>8} {'Seed42':>8} {'Seed137':>8} {'Seed2718':>8} {'Avg':>8}")
    print("-" * 56)
    for topo in TOPOLOGIES:
        gens_to_opt = []
        for seed in SEEDS:
            rows = all_data[topo][seed]
            gen_opt = None
            for r in rows:
                if r["bestFitness"] >= 1.0:
                    gen_opt = r["generation"]
                    break
            gens_to_opt.append(gen_opt)

        avg_gen = statistics.mean([g for g in gens_to_opt if g is not None]) if any(g is not None for g in gens_to_opt) else None
        seed_strs = [f"{g:>8}" if g is not None else f"{'N/A':>8}" for g in gens_to_opt]
        avg_str = f"{avg_gen:>8.1f}" if avg_gen is not None else f"{'N/A':>8}"
        print(f"{topo:<18} {LAMBDA2[topo]:>8.3f} {seed_strs[0]} {seed_strs[1]} {seed_strs[2]} {avg_str}")


if __name__ == "__main__":
    analyze()
