#!/usr/bin/env python3
"""
Direct test of Clio's functor hypothesis: does # significant PCs = beta_1?

Topologies tested:
  star (8 nodes):       beta_1 = 0  (tree, 7 edges, 8 nodes)
  ring (7 nodes):       beta_1 = 1  (6 edges + 1 closing edge = 7 edges, 7 nodes)
  ring-chord1 (8 nodes): beta_1 = 2  (ring + 1 chord)
  ring-chord2 (8 nodes): beta_1 = 3  (ring + 2 chords)
  ring-chord3 (8 nodes): beta_1 = 4  (ring + 3 chords)

Prediction: effective rank should scale with beta_1.
"""

import math
import os
import sys
from collections import defaultdict


def covariance_matrix(data_rows):
    n = len(data_rows)
    if n < 2:
        return []
    p = len(data_rows[0])
    col_means = [sum(row[j] for row in data_rows) / n for j in range(p)]
    centered = [[row[j] - col_means[j] for j in range(p)] for row in data_rows]
    cov = [[0.0] * p for _ in range(p)]
    for i in range(p):
        for j in range(i, p):
            s = sum(centered[k][i] * centered[k][j] for k in range(n))
            cov[i][j] = s / (n - 1)
            cov[j][i] = cov[i][j]
    return cov


def power_iteration(matrix, n_components, max_iter=200, tol=1e-10):
    n = len(matrix)
    M = [row[:] for row in matrix]
    eigenvalues = []
    for _ in range(min(n_components, n)):
        v = [1.0 / math.sqrt(n) if i % 2 == 0 else -1.0 / math.sqrt(n) for i in range(n)]
        norm = math.sqrt(sum(x*x for x in v))
        if norm > 0:
            v = [x / norm for x in v]
        eigenvalue = 0.0
        for _ in range(max_iter):
            Mv = [sum(M[i][j] * v[j] for j in range(n)) for i in range(n)]
            new_eigenvalue = sum(v[i] * Mv[i] for i in range(n))
            norm = math.sqrt(sum(x*x for x in Mv))
            if norm < tol:
                break
            new_v = [x / norm for x in Mv]
            diff = sum((a - b) ** 2 for a, b in zip(new_v, v))
            v = new_v
            eigenvalue = new_eigenvalue
            if diff < tol:
                break
        eigenvalues.append(max(eigenvalue, 0.0))
        for i in range(n):
            for j in range(n):
                M[i][j] -= eigenvalue * v[i] * v[j]
    return sorted(eigenvalues, reverse=True)


def count_significant_pcs(eigenvalues, threshold=0.05):
    """Count PCs that explain > threshold fraction of total variance."""
    total = sum(eigenvalues)
    if total == 0:
        return 0
    return sum(1 for ev in eigenvalues if ev / total > threshold)


def effective_rank(eigenvalues):
    total = sum(ev for ev in eigenvalues if ev > 0)
    if total == 0:
        return 1.0
    proportions = [ev / total for ev in eigenvalues if ev > 0]
    entropy = -sum(p * math.log(p) for p in proportions if p > 0)
    return math.exp(entropy)


GENOME_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "genomes")

BETA1_MAP = {
    "star": 0,
    "ring": 1,
    "ring-chord1": 2,
    "ring-chord2": 3,
    "ring-chord3": 4,
}

CHECKPOINTS = [100, 200, 300, 500]

# Discover runs
runs = defaultdict(lambda: defaultdict(list))  # topo -> seed -> path
for entry in sorted(os.listdir(GENOME_DIR)):
    entry_path = os.path.join(GENOME_DIR, entry)
    if not os.path.isdir(entry_path):
        continue
    parts = entry.rsplit("_seed", 1)
    if len(parts) != 2:
        continue
    prefix, seed = parts
    prefix_parts = prefix.rsplit("_", 1)
    if len(prefix_parts) != 2:
        continue
    topo, domain = prefix_parts
    if domain != "nk4":
        continue
    runs[topo][seed] = entry_path

print("=" * 80)
print("CLIO'S FUNCTOR HYPOTHESIS TEST: # significant PCs vs beta_1")
print("=" * 80)
print()
print("Topologies:")
for topo in sorted(BETA1_MAP.keys(), key=lambda t: BETA1_MAP[t]):
    seeds = len(runs[topo])
    print(f"  {topo:<15} beta_1={BETA1_MAP[topo]}  ({seeds} seeds)")
print()

for gen in CHECKPOINTS:
    gen_file = f"gen_{gen}.csv"
    print(f"\n{'='*80}")
    print(f"GENERATION {gen}")
    print(f"{'='*80}")
    print(f"{'Topology':<15} {'beta_1':>6} | {'EffRank':>8} {'#PC>5%':>7} {'#PC>10%':>8} "
          f"| {'PC1%':>6} {'PC2%':>6} {'PC3%':>6} {'PC4%':>6} {'PC5%':>6}")
    print("-" * 80)

    summary = []  # (beta_1, eff_ranks_per_seed)

    for topo in sorted(BETA1_MAP.keys(), key=lambda t: BETA1_MAP[t]):
        b1 = BETA1_MAP[topo]
        eff_ranks = []
        sig_pcs_5 = []
        sig_pcs_10 = []
        pc_pcts = []

        for seed, run_dir in sorted(runs[topo].items()):
            gen_path = os.path.join(run_dir, gen_file)
            if not os.path.exists(gen_path):
                continue
            genomes = []
            with open(gen_path) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    genomes.append([float(x) for x in line.split(",")])
            if len(genomes) < 3:
                continue

            cov = covariance_matrix(genomes)
            if not cov:
                continue
            evs = power_iteration(cov, min(len(cov), 20))
            er = effective_rank(evs)
            n5 = count_significant_pcs(evs, 0.05)
            n10 = count_significant_pcs(evs, 0.10)
            total = sum(evs)
            pcts = [(ev / total * 100) if total > 0 else 0 for ev in evs[:5]]
            while len(pcts) < 5:
                pcts.append(0.0)

            eff_ranks.append(er)
            sig_pcs_5.append(n5)
            sig_pcs_10.append(n10)
            pc_pcts.append(pcts)

        if not eff_ranks:
            continue

        avg_er = sum(eff_ranks) / len(eff_ranks)
        avg_n5 = sum(sig_pcs_5) / len(sig_pcs_5)
        avg_n10 = sum(sig_pcs_10) / len(sig_pcs_10)
        avg_pcts = [sum(p[i] for p in pc_pcts) / len(pc_pcts) for i in range(5)]

        print(f"{topo:<15} {b1:>6} | {avg_er:>8.2f} {avg_n5:>7.1f} {avg_n10:>8.1f} "
              f"| {avg_pcts[0]:>5.1f}% {avg_pcts[1]:>5.1f}% {avg_pcts[2]:>5.1f}% "
              f"{avg_pcts[3]:>5.1f}% {avg_pcts[4]:>5.1f}%")

        summary.append((b1, eff_ranks))

    # Correlation analysis
    if len(summary) >= 3:
        beta1s = []
        mean_ers = []
        for b1, ers in summary:
            beta1s.append(b1)
            mean_ers.append(sum(ers) / len(ers))

        # Pearson correlation
        n = len(beta1s)
        mean_b = sum(beta1s) / n
        mean_e = sum(mean_ers) / n
        cov_be = sum((beta1s[i] - mean_b) * (mean_ers[i] - mean_e) for i in range(n))
        var_b = sum((b - mean_b) ** 2 for b in beta1s)
        var_e = sum((e - mean_e) ** 2 for e in mean_ers)
        if var_b > 0 and var_e > 0:
            r = cov_be / math.sqrt(var_b * var_e)
            print(f"\n  Pearson r(beta_1, effective_rank) = {r:.4f}")
            print(f"  Direction: {'POSITIVE (more cycles = higher rank)' if r > 0 else 'NEGATIVE (more cycles = lower rank)'}")
        else:
            print(f"\n  Cannot compute correlation (zero variance)")

        print(f"\n  beta_1 -> mean effective rank:")
        for b1, ers in summary:
            bar = "#" * int(sum(ers) / len(ers) * 3)
            print(f"    beta_1={b1}: {sum(ers)/len(ers):.2f}  {bar}")

print("\n\n" + "=" * 80)
print("INTERPRETATION")
print("=" * 80)
print("""
Clio's hypothesis: # significant PCs = beta_1 (the first Betti number).

The functor prediction is that the number of independent directions of
genetic diversity in the population should equal the number of independent
cycles in the communication topology.

Key observations from the data above:
- If effective rank INCREASES with beta_1: topology creates new diversity dimensions
- If effective rank DECREASES with beta_1: cycles CONCENTRATE diversity (fewer but stronger modes)
- If # significant PCs matches beta_1: strong evidence for the functor hypothesis
""")
