#!/usr/bin/env python3
"""
PCA eigenvalue spectrum analysis for topology-diversity experiments.

Motivation (Clio Vega): The eigenvalue RATIOS of PCA on population genomes
may encode irreducible representation multiplicities of the cactus group,
giving a representation-theoretic fingerprint of how topology affects diversity.

Instead of collapsing diversity to a single scalar (mean pairwise distance),
track the FULL eigenvalue spectrum of the population covariance matrix.

Two modes of operation:
  1. GENOME MODE (--genome-dir): reads per-generation genome dumps (CSV with
     one row per individual, columns = loci). Computes real PCA on population.
  2. SUMMARY MODE (--results-dir): reads standard simulation output (only has
     scalar diversity). Constructs a PROXY spectrum from the diversity time
     series across seeds -- a between-run covariance decomposition.

Genome mode is preferred; summary mode is a useful stopgap.

Usage:
  # Genome mode (requires --dump-genomes output from Haskell sim)
  python experiments/analyze_pca_spectrum.py --genome-dir results/star_vs_cycle/genomes/

  # Summary mode (works with existing CSV output)
  python experiments/analyze_pca_spectrum.py --results-dir results/star_vs_cycle/
"""

import argparse
import csv
import math
import os
import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# Shared utilities
# ---------------------------------------------------------------------------

def mean(xs):
    return sum(xs) / len(xs) if xs else 0.0


def std(xs):
    if len(xs) < 2:
        return 0.0
    m = mean(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))


def load_csv(path):
    """Load CSV, skip non-CSV header lines."""
    with open(path) as f:
        lines = [l for l in f if l.strip() and
                 (l.strip()[0].isdigit() or l.strip().startswith("generation"))]
    reader = csv.DictReader(lines)
    return list(reader)


# ---------------------------------------------------------------------------
# Eigenvalue decomposition (pure Python, no numpy dependency)
# ---------------------------------------------------------------------------

def mat_mul(A, B):
    """Multiply two matrices (lists of lists)."""
    n = len(A)
    m = len(B[0])
    k = len(B)
    C = [[0.0] * m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            s = 0.0
            for p in range(k):
                s += A[i][p] * B[p][j]
            C[i][j] = s
    return C


def mat_transpose(A):
    """Transpose a matrix."""
    n = len(A)
    m = len(A[0])
    return [[A[i][j] for i in range(n)] for j in range(m)]


def mat_vec_mul(A, v):
    """Multiply matrix by vector."""
    n = len(A)
    return [sum(A[i][j] * v[j] for j in range(len(v))) for i in range(n)]


def vec_norm(v):
    return math.sqrt(sum(x * x for x in v))


def vec_dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def power_iteration(matrix, n_components, max_iter=200, tol=1e-10):
    """
    Extract top eigenvalues via deflated power iteration.

    For the covariance matrices we deal with (typically <= 100x100 for genome
    length, or <= 50x50 for seed count), this is fast and sufficient.
    """
    n = len(matrix)
    # Work on a copy
    M = [row[:] for row in matrix]
    eigenvalues = []

    for _ in range(min(n_components, n)):
        # Random-ish starting vector (deterministic)
        v = [1.0 / math.sqrt(n) if i % 2 == 0 else -1.0 / math.sqrt(n)
             for i in range(n)]
        norm = vec_norm(v)
        if norm > 0:
            v = [x / norm for x in v]

        eigenvalue = 0.0
        for _ in range(max_iter):
            Mv = mat_vec_mul(M, v)
            new_eigenvalue = vec_dot(v, Mv)
            norm = vec_norm(Mv)
            if norm < tol:
                break
            new_v = [x / norm for x in Mv]
            # Check convergence
            diff = sum((a - b) ** 2 for a, b in zip(new_v, v))
            v = new_v
            eigenvalue = new_eigenvalue
            if diff < tol:
                break

        eigenvalues.append(max(eigenvalue, 0.0))  # Clamp negative numerical noise

        # Deflate: M = M - eigenvalue * v * v^T
        for i in range(n):
            for j in range(n):
                M[i][j] -= eigenvalue * v[i] * v[j]

    return eigenvalues


def covariance_matrix(data_rows):
    """
    Compute covariance matrix from data rows (list of lists of floats).
    Each row is an observation, each column is a variable.
    Centers the data first.
    """
    n = len(data_rows)
    if n < 2:
        return []
    p = len(data_rows[0])

    # Compute column means
    col_means = [0.0] * p
    for row in data_rows:
        for j in range(p):
            col_means[j] += row[j]
    col_means = [m / n for m in col_means]

    # Center the data
    centered = [[row[j] - col_means[j] for j in range(p)] for row in data_rows]

    # Covariance = (1/(n-1)) * X^T X
    cov = [[0.0] * p for _ in range(p)]
    for i in range(p):
        for j in range(i, p):
            s = sum(centered[k][i] * centered[k][j] for k in range(n))
            cov[i][j] = s / (n - 1)
            cov[j][i] = cov[i][j]

    return cov


# ---------------------------------------------------------------------------
# Spectrum analysis
# ---------------------------------------------------------------------------

def compute_spectrum(cov, n_components=None):
    """
    Compute eigenvalue spectrum from a covariance matrix.
    Returns sorted eigenvalues (descending) and derived metrics.
    """
    p = len(cov)
    if p == 0:
        return {"eigenvalues": [], "ratios": [], "cumulative_variance": [],
                "effective_rank": 0, "spectral_entropy": 0.0}

    if n_components is None:
        n_components = min(p, 20)  # Top 20 components usually sufficient

    eigenvalues = power_iteration(cov, n_components)
    eigenvalues.sort(reverse=True)

    # Filter out near-zero eigenvalues
    eigenvalues = [max(ev, 0.0) for ev in eigenvalues]

    total_var = sum(eigenvalues)
    if total_var == 0:
        return {"eigenvalues": eigenvalues, "ratios": [0.0] * len(eigenvalues),
                "cumulative_variance": [0.0] * len(eigenvalues),
                "effective_rank": 0, "spectral_entropy": 0.0}

    # Ratios relative to lambda_1
    lambda_1 = eigenvalues[0] if eigenvalues[0] > 0 else 1.0
    ratios = [ev / lambda_1 for ev in eigenvalues]

    # Cumulative variance explained
    cumvar = []
    running = 0.0
    for ev in eigenvalues:
        running += ev / total_var
        cumvar.append(running)

    # Effective rank (exponential of spectral entropy)
    # H = -sum(p_i * log(p_i)) where p_i = lambda_i / sum(lambda)
    proportions = [ev / total_var for ev in eigenvalues if ev > 0]
    spectral_entropy = -sum(p * math.log(p) for p in proportions if p > 0)
    effective_rank = math.exp(spectral_entropy) if spectral_entropy > 0 else 1.0

    return {
        "eigenvalues": eigenvalues,
        "ratios": ratios,
        "cumulative_variance": cumvar,
        "effective_rank": effective_rank,
        "spectral_entropy": spectral_entropy,
    }


# ---------------------------------------------------------------------------
# Mode 1: Genome-based PCA (requires --dump-genomes output)
# ---------------------------------------------------------------------------

def analyze_genome_dir(genome_dir, checkpoints):
    """
    Analyze genome dumps. Expected file structure:
      genome_dir/<topology>_<domain>_seed<N>/gen_<G>.csv

    Each gen_<G>.csv has one row per individual (all islands concatenated),
    columns = locus values (0/1 for NK landscapes).
    """
    if not os.path.isdir(genome_dir):
        print(f"ERROR: Genome directory '{genome_dir}' not found.")
        print()
        print("To generate genome dumps, rebuild the Haskell simulation with")
        print("--dump-genomes flag:")
        print()
        print("  cabal run topology-sim -- --dump-genomes results/genomes/ \\")
        print("    --domain nk4 ring 8 50 10 5 500 42")
        print()
        print("This will write per-generation genome CSVs to the specified directory.")
        return

    # Discover topology/domain/seed combos from directory structure
    runs = defaultdict(lambda: defaultdict(dict))  # topo -> domain -> seed -> dir
    for entry in sorted(os.listdir(genome_dir)):
        entry_path = os.path.join(genome_dir, entry)
        if not os.path.isdir(entry_path):
            continue
        # Parse directory name: <topology>_<domain>_seed<N>
        parts = entry.rsplit("_seed", 1)
        if len(parts) != 2:
            continue
        prefix, seed = parts
        # Split prefix into topology and domain
        # Domain is the last underscore-separated token (nk0, nk2, nk4, etc.)
        prefix_parts = prefix.rsplit("_", 1)
        if len(prefix_parts) != 2:
            continue
        topo, domain = prefix_parts
        runs[topo][domain][seed] = entry_path

    if not runs:
        print(f"ERROR: No genome dump directories found in '{genome_dir}'.")
        print("Expected subdirectories like: ring_nk4_seed42/")
        return

    topologies = sorted(runs.keys())
    domains = sorted(set(d for topo in runs for d in runs[topo]))

    print("=" * 78)
    print("PCA EIGENVALUE SPECTRUM ANALYSIS — Genome Mode")
    print("=" * 78)
    print()
    print(f"Genome directory: {genome_dir}")
    print(f"Topologies: {', '.join(topologies)}")
    print(f"Domains: {', '.join(domains)}")
    total_runs = sum(len(runs[t][d]) for t in runs for d in runs[t])
    print(f"Total runs: {total_runs}")
    print()

    for domain in domains:
        print(f"\n{'='*78}")
        print(f"DOMAIN: {domain.upper()}")
        print(f"{'='*78}")

        for gen in checkpoints:
            gen_file = f"gen_{gen}.csv"
            print(f"\n--- Generation {gen} ---")
            print(f"{'Topology':<20} {'EffRank':>8} {'H(spec)':>8} "
                  f"{'lambda_1':>10} {'r2':>6} {'r3':>6} {'r4':>6} {'r5':>6}")
            print("-" * 78)

            for topo in topologies:
                if domain not in runs[topo]:
                    continue

                # Collect genomes across all seeds for this topology at this gen
                all_spectra = []
                for seed, run_dir in sorted(runs[topo][domain].items()):
                    gen_path = os.path.join(run_dir, gen_file)
                    if not os.path.exists(gen_path):
                        continue

                    # Load genome matrix
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

                    spec = compute_spectrum(cov, n_components=min(len(cov), 20))
                    all_spectra.append(spec)

                if not all_spectra:
                    print(f"{topo:<20} {'N/A':>8}")
                    continue

                # Average spectrum metrics across seeds
                avg_eff_rank = mean([s["effective_rank"] for s in all_spectra])
                avg_entropy = mean([s["spectral_entropy"] for s in all_spectra])
                avg_lambda1 = mean([s["eigenvalues"][0] for s in all_spectra
                                    if s["eigenvalues"]])

                # Average ratios (pad shorter spectra with 0)
                max_len = max(len(s["ratios"]) for s in all_spectra)
                avg_ratios = []
                for i in range(min(max_len, 5)):
                    vals = [s["ratios"][i] for s in all_spectra
                            if i < len(s["ratios"])]
                    avg_ratios.append(mean(vals) if vals else 0.0)

                # Pad to 5
                while len(avg_ratios) < 5:
                    avg_ratios.append(0.0)

                print(f"{topo:<20} {avg_eff_rank:>8.2f} {avg_entropy:>8.3f} "
                      f"{avg_lambda1:>10.4f} "
                      f"{avg_ratios[1]:>6.3f} {avg_ratios[2]:>6.3f} "
                      f"{avg_ratios[3]:>6.3f} {avg_ratios[4]:>6.3f}")

            # Eta-squared on effective rank across topologies
            groups = []
            for topo in topologies:
                if domain not in runs[topo]:
                    continue
                eff_ranks = []
                for seed, run_dir in sorted(runs[topo][domain].items()):
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
                    spec = compute_spectrum(cov, n_components=min(len(cov), 20))
                    eff_ranks.append(spec["effective_rank"])
                groups.append(eff_ranks)

            if len(groups) >= 2 and all(groups):
                eta2 = eta_squared_groups(groups)
                print(f"\n  eta^2(effective_rank) across topologies: {eta2:.4f}")


def eta_squared_groups(groups):
    """One-way eta-squared."""
    all_vals = [v for g in groups for v in g]
    if not all_vals:
        return 0.0
    grand = mean(all_vals)
    ss_total = sum((v - grand) ** 2 for v in all_vals)
    if ss_total == 0:
        return 0.0
    ss_between = sum(len(g) * (mean(g) - grand) ** 2 for g in groups if g)
    return ss_between / ss_total


# ---------------------------------------------------------------------------
# Mode 2: Summary-based proxy spectrum (works with existing CSV output)
# ---------------------------------------------------------------------------

def analyze_results_dir(results_dir, checkpoints):
    """
    Proxy spectrum analysis from standard simulation output.

    Since we only have scalar diversity per run, we construct a between-run
    covariance matrix: rows = seeds, columns = generations. PCA on this
    matrix reveals the dominant temporal MODES of diversity variation.

    This is NOT the same as PCA on genomes (which reveals the dominant
    allelic variation modes). But it provides:
      - How many independent temporal modes drive diversity differences
      - Whether topologies differ in their temporal eigenstructure
      - A proxy for the "complexity" of diversity dynamics

    For Clio's representation-theoretic question: if different topologies
    produce different eigenvalue ratio profiles in the TEMPORAL covariance,
    this still suggests different symmetry structures.
    """
    if not os.path.isdir(results_dir):
        print(f"ERROR: Results directory '{results_dir}' not found.")
        return

    # Discover topology/domain combos from filenames
    # Expected format: <topology>_<domain>_seed<N>.csv
    data = defaultdict(lambda: defaultdict(dict))
    for fname in sorted(os.listdir(results_dir)):
        if not fname.endswith(".csv"):
            continue
        # Parse: topology_domain_seedN.csv
        base = fname[:-4]  # strip .csv
        parts = base.rsplit("_seed", 1)
        if len(parts) != 2:
            continue
        prefix, seed = parts
        # Domain is last underscore-separated token
        prefix_parts = prefix.rsplit("_", 1)
        if len(prefix_parts) != 2:
            continue
        topo, domain = prefix_parts
        path = os.path.join(results_dir, fname)
        rows = load_csv(path)
        gen_data = {}
        for row in rows:
            gen = int(row["generation"])
            gen_data[gen] = {
                "diversity": float(row["diversity"]),
                "meanFitness": float(row["meanFitness"]),
                "bestFitness": float(row["bestFitness"]),
            }
        data[topo][domain][seed] = gen_data

    if not data:
        print(f"ERROR: No CSV files found in '{results_dir}'.")
        return

    topologies = sorted(data.keys())
    domains = sorted(set(d for topo in data for d in data[topo]))

    print("=" * 78)
    print("PCA EIGENVALUE SPECTRUM ANALYSIS — Summary Mode (Proxy)")
    print("=" * 78)
    print()
    print("NOTE: This analysis uses TEMPORAL covariance across seeds, not")
    print("genome-level PCA. It reveals dominant modes of diversity variation")
    print("over time. For true genome PCA, rebuild the Haskell sim with")
    print("--dump-genomes (see below).")
    print()
    print(f"Results directory: {results_dir}")
    print(f"Topologies: {', '.join(topologies)}")
    print(f"Domains: {', '.join(domains)}")
    total_runs = sum(len(data[t][d]) for t in data for d in data[t])
    print(f"Total runs: {total_runs}")
    print()

    for domain in domains:
        print(f"\n{'='*78}")
        print(f"DOMAIN: {domain.upper()}")
        print(f"{'='*78}")

        for topo in topologies:
            if domain not in data[topo]:
                continue

            seeds = sorted(data[topo][domain].keys())
            if len(seeds) < 3:
                print(f"\n{topo}: too few seeds ({len(seeds)}) for PCA")
                continue

            # Find common generations across all seeds
            common_gens = None
            for seed in seeds:
                gens = set(data[topo][domain][seed].keys())
                if common_gens is None:
                    common_gens = gens
                else:
                    common_gens &= gens
            common_gens = sorted(common_gens)

            if len(common_gens) < 3:
                print(f"\n{topo}: too few common generations for PCA")
                continue

            # Build matrix: rows = seeds, columns = generations
            # Each cell = diversity at that (seed, generation)
            diversity_matrix = []
            for seed in seeds:
                row = [data[topo][domain][seed][gen]["diversity"]
                       for gen in common_gens]
                diversity_matrix.append(row)

            cov = covariance_matrix(diversity_matrix)
            if not cov:
                continue

            n_comp = min(len(seeds), len(common_gens), 10)
            spec = compute_spectrum(cov, n_components=n_comp)

            print(f"\n--- {topo} ({len(seeds)} seeds, {len(common_gens)} generations) ---")
            print(f"  Effective rank: {spec['effective_rank']:.2f}")
            print(f"  Spectral entropy: {spec['spectral_entropy']:.3f}")
            print()

            # Print eigenvalue table
            print(f"  {'Component':>10} {'Eigenvalue':>12} {'Ratio':>8} "
                  f"{'Var%':>8} {'CumVar%':>8}")
            print(f"  {'-'*50}")
            for i, (ev, ratio, cumvar) in enumerate(
                    zip(spec["eigenvalues"], spec["ratios"],
                        spec["cumulative_variance"])):
                var_pct = (ev / sum(spec["eigenvalues"]) * 100
                           if sum(spec["eigenvalues"]) > 0 else 0)
                print(f"  {'PC'+str(i+1):>10} {ev:>12.6f} {ratio:>8.4f} "
                      f"{var_pct:>7.1f}% {cumvar*100:>7.1f}%")

        # Cross-topology comparison at this domain
        print(f"\n{'-'*78}")
        print(f"CROSS-TOPOLOGY COMPARISON — {domain.upper()}")
        print(f"{'-'*78}")
        print(f"{'Topology':<20} {'EffRank':>8} {'H(spec)':>8} "
              f"{'PC1%':>7} {'PC2%':>7} {'PC3%':>7} {'r2/r1':>7} {'r3/r1':>7}")
        print("-" * 78)

        all_eff_ranks = {}
        for topo in topologies:
            if domain not in data[topo]:
                continue
            seeds = sorted(data[topo][domain].keys())
            if len(seeds) < 3:
                continue

            common_gens = None
            for seed in seeds:
                gens = set(data[topo][domain][seed].keys())
                if common_gens is None:
                    common_gens = gens
                else:
                    common_gens &= gens
            common_gens = sorted(common_gens)
            if len(common_gens) < 3:
                continue

            diversity_matrix = []
            for seed in seeds:
                row = [data[topo][domain][seed][gen]["diversity"]
                       for gen in common_gens]
                diversity_matrix.append(row)

            cov = covariance_matrix(diversity_matrix)
            if not cov:
                continue

            n_comp = min(len(seeds), len(common_gens), 10)
            spec = compute_spectrum(cov, n_components=n_comp)

            total = sum(spec["eigenvalues"]) if spec["eigenvalues"] else 1.0
            pcs = [(ev / total * 100) if total > 0 else 0.0
                   for ev in spec["eigenvalues"][:3]]
            while len(pcs) < 3:
                pcs.append(0.0)
            ratios = spec["ratios"][:3]
            while len(ratios) < 3:
                ratios.append(0.0)

            print(f"{topo:<20} {spec['effective_rank']:>8.2f} "
                  f"{spec['spectral_entropy']:>8.3f} "
                  f"{pcs[0]:>6.1f}% {pcs[1]:>6.1f}% {pcs[2]:>6.1f}% "
                  f"{ratios[1]:>7.4f} {ratios[2]:>7.4f}")

            all_eff_ranks[topo] = spec["effective_rank"]

        if len(all_eff_ranks) >= 2:
            print()
            print("  Effective rank comparison:")
            sorted_topos = sorted(all_eff_ranks.items(), key=lambda x: -x[1])
            for topo, er in sorted_topos:
                bar = "#" * int(er * 5)
                print(f"    {topo:<20} {er:.2f}  {bar}")

    # Summary and instructions for genome mode
    print()
    print("=" * 78)
    print("NEXT STEPS: Full Genome-Level PCA")
    print("=" * 78)
    print()
    print("The temporal proxy analysis above shows between-run variation modes.")
    print("For Clio's representation-theoretic fingerprint, we need per-genome PCA.")
    print()
    print("To enable genome dumps in the Haskell simulation:")
    print()
    print("  cabal run topology-sim -- --dump-genomes results/genomes/ \\")
    print("    --domain nk4 ring 8 50 10 5 500 42")
    print()
    print("Then re-run this script in genome mode:")
    print()
    print("  python experiments/analyze_pca_spectrum.py \\")
    print("    --genome-dir results/genomes/")
    print()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="PCA eigenvalue spectrum analysis for topology experiments."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--genome-dir",
                       help="Directory with genome dumps (genome mode)")
    group.add_argument("--results-dir",
                       help="Directory with standard CSV results (summary mode)")
    parser.add_argument("--checkpoints", type=str, default="50,100,200,500",
                        help="Comma-separated generation checkpoints "
                             "(genome mode only, default: 50,100,200,500)")

    args = parser.parse_args()
    checkpoints = [int(x) for x in args.checkpoints.split(",")]

    if args.genome_dir:
        analyze_genome_dir(args.genome_dir, checkpoints)
    else:
        analyze_results_dir(args.results_dir, checkpoints)


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/..")
    main()
