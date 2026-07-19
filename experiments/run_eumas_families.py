#!/usr/bin/env python3
"""
EUMAS graph family experiments.

Runs Family 1 (5 graphs, same beta_1=2, varying arrangement) and
Family 2 (4 graphs, varying beta_1) through an island-model GA on
OneMax and NK4 domains. Matches the Haskell topology-sim parameters:

  - 12 islands, pop 50 per island, migration interval 10,
    5 migrants per event, 500 generations, 30 seeds.

The migration model matches IslandGA.hs:
  - For each island i, its neighbors (topo[i]) send their top-5 fittest
    individuals, which replace the worst individuals in island i.

Domains:
  - OneMax (L=100): smooth, negative control
  - NK4 (N=100, K=4): rugged epistatic landscape

Output: CSV files in results/eumas_families/{domain}/{graph_name}_seed{s}.csv
"""

import os
import sys
import time
import csv
import numpy as np

# Add graph_families to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'categorical-evolution', 'cais2026'))
from graph_families import get_family1, get_family2, betti_1, lambda_2

import networkx as nx


# ---------------------------------------------------------------------------
# Parameters (matching Haskell experiments exactly)
# ---------------------------------------------------------------------------

NUM_ISLANDS = 12
POP_SIZE = 50        # per island
MIG_INTERVAL = 10
NUM_MIGRANTS = 5
TOTAL_GENS = 500
GENOME_LENGTH = 100
NK_K = 4
TOURNAMENT_SIZE = 3

SEEDS = [42, 137, 2718, 314, 1618, 7, 99, 256, 512, 1024,
         2048, 4096, 8192, 16384, 32768, 11, 23, 37, 53, 71,
         89, 97, 113, 131, 151, 173, 191, 199, 211, 223]

RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', 'results', 'eumas_families')


# ---------------------------------------------------------------------------
# NK4 landscape (matching NKLandscape.hs)
# ---------------------------------------------------------------------------

def build_nk_table(n=100, k=4, landscape_seed=12345):
    """Build NK fitness table matching Haskell implementation.

    Each locus i has a lookup table of size 2^(K+1) mapping the (K+1)-bit
    pattern (locus i and its K neighbors in circular fashion) to a fitness
    contribution in [0,1].

    The Haskell uses: mkStdGen (landscapeSeed * 1000003 + i * 7919)
    We replicate the table structure but use numpy RNG seeded identically
    (the actual random values will differ from Haskell, but the landscape
    structure - N, K, circular adjacency - is the same).
    """
    table_size = 2 ** (k + 1)
    tables = []
    for i in range(n):
        rng = np.random.default_rng(landscape_seed * 1000003 + i * 7919)
        tables.append(rng.random(table_size))
    return tables


# Precompute NK4 table (used across all runs)
NK4_TABLE = build_nk_table(GENOME_LENGTH, NK_K)


def _precompute_nk4_indices():
    """Precompute the bit-gathering indices for NK4 evaluation."""
    n = GENOME_LENGTH
    k = NK_K
    # For each locus i, positions (i, i+1, ..., i+K) mod N
    indices = np.zeros((n, k + 1), dtype=np.int32)
    for i in range(n):
        for j in range(k + 1):
            indices[i, j] = (i + j) % n
    return indices


NK4_INDICES = _precompute_nk4_indices()
# Stack tables into 2D array for vectorized lookup
NK4_TABLE_ARRAY = np.array(NK4_TABLE)  # shape (100, 32)


def nk4_fitness_pop(pop):
    """Fully vectorized NK4 fitness for entire population.

    Uses precomputed index arrays and table lookups.
    """
    n_ind = pop.shape[0]
    n_loci = GENOME_LENGTH
    k = NK_K

    # Gather the K+1 bits for each locus for all individuals
    # pop shape: (n_ind, n_loci), NK4_INDICES shape: (n_loci, K+1)
    # gathered shape: (n_ind, n_loci, K+1)
    gathered = pop[:, NK4_INDICES]  # shape (n_ind, n_loci, K+1)

    # Convert K+1 bits to integer index using powers of 2
    powers = 2 ** np.arange(k, -1, -1, dtype=np.int32)  # [16, 8, 4, 2, 1]
    table_indices = (gathered * powers[np.newaxis, np.newaxis, :]).sum(axis=2)  # (n_ind, n_loci)

    # Look up fitness contributions: for each locus i, use NK4_TABLE_ARRAY[i, index]
    # Use advanced indexing
    locus_range = np.arange(n_loci)
    contributions = NK4_TABLE_ARRAY[locus_range[np.newaxis, :], table_indices]  # (n_ind, n_loci)

    return contributions.mean(axis=1)


# ---------------------------------------------------------------------------
# OneMax
# ---------------------------------------------------------------------------

def onemax_fitness_pop(pop):
    """OneMax fitness for entire population."""
    return pop.sum(axis=1).astype(float) / GENOME_LENGTH


# ---------------------------------------------------------------------------
# GA operators (matching IslandGA.hs)
# ---------------------------------------------------------------------------

def tournament_select(rng, pop, fitnesses):
    """Tournament selection, size 3, with replacement. Returns new population."""
    n = len(pop)
    selected = np.empty_like(pop)
    for i in range(n):
        contestants = rng.integers(0, n, size=TOURNAMENT_SIZE)
        best = contestants[np.argmax(fitnesses[contestants])]
        selected[i] = pop[best]
    return selected


def crossover(rng, pop):
    """Uniform crossover matching Haskell nkCrossover.

    The Haskell NK crossover generates a random coin per locus and picks
    from parent1 or parent2. We pair individuals sequentially.
    """
    result = pop.copy()
    n = len(pop)
    for i in range(0, n - 1, 2):
        mask = rng.integers(0, 2, size=GENOME_LENGTH, dtype=np.int8)
        child1 = np.where(mask == 0, pop[i], pop[i + 1])
        child2 = np.where(mask == 0, pop[i + 1], pop[i])
        result[i] = child1
        result[i + 1] = child2
    return result


def mutate(rng, pop):
    """Per-locus bit-flip mutation matching Haskell nkMutate.

    Haskell nkMutate: for each locus, generate random int in [0, N-1];
    if it equals 0, flip the bit. This gives rate 1/N per locus.
    """
    rolls = rng.integers(0, GENOME_LENGTH, size=pop.shape)
    mask = rolls == 0
    result = pop.copy()
    result[mask] = 1 - result[mask]
    return result


def compute_diversity(rng, pop, n_pairs=20):
    """Mean pairwise Hamming distance, sampled. Matches IslandGA.hs."""
    n = len(pop)
    if n <= 1:
        return 0.0
    total = 0.0
    for _ in range(n_pairs):
        a = rng.integers(0, n)
        b = rng.integers(0, n - 1)
        if b >= a:
            b += 1
        total += np.sum(pop[a] != pop[b]) / GENOME_LENGTH
    return total / n_pairs


# ---------------------------------------------------------------------------
# Migration (matching IslandGA.hs migrate function)
# ---------------------------------------------------------------------------

def migrate(islands, fitnesses_per_island, adj_list, n_migrants):
    """Migrate top-n_migrants from each neighbor into each island.

    adj_list[i] = list of neighbor indices that send TO island i.
    Replace worst n_incoming individuals in island i with best from neighbors.
    """
    n_islands = len(islands)
    new_islands = [isl.copy() for isl in islands]
    new_fitnesses = [f.copy() for f in fitnesses_per_island]

    for i in range(n_islands):
        neighbors = adj_list[i]
        if not neighbors:
            continue

        # Collect top n_migrants from each neighbor
        incoming_inds = []
        incoming_fits = []
        for nb in neighbors:
            sorted_idx = np.argsort(fitnesses_per_island[nb])[::-1]  # best first
            top_idx = sorted_idx[:n_migrants]
            incoming_inds.append(islands[nb][top_idx])
            incoming_fits.append(fitnesses_per_island[nb][top_idx])

        incoming_inds = np.vstack(incoming_inds)
        incoming_fits = np.concatenate(incoming_fits)
        n_incoming = len(incoming_fits)

        # Replace worst n_incoming in island i
        sorted_idx_i = np.argsort(new_fitnesses[i])  # worst first (ascending)
        replace_idx = sorted_idx_i[:n_incoming]

        new_islands[i][replace_idx] = incoming_inds
        new_fitnesses[i][replace_idx] = incoming_fits

    return new_islands, new_fitnesses


# ---------------------------------------------------------------------------
# Convert networkx graph to adjacency list
# ---------------------------------------------------------------------------

def graph_to_adj_list(G):
    """Convert undirected networkx graph to adjacency list for migration.

    For undirected graphs, neighbors of i are both senders and receivers.
    adj_list[i] = list of nodes that send migrants TO i.
    """
    n = G.number_of_nodes()
    adj_list = [[] for _ in range(n)]
    for i in range(n):
        adj_list[i] = list(G.neighbors(i))
    return adj_list


# ---------------------------------------------------------------------------
# Run one GA replicate
# ---------------------------------------------------------------------------

def run_one(adj_list, seed, fitness_fn, domain_name):
    """Run one GA replicate, return list of (gen, mean_fit, best_fit, diversity)."""
    rng = np.random.default_rng(seed)

    # Initialize islands
    islands = []
    for _ in range(NUM_ISLANDS):
        pop = rng.integers(0, 2, size=(POP_SIZE, GENOME_LENGTH), dtype=np.int8)
        islands.append(pop)

    # Evaluate initial populations
    fitnesses = [fitness_fn(isl) for isl in islands]

    # Initial stats
    all_pop = np.vstack(islands)
    all_fit = np.concatenate(fitnesses)
    div0 = compute_diversity(rng, all_pop)
    stats = [(0, float(all_fit.mean()), float(all_fit.max()), div0)]

    for gen in range(1, TOTAL_GENS + 1):
        # Evolve each island
        new_islands = []
        new_fitnesses = []
        for isl, fit in zip(islands, fitnesses):
            selected = tournament_select(rng, isl, fit)
            crossed = crossover(rng, selected)
            mutated = mutate(rng, crossed)
            new_fit = fitness_fn(mutated)
            new_islands.append(mutated)
            new_fitnesses.append(new_fit)

        islands = new_islands
        fitnesses = new_fitnesses

        # Migrate at interval
        if gen % MIG_INTERVAL == 0:
            islands, fitnesses = migrate(islands, fitnesses, adj_list, NUM_MIGRANTS)

            # Record stats at checkpoints
            all_pop = np.vstack(islands)
            all_fit = np.concatenate(fitnesses)
            div = compute_diversity(rng, all_pop)
            stats.append((gen, float(all_fit.mean()), float(all_fit.max()), div))

    return stats


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------

def run_experiment():
    """Run all graph families on all domains."""
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Build graph list
    family1 = get_family1()
    family2 = get_family2()
    all_graphs = family1 + family2

    domains = [
        ("onemax", onemax_fitness_pop),
        ("nk4", nk4_fitness_pop),
    ]

    total_runs = len(all_graphs) * len(domains) * len(SEEDS)
    count = 0
    t_start = time.time()

    for domain_name, fitness_fn in domains:
        domain_dir = os.path.join(RESULTS_DIR, domain_name)
        os.makedirs(domain_dir, exist_ok=True)

        for graph_name, graph_desc, G in all_graphs:
            adj_list = graph_to_adj_list(G)

            for seed in SEEDS:
                count += 1
                # Use filesystem-safe name
                safe_name = graph_name.replace("β", "b")
                outfile = os.path.join(domain_dir, f"{safe_name}_seed{seed}.csv")

                if os.path.exists(outfile):
                    elapsed = time.time() - t_start
                    print(f"[{count}/{total_runs}] SKIP (exists): {domain_name} {graph_name} seed={seed}")
                    continue

                elapsed = time.time() - t_start
                rate = count / elapsed if elapsed > 0 else 0
                eta = (total_runs - count) / rate if rate > 0 else 0
                print(f"[{count}/{total_runs}] Running: {domain_name} {graph_name} seed={seed} "
                      f"(elapsed={elapsed:.0f}s, ETA={eta:.0f}s)")

                stats = run_one(adj_list, seed, fitness_fn, domain_name)

                # Write CSV
                with open(outfile, 'w', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(['generation', 'meanFitness', 'bestFitness', 'diversity'])
                    for gen, mean_f, best_f, div in stats:
                        writer.writerow([gen, f"{mean_f:.6f}", f"{best_f:.6f}", f"{div:.6f}"])

    elapsed = time.time() - t_start
    print(f"\nAll {total_runs} runs complete in {elapsed:.1f}s.")
    print(f"Results in: {RESULTS_DIR}/")


# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------

def summarize():
    """Print summary statistics from completed results."""
    family1 = get_family1()
    family2 = get_family2()
    all_graphs = family1 + family2

    domains = ["onemax", "nk4"]

    print("\n" + "=" * 80)
    print("  EUMAS Graph Family Experiment — Summary Statistics")
    print("=" * 80)

    for domain_name in domains:
        domain_dir = os.path.join(RESULTS_DIR, domain_name)
        if not os.path.isdir(domain_dir):
            continue

        print(f"\n--- Domain: {domain_name} ---")
        print(f"{'Graph':<10} {'β₁':>3} {'λ₂':>8} {'MeanDiv':>10} {'StdDiv':>8} "
              f"{'MeanBest':>10} {'StdBest':>8} {'n_seeds':>7}")
        print("-" * 75)

        for graph_name, graph_desc, G in all_graphs:
            safe_name = graph_name.replace("β", "b")
            A = nx.to_numpy_array(G, dtype=float)
            b1 = betti_1(G)
            lam2 = lambda_2(A)

            # Read final-generation stats from each seed
            final_divs = []
            final_bests = []

            for seed in SEEDS:
                outfile = os.path.join(domain_dir, f"{safe_name}_seed{seed}.csv")
                if not os.path.exists(outfile):
                    continue
                with open(outfile, 'r') as f:
                    reader = csv.DictReader(f)
                    rows = list(reader)
                    if rows:
                        last = rows[-1]
                        final_divs.append(float(last['diversity']))
                        final_bests.append(float(last['bestFitness']))

            if final_divs:
                mean_div = np.mean(final_divs)
                std_div = np.std(final_divs)
                mean_best = np.mean(final_bests)
                std_best = np.std(final_bests)
                print(f"{graph_name:<10} {b1:>3} {lam2:>8.5f} {mean_div:>10.6f} {std_div:>8.6f} "
                      f"{mean_best:>10.6f} {std_best:>8.6f} {len(final_divs):>7}")
            else:
                print(f"{graph_name:<10} {b1:>3} {lam2:>8.5f} {'(no data)':>10}")

    print()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="EUMAS graph family experiments")
    parser.add_argument("--summary", action="store_true", help="Print summary of existing results")
    parser.add_argument("--seeds", type=int, default=30, help="Number of seeds to use (default: 30)")
    args = parser.parse_args()

    if args.seeds < 30:
        SEEDS = SEEDS[:args.seeds]

    if args.summary:
        summarize()
    else:
        run_experiment()
        summarize()
