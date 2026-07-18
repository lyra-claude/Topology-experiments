#!/usr/bin/env python3
"""
decompose_direction.py — Within/between-island decomposition of the POOLED
final-generation diversity from Direction Experiment #2.

METRIC (from run_eumas_families.compute_diversity, called on np.vstack(islands)):
    Pooled diversity = mean pairwise (normalised) Hamming distance over the
    pooled population of all 5 islands, SAMPLED with n_pairs=20 using the run's
    own rng stream. This is an unbiased sampled estimator of the EXACT pooled
    mean pairwise Hamming distance over all C(N,2) pairs (N=250).

DECOMPOSITION IDENTITY (AMOVA / pair-count identity for mean pairwise distance):
    Let D_ij = normalised Hamming distance between individuals i,j.
    Pooled exact mean = ( sum over all pairs D_ij ) / C(N,2).
    Partition pairs into same-island and different-island:
        within_sum   = sum of D over same-island pairs
        between_sum  = sum of D over different-island pairs
    within_mean  = within_sum  / (#same-island pairs)   [averaged over all within pairs]
    between_mean = between_sum / (#diff-island pairs)
    w = (#same-island pairs) / C(N,2)
    Pooled_exact = w * within_mean + (1-w) * between_mean   <-- reconstruction gate

    With 5 islands of 50: #same = 5*C(50,2)=6125, #diff = C(250,2)-6125=24925,
    total = C(250,2)=31125, w = 6125/31125 = 0.19678.

We ALSO reproduce the sampled (n_pairs=20) pooled number bit-for-bit by replaying
the exact rng call sequence, to prove the re-run matches the published CSVs, then
report the EXACT pooled (recommended, noise-free) alongside.

This re-runs each (landscape, topology, seed) with the SAME seed/config as the
published run and captures per-island final populations. It writes:
    experiments/results/direction/within_between_decomposition.csv
"""

import csv
import itertools
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import run_eumas_families as R  # noqa: E402
import run_direction as D       # noqa: E402

RESULTS_DIR = os.path.join(HERE, "results", "direction")
OUT_CSV = os.path.join(RESULTS_DIR, "within_between_decomposition.csv")

SEEDS = list(range(30))
GENS = 500
DOMAINS = ["nk4", "onemax"]
TOPOLOGIES = ["star_hub_out", "star_hub_in", "star_undirected"]


def exact_mean_pairwise(pop):
    """Exact mean normalised Hamming distance over all C(n,2) pairs of pop."""
    n = len(pop)
    if n <= 1:
        return 0.0
    total = 0.0
    npairs = 0
    for a in range(n):
        # vectorised: distances from a to all b>a
        diff = np.sum(pop[a + 1:] != pop[a], axis=1)
        total += diff.sum()
        npairs += len(diff)
    return (total / R.GENOME_LENGTH) / npairs


def exact_between_islands(islands):
    """Exact mean normalised Hamming distance over all different-island pairs."""
    total = 0.0
    npairs = 0
    for ia, ib in itertools.combinations(range(len(islands)), 2):
        A = islands[ia]
        B = islands[ib]
        # all cross pairs: |A| x |B|
        for a in range(len(A)):
            diff = np.sum(B != A[a], axis=1)
            total += diff.sum()
            npairs += len(diff)
    return (total / R.GENOME_LENGTH) / npairs


def run_one_capture(adj_list, seed, fitness_fn, gens=GENS):
    """
    Replicates run_direction.run_one_directed EXACTLY (same rng calls),
    but returns the per-island final populations AND the final sampled
    diversity (to prove bit-for-bit reproduction of the published CSV).

    Returns: (islands_final, final_sampled_div, final_mean_fit, final_best_fit)
    """
    n_islands = len(adj_list)
    pop_size = R.POP_SIZE
    mig_interval = R.MIG_INTERVAL
    n_migrants = R.NUM_MIGRANTS
    genome_length = R.GENOME_LENGTH

    rng = np.random.default_rng(seed)

    islands = []
    for _ in range(n_islands):
        pop = rng.integers(0, 2, size=(pop_size, genome_length), dtype=np.int8)
        islands.append(pop)

    fitnesses = [fitness_fn(isl) for isl in islands]

    # Gen-0 stats (rng consumed identically to run_one_directed)
    all_pop = np.vstack(islands)
    all_fit = np.concatenate(fitnesses)
    _ = R.compute_diversity(rng, all_pop)

    last_sampled_div = None
    last_mean_fit = None
    last_best_fit = None

    for gen in range(1, gens + 1):
        new_islands = []
        new_fitnesses = []
        for isl, fit in zip(islands, fitnesses):
            selected = R.tournament_select(rng, isl, fit)
            crossed = R.crossover(rng, selected)
            mutated = R.mutate(rng, crossed)
            new_fit = fitness_fn(mutated)
            new_islands.append(mutated)
            new_fitnesses.append(new_fit)

        islands = new_islands
        fitnesses = new_fitnesses

        if gen % mig_interval == 0:
            islands, fitnesses = R.migrate(islands, fitnesses, adj_list, n_migrants)
            all_pop = np.vstack(islands)
            all_fit = np.concatenate(fitnesses)
            div = R.compute_diversity(rng, all_pop)
            last_sampled_div = div
            last_mean_fit = float(all_fit.mean())
            last_best_fit = float(all_fit.max())

    # islands here are the FINAL post-migration populations at gen=500,
    # matching the state whose sampled diversity was recorded.
    return islands, last_sampled_div, last_mean_fit, last_best_fit


def load_published_final(domain, topo, seed):
    """Return published (diversity, bestFitness) from stored CSV last row."""
    fp = os.path.join(RESULTS_DIR, domain, topo, f"seed{seed}.csv")
    with open(fp, newline="") as f:
        rows = list(csv.DictReader(f))
    last = rows[-1]
    return float(last["diversity"]), float(last["bestFitness"])


def main():
    rows = []
    n_islands = D.NUM_ISLANDS
    isl_size = R.POP_SIZE
    N = n_islands * isl_size
    same_pairs = n_islands * (isl_size * (isl_size - 1) // 2)
    total_pairs = N * (N - 1) // 2
    diff_pairs = total_pairs - same_pairs
    w = same_pairs / total_pairs
    print(f"N={N}  same_pairs={same_pairs}  diff_pairs={diff_pairs}  "
          f"total_pairs={total_pairs}  w={w:.6f}")

    max_repro_err = 0.0
    max_recon_err = 0.0
    n_repro_mismatch = 0

    for domain in DOMAINS:
        fitness_fn = D.DOMAINS[domain]
        for topo in TOPOLOGIES:
            make_adj = D.TOPOLOGIES[topo]
            adj = make_adj()
            for seed in SEEDS:
                islands, samp_div, mean_fit, best_fit = run_one_capture(
                    adj, seed, fitness_fn, gens=GENS
                )
                pub_div, pub_fit = load_published_final(domain, topo, seed)

                # bit-for-bit reproduction check (sampled diversity, 6 dp in CSV)
                repro_err = abs(samp_div - pub_div)
                fit_err = abs(best_fit - pub_fit)
                if repro_err > 5e-7 or fit_err > 5e-7:
                    n_repro_mismatch += 1
                    print(f"  REPRO MISMATCH {domain}/{topo}/seed{seed}: "
                          f"samp={samp_div:.6f} pub={pub_div:.6f} "
                          f"fit={best_fit:.6f} pubfit={pub_fit:.6f}")
                max_repro_err = max(max_repro_err, repro_err, fit_err)

                # exact decomposition
                within_each = [exact_mean_pairwise(isl) for isl in islands]
                within_mean = float(np.mean(within_each))
                between_mean = exact_between_islands(islands)
                pooled_exact = exact_mean_pairwise(np.vstack(islands))
                recon = w * within_mean + (1.0 - w) * between_mean
                recon_err = abs(recon - pooled_exact)
                max_recon_err = max(max_recon_err, recon_err)

                rows.append({
                    "domain": domain,
                    "topology": topo,
                    "seed": seed,
                    "within_island": within_mean,
                    "between_island": between_mean,
                    "pooled_exact": pooled_exact,
                    "pooled_reconstructed": recon,
                    "pooled_sampled_repro": samp_div,
                    "pooled_published": pub_div,
                    "final_best_fitness": best_fit,
                    "final_mean_fitness": mean_fit,
                })
            print(f"done {domain}/{topo}")

    os.makedirs(RESULTS_DIR, exist_ok=True)
    fieldnames = ["domain", "topology", "seed", "within_island", "between_island",
                  "pooled_exact", "pooled_reconstructed", "pooled_sampled_repro",
                  "pooled_published", "final_best_fitness", "final_mean_fitness"]
    with open(OUT_CSV, "w", newline="") as f:
        wri = csv.DictWriter(f, fieldnames=fieldnames)
        wri.writeheader()
        for r in rows:
            wri.writerow({k: (f"{v:.8f}" if isinstance(v, float) else v)
                          for k, v in r.items()})

    print(f"\nWrote {len(rows)} rows to {OUT_CSV}")
    print(f"Max |sampled_repro - published| (diversity & fitness): {max_repro_err:.2e}")
    print(f"  reproduction mismatches (>5e-7): {n_repro_mismatch}")
    print(f"Max |reconstructed - pooled_exact| (AMOVA identity): {max_recon_err:.2e}")


if __name__ == "__main__":
    main()
