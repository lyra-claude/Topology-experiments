#!/usr/bin/env python3
"""
run_diverse_seed.py — Diverse-Seed Hub-Out Experiment.

Follow-up to the Direction Experiment (run_direction.py). Breaks the confound
between "topology generates diverse hubs" and "diverse hubs transmit diversity
to leaves."

MANIPULATION (the ONLY change vs run_direction.run_one_directed):
    At generation 0, seed the HUB island random-uniform (diverse) but seed each
    LEAF island as a homogeneous clone — 50 copies of a single random genotype,
    so each leaf starts with ZERO within-island diversity. Everything downstream
    (selection, crossover, mutation, migration, RNG stream ordering) is identical
    to the direction experiment.

ESTIMAND (NEW — differs from direction experiment):
    We log LEAF WITHIN-ISLAND diversity per generation: the EXACT mean pairwise
    normalised Hamming distance computed WITHIN each leaf island separately, then
    averaged over the 4 leaves. This is NOT the pooled diversity of the direction
    experiment (which was between + within combined). It is a per-island (within)
    quantity, computed exactly (all C(50,2) pairs), not sampled.

    We ALSO log the pooled sampled diversity (compute_diversity on vstack) to keep
    the RNG stream identical to run_direction and to allow a cross-check. The
    per-generation hub within-diversity and mean between-island diversity are also
    logged for context. NONE of these extra logs consume RNG (they are computed
    exactly / deterministically), so the run reproduces run_direction's RNG stream
    bit-for-bit up to the seeding difference.

Conditions: star_hub_out, star_hub_in, star_undirected (reference).
Domains: nk4 (primary), onemax (reference).
Seeds: 0-29 (same as direction experiment). Gens: 500. Checkpoints every 10.

Output: experiments/results/diverse_seed/{domain}/{topo}/seed{s}.csv
    columns: generation, meanFitness, bestFitness, pooledSampledDiv,
             leafWithinDiv, hubWithinDiv, betweenDiv
"""

import argparse
import csv
import itertools
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import run_eumas_families as R  # noqa: E402
import run_direction as D       # noqa: E402

NUM_ISLANDS = D.NUM_ISLANDS      # 5
HUB = D.HUB                      # 0
LEAVES = D.LEAVES               # [1,2,3,4]
DOMAINS = D.DOMAINS
TOPOLOGIES = D.TOPOLOGIES

RESULTS_DIR = os.path.join(HERE, "results", "diverse_seed")


# ---------------------------------------------------------------------------
# Exact within-island diversity (all C(n,2) pairs, NOT sampled, consumes no RNG)
# ---------------------------------------------------------------------------

def exact_mean_pairwise(pop):
    """Exact mean normalised Hamming distance over all C(n,2) pairs of pop."""
    n = len(pop)
    if n <= 1:
        return 0.0
    total = 0.0
    npairs = 0
    for a in range(n):
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
        for a in range(len(A)):
            diff = np.sum(B != A[a], axis=1)
            total += diff.sum()
            npairs += len(diff)
    return (total / R.GENOME_LENGTH) / npairs


def record(gen, islands, fitnesses, samp_div):
    """Build one stats row. Exact within/between computed deterministically."""
    all_fit = np.concatenate(fitnesses)
    leaf_within = float(np.mean([exact_mean_pairwise(islands[l]) for l in LEAVES]))
    hub_within = float(exact_mean_pairwise(islands[HUB]))
    between = exact_between_islands(islands)
    return (
        gen,
        float(all_fit.mean()),
        float(all_fit.max()),
        float(samp_div),
        leaf_within,
        hub_within,
        float(between),
    )


# ---------------------------------------------------------------------------
# Diverse-seed initialisation
# ---------------------------------------------------------------------------

def init_islands_diverse_seed(rng):
    """
    Gen-0 seeding: hub random-uniform (diverse), each leaf a homogeneous clone
    (50 copies of one random genotype).

    RNG ORDERING NOTE: run_direction draws islands in order 0..4, each a full
    (50,100) integers draw. We preserve that exact draw ORDER and SHAPE: for
    each island index we draw a (50,100) block, then for leaves we OVERWRITE the
    block with 50 copies of its own first row. This keeps the RNG stream aligned
    with run_direction through gen 0 (the hub, drawn first at index 0, is
    bit-identical to run_direction's hub; leaves consume the same number of RNG
    draws but are collapsed to clones).
    """
    islands = []
    for node in range(NUM_ISLANDS):
        pop = rng.integers(0, 2, size=(R.POP_SIZE, R.GENOME_LENGTH), dtype=np.int8)
        if node != HUB:
            # collapse leaf to homogeneous clone of its own first genotype
            pop = np.tile(pop[0], (R.POP_SIZE, 1)).astype(np.int8)
        islands.append(pop)
    return islands


# ---------------------------------------------------------------------------
# Run one replicate (mirrors run_direction.run_one_directed, diverse seeding)
# ---------------------------------------------------------------------------

def run_one(adj_list, seed, fitness_fn, gens=500):
    n_islands = len(adj_list)
    assert n_islands == NUM_ISLANDS
    mig_interval = R.MIG_INTERVAL
    n_migrants = R.NUM_MIGRANTS

    rng = np.random.default_rng(seed)

    islands = init_islands_diverse_seed(rng)
    fitnesses = [fitness_fn(isl) for isl in islands]

    # Gen-0 stats. compute_diversity consumes RNG identically to run_direction.
    all_pop = np.vstack(islands)
    samp_div0 = R.compute_diversity(rng, all_pop)
    stats = [record(0, islands, fitnesses, samp_div0)]

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
            samp_div = R.compute_diversity(rng, all_pop)
            stats.append(record(gen, islands, fitnesses, samp_div))

    return stats


CSV_HEADER = ["generation", "meanFitness", "bestFitness", "pooledSampledDiv",
              "leafWithinDiv", "hubWithinDiv", "betweenDiv"]


def write_csv(outfile, stats):
    with open(outfile, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(CSV_HEADER)
        for row in stats:
            gen = row[0]
            vals = [f"{v:.6f}" for v in row[1:]]
            w.writerow([gen] + vals)


# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------

def run_sweep(domain, seeds, gens, topos):
    fitness_fn = DOMAINS[domain]
    total = len(topos) * len(seeds)
    n = 0
    t0 = time.time()
    for topo_name in topos:
        adj = TOPOLOGIES[topo_name]()
        topo_dir = os.path.join(RESULTS_DIR, domain, topo_name)
        os.makedirs(topo_dir, exist_ok=True)
        for seed in seeds:
            n += 1
            outfile = os.path.join(topo_dir, f"seed{seed}.csv")
            if os.path.exists(outfile):
                print(f"[{n}/{total}] SKIP (exists): {topo_name} {domain} seed={seed}")
                continue
            print(f"[{n}/{total}] {topo_name} {domain} seed={seed} "
                  f"(elapsed={time.time()-t0:.0f}s)", flush=True)
            stats = run_one(adj, seed, fitness_fn, gens=gens)
            write_csv(outfile, stats)
    print(f"\nSweep complete: {n} runs in {time.time()-t0:.0f}s -> {RESULTS_DIR}/")


def main():
    p = argparse.ArgumentParser(description="Diverse-seed hub_out experiment")
    p.add_argument("--sweep", action="store_true")
    p.add_argument("--domain", default="nk4", choices=list(DOMAINS))
    p.add_argument("--seeds", type=int, default=30)
    p.add_argument("--gens", type=int, default=500)
    p.add_argument("--topos", default="star_hub_out,star_hub_in,star_undirected")
    args = p.parse_args()

    topos = [t.strip() for t in args.topos.split(",") if t.strip()]
    seeds = list(range(args.seeds))

    if args.sweep:
        run_sweep(args.domain, seeds, args.gens, topos)
    else:
        p.print_help()


if __name__ == "__main__":
    main()
