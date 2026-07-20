#!/usr/bin/env python3
"""
Experiment #2 — Direction.

Tests whether migration *direction* affects standing diversity and fitness
independently of the undirected spectral summary λ₂.

Key insight: a hub-in star and a hub-out star have identical underlying
undirected graphs → identical λ₂.  If their diversity outcomes differ,
direction is a real factor that λ₂ cannot see.

Three topologies on n=5 islands (node 0 = hub, nodes 1-4 = leaves):

  star_hub_out   — hub → leaves (hub sends its best INTO every leaf; hub
                   receives nothing).  adj_list[leaf] = [0], adj_list[hub] = []
  star_hub_in    — leaves → hub  (each leaf sends its best INTO the hub; leaves
                   receive nothing).  adj_list[hub] = [1,2,3,4], adj_list[leaf] = []
  star_undirected — symmetric exchange (both directions).  Each leaf and the
                   hub both send to and receive from each other.
                   adj_list[hub] = [1,2,3,4], adj_list[leaf] = [0]

NOTE: do NOT symmetrise adj_list — directedness IS the experiment.

Full sweep (TODO — next session):
    10 seeds × 3 modes × NK4, 500 gens, matching Exp #1 conventions.
    Command:
        python experiments/run_direction.py --sweep --seeds 10 --gens 500 --domain nk4

Future extension (do NOT implement now):
    Richter-Neumann iso-spectral graph families (arXiv 2603.28151) to vary
    hub structure at fixed λ₂.  These let us sweep λ₂-matched directed
    alternatives and cleanly attribute any diversity difference to direction
    vs. volume.
"""

import argparse
import csv
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import run_eumas_families as R  # noqa: E402

# ---------------------------------------------------------------------------
# Parameters (defaults match Exp #1 / EUMAS families — override for smoke)
# ---------------------------------------------------------------------------

NUM_ISLANDS = 5     # 1 hub + 4 leaves
HUB = 0
LEAVES = [1, 2, 3, 4]

DOMAINS = {
    "onemax": R.onemax_fitness_pop,
    "nk0":   R.nk0_fitness_pop,
    "nk1":   R.nk1_fitness_pop,
    "nk2":   R.nk2_fitness_pop,
    "nk4":   R.nk4_fitness_pop,
}

RESULTS_DIR = os.path.join(HERE, "results", "direction")


# ---------------------------------------------------------------------------
# Directed adjacency lists for the three star modes
# ---------------------------------------------------------------------------

def make_adj_hub_out(n_islands=NUM_ISLANDS, hub=HUB):
    """
    star_hub_out: hub → every leaf.
    adj_list[leaf] = [hub]  (leaf *receives* from hub)
    adj_list[hub]  = []     (hub receives from nobody)
    """
    adj = [[] for _ in range(n_islands)]
    for node in range(n_islands):
        if node != hub:
            adj[node] = [hub]   # leaf receives from hub
        # hub receives from nobody
    return adj


def make_adj_hub_in(n_islands=NUM_ISLANDS, hub=HUB):
    """
    star_hub_in: every leaf → hub.
    adj_list[hub]  = [all leaves]  (hub receives from every leaf)
    adj_list[leaf] = []            (leaves receive from nobody)
    """
    adj = [[] for _ in range(n_islands)]
    leaves = [i for i in range(n_islands) if i != hub]
    adj[hub] = leaves   # hub receives from all leaves
    return adj


def make_adj_undirected(n_islands=NUM_ISLANDS, hub=HUB):
    """
    star_undirected: symmetric hub ↔ every leaf.
    adj_list[hub]  = [all leaves]
    adj_list[leaf] = [hub]
    This is the undirected star — same λ₂ as hub_out and hub_in.
    """
    adj = [[] for _ in range(n_islands)]
    leaves = [i for i in range(n_islands) if i != hub]
    adj[hub] = leaves
    for leaf in leaves:
        adj[leaf] = [hub]
    return adj


TOPOLOGIES = {
    "star_hub_out":    make_adj_hub_out,
    "star_hub_in":     make_adj_hub_in,
    "star_undirected": make_adj_undirected,
}


# ---------------------------------------------------------------------------
# Patched run_one for n=5 islands
# ---------------------------------------------------------------------------

def run_one_directed(adj_list, seed, fitness_fn, domain_name, gens=500):
    """
    Island-model GA with directed migration on n=5 islands.

    Reuses R.tournament_select, R.crossover, R.mutate, R.migrate,
    R.compute_diversity verbatim — only the parameter n_islands differs
    from run_eumas_families.run_one (which hard-codes 12).

    Returns list of (gen, mean_fit, best_fit, diversity) at every
    MIG_INTERVAL checkpoint.
    """
    n_islands = len(adj_list)
    pop_size = R.POP_SIZE
    mig_interval = R.MIG_INTERVAL
    n_migrants = R.NUM_MIGRANTS
    genome_length = R.GENOME_LENGTH

    rng = np.random.default_rng(seed)

    # Initialise islands
    islands = []
    for _ in range(n_islands):
        pop = rng.integers(0, 2, size=(pop_size, genome_length), dtype=np.int8)
        islands.append(pop)

    fitnesses = [fitness_fn(isl) for isl in islands]

    # Gen-0 stats
    all_pop = np.vstack(islands)
    all_fit = np.concatenate(fitnesses)
    div0 = R.compute_diversity(rng, all_pop)
    stats = [(0, float(all_fit.mean()), float(all_fit.max()), div0)]

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
            stats.append((gen, float(all_fit.mean()), float(all_fit.max()), div))

    return stats


# ---------------------------------------------------------------------------
# Smoke test — one run per mode, print per-island state after first migration
# ---------------------------------------------------------------------------

def smoke_test(domain="nk4", gens=20, seed=0):
    """
    Run a very short (gens=20) single-seed smoke test for each mode.

    Key behavioural check: after the FIRST migration event (gen=10) —
      hub_out: every leaf should share the hub's best genome (leaves converge)
      hub_in:  only the hub should have changed; leaves should be untouched

    Prints per-island best fitness before and after migration for both modes.
    """
    fitness_fn = DOMAINS[domain]
    print(f"\n=== Smoke test: domain={domain}, seed={seed}, gens={gens} ===\n")

    for topo_name, make_adj in TOPOLOGIES.items():
        adj = make_adj()
        print(f"--- {topo_name} ---")
        print(f"    adj_list: {adj}")

        # Replicate the inner loop up to gen 10 so we can inspect state
        rng = np.random.default_rng(seed)
        islands = []
        for _ in range(NUM_ISLANDS):
            pop = rng.integers(0, 2, size=(R.POP_SIZE, R.GENOME_LENGTH), dtype=np.int8)
            islands.append(pop)

        fitnesses = [fitness_fn(isl) for isl in islands]

        # Pre-migration: record best fitness per island
        pre_best = [float(f.max()) for f in fitnesses]
        print(f"    pre-migration best per island: {[f'{v:.4f}' for v in pre_best]}")

        # Evolve to gen 10 (first migration)
        for gen in range(1, R.MIG_INTERVAL + 1):
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

        # Record state just before migration
        pre_mig_best = [float(f.max()) for f in fitnesses]
        # For hub_out check: hub's best individual
        hub_best_genome = islands[HUB][np.argmax(fitnesses[HUB])].copy()

        # Apply migration
        islands_post, fitnesses_post = R.migrate(
            islands, fitnesses, adj, R.NUM_MIGRANTS
        )

        post_mig_best = [float(f.max()) for f in fitnesses_post]
        print(f"    pre-mig  best per island: {[f'{v:.4f}' for v in pre_mig_best]}")
        print(f"    post-mig best per island: {[f'{v:.4f}' for v in post_mig_best]}")

        # Directed behavioural checks
        if topo_name == "star_hub_out":
            # Every leaf should now contain the hub's best individual
            hub_best_fit = float(fitnesses[HUB].max())
            leaf_best_fits_post = [float(fitnesses_post[l].max()) for l in LEAVES]
            hub_best_among_leaves = all(
                f >= hub_best_fit - 1e-9 for f in leaf_best_fits_post
            )
            print(f"    CHECK hub_out: hub best={hub_best_fit:.4f}, "
                  f"all leaves now >= hub best? {hub_best_among_leaves}")
            # Verify hub is unchanged
            hub_changed = not np.array_equal(islands[HUB], islands_post[HUB])
            print(f"    CHECK hub_out: hub population changed? {hub_changed} "
                  f"(expected: False)")

        elif topo_name == "star_hub_in":
            # Hub should hold migrants from leaves; leaves should be UNCHANGED
            leaves_unchanged = all(
                np.array_equal(islands[l], islands_post[l]) for l in LEAVES
            )
            hub_changed = not np.array_equal(islands[HUB], islands_post[HUB])
            print(f"    CHECK hub_in: leaves unchanged? {leaves_unchanged} "
                  f"(expected: True)")
            print(f"    CHECK hub_in: hub changed? {hub_changed} "
                  f"(expected: True)")

        elif topo_name == "star_undirected":
            # Both hub and leaves should change
            hub_changed = not np.array_equal(islands[HUB], islands_post[HUB])
            any_leaf_changed = any(
                not np.array_equal(islands[l], islands_post[l]) for l in LEAVES
            )
            print(f"    CHECK undirected: hub changed? {hub_changed} "
                  f"(expected: True)")
            print(f"    CHECK undirected: any leaf changed? {any_leaf_changed} "
                  f"(expected: True)")

        # Also run the full short experiment to confirm no errors
        stats = run_one_directed(adj, seed, fitness_fn, domain_name=domain, gens=gens)
        final = stats[-1]
        print(f"    Full run OK: gen={final[0]}, mean_fit={final[1]:.4f}, "
              f"best_fit={final[2]:.4f}, diversity={final[3]:.4f}")
        print()

    print("=== Smoke test complete — all three modes executed without error ===")


# ---------------------------------------------------------------------------
# Full sweep runner
# ---------------------------------------------------------------------------

def run_sweep(domain="nk4", seeds=None, gens=500):
    """
    Full sweep: all 3 modes × n_seeds × domain, save CSVs.

    TODO (next session): run with --sweep --seeds 10 --gens 500 --domain nk4
    to mirror Exp #1's cubic12 / pilot_with_fitness structure.
    """
    if seeds is None:
        seeds = list(range(10))

    os.makedirs(RESULTS_DIR, exist_ok=True)
    fitness_fn = DOMAINS[domain]
    total = len(TOPOLOGIES) * len(seeds)
    n = 0
    t_start = time.time()
    all_rows = []

    for topo_name, make_adj in TOPOLOGIES.items():
        adj = make_adj()
        topo_dir = os.path.join(RESULTS_DIR, domain, topo_name)
        os.makedirs(topo_dir, exist_ok=True)

        for seed in seeds:
            n += 1
            outfile = os.path.join(topo_dir, f"seed{seed}.csv")

            if os.path.exists(outfile):
                print(f"[{n}/{total}] SKIP (exists): {topo_name} {domain} seed={seed}")
                continue

            elapsed = time.time() - t_start
            print(f"[{n}/{total}] {topo_name} {domain} seed={seed} "
                  f"(elapsed={elapsed:.0f}s)", flush=True)

            stats = run_one_directed(adj, seed, fitness_fn, domain_name=domain, gens=gens)

            with open(outfile, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["generation", "meanFitness", "bestFitness", "diversity"])
                for gen, mean_f, best_f, div in stats:
                    writer.writerow([gen, f"{mean_f:.6f}", f"{best_f:.6f}", f"{div:.6f}"])

            final = stats[-1]
            all_rows.append({
                "topology": topo_name,
                "domain": domain,
                "seed": seed,
                "final_gen": final[0],
                "final_mean_fitness": final[1],
                "final_best_fitness": final[2],
                "final_diversity": final[3],
            })

    print(f"\nSweep complete: {n} runs in {time.time() - t_start:.0f}s")
    print(f"Results in: {RESULTS_DIR}/")

    # Print summary table
    if all_rows:
        print("\n--- Summary (final diversity by topology) ---")
        for topo_name in TOPOLOGIES:
            rows = [r for r in all_rows if r["topology"] == topo_name]
            if rows:
                divs = [r["final_diversity"] for r in rows]
                fits = [r["final_best_fitness"] for r in rows]
                print(f"  {topo_name:<20} div={np.mean(divs):.4f}±{np.std(divs):.4f}  "
                      f"best_fit={np.mean(fits):.4f}±{np.std(fits):.4f}  (n={len(rows)})")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description="Experiment #2: directed star migration on n=5 islands"
    )
    p.add_argument("--smoke", action="store_true",
                   help="Run smoke test (one short run per mode, print checks)")
    p.add_argument("--sweep", action="store_true",
                   help="Run full sweep (TODO: next session)")
    p.add_argument("--domain", default="nk4", choices=list(DOMAINS),
                   help="Fitness landscape (default: nk4)")
    p.add_argument("--seeds", type=int, default=10,
                   help="Number of seeds for sweep (default: 10)")
    p.add_argument("--gens", type=int, default=500,
                   help="Generations per run (default: 500; smoke uses 20)")
    p.add_argument("--seed", type=int, default=0,
                   help="Single seed for smoke test (default: 0)")
    args = p.parse_args()

    if args.smoke:
        smoke_test(domain=args.domain, seed=args.seed)
    elif args.sweep:
        run_sweep(domain=args.domain, seeds=list(range(args.seeds)), gens=args.gens)
    else:
        p.print_help()


if __name__ == "__main__":
    main()
