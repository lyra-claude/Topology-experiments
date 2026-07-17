#!/usr/bin/env python3
"""
Thin wrapper over run_pilot.py that also records final fitness.

Identical to run_pilot.py in every respect (graph set, domains, seeds, gens,
RNG seeding) — the only additions are two extra columns in the output:

    final_best_fitness   — best individual fitness across all islands at gen 500
    final_mean_fitness   — mean individual fitness across all islands at gen 500

These are already computed by run_cubic_spectrum.run() / run_eumas_families.run_one;
run_pilot.py simply discarded them.  This script captures them.

Output: experiments/results/volume_vs_spectrum/pilot_with_fitness.csv
DO NOT overwrite pilot.csv — this is a separate file.

After the run, compare final_diversity values against pilot.csv to verify
bit-for-bit reproduction (determinism check).
"""

import csv
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import run_cubic_spectrum as RCS  # noqa: E402

DOMAINS = ["onemax", "nk4"]
SEEDS = list(range(20))
OUT_DIR = os.path.join(HERE, "results", "volume_vs_spectrum")
OUT_PATH = os.path.join(OUT_DIR, "pilot_with_fitness.csv")
ORIG_PATH = os.path.join(OUT_DIR, "pilot.csv")

FIELDNAMES = [
    "graph_id",
    "lambda2",
    "domain",
    "seed",
    "final_diversity",
    "final_best_fitness",
    "final_mean_fitness",
]


def load_orig_diversity(path):
    """Load pilot.csv into a dict keyed by (graph_id, domain, seed) -> final_diversity."""
    out = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            key = (row["graph_id"], row["domain"], int(row["seed"]))
            out[key] = float(row["final_diversity"])
    return out


def main():
    graphset = RCS.load_graphset()
    graph_ids = [g["id"] for g in graphset["graphs"]]
    total = len(graph_ids) * len(DOMAINS) * len(SEEDS)

    os.makedirs(OUT_DIR, exist_ok=True)

    t_start = time.time()
    n = 0
    rows = []
    for gid in graph_ids:
        for domain in DOMAINS:
            for seed in SEEDS:
                res = RCS.run(gid, domain, seed, gens=500)
                rows.append({
                    "graph_id": gid,
                    "lambda2": round(res["lambda2"], 4),
                    "domain": domain,
                    "seed": seed,
                    "final_diversity": res["final_diversity"],
                    "final_best_fitness": res["final_best_fitness"],
                    "final_mean_fitness": res["final_mean_fitness"],
                })
                n += 1
                if n % 50 == 0 or n == total:
                    elapsed = time.time() - t_start
                    print(
                        f"  [{n}/{total}] {elapsed:.0f}s elapsed "
                        f"({elapsed/n:.2f}s/run) -- last: {gid} {domain} s{seed} "
                        f"div={res['final_diversity']:.4f} "
                        f"best={res['final_best_fitness']:.4f}",
                        flush=True,
                    )

    with open(OUT_PATH, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES)
        w.writeheader()
        w.writerows(rows)

    print(f"\nWrote {len(rows)} rows to {OUT_PATH}")
    print(f"Total wall-clock: {time.time() - t_start:.0f}s")

    # --- Determinism check ---
    if not os.path.exists(ORIG_PATH):
        print(f"\nWARNING: {ORIG_PATH} not found — cannot run determinism check.")
        return

    orig = load_orig_diversity(ORIG_PATH)
    max_diff = 0.0
    mismatches = 0
    for row in rows:
        key = (row["graph_id"], row["domain"], row["seed"])
        if key not in orig:
            print(f"  MISSING in pilot.csv: {key}")
            mismatches += 1
            continue
        diff = abs(row["final_diversity"] - orig[key])
        if diff > max_diff:
            max_diff = diff

    print(f"\n--- Determinism check ---")
    print(f"  Rows in new CSV:   {len(rows)}")
    print(f"  Rows in pilot.csv: {len(orig)}")
    print(f"  Missing keys:      {mismatches}")
    print(f"  Max |Δ final_diversity|: {max_diff:.2e}")
    if max_diff <= 1e-9:
        print("  PASS — diversity reproduced exactly (max diff <= 1e-9).")
    else:
        print("  FAIL — diversity does NOT reproduce exactly. Fitness columns are NOT valid.")
        print("  Do NOT use pilot_with_fitness.csv. Investigate RNG seeding.")


if __name__ == "__main__":
    main()
