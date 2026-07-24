#!/usr/bin/env python3
"""
stats_diverse_seed.py — Per-seed inferential stats for the diverse-seed experiment.

Computes per-seed floor and decay slope for leafWithinDiv (and hubWithinDiv at gen 500),
then runs Mann-Whitney U + Cliff's delta comparisons with Bonferroni correction.

Saves results/diverse_seed/STATS.md and prints a summary.
"""
import csv
import os

import numpy as np
from scipy.stats import mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results", "diverse_seed")
TOPOLOGIES = ["star_hub_out", "star_hub_in", "star_undirected"]
DOMAINS = ["nk4", "onemax"]


def load_seed_data(domain, topo):
    """Return list of {seed, gens, leaf_div, hub_div} dicts for all 30 seeds."""
    topo_dir = os.path.join(RESULTS_DIR, domain, topo)
    seeds = []
    for fn in sorted(os.listdir(topo_dir),
                     key=lambda s: int(s.replace("seed", "").replace(".csv", ""))):
        if not fn.endswith(".csv"):
            continue
        with open(os.path.join(topo_dir, fn), newline="") as f:
            data = list(csv.DictReader(f))
        gens = np.array([int(r["generation"]) for r in data])
        leaf_div = np.array([float(r["leafWithinDiv"]) for r in data])
        hub_div = np.array([float(r["hubWithinDiv"]) for r in data])
        seed_num = int(fn.replace("seed", "").replace(".csv", ""))
        seeds.append({"seed": seed_num, "gens": gens, "leaf_div": leaf_div, "hub_div": hub_div})
    return seeds


def per_seed_metrics(seeds):
    """
    For each seed, compute:
      - floor: leafWithinDiv at gen 500 (final checkpoint)
      - peak_val, peak_gen: max leafWithinDiv and its generation
      - decay_slope: linear regression slope from peak_gen to gen 500
      - frac_drop: (peak - floor) / peak
      - hub_floor: hubWithinDiv at gen 500
    Returns list of dicts, one per seed.
    """
    metrics = []
    for s in seeds:
        gens = s["gens"]
        ld = s["leaf_div"]
        hd = s["hub_div"]

        # floor at gen 500
        final_idx = np.where(gens == 500)[0]
        if len(final_idx) == 0:
            raise ValueError(f"Seed {s['seed']}: no gen 500 checkpoint")
        floor_val = float(ld[final_idx[0]])
        hub_floor = float(hd[final_idx[0]])

        # per-seed peak
        peak_idx = int(np.argmax(ld))
        peak_gen = int(gens[peak_idx])
        peak_val = float(ld[peak_idx])

        # decay slope: OLS on [peak_gen .. 500]
        mask = gens >= peak_gen
        decay_gens = gens[mask].astype(float)
        decay_ld = ld[mask]
        if len(decay_gens) >= 2:
            slope = float(np.polyfit(decay_gens, decay_ld, 1)[0])
        else:
            slope = float("nan")

        frac_drop = (peak_val - floor_val) / peak_val if peak_val > 0 else float("nan")

        metrics.append({
            "seed": s["seed"],
            "floor": floor_val,
            "peak_val": peak_val,
            "peak_gen": peak_gen,
            "decay_slope": slope,
            "frac_drop": frac_drop,
            "hub_floor": hub_floor,
        })
    return metrics


def cliffs_delta(a, b):
    """
    Cliff's delta: (# pairs where a > b  -  # pairs where a < b) / (n_a * n_b).
    Positive means a tends to be larger than b.
    """
    a = np.asarray(a)
    b = np.asarray(b)
    n_a, n_b = len(a), len(b)
    greater = np.sum(a[:, None] > b[None, :])
    less = np.sum(a[:, None] < b[None, :])
    return (greater - less) / (n_a * n_b)


def mann_whitney(a, b):
    """Two-sided Mann-Whitney U; returns (U, p)."""
    return mannwhitneyu(a, b, alternative="two-sided")


def compare(name_a, vals_a, name_b, vals_b, alpha_bonf):
    """
    Run Mann-Whitney + Cliff's delta for vals_a vs vals_b.
    Returns a result dict.
    """
    a = np.array(vals_a)
    b = np.array(vals_b)
    U, p = mann_whitney(a, b)
    d = cliffs_delta(a, b)
    med_a = float(np.median(a))
    med_b = float(np.median(b))
    sig = p < alpha_bonf
    # direction: positive d means a > b
    direction = f"{name_a} > {name_b}" if d > 0 else f"{name_a} < {name_b}"
    return {
        "comparison": f"{name_a} vs {name_b}",
        "med_a": med_a,
        "med_b": med_b,
        "U": U,
        "p_raw": p,
        "p_bonf_alpha": alpha_bonf,
        "significant": sig,
        "cliffs_d": d,
        "direction": direction,
    }


def format_row(r):
    sig_str = "YES" if r["significant"] else "no"
    return (f"| {r['comparison']} | {r['med_a']:.5f} | {r['med_b']:.5f} | "
            f"{r['p_raw']:.4e} | {r['p_bonf_alpha']:.4e} | "
            f"{r['cliffs_d']:+.3f} | {sig_str} | {r['direction']} |")


def analyze_domain(domain, md_lines, n_primary_comparisons=5):
    """Full analysis for one domain."""
    md_lines.append(f"\n## Domain: {domain}\n")

    # Load per-seed metrics for all 3 topologies
    topo_metrics = {}
    for topo in TOPOLOGIES:
        seeds_data = load_seed_data(domain, topo)
        topo_metrics[topo] = per_seed_metrics(seeds_data)

    # Summaries
    md_lines.append("### Per-seed metric summaries\n")
    md_lines.append("| Condition | floor median | floor IQR | decay_slope median | frac_drop median |")
    md_lines.append("|---|---|---|---|---|")
    for topo in TOPOLOGIES:
        m = topo_metrics[topo]
        floors = [x["floor"] for x in m]
        slopes = [x["decay_slope"] for x in m]
        fracs = [x["frac_drop"] for x in m]
        q25_f, q75_f = np.percentile(floors, [25, 75])
        md_lines.append(
            f"| {topo} | {np.median(floors):.5f} | [{q25_f:.5f}, {q75_f:.5f}] | "
            f"{np.median(slopes):.4e} | {np.median(fracs):.3f} |"
        )

    # --- Floor comparisons (3 pairs) ---
    floor_comparisons = [
        ("hub_out", "star_hub_out", "hub_in", "star_hub_in"),
        ("hub_out", "star_hub_out", "undirected", "star_undirected"),
        ("hub_in", "star_hub_in", "undirected", "star_undirected"),
    ]
    # --- Decay slope comparisons (2 pairs) ---
    slope_comparisons = [
        ("hub_out", "star_hub_out", "hub_in", "star_hub_in"),
        ("hub_out", "star_hub_out", "undirected", "star_undirected"),
    ]
    # Bonferroni over 5 comparisons total within each domain
    n_comp = n_primary_comparisons
    alpha_bonf = 0.05 / n_comp

    md_lines.append(f"\n### Floor (leafWithinDiv at gen 500) comparisons")
    md_lines.append(f"Bonferroni alpha = 0.05 / {n_comp} = {alpha_bonf:.4e}\n")
    md_lines.append("| Comparison | med_A | med_B | p_raw | p_bonf_alpha | Cliff's d | sig? | direction |")
    md_lines.append("|---|---|---|---|---|---|---|---|")
    floor_results = []
    for na, ta, nb, tb in floor_comparisons:
        fa = [x["floor"] for x in topo_metrics[ta]]
        fb = [x["floor"] for x in topo_metrics[tb]]
        r = compare(na, fa, nb, fb, alpha_bonf)
        floor_results.append(r)
        md_lines.append(format_row(r))

    md_lines.append(f"\n### Decay slope (peak→500) comparisons")
    md_lines.append(f"Bonferroni alpha = 0.05 / {n_comp} = {alpha_bonf:.4e}\n")
    md_lines.append("| Comparison | med_A | med_B | p_raw | p_bonf_alpha | Cliff's d | sig? | direction |")
    md_lines.append("|---|---|---|---|---|---|---|---|")
    slope_results = []
    for na, ta, nb, tb in slope_comparisons:
        sa = [x["decay_slope"] for x in topo_metrics[ta]]
        sb = [x["decay_slope"] for x in topo_metrics[tb]]
        r = compare(na, sa, nb, sb, alpha_bonf)
        slope_results.append(r)
        md_lines.append(format_row(r))

    # Hub within-diversity at gen 500
    md_lines.append(f"\n### Hub within-island diversity at gen 500 (hub_out vs hub_in)")
    md_lines.append("| Comparison | med_A | med_B | p_raw | Cliff's d | direction |")
    md_lines.append("|---|---|---|---|---|---|")
    hub_a = [x["hub_floor"] for x in topo_metrics["star_hub_out"]]
    hub_b = [x["hub_floor"] for x in topo_metrics["star_hub_in"]]
    Uh, ph = mann_whitney(hub_a, hub_b)
    dh = cliffs_delta(hub_a, hub_b)
    mh_a, mh_b = float(np.median(hub_a)), float(np.median(hub_b))
    dir_h = "hub_out > hub_in" if dh > 0 else "hub_out < hub_in"
    md_lines.append(f"| hub_out vs hub_in (hub div) | {mh_a:.5f} | {mh_b:.5f} | "
                    f"{ph:.4e} | {dh:+.3f} | {dir_h} |")

    return {
        "floor_results": floor_results,
        "slope_results": slope_results,
        "hub_hub": {"med_a": mh_a, "med_b": mh_b, "p": ph, "d": dh, "dir": dir_h},
        "topo_metrics": topo_metrics,
        "alpha_bonf": alpha_bonf,
    }


def main():
    md = ["# Diverse-Seed Experiment — Inferential Statistics\n",
          "Generated by `experiments/stats_diverse_seed.py`.\n",
          "Estimand: `leafWithinDiv` = mean pairwise Hamming within each leaf, "
          "averaged over 4 leaves. Stats computed PER SEED (n=30), not from "
          "seed-averaged trajectories.\n",
          "Bonferroni correction applied across 5 comparisons per domain "
          "(3 floor + 2 decay-slope). Hub div comparison is exploratory "
          "(no Bonferroni applied).\n",
          "Cliff's delta sign convention: positive = group A > group B.\n"]

    domain_results = {}
    for domain in DOMAINS:
        res = analyze_domain(domain, md)
        domain_results[domain] = res

    # Plain-English crux summaries
    md.append("\n## Crux Comparisons — Plain-English Summary\n")
    for domain in DOMAINS:
        res = domain_results[domain]
        fr = res["floor_results"]
        sr = res["slope_results"]
        ab = res["alpha_bonf"]

        # floor hub_out vs hub_in
        r = fr[0]
        sig_str = f"**significant after Bonferroni** (p={r['p_raw']:.3e} < {ab:.3e})" if r["significant"] \
            else f"NOT significant after Bonferroni (p={r['p_raw']:.3e} >= {ab:.3e})"
        md.append(f"**[{domain}] Floor hub_out vs hub_in:** {r['direction']} "
                  f"(medians {r['med_a']:.5f} vs {r['med_b']:.5f}, "
                  f"Cliff's d={r['cliffs_d']:+.3f}); {sig_str}.\n")

        # decay slope hub_out vs hub_in
        r = sr[0]
        sig_str = f"**significant after Bonferroni** (p={r['p_raw']:.3e} < {ab:.3e})" if r["significant"] \
            else f"NOT significant after Bonferroni (p={r['p_raw']:.3e} >= {ab:.3e})"
        md.append(f"**[{domain}] Decay slope hub_out vs hub_in:** {r['direction']} "
                  f"(medians {r['med_a']:.5f} vs {r['med_b']:.5f}, "
                  f"Cliff's d={r['cliffs_d']:+.3f}); {sig_str}.\n")

    out_path = os.path.join(RESULTS_DIR, "STATS.md")
    with open(out_path, "w") as f:
        f.write("\n".join(md))
    print(f"Wrote {out_path}")

    # Also print summary table to stdout
    print("\n=== CRUX TABLE ===")
    print(f"{'Comparison':<45} {'med_A':>10} {'med_B':>10} {'p_raw':>12} {'Cliff_d':>9} {'sig?':>5}")
    for domain in DOMAINS:
        res = domain_results[domain]
        for label, r in [("floor", res["floor_results"][0]),
                         ("slope", res["slope_results"][0])]:
            print(f"[{domain}] {label} hub_out vs hub_in: "
                  f"med_A={r['med_a']:.5f} med_B={r['med_b']:.5f} "
                  f"p={r['p_raw']:.4e} d={r['cliffs_d']:+.3f} "
                  f"sig={'YES' if r['significant'] else 'no'}")
        hh = res["hub_hub"]
        print(f"[{domain}] hub_div floor hub_out vs hub_in: "
              f"med_A={hh['med_a']:.5f} med_B={hh['med_b']:.5f} "
              f"p={hh['p']:.4e} d={hh['d']:+.3f}")


if __name__ == "__main__":
    main()
