#!/usr/bin/env python3
"""
nk_grading_analysis.py — Per-K hub_out-vs-hub_in Cliff's delta analysis.

Reuses EXACTLY the per_seed_metrics, cliffs_delta, and mann_whitney definitions
from stats_diverse_seed.py. Do not alter those functions.

For each K in {0, 1, 2, 4}, computes hub_out vs hub_in (a=hub_out, b=hub_in):
  PRIMARY:   Cliff's delta on leaf within-diversity DECAY-SLOPE
  SECONDARY: Cliff's delta on hub within-diversity at gen 500
  FLOOR:     Cliff's delta on leaf within-diversity at gen 500

Sign convention: delta = (#(a>b) - #(a<b)) / (n_a * n_b)
  positive delta = hub_out value LARGER than hub_in value.
  For decay_slope: more negative slope = faster decay.
  Negative delta on decay_slope => hub_out decays faster than hub_in.
"""
import csv
import os
import numpy as np
from scipy.stats import mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results", "diverse_seed")
NK_GRADING_DIR = os.path.join(HERE, "results", "nk_grading")
NK_VALUES = [0, 1, 2, 4]
TOPOLOGIES = ["star_hub_out", "star_hub_in"]


# ── Exact copies of stats_diverse_seed.py definitions ──────────────────────

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


# ── NK-grading analysis ─────────────────────────────────────────────────────

def analyze_nk(k):
    """Compute all three metrics for a single K value."""
    domain = f"nk{k}"
    out_seeds = load_seed_data(domain, "star_hub_out")
    in_seeds  = load_seed_data(domain, "star_hub_in")

    out_m = per_seed_metrics(out_seeds)
    in_m  = per_seed_metrics(in_seeds)

    # Extract scalar lists
    out_decay  = np.array([m["decay_slope"] for m in out_m])
    in_decay   = np.array([m["decay_slope"]  for m in in_m])
    out_hub500 = np.array([m["hub_floor"]    for m in out_m])
    in_hub500  = np.array([m["hub_floor"]    for m in in_m])
    out_floor  = np.array([m["floor"]        for m in out_m])
    in_floor   = np.array([m["floor"]        for m in in_m])

    # Peak gen context
    out_peak_gens = np.array([m["peak_gen"] for m in out_m])
    in_peak_gens  = np.array([m["peak_gen"] for m in in_m])

    # PRIMARY: decay slope  (a=hub_out, b=hub_in)
    d_decay = cliffs_delta(out_decay, in_decay)
    _, p_decay = mann_whitney(out_decay, in_decay)

    # SECONDARY: hub_floor at gen 500
    d_hub = cliffs_delta(out_hub500, in_hub500)
    _, p_hub = mann_whitney(out_hub500, in_hub500)

    # FLOOR: leaf floor at gen 500
    d_floor = cliffs_delta(out_floor, in_floor)
    _, p_floor = mann_whitney(out_floor, in_floor)

    return {
        "k": k,
        "n_out": len(out_m),
        "n_in":  len(in_m),
        # PRIMARY
        "d_decay":      d_decay,
        "p_decay":      p_decay,
        "mean_out_decay": float(np.mean(out_decay)),
        "mean_in_decay":  float(np.mean(in_decay)),
        # SECONDARY
        "d_hub":        d_hub,
        "p_hub":        p_hub,
        "mean_out_hub": float(np.mean(out_hub500)),
        "mean_in_hub":  float(np.mean(in_hub500)),
        # FLOOR
        "d_floor":      d_floor,
        "p_floor":      p_floor,
        "mean_out_floor": float(np.mean(out_floor)),
        "mean_in_floor":  float(np.mean(in_floor)),
        # Context
        "mean_peak_gen_out": float(np.mean(out_peak_gens)),
        "mean_peak_gen_in":  float(np.mean(in_peak_gens)),
    }


def fmt_p(p):
    """Format p to 2 sig figs."""
    if p >= 0.01:
        return f"{p:.2f}"
    elif p >= 0.001:
        return f"{p:.4f}"
    else:
        return f"{p:.2e}"


def build_analysis_md(results):
    lines = []
    lines.append("# NK-Grading Hub_Out vs Hub_In: Blind Analysis")
    lines.append("")
    lines.append("Generated by `experiments/nk_grading_analysis.py`.")
    lines.append("")
    lines.append("## Method")
    lines.append("")
    lines.append("Per-seed metrics computed using EXACT definitions from `stats_diverse_seed.py`:")
    lines.append("- **floor**: `leafWithinDiv` at gen 500.")
    lines.append("- **peak_gen**: generation of maximum `leafWithinDiv` per seed.")
    lines.append("- **decay_slope**: OLS slope of `leafWithinDiv` from peak_gen to gen 500.")
    lines.append("- **hub_floor**: `hubWithinDiv` at gen 500.")
    lines.append("")
    lines.append("Mann-Whitney two-sided, n=30 per condition.")
    lines.append("")
    lines.append("**Sign convention**: δ = (#(a>b) − #(a<b)) / (n_a × n_b), a=hub_out, b=hub_in.")
    lines.append("Positive δ = hub_out value LARGER than hub_in. For decay_slope (negative = decay),")
    lines.append("a NEGATIVE δ means hub_out decays FASTER than hub_in.")
    lines.append("")

    # ── Main table ──
    lines.append("## Main Results Table")
    lines.append("")
    lines.append("| K | δ_leaf_decay | p_leaf_decay | δ_hub_gen500 | p_hub_gen500 | δ_leaf_gen500 | p_leaf_gen500 |")
    lines.append("|---|---|---|---|---|---|---|")
    for r in results:
        lines.append(
            f"| {r['k']} "
            f"| {r['d_decay']:+.2f} | {fmt_p(r['p_decay'])} "
            f"| {r['d_hub']:+.2f} | {fmt_p(r['p_hub'])} "
            f"| {r['d_floor']:+.2f} | {fmt_p(r['p_floor'])} |"
        )
    lines.append("")

    # ── Per-condition means ──
    lines.append("## Per-Condition Means (for sign interpretation)")
    lines.append("")
    lines.append("### PRIMARY: Leaf Decay Slope (hub_out mean / hub_in mean)")
    lines.append("")
    lines.append("| K | hub_out mean slope | hub_in mean slope | mean peak_gen (hub_out) | mean peak_gen (hub_in) |")
    lines.append("|---|---|---|---|---|")
    for r in results:
        lines.append(
            f"| {r['k']} "
            f"| {r['mean_out_decay']:+.6f} "
            f"| {r['mean_in_decay']:+.6f} "
            f"| {r['mean_peak_gen_out']:.1f} "
            f"| {r['mean_peak_gen_in']:.1f} |"
        )
    lines.append("")

    lines.append("### SECONDARY: Hub Within-Div at Gen 500 (hub_out mean / hub_in mean)")
    lines.append("")
    lines.append("| K | hub_out mean | hub_in mean |")
    lines.append("|---|---|---|")
    for r in results:
        lines.append(
            f"| {r['k']} "
            f"| {r['mean_out_hub']:.6f} "
            f"| {r['mean_in_hub']:.6f} |"
        )
    lines.append("")

    lines.append("### FLOOR: Leaf Within-Div at Gen 500 (hub_out mean / hub_in mean)")
    lines.append("")
    lines.append("| K | hub_out mean | hub_in mean |")
    lines.append("|---|---|---|")
    for r in results:
        lines.append(
            f"| {r['k']} "
            f"| {r['mean_out_floor']:.6f} "
            f"| {r['mean_in_floor']:.6f} |"
        )
    lines.append("")

    # ── Narrative ──
    lines.append("## Narrative: sign(δ) behavior across K = 0 → 1 → 2 → 4")
    lines.append("")

    # Build narrative from actual values
    def sign_str(d, p, threshold=0.05):
        if p >= threshold:
            return "null (p≥0.05)"
        return "positive" if d > 0 else "negative"

    lines.append("### PRIMARY — Leaf decay slope δ (hub_out vs hub_in)")
    lines.append("")
    sign_seq = [(r['k'], r['d_decay'], r['p_decay']) for r in results]
    decay_lines = []
    for k, d, p in sign_seq:
        sig = "p<0.05" if p < 0.05 else "p≥0.05 (not significant)"
        decay_lines.append(
            f"- **K={k}**: δ={d:+.2f}, p={fmt_p(p)} ({sig}). "
            f"hub_out mean slope {[r['mean_out_decay'] for r in results if r['k']==k][0]:+.6f} "
            f"vs hub_in {[r['mean_in_decay'] for r in results if r['k']==k][0]:+.6f}."
        )
    lines.extend(decay_lines)

    # Detect sign flips
    signs = [(k, d, p) for k, d, p in sign_seq]
    flips = []
    for i in range(1, len(signs)):
        k_prev, d_prev, p_prev = signs[i-1]
        k_curr, d_curr, p_curr = signs[i]
        if (d_prev > 0) != (d_curr > 0):
            flips.append(f"K={k_prev}→K={k_curr}")
    if flips:
        lines.append(f"\nSign flip(s) in decay slope δ: {', '.join(flips)}.")
    else:
        lines.append("\nNo sign flip in decay slope δ across K = 0, 1, 2, 4.")

    lines.append("")

    lines.append("### SECONDARY — Hub within-div at gen 500 δ")
    lines.append("")
    hub_seq = [(r['k'], r['d_hub'], r['p_hub']) for r in results]
    for k, d, p in hub_seq:
        sig = "p<0.05" if p < 0.05 else "p≥0.05 (not significant)"
        lines.append(
            f"- **K={k}**: δ={d:+.2f}, p={fmt_p(p)} ({sig}). "
            f"hub_out mean {[r['mean_out_hub'] for r in results if r['k']==k][0]:.6f} "
            f"vs hub_in {[r['mean_in_hub'] for r in results if r['k']==k][0]:.6f}."
        )

    hub_signs = [(k, d, p) for k, d, p in hub_seq]
    hub_flips = []
    for i in range(1, len(hub_signs)):
        k_prev, d_prev, _ = hub_signs[i-1]
        k_curr, d_curr, _ = hub_signs[i]
        if (d_prev > 0) != (d_curr > 0):
            hub_flips.append(f"K={k_prev}→K={k_curr}")
    if hub_flips:
        lines.append(f"\nSign flip(s) in hub div δ: {', '.join(hub_flips)}.")
    else:
        lines.append("\nNo sign flip in hub div δ across K = 0, 1, 2, 4.")

    lines.append("")

    lines.append("### FLOOR — Leaf within-div at gen 500 δ")
    lines.append("")
    floor_seq = [(r['k'], r['d_floor'], r['p_floor']) for r in results]
    for k, d, p in floor_seq:
        sig = "p<0.05" if p < 0.05 else "p≥0.05 (not significant)"
        lines.append(
            f"- **K={k}**: δ={d:+.2f}, p={fmt_p(p)} ({sig}). "
            f"hub_out mean {[r['mean_out_floor'] for r in results if r['k']==k][0]:.6f} "
            f"vs hub_in {[r['mean_in_floor'] for r in results if r['k']==k][0]:.6f}."
        )

    floor_signs = [(k, d, p) for k, d, p in floor_seq]
    floor_flips = []
    for i in range(1, len(floor_signs)):
        k_prev, d_prev, _ = floor_signs[i-1]
        k_curr, d_curr, _ = floor_signs[i]
        if (d_prev > 0) != (d_curr > 0):
            floor_flips.append(f"K={k_prev}→K={k_curr}")
    if floor_flips:
        lines.append(f"\nSign flip(s) in leaf floor δ: {', '.join(floor_flips)}.")
    else:
        lines.append("\nNo sign flip in leaf floor δ across K = 0, 1, 2, 4.")

    lines.append("")
    lines.append("---")
    lines.append("*Blind analysis — no pre-registered prediction consulted during computation.*")
    lines.append("")

    return "\n".join(lines)


def main():
    os.makedirs(NK_GRADING_DIR, exist_ok=True)

    results = []
    for k in NK_VALUES:
        r = analyze_nk(k)
        results.append(r)
        print(f"K={k}: decay δ={r['d_decay']:+.4f} p={r['p_decay']:.4e} | "
              f"hub δ={r['d_hub']:+.4f} p={r['p_hub']:.4e} | "
              f"floor δ={r['d_floor']:+.4f} p={r['p_floor']:.4e}")

    md = build_analysis_md(results)
    out_path = os.path.join(NK_GRADING_DIR, "ANALYSIS.md")
    with open(out_path, "w") as f:
        f.write(md)
    print(f"\nWrote {out_path}")

    # Print the 4-row summary table to stdout
    print("\n=== 4-ROW SUMMARY TABLE ===")
    print(f"{'K':<4} {'d_decay':>10} {'p_decay':>12} {'d_hub':>8} {'p_hub':>12} {'d_floor':>8} {'p_floor':>12}")
    for r in results:
        print(f"{r['k']:<4} {r['d_decay']:+10.2f} {r['p_decay']:12.4e} "
              f"{r['d_hub']:+8.2f} {r['p_hub']:12.4e} "
              f"{r['d_floor']:+8.2f} {r['p_floor']:12.4e}")


if __name__ == "__main__":
    main()
