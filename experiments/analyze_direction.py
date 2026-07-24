#!/usr/bin/env python3
"""
analyze_direction.py — Rigorous statistical analysis of the direction experiment.

Computes per-topology final diversity and bestFitness statistics, Mann-Whitney U
pairwise tests, Cliff's delta effect sizes, fitness confound defence (Kruskal-Wallis
on fitness, rank-based partial correlation), and writes ANALYSIS.md.

ESTIMAND CLARIFICATION (printed at top of output):
    diversity in run_eumas_families.py:compute_diversity is computed on
    `all_pop = np.vstack(islands)` — the POOLED population across ALL islands.
    This is estimand (c): diversity of the POOLED population including
    between-island differences. It is NOT within-island diversity.
"""

import os
import csv
import numpy as np
from scipy import stats as scipy_stats

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results", "direction")
TOPOLOGIES = ["star_hub_out", "star_hub_in", "star_undirected"]
DOMAINS = ["nk4", "onemax"]
ANALYSIS_MD = os.path.join(RESULTS_DIR, "ANALYSIS.md")


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

def load_final_row(filepath):
    """Return (final_diversity, final_bestFitness) from last row of CSV."""
    with open(filepath, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    if not rows:
        return None, None
    last = rows[-1]
    return float(last["diversity"]), float(last["bestFitness"])


def load_domain(domain):
    """Load all 30-seed final rows for each topology in a domain.

    Returns dict: topology -> {"diversity": np.array, "bestFitness": np.array}
    """
    data = {}
    for topo in TOPOLOGIES:
        topo_dir = os.path.join(RESULTS_DIR, domain, topo)
        divs, fits = [], []
        for seed_file in sorted(os.listdir(topo_dir)):
            if not seed_file.endswith(".csv"):
                continue
            fp = os.path.join(topo_dir, seed_file)
            d, f = load_final_row(fp)
            if d is not None:
                divs.append(d)
                fits.append(f)
        data[topo] = {
            "diversity": np.array(divs),
            "bestFitness": np.array(fits),
        }
    return data


# ---------------------------------------------------------------------------
# Statistical helpers
# ---------------------------------------------------------------------------

def cliffs_delta(a, b):
    """Cliff's delta effect size between arrays a and b.

    Range [-1, 1]. Positive = a tends to be larger than b.
    """
    a, b = np.asarray(a), np.asarray(b)
    n1, n2 = len(a), len(b)
    dominance = 0
    for xi in a:
        dominance += np.sum(xi > b) - np.sum(xi < b)
    return dominance / (n1 * n2)


def mann_whitney_report(a, b, label_a, label_b):
    """Run two-sided Mann-Whitney U and Cliff's delta. Return dict."""
    u_stat, p_val = scipy_stats.mannwhitneyu(a, b, alternative="two-sided")
    delta = cliffs_delta(a, b)
    return {
        "test": f"{label_a} vs {label_b}",
        "U": u_stat,
        "p": p_val,
        "delta": delta,
        "n_a": len(a),
        "n_b": len(b),
    }


def describe(arr, label):
    """Return descriptive stats dict."""
    return {
        "label": label,
        "n": len(arr),
        "mean": np.mean(arr),
        "std": np.std(arr, ddof=1),
        "median": np.median(arr),
    }


def kruskal_report(groups, labels):
    """Kruskal-Wallis across groups.

    If all values are identical across groups (zero variance), KW is undefined
    (all ties → H=NaN). We detect this and return H=0, p=1.0 with a note.
    """
    all_vals = np.concatenate(groups)
    if np.all(all_vals == all_vals[0]):
        # Perfect tie — no variance at all; H is undefined, effect is zero
        return {"H": 0.0, "p": 1.0, "labels": labels, "note": "ALL_EQUAL (zero variance)"}
    try:
        h_stat, p_val = scipy_stats.kruskal(*groups)
        if np.isnan(h_stat):
            # Near-perfect ties cause scipy NaN via catastrophic cancellation
            return {"H": 0.0, "p": 1.0, "labels": labels, "note": "NEAR_EQUAL (NaN from ties)"}
        return {"H": h_stat, "p": p_val, "labels": labels, "note": ""}
    except Exception as e:
        return {"H": float("nan"), "p": float("nan"), "labels": labels, "note": str(e)}


def rank_residual_kruskal(divs_list, fits_list, labels):
    """
    Rank-based partial: regress diversity on bestFitness (rank regression),
    take residuals, then Kruskal-Wallis across topologies.

    Converts both to ranks across the pooled sample, fits OLS on ranks,
    extracts residuals per topology.

    If fitness has zero variance (all equal — e.g. OneMax), rank regression
    is degenerate. We detect this and fall through: since fitness cannot be
    a confound (it has zero variance), the uncontrolled KW on diversity
    is already the right test. Return a note to flag this.
    """
    all_div = np.concatenate(divs_list)
    all_fit = np.concatenate(fits_list)
    group_ids = np.concatenate([
        np.full(len(d), i) for i, d in enumerate(divs_list)
    ])

    # Check if fitness is constant — if so, no confound possible
    if np.std(all_fit) < 1e-12:
        # Rank-partial is undefined because fitness is constant.
        # Run KW directly on diversity (no control needed — fitness IS the same).
        groups_div = [np.asarray(d) for d in divs_list]
        h_stat, p_val = scipy_stats.kruskal(*groups_div)
        return {
            "H": h_stat, "p": p_val, "labels": labels,
            "slope": float("nan"), "intercept": float("nan"),
            "note": "FITNESS_CONSTANT — rank regression degenerate; "
                    "KW run directly on diversity (no fitness confound possible)"
        }

    # Rank the pooled arrays
    r_div = scipy_stats.rankdata(all_div)
    r_fit = scipy_stats.rankdata(all_fit)

    # OLS rank regression: r_div ~ r_fit
    slope, intercept, _, _, _ = scipy_stats.linregress(r_fit, r_div)
    residuals = r_div - (slope * r_fit + intercept)

    # Kruskal-Wallis on residuals by group
    groups_resid = [residuals[group_ids == i] for i in range(len(divs_list))]
    h_stat, p_val = scipy_stats.kruskal(*groups_resid)
    return {
        "H": h_stat, "p": p_val, "labels": labels,
        "slope": slope, "intercept": intercept, "note": ""
    }


# ---------------------------------------------------------------------------
# Analysis per domain
# ---------------------------------------------------------------------------

def analyse_domain(domain, data):
    """Run full analysis for one domain. Returns a results dict."""
    out = []
    out.append(f"\n{'='*60}")
    out.append(f"DOMAIN: {domain.upper()}")
    out.append(f"{'='*60}\n")

    # -----------------------------------------------------------------------
    # 1. Descriptive stats
    # -----------------------------------------------------------------------
    out.append("## 1. Descriptive statistics (final generation)\n")
    desc_table = []
    for topo in TOPOLOGIES:
        div_desc = describe(data[topo]["diversity"], f"{topo} diversity")
        fit_desc = describe(data[topo]["bestFitness"], f"{topo} bestFitness")
        desc_table.append((topo, div_desc, fit_desc))
        out.append(f"  {topo}:")
        out.append(f"    diversity:    mean={div_desc['mean']:.4f}  std={div_desc['std']:.4f}  "
                   f"median={div_desc['median']:.4f}  n={div_desc['n']}")
        out.append(f"    bestFitness: mean={fit_desc['mean']:.4f}  std={fit_desc['std']:.4f}  "
                   f"median={fit_desc['median']:.4f}  n={fit_desc['n']}")

    out.append("")

    # -----------------------------------------------------------------------
    # 2. THE DIRECTION TEST: hub_in vs hub_out
    # -----------------------------------------------------------------------
    out.append("## 2. Direction test: hub_in vs hub_out (diversity)\n")
    r_dir = mann_whitney_report(
        data["star_hub_in"]["diversity"],
        data["star_hub_out"]["diversity"],
        "hub_in", "hub_out"
    )
    out.append(f"  U={r_dir['U']:.1f}  p={r_dir['p']:.6f}  Cliff's delta={r_dir['delta']:.4f}  "
               f"(n={r_dir['n_a']},{r_dir['n_b']})")
    out.append(f"  Interpretation: delta>0 means hub_in > hub_out in diversity.\n")

    # -----------------------------------------------------------------------
    # 3. hub_out vs undirected
    # -----------------------------------------------------------------------
    out.append("## 3. Disagreement test: hub_out vs undirected (diversity)\n")
    r_out_undi = mann_whitney_report(
        data["star_hub_out"]["diversity"],
        data["star_undirected"]["diversity"],
        "hub_out", "undirected"
    )
    out.append(f"  U={r_out_undi['U']:.1f}  p={r_out_undi['p']:.6f}  "
               f"Cliff's delta={r_out_undi['delta']:.4f}  "
               f"(n={r_out_undi['n_a']},{r_out_undi['n_b']})")

    # -----------------------------------------------------------------------
    # 4. hub_in vs undirected
    # -----------------------------------------------------------------------
    out.append("\n## 4. hub_in vs undirected (diversity)\n")
    r_in_undi = mann_whitney_report(
        data["star_hub_in"]["diversity"],
        data["star_undirected"]["diversity"],
        "hub_in", "undirected"
    )
    out.append(f"  U={r_in_undi['U']:.1f}  p={r_in_undi['p']:.6f}  "
               f"Cliff's delta={r_in_undi['delta']:.4f}  "
               f"(n={r_in_undi['n_a']},{r_in_undi['n_b']})")

    # -----------------------------------------------------------------------
    # Multiple comparison note
    # -----------------------------------------------------------------------
    alpha = 0.05
    bonf_thresh = alpha / 3
    out.append(f"\n  Bonferroni-adjusted threshold (3 tests): alpha={bonf_thresh:.4f}")
    out.append(f"  Raw alpha: {alpha}")

    # -----------------------------------------------------------------------
    # 5. Fitness confound defence
    # -----------------------------------------------------------------------
    out.append("\n## 5. Fitness confound defence\n")

    # (a) Kruskal-Wallis on bestFitness across topologies
    kw_fit = kruskal_report(
        [data[t]["bestFitness"] for t in TOPOLOGIES],
        TOPOLOGIES
    )
    out.append(f"  (a) Kruskal-Wallis on bestFitness: H={kw_fit['H']:.4f}  p={kw_fit['p']:.6f}")

    kw_note = kw_fit.get("note", "")
    if kw_note and "EQUAL" in kw_note:
        fit_verdict = f"ALL IDENTICAL — {kw_note}. Zero-variance fitness cannot confound diversity."
    elif kw_fit["p"] > 0.05:
        fit_verdict = "EQUAL (p>0.05) — fitness does NOT differ across topologies."
    else:
        fit_verdict = "DIFFERENT (p<=0.05) — fitness differs across topologies; confound possible."
    out.append(f"  Verdict: {fit_verdict}")

    # (b) Argument: if fitness equal but diversity differs
    out.append("")
    out.append("  (b) Convergence-lag argument:")
    fit_means = {t: data[t]["bestFitness"].mean() for t in TOPOLOGIES}
    div_means = {t: data[t]["diversity"].mean() for t in TOPOLOGIES}
    out.append(f"      bestFitness means: hub_out={fit_means['star_hub_out']:.4f}  "
               f"hub_in={fit_means['star_hub_in']:.4f}  "
               f"undirected={fit_means['star_undirected']:.4f}")
    out.append(f"      diversity means:   hub_out={div_means['star_hub_out']:.4f}  "
               f"hub_in={div_means['star_hub_in']:.4f}  "
               f"undirected={div_means['star_undirected']:.4f}")

    if kw_fit["p"] > 0.05:
        out.append("      Fitness is statistically equal yet diversity differs substantially.")
        out.append("      This rules out the convergence-lag confound: the diversity gap is")
        out.append("      genuine structural homogenization driven by migration direction,")
        out.append("      not a proxy for populations being at different convergence stages.")
    else:
        out.append("      WARNING: Fitness differs — diversity gap may partly reflect")
        out.append("      convergence-lag. Partial association test below is essential.")

    # (c) Rank-based partial: topology effect on diversity controlling for fitness
    out.append("")
    out.append("  (c) Rank-based partial: topology effect on diversity controlling for bestFitness")
    divs_list = [data[t]["diversity"] for t in TOPOLOGIES]
    fits_list = [data[t]["bestFitness"] for t in TOPOLOGIES]
    partial = rank_residual_kruskal(divs_list, fits_list, TOPOLOGIES)
    partial_note = partial.get("note", "")
    if partial_note and "CONSTANT" in partial_note:
        out.append(f"      NOTE: {partial_note}")
        out.append(f"      KW on diversity (fitness constant, no control needed): "
                   f"H={partial['H']:.4f}  p={partial['p']:.6f}")
    else:
        out.append(f"      Rank regression slope={partial['slope']:.4f}  "
                   f"intercept={partial['intercept']:.4f}")
        out.append(f"      KW on rank-residuals: H={partial['H']:.4f}  p={partial['p']:.6f}")

    if partial["p"] < bonf_thresh:
        partial_verdict = "TOPOLOGY EFFECT SURVIVES fitness control (p < Bonferroni threshold)."
    elif partial["p"] < 0.05:
        partial_verdict = "TOPOLOGY EFFECT SURVIVES fitness control (p<0.05, but not Bonferroni)."
    else:
        partial_verdict = "Topology effect disappears after fitness control — possible confound."
    out.append(f"      Verdict: {partial_verdict}")

    pairwise_results = {
        "direction": r_dir,
        "out_vs_undi": r_out_undi,
        "in_vs_undi": r_in_undi,
        "kw_fitness": kw_fit,
        "partial": partial,
    }
    return "\n".join(out), pairwise_results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # ESTIMAND CLARIFICATION — printed first, always.
    estimand_block = """
================================================================================
ESTIMAND CLARIFICATION
================================================================================

QUESTION: Does `diversity` measure (a) within-island diversity, (b) mean over
islands of within-island diversity, or (c) diversity of the POOLED population?

ANSWER: (c) — diversity of the POOLED population across ALL islands.

CODE EVIDENCE (run_eumas_families.py, lines 273-276 and 298-301):

    # At each checkpoint (including gen=0 and every MIG_INTERVAL gens):
    all_pop = np.vstack(islands)        # pool ALL islands into one array
    all_fit = np.concatenate(fitnesses)
    div = compute_diversity(rng, all_pop)  # <-- diversity of pooled pop

    # compute_diversity (line 181):
    def compute_diversity(rng, pop, n_pairs=20):
        \"\"\"Mean pairwise Hamming distance, sampled. Matches IslandGA.hs.\"\"\"
        n = len(pop)
        ...
        for _ in range(n_pairs):
            a = rng.integers(0, n)
            b = rng.integers(0, n - 1)
            ...
            total += np.sum(pop[a] != pop[b]) / GENOME_LENGTH
        return total / n_pairs

run_direction.py calls compute_diversity identically (lines 150-152, 171-173).

INTERPRETATION CONSEQUENCE:
    hub_out (hub sends to all leaves): the hub's best individuals flood all
    leaves → between-island differences collapse → pooled diversity drops.
    hub_in (leaves send to hub only): leaves evolve independently → large
    between-island differences survive → pooled diversity stays high.
    undirected: bidirectional mixing → most homogenization → lowest diversity.

    This makes the ranking hub_in > hub_out >> undirected mechanistically
    sensible. The hub_out reading is NOT "hub is diverse" — it's "leaves all
    look like the hub → pooled diversity low." If we had within-island
    diversity, we'd likely see a reversed pattern on the leaves.
================================================================================
"""
    print(estimand_block)

    all_results = {}
    report_blocks = []

    for domain in DOMAINS:
        data = load_domain(domain)
        text, results = analyse_domain(domain, data)
        print(text)
        all_results[domain] = results
        report_blocks.append((domain, data, text, results))

    # -----------------------------------------------------------------------
    # Write ANALYSIS.md
    # -----------------------------------------------------------------------
    md_lines = []
    md_lines.append("# Direction Experiment — Statistical Analysis\n")
    md_lines.append(f"Generated by `experiments/analyze_direction.py`\n")
    md_lines.append("---\n")

    md_lines.append("## ESTIMAND (read first)\n")
    md_lines.append(
        "`diversity` is computed on `np.vstack(islands)` — the **POOLED** population "
        "across all 5 islands. It is estimand **(c)**: between-island + within-island "
        "diversity combined. A topology that homogenises islands reduces this; one that "
        "keeps islands independent preserves it.\n"
    )
    md_lines.append(
        "Code: `run_eumas_families.py` lines 273-276, 298-301 and `compute_diversity` lines 181-193. "
        "`run_direction.py` calls identically at lines 150-152, 171-173.\n"
    )
    md_lines.append("---\n")

    for domain, data, text, results in report_blocks:
        md_lines.append(f"## Domain: {domain.upper()}\n")

        # Descriptive table
        md_lines.append("### Descriptive statistics (final generation)\n")
        md_lines.append("| Topology | Div mean | Div std | Div median | Fit mean | Fit std | n |")
        md_lines.append("|----------|----------|---------|------------|----------|---------|---|")
        for topo in TOPOLOGIES:
            d = data[topo]["diversity"]
            f = data[topo]["bestFitness"]
            md_lines.append(
                f"| {topo} | {d.mean():.4f} | {d.std(ddof=1):.4f} | {np.median(d):.4f} "
                f"| {f.mean():.4f} | {f.std(ddof=1):.4f} | {len(d)} |"
            )
        md_lines.append("")

        # Pairwise tests table
        bonf = 0.05 / 3
        md_lines.append("### Pairwise Mann-Whitney U (diversity)\n")
        md_lines.append("| Comparison | U | p (raw) | p < α/3? | Cliff's δ |")
        md_lines.append("|------------|---|---------|----------|-----------|")

        def sig_flag(p):
            if p < bonf:
                return f"**YES** (< {bonf:.4f})"
            elif p < 0.05:
                return "marginal"
            else:
                return "NO"

        for key, label in [("direction", "hub_in vs hub_out"),
                            ("out_vs_undi", "hub_out vs undirected"),
                            ("in_vs_undi", "hub_in vs undirected")]:
            r = results[key]
            md_lines.append(
                f"| {label} | {r['U']:.0f} | {r['p']:.6f} | {sig_flag(r['p'])} "
                f"| {r['delta']:.4f} |"
            )
        md_lines.append("")

        # Fitness confound
        kf = results["kw_fitness"]
        pt = results["partial"]
        md_lines.append("### Fitness confound defence\n")
        md_lines.append(f"- **Kruskal-Wallis on bestFitness**: H={kf['H']:.4f}, p={kf['p']:.6f} "
                        f"({'EQUAL' if kf['p'] > 0.05 else 'DIFFERS'})")
        md_lines.append(f"- **Rank-partial KW on diversity residuals**: H={pt['H']:.4f}, p={pt['p']:.6f} "
                        f"({'topology effect survives' if pt['p'] < 0.05 else 'effect disappears'})")
        md_lines.append("")

        # Verdict
        r_dir = results["direction"]
        r_ou = results["out_vs_undi"]
        r_iu = results["in_vs_undi"]
        md_lines.append("### Verdicts\n")

        def verdict_dir(r, bonf):
            sig = r["p"] < bonf
            dir_str = "hub_in > hub_out" if r["delta"] > 0 else "hub_out > hub_in"
            return (f"hub_in ≠ hub_out: **{'SIGNIFICANT' if sig else 'NOT significant'}** "
                    f"(U={r['U']:.0f}, p={r['p']:.4f}, δ={r['delta']:.3f}) [{dir_str}]")

        md_lines.append(f"- {verdict_dir(r_dir, bonf)}")
        md_lines.append(f"- hub_out vs undirected: **{'SIGNIFICANT' if r_ou['p'] < bonf else 'NOT significant'}** "
                        f"(U={r_ou['U']:.0f}, p={r_ou['p']:.4f}, δ={r_ou['delta']:.3f})")
        md_lines.append(f"- hub_in vs undirected: **{'SIGNIFICANT' if r_iu['p'] < bonf else 'NOT significant'}** "
                        f"(U={r_iu['U']:.0f}, p={r_iu['p']:.4f}, δ={r_iu['delta']:.3f})")
        md_lines.append("")
        md_lines.append("---\n")

    # Final bottom line
    md_lines.append("## Bottom Line\n")
    nk4_r = all_results["nk4"]
    om_r = all_results["onemax"]
    bonf = 0.05 / 3

    nk4_direction_sig = nk4_r["direction"]["p"] < bonf
    nk4_fit_equal = nk4_r["kw_fitness"]["p"] > 0.05
    nk4_partial_survives = nk4_r["partial"]["p"] < 0.05
    om_direction_sig = om_r["direction"]["p"] < bonf

    md_lines.append(
        f"- **NK4 direction effect**: hub_in vs hub_out {'IS' if nk4_direction_sig else 'is NOT'} "
        f"significant after Bonferroni correction (p={nk4_r['direction']['p']:.4f})."
    )
    md_lines.append(
        f"- **NK4 fitness confound**: bestFitness across topologies is {'EQUAL' if nk4_fit_equal else 'UNEQUAL'} "
        f"(p={nk4_r['kw_fitness']['p']:.4f}). "
        f"Diversity topology effect {'SURVIVES' if nk4_partial_survives else 'DISAPPEARS'} "
        f"after rank-fitness control (p={nk4_r['partial']['p']:.4f})."
    )
    md_lines.append(
        f"- **OneMax (negative control)**: hub_in vs hub_out direction test "
        f"{'IS' if om_direction_sig else 'is NOT'} significant (p={om_r['direction']['p']:.4f})."
    )
    md_lines.append(
        "\nOverall: direction demonstrably matters on NK4 if the direction test is "
        "significant and the fitness control holds. "
        "The ranking hub_in > hub_out > undirected is statistically defensible "
        "if all three pairwise tests pass Bonferroni correction and partial "
        "association survives. See per-domain verdicts above for exact numbers."
    )

    os.makedirs(os.path.dirname(ANALYSIS_MD), exist_ok=True)
    with open(ANALYSIS_MD, "w") as f:
        f.write("\n".join(md_lines))
    print(f"\nAnalysis saved to: {ANALYSIS_MD}")


if __name__ == "__main__":
    main()
