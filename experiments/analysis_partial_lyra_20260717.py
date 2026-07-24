"""
Partial correlation analysis: does fitness explain the lambda2->diversity relationship?
Generated: 2026-07-17

Computes for each domain (NK4, OneMax):
  1. Graph-level (n=10, mean over seeds):
     - Spearman(lambda2, final_diversity)
     - Spearman(lambda2, final_best_fitness), Spearman(lambda2, final_mean_fitness)
     - Spearman(final_mean_fitness, final_diversity)
     - Partial Spearman(lambda2, diversity | mean_fitness) via rank-then-residualize
     - Partial Spearman(lambda2, diversity | best_fitness) same method
     - Both also via closed-form r_xy.z formula as cross-check
  2. Run-level (n=200, all individual rows):
     - Same correlations
     - Partial rho + p-value (t-test, df=n-3=197)
  3. Plain-language NK4 verdict.
"""

import csv
import math
import sys
from collections import defaultdict

import numpy as np
from scipy import stats


# ---------------------------------------------------------------------------
# Load CSV
# ---------------------------------------------------------------------------
DATA_PATH = (
    "/home/lyra/projects/Topology-experiments/experiments/results/"
    "volume_vs_spectrum/pilot_with_fitness.csv"
)

rows = []
with open(DATA_PATH, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        rows.append({
            "graph_id": row["graph_id"],
            "lambda2": float(row["lambda2"]),
            "domain": row["domain"],
            "seed": int(row["seed"]),
            "final_diversity": float(row["final_diversity"]),
            "final_best_fitness": float(row["final_best_fitness"]),
            "final_mean_fitness": float(row["final_mean_fitness"]),
        })

print(f"Loaded {len(rows)} rows.")
domains = sorted(set(r["domain"] for r in rows))
print(f"Domains: {domains}")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def spearman(x, y):
    """Return (rho, p-value) via scipy."""
    return stats.spearmanr(x, y)


def partial_spearman_residual(x, y, z):
    """
    Partial Spearman correlation of x and y controlling for z.
    Method: rank-transform all three, then residualize x-ranks and y-ranks
    on z-ranks via OLS, correlate residuals with Pearson.
    Returns (rho, p-value) using t-test with df = n-3.
    """
    n = len(x)
    rx = stats.rankdata(x).astype(float)
    ry = stats.rankdata(y).astype(float)
    rz = stats.rankdata(z).astype(float)

    # Residualize rx on rz
    slope_xz, intercept_xz, _, _, _ = stats.linregress(rz, rx)
    resid_x = rx - (slope_xz * rz + intercept_xz)

    # Residualize ry on rz
    slope_yz, intercept_yz, _, _, _ = stats.linregress(rz, ry)
    resid_y = ry - (slope_yz * rz + intercept_yz)

    r, _ = stats.pearsonr(resid_x, resid_y)
    # t-test for partial correlation, df = n - 3
    df = n - 3
    if abs(r) >= 1.0:
        p = 0.0
    else:
        t = r * math.sqrt(df / (1.0 - r * r))
        p = 2.0 * stats.t.sf(abs(t), df)
    return r, p


def partial_spearman_formula(x, y, z):
    """
    Partial Spearman via closed-form: r_xy.z = (r_xy - r_xz*r_yz) / sqrt((1-r_xz^2)(1-r_yz^2))
    Returns (rho, p-value) with df = n-3.
    """
    n = len(x)
    r_xy, _ = stats.spearmanr(x, y)
    r_xz, _ = stats.spearmanr(x, z)
    r_yz, _ = stats.spearmanr(y, z)
    denom = math.sqrt((1 - r_xz**2) * (1 - r_yz**2))
    if denom == 0:
        return float("nan"), float("nan")
    r_partial = (r_xy - r_xz * r_yz) / denom
    df = n - 3
    if abs(r_partial) >= 1.0:
        p = 0.0
    else:
        t = r_partial * math.sqrt(df / (1.0 - r_partial**2))
        p = 2.0 * stats.t.sf(abs(t), df)
    return r_partial, p


def fmt(rho, p):
    return f"rho={rho:+.4f}, p={p:.4f}"


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

SEP = "=" * 70

for domain in domains:
    print(f"\n{SEP}")
    print(f"DOMAIN: {domain.upper()}")
    print(SEP)

    domain_rows = [r for r in rows if r["domain"] == domain]
    n_runs = len(domain_rows)
    print(f"  Run-level n = {n_runs}")

    # ------------------------------------------------------------------
    # 1. Graph-level (aggregate mean over seeds per graph)
    # ------------------------------------------------------------------
    graph_data = defaultdict(lambda: {
        "lambda2": None,
        "final_diversity": [],
        "final_best_fitness": [],
        "final_mean_fitness": [],
    })
    for r in domain_rows:
        g = r["graph_id"]
        graph_data[g]["lambda2"] = r["lambda2"]
        graph_data[g]["final_diversity"].append(r["final_diversity"])
        graph_data[g]["final_best_fitness"].append(r["final_best_fitness"])
        graph_data[g]["final_mean_fitness"].append(r["final_mean_fitness"])

    graphs = sorted(graph_data.keys())
    n_graphs = len(graphs)
    print(f"  Graph-level n = {n_graphs}")

    lam2 = np.array([graph_data[g]["lambda2"] for g in graphs])
    div_g = np.array([np.mean(graph_data[g]["final_diversity"]) for g in graphs])
    best_g = np.array([np.mean(graph_data[g]["final_best_fitness"]) for g in graphs])
    mean_g = np.array([np.mean(graph_data[g]["final_mean_fitness"]) for g in graphs])

    print(f"\n  --- Graph-level Spearman correlations (n={n_graphs}) ---")

    rho_ld, p_ld = spearman(lam2, div_g)
    print(f"  Spearman(lambda2, diversity):          {fmt(rho_ld, p_ld)}")

    rho_lb, p_lb = spearman(lam2, best_g)
    sign_lb = "HIGHER" if rho_lb > 0 else "LOWER"
    print(f"  Spearman(lambda2, best_fitness):       {fmt(rho_lb, p_lb)}  → higher lambda2 → {sign_lb} best_fitness")

    rho_lm, p_lm = spearman(lam2, mean_g)
    sign_lm = "HIGHER" if rho_lm > 0 else "LOWER"
    print(f"  Spearman(lambda2, mean_fitness):       {fmt(rho_lm, p_lm)}  → higher lambda2 → {sign_lm} mean_fitness")

    rho_md, p_md = spearman(mean_g, div_g)
    print(f"  Spearman(mean_fitness, diversity):     {fmt(rho_md, p_md)}")

    # Partial: lambda2 vs diversity | mean_fitness
    pr_resid_m, pp_resid_m = partial_spearman_residual(lam2, div_g, mean_g)
    pr_form_m, pp_form_m = partial_spearman_formula(lam2, div_g, mean_g)
    print(f"\n  Partial Spearman(lambda2, diversity | mean_fitness):")
    print(f"    Residual method:  {fmt(pr_resid_m, pp_resid_m)}")
    print(f"    Formula method:   {fmt(pr_form_m, pp_form_m)}")

    # Partial: lambda2 vs diversity | best_fitness
    pr_resid_b, pp_resid_b = partial_spearman_residual(lam2, div_g, best_g)
    pr_form_b, pp_form_b = partial_spearman_formula(lam2, div_g, best_g)
    print(f"\n  Partial Spearman(lambda2, diversity | best_fitness):")
    print(f"    Residual method:  {fmt(pr_resid_b, pp_resid_b)}")
    print(f"    Formula method:   {fmt(pr_form_b, pp_form_b)}")

    # ------------------------------------------------------------------
    # 2. Run-level (n=200 per domain)
    # ------------------------------------------------------------------
    lam2_r = np.array([r["lambda2"] for r in domain_rows])
    div_r = np.array([r["final_diversity"] for r in domain_rows])
    best_r = np.array([r["final_best_fitness"] for r in domain_rows])
    mean_r = np.array([r["final_mean_fitness"] for r in domain_rows])

    print(f"\n  --- Run-level Spearman correlations (n={n_runs}) ---")

    rho_ld_r, p_ld_r = spearman(lam2_r, div_r)
    print(f"  Spearman(lambda2, diversity):          {fmt(rho_ld_r, p_ld_r)}")

    rho_lm_r, p_lm_r = spearman(lam2_r, mean_r)
    sign_lm_r = "HIGHER" if rho_lm_r > 0 else "LOWER"
    print(f"  Spearman(lambda2, mean_fitness):       {fmt(rho_lm_r, p_lm_r)}  → higher lambda2 → {sign_lm_r} mean_fitness")

    rho_md_r, p_md_r = spearman(mean_r, div_r)
    print(f"  Spearman(mean_fitness, diversity):     {fmt(rho_md_r, p_md_r)}")

    pr_r, pp_r = partial_spearman_residual(lam2_r, div_r, mean_r)
    pr_f, pp_f = partial_spearman_formula(lam2_r, div_r, mean_r)
    print(f"\n  Partial Spearman(lambda2, diversity | mean_fitness) run-level:")
    print(f"    Residual method:  {fmt(pr_r, pp_r)}  (df={n_runs - 3})")
    print(f"    Formula method:   {fmt(pr_f, pp_f)}  (df={n_runs - 3})")

    # ------------------------------------------------------------------
    # 3. NK4 plain-language verdict
    # ------------------------------------------------------------------
    if domain == "nk4":
        print(f"\n  --- NK4 VERDICT ---")
        print(f"  Raw graph-level Spearman(lambda2, diversity):           rho={rho_ld:+.4f}")
        print(f"  Partial graph-level (controlling mean_fitness):         rho={pr_resid_m:+.4f}  [residual], {pr_form_m:+.4f}  [formula]")
        print(f"  Partial graph-level (controlling best_fitness):         rho={pr_resid_b:+.4f}  [residual], {pr_form_b:+.4f}  [formula]")
        print(f"  Raw run-level Spearman(lambda2, diversity):             rho={rho_ld_r:+.4f}")
        print(f"  Partial run-level (controlling mean_fitness):           rho={pr_r:+.4f}  [residual], {pr_f:+.4f}  [formula]")
        print(f"  Sign: higher lambda2 → {sign_lm} mean_fitness (graph-level)")

        # Assess survival vs collapse
        raw_abs = abs(rho_ld)
        partial_abs_m = abs(pr_resid_m)
        attenuation_pct = 100.0 * (raw_abs - partial_abs_m) / raw_abs if raw_abs > 0 else float("nan")
        print(f"  Attenuation from raw to partial (mean_fitness control): {attenuation_pct:.1f}%")

        raw_abs_r = abs(rho_ld_r)
        partial_abs_r = abs(pr_r)
        attenuation_r = 100.0 * (raw_abs_r - partial_abs_r) / raw_abs_r if raw_abs_r > 0 else float("nan")
        print(f"  Attenuation run-level:                                  {attenuation_r:.1f}%")

        if partial_abs_m >= 0.7 * raw_abs:
            verdict = "SURVIVES — partial rho stays within 30% of raw; fitness does not explain away the effect."
        elif partial_abs_m >= 0.4 * raw_abs:
            verdict = "PARTIAL MEDIATION — partial rho attenuated but still substantial; fitness partially explains the effect."
        else:
            verdict = "COLLAPSES — partial rho drops sharply; fitness largely explains the lambda2->diversity relationship."
        print(f"\n  ONE-LINE VERDICT: {verdict}")

print(f"\n{SEP}")
print("Done.")
