"""
Adversarial analysis of the claim:
  "In the NK4 domain, graph connectivity lambda2 is negatively associated
   with final within-island diversity." (Spearman ~ -0.84, n=10 graphs)

stdlib + numpy + scipy only.
"""

import csv
import sys
from itertools import combinations
import numpy as np
from scipy import stats

CSV_PATH = "experiments/results/volume_vs_spectrum/pilot.csv"

# ── Load data ──────────────────────────────────────────────────────────────

rows = []
with open(CSV_PATH, newline="") as f:
    reader = csv.DictReader(f)
    header = reader.fieldnames
    print(f"Columns: {header}")
    for row in reader:
        rows.append(row)

print(f"Total rows: {len(rows)}")

# Separate domains
nk4_rows   = [r for r in rows if r["domain"] == "nk4"]
om_rows    = [r for r in rows if r["domain"] == "onemax"]
print(f"NK4 rows: {len(nk4_rows)},  OneMax rows: {len(om_rows)}")

# ── Fitness column check ────────────────────────────────────────────────────
fitness_col = None
for candidate in ["final_fitness", "best_fitness", "mean_fitness", "fitness"]:
    if candidate in header:
        fitness_col = candidate
        break

print(f"\nFitness column detected: {fitness_col!r}")
if fitness_col is None:
    print("  WARNING: No fitness column in data. Task 2 (fitness mediation) "
          "CANNOT BE COMPUTED from this dataset.")

# ── Build graph-level aggregates ────────────────────────────────────────────

def graph_agg(row_subset, use_median=False):
    """Return dict: graph_id -> (lambda2, agg_diversity [, agg_fitness])"""
    from collections import defaultdict
    buckets = defaultdict(list)
    for r in row_subset:
        gid = r["graph_id"]
        buckets[gid].append(r)
    result = {}
    for gid, rs in buckets.items():
        lam  = float(rs[0]["lambda2"])
        divs = [float(r["final_diversity"]) for r in rs]
        agg_div = float(np.median(divs)) if use_median else float(np.mean(divs))
        entry = {"lambda2": lam, "diversity": agg_div, "div_vals": divs}
        if fitness_col:
            fits = [float(r[fitness_col]) for r in rs]
            entry["fitness"] = float(np.median(fits)) if use_median else float(np.mean(fits))
            entry["fit_vals"] = fits
        result[gid] = entry
    return result

nk4_agg  = graph_agg(nk4_rows, use_median=False)
nk4_agg_med = graph_agg(nk4_rows, use_median=True)

graphs = sorted(nk4_agg.keys())
lam2   = np.array([nk4_agg[g]["lambda2"]   for g in graphs])
div    = np.array([nk4_agg[g]["diversity"]  for g in graphs])

print(f"\nNumber of graphs (NK4): {len(graphs)}")
print(f"Graph IDs: {graphs}")
print(f"\nlambda2 values: {np.round(lam2, 4)}")
print(f"mean final_diversity: {np.round(div, 4)}")

# Baseline Spearman
rho_base, p_base = stats.spearmanr(lam2, div)
print(f"\nBaseline Spearman (lambda2 vs mean diversity, NK4): rho={rho_base:.4f}, p={p_base:.4f}")

# ══════════════════════════════════════════════════════════════════════════
# TASK 1 — Leave-one-out robustness
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("TASK 1: Leave-one-out Spearman")
print("="*60)

loo_rhos = []
for i, g in enumerate(graphs):
    mask = [j for j in range(len(graphs)) if j != i]
    r, _ = stats.spearmanr(lam2[mask], div[mask])
    loo_rhos.append(r)
    print(f"  drop {g:30s}  rho={r:.4f}")

loo_min = min(loo_rhos)
loo_max = max(loo_rhos)
print(f"\nLOO min rho: {loo_min:.4f}")
print(f"LOO max rho: {loo_max:.4f}")
sign_flips = sum(1 for r in loo_rhos if r > 0)
drops_below_07 = sum(1 for r in loo_rhos if abs(r) < 0.7)
print(f"Sign flips: {sign_flips}")
print(f"Drops |rho|<0.7: {drops_below_07}")

# ══════════════════════════════════════════════════════════════════════════
# TASK 2 — Fitness mediation
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("TASK 2: Fitness mediation")
print("="*60)

if fitness_col is None:
    print("  CANNOT COMPUTE: no fitness column in pilot.csv.")
    print("  The dataset records only final_diversity and lambda2.")
    print("  Mediation test requires a fitness column that is absent.")
    print("  VERDICT: UNTESTABLE from this data — this is itself a finding.")
    print("  The claim is neither confirmed nor broken on fitness grounds;")
    print("  the data simply does not contain fitness measurements.")
else:
    fit = np.array([nk4_agg[g]["fitness"] for g in graphs])
    rho_lam_fit, p_lam_fit   = stats.spearmanr(lam2, fit)
    rho_fit_div, p_fit_div   = stats.spearmanr(fit, div)
    print(f"  (a) Spearman lambda2 vs mean fitness:    rho={rho_lam_fit:.4f}, p={p_lam_fit:.4f}")
    print(f"  (b) Spearman mean fitness vs diversity:  rho={rho_fit_div:.4f}, p={p_fit_div:.4f}")

    # Partial correlation: lambda2 ~ diversity | fitness
    # Residualise both on fitness, then correlate residuals
    def residuals(x, z):
        slope, intercept, *_ = stats.linregress(z, x)
        return x - (slope * z + intercept)

    res_lam = residuals(lam2, fit)
    res_div = residuals(div, fit)
    rho_partial, p_partial = stats.pearsonr(res_lam, res_div)
    print(f"  (c) Partial corr lambda2~diversity | fitness (Pearson of residuals):")
    print(f"      r={rho_partial:.4f}, p={p_partial:.4f}")

    if abs(rho_partial) < 0.3 or p_partial > 0.05:
        print("  VERDICT: FITNESS CONFOUND — lambda2–diversity relation killed by controlling for fitness.")
    else:
        print("  VERDICT: SURVIVES partial control — lambda2 predicts diversity even after fitness.")

# ══════════════════════════════════════════════════════════════════════════
# TASK 3 — Aggregation sensitivity (median vs mean)
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("TASK 3: Aggregation sensitivity (median vs mean)")
print("="*60)

div_med = np.array([nk4_agg_med[g]["diversity"] for g in graphs])
rho_med, p_med = stats.spearmanr(lam2, div_med)
print(f"  Spearman with MEDIAN-over-seeds:  rho={rho_med:.4f}, p={p_med:.4f}")
print(f"  Spearman with MEAN-over-seeds:    rho={rho_base:.4f}, p={p_base:.4f}")
delta = abs(rho_base) - abs(rho_med)
print(f"  |rho| change (mean-median): {delta:+.4f}")
if abs(delta) > 0.1:
    print("  VERDICT: SENSITIVE — median aggregation materially changes magnitude.")
else:
    print("  VERDICT: ROBUST — mean and median give similar magnitudes.")

# ══════════════════════════════════════════════════════════════════════════
# TASK 4 — Is diversity just measuring convergence (not-yet-converged)?
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("TASK 4: Diversity as convergence proxy (variance check)")
print("="*60)

print("\n  Per-graph NK4 diversity statistics:")
print(f"  {'graph_id':30s}  lambda2  mean_div  std_div  cv")
div_means = []
div_stds  = []
for g in graphs:
    d = nk4_agg[g]["div_vals"]
    mn  = np.mean(d)
    sd  = np.std(d, ddof=1)
    cv  = sd / mn if mn > 0 else float("nan")
    div_means.append(mn)
    div_stds.append(sd)
    print(f"  {g:30s}  {nk4_agg[g]['lambda2']:.4f}   {mn:.4f}    {sd:.4f}   {cv:.3f}")

# Correlation: lambda2 vs within-graph std of diversity
div_stds_arr = np.array(div_stds)
rho_var, p_var = stats.spearmanr(lam2, div_stds_arr)
print(f"\n  Spearman(lambda2, within-graph std of diversity): rho={rho_var:.4f}, p={p_var:.4f}")

# Median absolute deviation per graph — is high-diversity driven by outlier seeds?
mads = []
for g in graphs:
    d = np.array(nk4_agg[g]["div_vals"])
    mads.append(np.median(np.abs(d - np.median(d))))
mads_arr = np.array(mads)
rho_mad, p_mad = stats.spearmanr(lam2, mads_arr)
print(f"  Spearman(lambda2, MAD of diversity):              rho={rho_mad:.4f}, p={p_mad:.4f}")

if fitness_col:
    # Are high-diversity graphs also lower fitness?
    fit_arr = np.array([nk4_agg[g]["fitness"] for g in graphs])
    rho_df, p_df = stats.spearmanr(div, fit_arr)
    print(f"\n  Spearman(mean diversity, mean fitness):  rho={rho_df:.4f}, p={p_df:.4f}")
    if rho_df < -0.5:
        print("  High-diversity graphs have systematically LOWER fitness → convergence story plausible.")
    else:
        print("  No strong fitness-diversity anti-correlation → convergence story weak.")
else:
    print("\n  Fitness column absent — cannot directly check whether high-diversity = low-fitness.")
    print("  Using variance of diversity as indirect convergence proxy:")
    # If low-lambda2 graphs are simply not yet converged, we'd expect BOTH higher mean diversity
    # AND higher variance (unsettled trajectories). Test that.
    rho_mean_std, _ = stats.spearmanr(div_means, div_stds)
    print(f"  Spearman(mean diversity, std diversity) across graphs: rho={rho_mean_std:.4f}")
    if rho_mean_std > 0.6:
        print("  High-diversity graphs also have high variance → consistent with 'not yet converged'.")
        print("  VERDICT: PARTIAL — cannot rule out convergence story without fitness data.")
    else:
        print("  High-diversity graphs do NOT have disproportionate variance → variance pattern")
        print("  does not support a simple 'not yet converged' explanation.")
        print("  VERDICT: SURVIVES — variance pattern does not confirm convergence artifact.")

# ══════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"Baseline Spearman (NK4, lambda2 vs mean diversity): rho={rho_base:.4f}, p={p_base:.4f}")
print(f"LOO range: [{loo_min:.4f}, {loo_max:.4f}]  sign flips={sign_flips}  |rho|<0.7: {drops_below_07}")
print(f"Median agg Spearman: rho={rho_med:.4f}  (delta |rho| from mean: {delta:+.4f})")
