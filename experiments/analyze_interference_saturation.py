#!/usr/bin/env python3
"""
Interference Saturation Analysis — Testing Clio's Predictions

Clio predicts that different interference configurations should have
different saturation CEILINGS for eta-squared, based on effective-n:
  - separated: effective-n = 2L (highest ceiling)
  - adjacent:  effective-n = 2L-1
  - nested:    effective-n = 2L-2 (lowest ceiling)

She also asks: does the Phase 1 power-law exponent (alpha) differ
across configurations?

Analysis:
  1. Compute eta-squared time series for each pairwise comparison
     AND for each config vs the other two combined
  2. Fit logistic saturation: eta2(t) = C / (1 + exp(-r*(t - t0)))
  3. Fit Phase 1 power law: log(eta2) = alpha * log(t) + c  (gen 10-100)
  4. Compare ceilings and alphas across datasets (OneMax, maze, sudoku)
"""

import os
import csv
import math
import json
from collections import defaultdict

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import pearsonr

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE_DIR = "/home/lyra/projects/Topology-experiments"
OUTPUT_DIR = os.path.join(BASE_DIR, "results", "interference-analysis")
os.makedirs(OUTPUT_DIR, exist_ok=True)

TOPOLOGIES = ["interference-adjacent", "interference-separated", "interference-nested"]
TOPO_SHORT = {
    "interference-adjacent": "adjacent",
    "interference-separated": "separated",
    "interference-nested": "nested",
}
TOPO_COLORS = {
    "interference-adjacent": "#e74c3c",
    "interference-separated": "#2ecc71",
    "interference-nested": "#3498db",
}

DATASETS = {
    "onemax": os.path.join(BASE_DIR, "results", "interference"),
    "maze": os.path.join(BASE_DIR, "results", "interference-maze"),
    "sudoku": os.path.join(BASE_DIR, "results", "interference-sudoku"),
}


def load_csv(path):
    with open(path) as f:
        lines = [l for l in f if l.strip() and
                 (l.strip()[0].isdigit() or l.strip().startswith("generation"))]
    reader = csv.DictReader(lines)
    return list(reader)


def load_dataset(results_dir):
    """Load all interference data. Returns {topo: {gen: [diversity_values]}}."""
    data = defaultdict(lambda: defaultdict(list))
    file_counts = defaultdict(int)
    for topo in TOPOLOGIES:
        for fname in os.listdir(results_dir):
            if not fname.startswith(topo + "_seed"):
                continue
            path = os.path.join(results_dir, fname)
            rows = load_csv(path)
            file_counts[topo] += 1
            for row in rows:
                gen = int(row["generation"])
                data[topo][gen].append(float(row["diversity"]))
    return data, file_counts


def eta_squared(groups):
    """One-way ANOVA eta-squared."""
    all_vals = [x for g in groups for x in g]
    if not all_vals:
        return 0.0
    grand_mean = np.mean(all_vals)
    ss_total = sum((x - grand_mean) ** 2 for x in all_vals)
    if ss_total == 0:
        return 0.0
    ss_between = sum(len(g) * (np.mean(g) - grand_mean) ** 2 for g in groups if g)
    return ss_between / ss_total


def compute_eta2_timeseries(data, topologies):
    """Compute eta^2(config -> diversity) at every generation."""
    all_gens = set()
    for topo in topologies:
        all_gens.update(data[topo].keys())
    all_gens = sorted(all_gens)

    eta2_by_gen = {}
    for gen in all_gens:
        groups = [data[topo].get(gen, []) for topo in topologies]
        if all(len(g) > 0 for g in groups):
            eta2_by_gen[gen] = eta_squared(groups)
    return eta2_by_gen


def compute_pairwise_eta2(data, topo1, topo2):
    """Compute eta^2 between two specific topologies over time."""
    all_gens = sorted(set(data[topo1].keys()) & set(data[topo2].keys()))
    eta2_by_gen = {}
    for gen in all_gens:
        g1 = data[topo1].get(gen, [])
        g2 = data[topo2].get(gen, [])
        if g1 and g2:
            eta2_by_gen[gen] = eta_squared([g1, g2])
    return eta2_by_gen


def logistic(t, C, r, t0):
    """Logistic saturation model: eta2(t) = C / (1 + exp(-r*(t - t0)))."""
    return C / (1.0 + np.exp(-r * (t - t0)))


def fit_logistic(gens, eta2s):
    """Fit logistic saturation model. Returns (C, r, t0, R2) or None."""
    gens = np.array(gens, dtype=float)
    eta2s = np.array(eta2s, dtype=float)

    # Filter out gen 0 and very small values
    mask = (gens >= 10) & (eta2s > 0.001)
    if mask.sum() < 5:
        return None

    g_fit = gens[mask]
    e_fit = eta2s[mask]

    # Initial guesses
    C0 = max(e_fit) * 1.1
    r0 = 0.02
    t0_0 = g_fit[len(g_fit)//3]

    try:
        popt, pcov = curve_fit(logistic, g_fit, e_fit,
                               p0=[C0, r0, t0_0],
                               bounds=([0.001, 0.0001, 0], [1.0, 1.0, 500]),
                               maxfev=10000)
        C, r, t0 = popt

        # R-squared
        predicted = logistic(g_fit, *popt)
        ss_res = np.sum((e_fit - predicted) ** 2)
        ss_tot = np.sum((e_fit - np.mean(e_fit)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

        return C, r, t0, r2
    except (RuntimeError, ValueError):
        return None


def fit_power_law(gens, eta2s, min_gen=10, max_gen=100):
    """Fit log(eta2) = alpha * log(gen) + c in [min_gen, max_gen]."""
    log_g, log_e = [], []
    for g, e in zip(gens, eta2s):
        if min_gen <= g <= max_gen and e > 0.001:
            log_g.append(math.log(g))
            log_e.append(math.log(e))

    if len(log_g) < 3:
        return None, None, None

    log_g = np.array(log_g)
    log_e = np.array(log_e)
    n = len(log_g)

    # OLS in log-log space
    sum_x = log_g.sum()
    sum_y = log_e.sum()
    sum_xy = (log_g * log_e).sum()
    sum_x2 = (log_g ** 2).sum()

    denom = n * sum_x2 - sum_x ** 2
    if denom == 0:
        return None, None, None

    alpha = (n * sum_xy - sum_x * sum_y) / denom
    c = (sum_y - alpha * sum_x) / n

    mean_y = sum_y / n
    ss_tot = ((log_e - mean_y) ** 2).sum()
    ss_res = ((log_e - (alpha * log_g + c)) ** 2).sum()
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    return alpha, c, r2


def analyze_dataset(name, results_dir):
    """Full analysis for one dataset. Returns dict of results."""
    print(f"\n{'='*70}")
    print(f"  DATASET: {name.upper()}")
    print(f"{'='*70}")

    data, counts = load_dataset(results_dir)
    for topo in TOPOLOGIES:
        print(f"  {TOPO_SHORT[topo]}: {counts[topo]} runs")

    results = {"name": name, "configs": {}}

    # --- 3-way eta-squared ---
    eta2_3way = compute_eta2_timeseries(data, TOPOLOGIES)
    gens_3way = sorted(eta2_3way.keys())
    vals_3way = [eta2_3way[g] for g in gens_3way]

    print(f"\n  3-way eta^2 at key generations:")
    for g in [0, 50, 100, 200, 300, 400, 500]:
        if g in eta2_3way:
            print(f"    Gen {g:>4}: eta^2 = {eta2_3way[g]:.4f}")

    # Logistic fit for 3-way
    fit_3way = fit_logistic(gens_3way, vals_3way)
    if fit_3way:
        C, r, t0, r2 = fit_3way
        print(f"\n  3-way logistic fit: C={C:.4f}, r={r:.4f}, t0={t0:.1f}, R2={r2:.4f}")
        results["3way"] = {"ceiling": C, "rate": r, "t0": t0, "R2": r2}
    else:
        print(f"\n  3-way logistic fit: FAILED")
        results["3way"] = None

    # Power law for 3-way (Phase 1: gen 10-100)
    alpha_3, c_3, r2_3 = fit_power_law(gens_3way, vals_3way, 10, 100)
    if alpha_3 is not None:
        print(f"  3-way Phase 1 power law: alpha={alpha_3:.4f}, R2={r2_3:.4f}")
        if results["3way"]:
            results["3way"]["alpha_phase1"] = alpha_3
            results["3way"]["alpha_r2"] = r2_3

    # --- Per-config analysis: each config vs the other two ---
    # This measures how much each specific config deviates
    print(f"\n  Per-config eta^2 (config vs others combined):")
    for topo in TOPOLOGIES:
        others = [t for t in TOPOLOGIES if t != topo]
        # Build pairwise eta^2 between this config and all others combined
        all_gens = sorted(set(data[topo].keys()))
        eta2_ts = {}
        for gen in all_gens:
            this_vals = data[topo].get(gen, [])
            other_vals = []
            for ot in others:
                other_vals.extend(data[ot].get(gen, []))
            if this_vals and other_vals:
                eta2_ts[gen] = eta_squared([this_vals, other_vals])

        gens = sorted(eta2_ts.keys())
        vals = [eta2_ts[g] for g in gens]

        config_result = {}

        # Max eta^2
        if vals:
            max_e = max(vals)
            max_g = gens[vals.index(max_e)]
            print(f"    {TOPO_SHORT[topo]} vs others: max eta^2 = {max_e:.4f} at gen {max_g}")
            config_result["max_eta2"] = max_e
            config_result["max_gen"] = max_g

        # Logistic fit
        fit = fit_logistic(gens, vals)
        if fit:
            C, r, t0, r2 = fit
            print(f"      logistic: C={C:.4f}, r={r:.4f}, t0={t0:.1f}, R2={r2:.4f}")
            config_result["ceiling"] = C
            config_result["rate"] = r
            config_result["t0"] = t0
            config_result["logistic_R2"] = r2

        # Power law
        alpha, c, r2 = fit_power_law(gens, vals, 10, 100)
        if alpha is not None:
            print(f"      Phase 1 power law: alpha={alpha:.4f}, R2={r2:.4f}")
            config_result["alpha_phase1"] = alpha
            config_result["alpha_r2"] = r2

        results["configs"][TOPO_SHORT[topo]] = config_result

    # --- Pairwise analysis ---
    print(f"\n  Pairwise eta^2 saturation ceilings:")
    pairs = [
        ("interference-separated", "interference-adjacent"),
        ("interference-separated", "interference-nested"),
        ("interference-adjacent", "interference-nested"),
    ]
    results["pairwise"] = {}
    for t1, t2 in pairs:
        eta2_pw = compute_pairwise_eta2(data, t1, t2)
        gens = sorted(eta2_pw.keys())
        vals = [eta2_pw[g] for g in gens]

        pair_key = f"{TOPO_SHORT[t1]}_vs_{TOPO_SHORT[t2]}"
        pair_result = {}

        if vals:
            max_e = max(vals)
            max_g = gens[vals.index(max_e)]
            pair_result["max_eta2"] = max_e

        fit = fit_logistic(gens, vals)
        if fit:
            C, r, t0, r2 = fit
            print(f"    {TOPO_SHORT[t1]} vs {TOPO_SHORT[t2]}: C={C:.4f}, r={r:.4f}, t0={t0:.1f}, R2={r2:.4f}")
            pair_result["ceiling"] = C
            pair_result["rate"] = r
            pair_result["t0"] = t0
            pair_result["logistic_R2"] = r2
        else:
            print(f"    {TOPO_SHORT[t1]} vs {TOPO_SHORT[t2]}: logistic fit FAILED (max eta2={max_e:.4f})")

        alpha, c, r2 = fit_power_law(gens, vals, 10, 100)
        if alpha is not None:
            pair_result["alpha_phase1"] = alpha
            pair_result["alpha_r2"] = r2

        results["pairwise"][pair_key] = pair_result

    # --- Mean diversity trajectories (for context) ---
    print(f"\n  Mean diversity at gen 500:")
    for topo in TOPOLOGIES:
        vals = data[topo].get(500, [])
        if vals:
            print(f"    {TOPO_SHORT[topo]}: {np.mean(vals):.4f} +/- {np.std(vals)/np.sqrt(len(vals)):.4f}")

    return results, data, eta2_3way


def plot_results(all_results, all_data, all_eta2):
    """Generate comprehensive plots."""

    n_datasets = len(all_results)
    fig, axes = plt.subplots(n_datasets, 3, figsize=(18, 5 * n_datasets))
    if n_datasets == 1:
        axes = axes.reshape(1, -1)

    for i, (name, results) in enumerate(all_results.items()):
        data = all_data[name]
        eta2_3way = all_eta2[name]

        # Panel 1: Diversity trajectories
        ax = axes[i, 0]
        for topo in TOPOLOGIES:
            gens = sorted(data[topo].keys())
            means = [np.mean(data[topo][g]) for g in gens]
            sems = [np.std(data[topo][g]) / np.sqrt(len(data[topo][g])) for g in gens]
            color = TOPO_COLORS[topo]
            ax.plot(gens, means, color=color, linewidth=2, label=TOPO_SHORT[topo])
            ax.fill_between(gens,
                            [m - s for m, s in zip(means, sems)],
                            [m + s for m, s in zip(means, sems)],
                            color=color, alpha=0.15)
        ax.set_ylabel("Diversity", fontsize=11)
        ax.set_title(f"{name.upper()} — Diversity Trajectories", fontsize=12)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

        # Panel 2: 3-way eta^2 with logistic fit
        ax = axes[i, 1]
        gens = sorted(eta2_3way.keys())
        vals = [eta2_3way[g] for g in gens]
        ax.plot(gens, vals, 'ko-', markersize=3, linewidth=1.5, label='Observed')

        if results["3way"]:
            C = results["3way"]["ceiling"]
            r = results["3way"]["rate"]
            t0 = results["3way"]["t0"]
            t_fit = np.linspace(10, 500, 200)
            ax.plot(t_fit, logistic(t_fit, C, r, t0), 'r--', linewidth=2,
                    label=f'Logistic: C={C:.3f}')
            ax.axhline(y=C, color='r', linestyle=':', alpha=0.5, label=f'Ceiling={C:.3f}')

        ax.set_ylabel("$\\eta^2$", fontsize=11)
        ax.set_title(f"{name.upper()} — 3-way $\\eta^2$ + Logistic Fit", fontsize=12)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(bottom=-0.005)

        # Panel 3: Log-log Phase 1 with power law fits
        ax = axes[i, 2]
        # Plot the 3-way eta^2 in log-log
        log_gens = [g for g in gens if g >= 10 and eta2_3way[g] > 0.001]
        log_vals = [eta2_3way[g] for g in log_gens]
        if log_gens:
            ax.plot(log_gens, log_vals, 'ko-', markersize=3, linewidth=1, label='3-way')

        # Per-config in log-log
        for topo in TOPOLOGIES:
            others = [t for t in TOPOLOGIES if t != topo]
            eta2_config = {}
            for gen in sorted(data[topo].keys()):
                this_vals = data[topo].get(gen, [])
                other_vals = []
                for ot in others:
                    other_vals.extend(data[ot].get(gen, []))
                if this_vals and other_vals:
                    eta2_config[gen] = eta_squared([this_vals, other_vals])

            cg = sorted(g for g in eta2_config if g >= 10 and eta2_config[g] > 0.001)
            cv = [eta2_config[g] for g in cg]
            if cg:
                ax.plot(cg, cv, color=TOPO_COLORS[topo], marker='o', markersize=2,
                        linewidth=1, alpha=0.7, label=f'{TOPO_SHORT[topo]} vs others')

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel("Generation (log)", fontsize=11)
        ax.set_ylabel("$\\eta^2$ (log)", fontsize=11)
        ax.set_title(f"{name.upper()} — Log-Log (Power Law Test)", fontsize=12)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3, which='both')

    for ax in axes[-1]:
        ax.set_xlabel("Generation", fontsize=11)

    plt.tight_layout()
    path = os.path.join(OUTPUT_DIR, "interference_saturation_analysis.png")
    plt.savefig(path, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved: {path}")

    # --- Ceiling comparison plot ---
    fig2, ax = plt.subplots(1, 1, figsize=(10, 6))

    for name, results in all_results.items():
        ceilings = {}
        for config_name, config_data in results["configs"].items():
            if "ceiling" in config_data:
                ceilings[config_name] = config_data["ceiling"]

        if ceilings:
            configs = ["separated", "adjacent", "nested"]
            x_pos = {"separated": 0, "adjacent": 1, "nested": 2}
            dataset_offset = {"onemax": -0.2, "maze": 0.0, "sudoku": 0.2}
            dataset_colors = {"onemax": "#2ecc71", "maze": "#3498db", "sudoku": "#e74c3c"}

            for config in configs:
                if config in ceilings:
                    ax.bar(x_pos[config] + dataset_offset.get(name, 0),
                           ceilings[config], width=0.18,
                           color=dataset_colors.get(name, 'gray'),
                           label=name if config == configs[0] else None,
                           alpha=0.8)

    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["Separated\n(eff-n=2L)", "Adjacent\n(eff-n=2L-1)", "Nested\n(eff-n=2L-2)"])
    ax.set_ylabel("Logistic Ceiling ($C$)", fontsize=12)
    ax.set_title("Saturation Ceilings by Configuration\n(Clio predicts: separated > adjacent > nested)", fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')

    path2 = os.path.join(OUTPUT_DIR, "ceiling_comparison.png")
    plt.savefig(path2, dpi=150, bbox_inches='tight')
    print(f"Ceiling comparison saved: {path2}")


def main():
    print("=" * 70)
    print("  INTERFERENCE SATURATION ANALYSIS — Clio's Predictions")
    print("=" * 70)
    print()
    print("  Clio's predictions:")
    print("    1. Separated has HIGHEST saturation ceiling (effective-n = 2L)")
    print("    2. Adjacent has MIDDLE ceiling (effective-n = 2L-1)")
    print("    3. Nested has LOWEST ceiling (effective-n = 2L-2)")
    print("    4. Phase 1 power-law exponent (alpha) may differ across configs")

    all_results = {}
    all_data = {}
    all_eta2 = {}

    for name, results_dir in DATASETS.items():
        if not os.path.isdir(results_dir):
            print(f"\n  SKIPPING {name}: directory not found ({results_dir})")
            continue

        results, data, eta2_3way = analyze_dataset(name, results_dir)
        all_results[name] = results
        all_data[name] = data
        all_eta2[name] = eta2_3way

    # --- Cross-dataset summary ---
    print(f"\n{'='*70}")
    print(f"  CROSS-DATASET SUMMARY")
    print(f"{'='*70}")

    print(f"\n  {'Dataset':<10} {'Config':<12} {'Ceiling':>10} {'Alpha_P1':>10} {'Max_eta2':>10}")
    print(f"  {'-'*10} {'-'*12} {'-'*10} {'-'*10} {'-'*10}")

    for name, results in all_results.items():
        for config in ["separated", "adjacent", "nested"]:
            cr = results["configs"].get(config, {})
            ceiling = f"{cr['ceiling']:.4f}" if "ceiling" in cr else "N/A"
            alpha = f"{cr['alpha_phase1']:.4f}" if "alpha_phase1" in cr else "N/A"
            max_e = f"{cr['max_eta2']:.4f}" if "max_eta2" in cr else "N/A"
            print(f"  {name:<10} {config:<12} {ceiling:>10} {alpha:>10} {max_e:>10}")

    # --- Test Clio's ordering prediction ---
    print(f"\n  CLIO'S PREDICTION TEST: separated > adjacent > nested (ceiling)")
    for name, results in all_results.items():
        ceilings = {}
        for config in ["separated", "adjacent", "nested"]:
            cr = results["configs"].get(config, {})
            if "ceiling" in cr:
                ceilings[config] = cr["ceiling"]

        if len(ceilings) == 3:
            ordering = sorted(ceilings.keys(), key=lambda k: ceilings[k], reverse=True)
            predicted = ["separated", "adjacent", "nested"]
            match = ordering == predicted
            print(f"    {name}: {' > '.join(ordering)} (ceilings: {', '.join(f'{ceilings[k]:.4f}' for k in ordering)})")
            print(f"      Predicted ordering: {match}")
        else:
            print(f"    {name}: insufficient logistic fits")

    # --- 3-way ceilings ---
    print(f"\n  3-WAY CEILINGS (all configs together):")
    for name, results in all_results.items():
        if results["3way"]:
            print(f"    {name}: C={results['3way']['ceiling']:.4f}, "
                  f"alpha_P1={results['3way'].get('alpha_phase1', 'N/A')}")

    # --- Alpha comparison ---
    print(f"\n  PHASE 1 POWER LAW EXPONENTS (alpha):")
    for name, results in all_results.items():
        alphas = {}
        for config in ["separated", "adjacent", "nested"]:
            cr = results["configs"].get(config, {})
            if "alpha_phase1" in cr:
                alphas[config] = cr["alpha_phase1"]
        if alphas:
            print(f"    {name}: {', '.join(f'{k}={v:.3f}' for k, v in alphas.items())}")
            vals = list(alphas.values())
            if len(vals) >= 2:
                spread = max(vals) - min(vals)
                print(f"      Spread: {spread:.3f} (same if < 0.3)")

    # --- Generate plots ---
    if all_results:
        plot_results(all_results, all_data, all_eta2)

    # --- Save JSON ---
    # Convert to serializable
    json_results = {}
    for name, results in all_results.items():
        jr = {"name": name, "configs": {}, "pairwise": {}}
        if results["3way"]:
            jr["3way"] = {k: float(v) for k, v in results["3way"].items()}
        for k, v in results["configs"].items():
            jr["configs"][k] = {kk: float(vv) for kk, vv in v.items()}
        for k, v in results["pairwise"].items():
            jr["pairwise"][k] = {kk: float(vv) for kk, vv in v.items()}
        json_results[name] = jr

    json_path = os.path.join(OUTPUT_DIR, "saturation_results.json")
    with open(json_path, 'w') as f:
        json.dump(json_results, f, indent=2)
    print(f"\nJSON results saved: {json_path}")

    print(f"\n{'='*70}")
    print(f"  DONE")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
