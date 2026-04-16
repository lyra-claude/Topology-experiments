#!/usr/bin/env python3
"""
Temporal dynamics of the topology effect: power-law analysis.

Computes eta-squared(topology -> diversity) at every generation checkpoint
across NK domains, then fits log(eta^2) = alpha * log(gen) + c to test
Clio's prediction that the topology effect grows as a power law.

Also detects temporal jumps (sudden increases in eta^2 that might indicate
phase transitions rather than smooth growth).
"""

import os
import csv
import math
import sys
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

RESULTS_DIR = "/tmp/directed_multidomain"

DOMAINS = ["nk0", "nk2", "nk4", "nk6"]
DOMAIN_LABELS = {"nk0": "NK0 (smooth)", "nk2": "NK2", "nk4": "NK4", "nk6": "NK6 (rugged)"}

CYCLE_COUNTS = {
    "dag-layer": 0,
    "dag-wide": 0,
    "lowcyc-1": 3,
    "bidir-ring": 10,
    "two-cliques": 14,
    "mesh-cyclic": 20,
    "dense-triangles": 29,
    "ring-skip2": 47,
}


def load_csv(path):
    with open(path) as f:
        lines = [l for l in f if l.strip() and
                 (l.strip()[0].isdigit() or l.strip().startswith("generation"))]
    reader = csv.DictReader(lines)
    return list(reader)


def eta_squared(groups):
    """One-way ANOVA eta-squared."""
    all_vals = []
    for g in groups:
        all_vals.extend(g)
    if not all_vals:
        return 0.0
    grand_mean = sum(all_vals) / len(all_vals)
    ss_total = sum((x - grand_mean)**2 for x in all_vals)
    if ss_total == 0:
        return 0.0
    ss_between = sum(len(g) * (sum(g)/len(g) - grand_mean)**2 for g in groups if g)
    return ss_between / ss_total


def load_domain_data(domain):
    """Load all results for a domain. Returns {topology: {seed: [rows]}}."""
    domain_dir = os.path.join(RESULTS_DIR, domain)
    if not os.path.isdir(domain_dir):
        return {}
    data = defaultdict(dict)
    for fname in os.listdir(domain_dir):
        if not fname.endswith(".csv"):
            continue
        parts = fname.replace(".csv", "").rsplit("_seed", 1)
        if len(parts) != 2:
            continue
        topo, seed = parts[0], int(parts[1])
        rows = load_csv(os.path.join(domain_dir, fname))
        data[topo][seed] = rows
    return data


def compute_eta2_timeseries(data):
    """Compute eta^2(topology -> diversity) at every generation in the data.

    Returns: dict {generation: eta2_value}
    """
    topos = sorted(data.keys(), key=lambda t: CYCLE_COUNTS.get(t, 0))

    # Determine all generations available
    all_gens = set()
    for topo in topos:
        for seed, rows in data[topo].items():
            for row in rows:
                all_gens.add(int(row["generation"]))
    all_gens = sorted(all_gens)

    eta2_by_gen = {}
    for gen in all_gens:
        groups = []
        for topo in topos:
            divs = []
            for seed, rows in data[topo].items():
                for row in rows:
                    if int(row["generation"]) == gen:
                        divs.append(float(row["diversity"]))
                        break
            groups.append(divs)
        eta2_by_gen[gen] = eta_squared(groups)

    return eta2_by_gen


def detect_jumps(gens, eta2s, threshold=2.0):
    """Detect sudden jumps in eta^2 timeseries.

    A jump is where the derivative (delta eta^2 / delta gen) exceeds
    threshold * median_derivative.
    """
    if len(gens) < 3:
        return []

    deltas = []
    for i in range(1, len(gens)):
        d = (eta2s[i] - eta2s[i-1]) / max(gens[i] - gens[i-1], 1)
        deltas.append(d)

    median_delta = sorted(deltas)[len(deltas)//2]
    mad = sorted([abs(d - median_delta) for d in deltas])[len(deltas)//2]
    if mad == 0:
        mad = max(abs(d) for d in deltas) / 10.0

    jumps = []
    for i, d in enumerate(deltas):
        if mad > 0 and abs(d - median_delta) > threshold * mad:
            jumps.append((gens[i+1], eta2s[i+1], d))

    return jumps


def fit_power_law(gens, eta2s, min_gen=10):
    """Fit log(eta^2) = alpha * log(gen) + c for gen >= min_gen and eta^2 > 0.

    Returns: (alpha, c, r_squared)
    """
    log_g = []
    log_e = []
    for g, e in zip(gens, eta2s):
        if g >= min_gen and e > 0.001:  # filter near-zero values
            log_g.append(math.log(g))
            log_e.append(math.log(e))

    if len(log_g) < 3:
        return None, None, None

    # OLS
    n = len(log_g)
    sum_x = sum(log_g)
    sum_y = sum(log_e)
    sum_xy = sum(x*y for x, y in zip(log_g, log_e))
    sum_x2 = sum(x**2 for x in log_g)

    denom = n * sum_x2 - sum_x**2
    if denom == 0:
        return None, None, None

    alpha = (n * sum_xy - sum_x * sum_y) / denom
    c = (sum_y - alpha * sum_x) / n

    # R^2
    mean_y = sum_y / n
    ss_tot = sum((y - mean_y)**2 for y in log_e)
    ss_res = sum((y - (alpha*x + c))**2 for x, y in zip(log_g, log_e))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    return alpha, c, r2


def main():
    output_dir = "/home/lyra/projects/Topology-experiments/experiments"

    print("=" * 70)
    print("  TEMPORAL DYNAMICS OF TOPOLOGY EFFECT — POWER LAW ANALYSIS")
    print("=" * 70)

    all_eta2 = {}

    for domain in DOMAINS:
        print(f"\nLoading {domain}...")
        data = load_domain_data(domain)
        if not data:
            print(f"  WARNING: No data for {domain}")
            continue

        n_runs = sum(len(seeds) for seeds in data.values())
        print(f"  {len(data)} topologies, {n_runs} runs")

        eta2_ts = compute_eta2_timeseries(data)
        all_eta2[domain] = eta2_ts

        gens = sorted(eta2_ts.keys())
        eta2s = [eta2_ts[g] for g in gens]

        # Print timeseries
        print(f"\n  Gen   eta^2")
        print(f"  ---   -----")
        for g in gens:
            bar = "#" * int(eta2_ts[g] * 50)
            print(f"  {g:>3}   {eta2_ts[g]:.4f}  {bar}")

        # Power law fit (excluding gen 0)
        alpha, c, r2 = fit_power_law(gens, eta2s, min_gen=10)
        if alpha is not None:
            print(f"\n  Power law fit (gen >= 10, eta^2 > 0.001):")
            print(f"    log(eta^2) = {alpha:.4f} * log(gen) + {c:.4f}")
            print(f"    alpha (exponent) = {alpha:.4f}")
            print(f"    R^2 = {r2:.4f}")
        else:
            print(f"\n  Power law fit: insufficient data points")

        # Jump detection
        jumps = detect_jumps(gens, eta2s)
        if jumps:
            print(f"\n  Temporal jumps detected (> 2 MAD from median derivative):")
            for g, e, d in jumps:
                print(f"    Gen {g}: eta^2 = {e:.4f}, derivative = {d:.6f}")
        else:
            print(f"\n  No significant temporal jumps detected.")

    # =====================================================================
    # PLOT: log-log eta^2 vs generation
    # =====================================================================
    print(f"\n{'='*70}")
    print("  GENERATING PLOTS")
    print(f"{'='*70}")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    colors = {"nk0": "#4A90D9", "nk2": "#E8A838", "nk4": "#D94A4A", "nk6": "#6B3FA0"}
    markers = {"nk0": "o", "nk2": "s", "nk4": "^", "nk6": "D"}

    # Left: linear scale
    ax = axes[0]
    for domain in DOMAINS:
        if domain not in all_eta2:
            continue
        eta2_ts = all_eta2[domain]
        gens = sorted(eta2_ts.keys())
        eta2s = [eta2_ts[g] for g in gens]
        ax.plot(gens, eta2s, color=colors[domain], marker=markers[domain],
                markersize=4, linewidth=1.5, label=DOMAIN_LABELS[domain],
                markevery=5, alpha=0.9)

    ax.set_xlabel("Generation", fontsize=12)
    ax.set_ylabel("η² (topology → diversity)", fontsize=12)
    ax.set_title("Topology Effect Over Time", fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=-0.02)

    # Right: log-log scale
    ax = axes[1]
    for domain in DOMAINS:
        if domain not in all_eta2:
            continue
        eta2_ts = all_eta2[domain]
        gens = sorted(eta2_ts.keys())
        eta2s = [eta2_ts[g] for g in gens]

        # Filter for log-log: gen > 0, eta^2 > 0
        log_gens = [g for g, e in zip(gens, eta2s) if g > 0 and e > 0.001]
        log_eta2s = [e for g, e in zip(gens, eta2s) if g > 0 and e > 0.001]

        if log_gens:
            ax.plot(log_gens, log_eta2s, color=colors[domain], marker=markers[domain],
                    markersize=4, linewidth=1.5, label=DOMAIN_LABELS[domain],
                    markevery=3, alpha=0.9)

            # Power law fit line
            alpha, c, r2 = fit_power_law(gens, eta2s, min_gen=10)
            if alpha is not None and r2 is not None and r2 > 0.3:
                fit_gens = np.linspace(min(log_gens), max(log_gens), 100)
                fit_eta2 = np.exp(c) * fit_gens**alpha
                ax.plot(fit_gens, fit_eta2, color=colors[domain], linestyle='--',
                        linewidth=1, alpha=0.5)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("Generation (log)", fontsize=12)
    ax.set_ylabel("η² (log)", fontsize=12)
    ax.set_title("Power Law Test (log-log)", fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    plot_path = os.path.join(output_dir, "temporal_dynamics_powerlaw.png")
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    print(f"\n  Plot saved: {plot_path}")

    # =====================================================================
    # PLOT 2: Derivative (rate of change) to spot phase transitions
    # =====================================================================
    fig2, ax2 = plt.subplots(1, 1, figsize=(10, 5))

    for domain in DOMAINS:
        if domain not in all_eta2:
            continue
        eta2_ts = all_eta2[domain]
        gens = sorted(eta2_ts.keys())
        eta2s = [eta2_ts[g] for g in gens]

        if len(gens) < 2:
            continue

        # Compute centered derivative (finite differences)
        deriv_gens = []
        derivs = []
        for i in range(1, len(gens)-1):
            dg = gens[i+1] - gens[i-1]
            if dg > 0:
                deriv_gens.append(gens[i])
                derivs.append((eta2s[i+1] - eta2s[i-1]) / dg)

        ax2.plot(deriv_gens, derivs, color=colors[domain], marker=markers[domain],
                markersize=3, linewidth=1.2, label=DOMAIN_LABELS[domain],
                markevery=3, alpha=0.8)

    ax2.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
    ax2.set_xlabel("Generation", fontsize=12)
    ax2.set_ylabel("d(η²)/d(gen)", fontsize=12)
    ax2.set_title("Rate of Change of Topology Effect", fontsize=13)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)

    deriv_path = os.path.join(output_dir, "temporal_dynamics_derivative.png")
    plt.savefig(deriv_path, dpi=150, bbox_inches='tight')
    print(f"  Derivative plot saved: {deriv_path}")

    # =====================================================================
    # SUMMARY TABLE
    # =====================================================================
    print(f"\n{'='*70}")
    print("  POWER LAW SUMMARY")
    print(f"{'='*70}")
    print(f"  {'Domain':<20} {'alpha':>8} {'R^2':>8} {'eta2_max':>10} {'gen_max':>8}")
    print(f"  {'-'*20} {'-'*8} {'-'*8} {'-'*10} {'-'*8}")

    for domain in DOMAINS:
        if domain not in all_eta2:
            continue
        eta2_ts = all_eta2[domain]
        gens = sorted(eta2_ts.keys())
        eta2s = [eta2_ts[g] for g in gens]

        alpha, c, r2 = fit_power_law(gens, eta2s, min_gen=10)
        max_eta2 = max(eta2s)
        max_gen = gens[eta2s.index(max_eta2)]

        alpha_str = f"{alpha:.4f}" if alpha is not None else "N/A"
        r2_str = f"{r2:.4f}" if r2 is not None else "N/A"

        print(f"  {DOMAIN_LABELS[domain]:<20} {alpha_str:>8} {r2_str:>8} {max_eta2:>10.4f} {max_gen:>8}")

    # Also try piecewise fit (early vs late)
    print(f"\n  PIECEWISE POWER LAW (gen 10-100 vs gen 100-500)")
    print(f"  {'Domain':<20} {'alpha_early':>12} {'R2_early':>10} {'alpha_late':>12} {'R2_late':>10}")
    print(f"  {'-'*20} {'-'*12} {'-'*10} {'-'*12} {'-'*10}")

    for domain in DOMAINS:
        if domain not in all_eta2:
            continue
        eta2_ts = all_eta2[domain]
        gens = sorted(eta2_ts.keys())
        eta2s = [eta2_ts[g] for g in gens]

        # Early: gen 10-100
        early_g = [g for g in gens if 10 <= g <= 100]
        early_e = [eta2_ts[g] for g in early_g]
        a_early, c_early, r2_early = fit_power_law(early_g, early_e, min_gen=10)

        # Late: gen 100-500
        late_g = [g for g in gens if 100 <= g <= 500]
        late_e = [eta2_ts[g] for g in late_g]
        a_late, c_late, r2_late = fit_power_law(late_g, late_e, min_gen=100)

        a_e_str = f"{a_early:.4f}" if a_early is not None else "N/A"
        r2_e_str = f"{r2_early:.4f}" if r2_early is not None else "N/A"
        a_l_str = f"{a_late:.4f}" if a_late is not None else "N/A"
        r2_l_str = f"{r2_late:.4f}" if r2_late is not None else "N/A"

        print(f"  {DOMAIN_LABELS[domain]:<20} {a_e_str:>12} {r2_e_str:>10} {a_l_str:>12} {r2_l_str:>10}")

    print(f"\n{'='*70}")
    print("  DONE")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
