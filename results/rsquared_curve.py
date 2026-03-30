#!/usr/bin/env python3
"""
Generate publication-quality η² (eta-squared) curve over generations for GECCO paper.

Combines two datasets:
  - 10-run transient data (gens 0-100): higher statistical power
  - 3-run pilot data (gens 110-500): extended view

Computes one-way ANOVA effect size (η²) at each generation checkpoint,
with bootstrap 95% confidence intervals.

Output:
  - results/rsquared_over_time.png (300 DPI raster)
  - results/rsquared_over_time.pdf (vector)
  - Summary table printed to stdout
"""

import os
import csv
import io
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from collections import defaultdict

# ─── Configuration ────────────────────────────────────────────────────────────

TRANSIENT_DIR = os.path.join(os.path.dirname(__file__), "transient_onemax_10runs")
PILOT_DIR = os.path.join(os.path.dirname(__file__), "pilot_onemax")

TOPOLOGIES = [
    "disconnected", "ring", "star", "complete",
    "hypercube", "barbell", "watts-strogatz", "random-regular",
]

TRANSIENT_SEEDS = [42, 137, 2718, 314, 1618, 271, 577, 997, 1729, 4242]
PILOT_SEEDS = [42, 137, 2718]

N_BOOTSTRAP = 2000
RNG_SEED = 42

# ─── Data Loading ─────────────────────────────────────────────────────────────

def load_csv(filepath):
    """Load a CSV, handling potential cabal output prefix."""
    with open(filepath) as f:
        lines = f.readlines()
    # Find header line
    header_idx = 0
    for i, line in enumerate(lines):
        if line.strip().startswith("generation,"):
            header_idx = i
            break
    csv_text = "".join(lines[header_idx:])
    reader = csv.DictReader(io.StringIO(csv_text))
    data = {}
    for row in reader:
        gen = int(row['generation'])
        data[gen] = {
            'meanFitness': float(row['meanFitness']),
            'bestFitness': float(row['bestFitness']),
            'diversity': float(row['diversity']),
        }
    return data


def load_dataset(results_dir, seeds):
    """Load all CSVs for a dataset -> dict[topology][seed_idx] = {gen: {metrics}}"""
    all_data = defaultdict(list)
    for topo in TOPOLOGIES:
        for seed in seeds:
            filepath = os.path.join(results_dir, f"{topo}_seed{seed}.csv")
            if not os.path.exists(filepath):
                print(f"WARNING: Missing {filepath}")
                continue
            all_data[topo].append(load_csv(filepath))
    return all_data


# ─── Eta-squared Computation ─────────────────────────────────────────────────

def compute_eta_squared(groups):
    """
    Compute η² from a list of groups (each group is a list of values).
    η² = SS_between / SS_total
    Returns η² value, or NaN if undefined.
    """
    all_vals = []
    for g in groups:
        all_vals.extend(g)

    if len(all_vals) < 2:
        return np.nan

    grand_mean = np.mean(all_vals)
    ss_total = np.sum((np.array(all_vals) - grand_mean) ** 2)

    if ss_total == 0:
        return 0.0

    ss_between = 0.0
    for g in groups:
        if len(g) == 0:
            continue
        group_mean = np.mean(g)
        ss_between += len(g) * (group_mean - grand_mean) ** 2

    return ss_between / ss_total


def compute_eta_squared_at_gen(data, gen, metric):
    """
    Compute η² for a given metric at a given generation.
    data: dict[topology] -> list of run dicts (each run: {gen: {metrics}})
    """
    groups = []
    for topo in TOPOLOGIES:
        vals = []
        for run in data[topo]:
            if gen in run:
                vals.append(run[gen][metric])
        groups.append(vals)

    # Only compute if all groups have data
    if any(len(g) == 0 for g in groups):
        return np.nan

    return compute_eta_squared(groups)


def bootstrap_eta_squared(data, gen, metric, n_bootstrap=N_BOOTSTRAP, rng=None):
    """
    Bootstrap η² confidence interval.
    Resamples runs within each topology group.
    Returns (eta_sq, ci_low, ci_high).
    """
    if rng is None:
        rng = np.random.default_rng(RNG_SEED)

    # Collect groups
    groups = {}
    for topo in TOPOLOGIES:
        vals = []
        for run in data[topo]:
            if gen in run:
                vals.append(run[gen][metric])
        groups[topo] = np.array(vals)

    if any(len(v) == 0 for v in groups.values()):
        return np.nan, np.nan, np.nan

    # Observed η²
    observed = compute_eta_squared([groups[t].tolist() for t in TOPOLOGIES])

    # Bootstrap
    boot_vals = np.empty(n_bootstrap)
    for b in range(n_bootstrap):
        boot_groups = []
        for topo in TOPOLOGIES:
            g = groups[topo]
            idx = rng.integers(0, len(g), size=len(g))
            boot_groups.append(g[idx].tolist())
        boot_vals[b] = compute_eta_squared(boot_groups)

    ci_low = np.nanpercentile(boot_vals, 2.5)
    ci_high = np.nanpercentile(boot_vals, 97.5)

    return observed, ci_low, ci_high


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    print("Loading datasets...")
    transient_data = load_dataset(TRANSIENT_DIR, TRANSIENT_SEEDS)
    pilot_data = load_dataset(PILOT_DIR, PILOT_SEEDS)

    # Verify
    for topo in TOPOLOGIES:
        nt = len(transient_data[topo])
        np_ = len(pilot_data[topo])
        if nt != 10:
            print(f"  WARNING: transient {topo} has {nt} runs (expected 10)")
        if np_ != 3:
            print(f"  WARNING: pilot {topo} has {np_} runs (expected 3)")

    # Generation checkpoints
    transient_gens = list(range(0, 101, 10))   # 0, 10, ..., 100
    pilot_gens = list(range(110, 501, 10))      # 110, 120, ..., 500

    rng = np.random.default_rng(RNG_SEED)

    # ─── Compute η² for transient (10-run, gens 0-100) ───────────────────
    print("Computing η² for transient data (10 runs, gens 0-100)...")
    transient_results = []
    for gen in transient_gens:
        fit_eta, fit_lo, fit_hi = bootstrap_eta_squared(
            transient_data, gen, 'meanFitness', rng=rng)
        div_eta, div_lo, div_hi = bootstrap_eta_squared(
            transient_data, gen, 'diversity', rng=rng)
        transient_results.append({
            'gen': gen,
            'fit_eta': fit_eta, 'fit_lo': fit_lo, 'fit_hi': fit_hi,
            'div_eta': div_eta, 'div_lo': div_lo, 'div_hi': div_hi,
            'source': '10-run',
        })
        print(f"  Gen {gen:>3d}: η²(fitness)={fit_eta:.4f} [{fit_lo:.4f}, {fit_hi:.4f}]  "
              f"η²(diversity)={div_eta:.4f} [{div_lo:.4f}, {div_hi:.4f}]")

    # ─── Compute η² for pilot (3-run, gens 110-500) ──────────────────────
    print("\nComputing η² for pilot data (3 runs, gens 110-500)...")
    pilot_results = []
    for gen in pilot_gens:
        fit_eta, fit_lo, fit_hi = bootstrap_eta_squared(
            pilot_data, gen, 'meanFitness', rng=rng)
        div_eta, div_lo, div_hi = bootstrap_eta_squared(
            pilot_data, gen, 'diversity', rng=rng)
        pilot_results.append({
            'gen': gen,
            'fit_eta': fit_eta, 'fit_lo': fit_lo, 'fit_hi': fit_hi,
            'div_eta': div_eta, 'div_lo': div_lo, 'div_hi': div_hi,
            'source': '3-run',
        })

    # Also compute overlapping gens 0-100 from pilot for comparison
    print("\nComputing η² for pilot data at overlapping gens 0-100 (for comparison)...")
    pilot_overlap = []
    for gen in transient_gens:
        fit_eta = compute_eta_squared_at_gen(pilot_data, gen, 'meanFitness')
        div_eta = compute_eta_squared_at_gen(pilot_data, gen, 'diversity')
        pilot_overlap.append({'gen': gen, 'fit_eta': fit_eta, 'div_eta': div_eta})
        print(f"  Gen {gen:>3d}: η²(fitness)={fit_eta:.4f}  η²(diversity)={div_eta:.4f}")

    # ─── Combine results ─────────────────────────────────────────────────
    all_results = transient_results + pilot_results

    # ─── Find peaks ───────────────────────────────────────────────────────
    fit_peak = max(all_results, key=lambda r: r['fit_eta'] if not np.isnan(r['fit_eta']) else -1)
    div_peak = max(all_results, key=lambda r: r['div_eta'] if not np.isnan(r['div_eta']) else -1)

    print(f"\n{'='*70}")
    print(f"PEAK η² (fitness):   {fit_peak['fit_eta']:.4f} at gen {fit_peak['gen']} ({fit_peak['source']})")
    print(f"PEAK η² (diversity): {div_peak['div_eta']:.4f} at gen {div_peak['gen']} ({div_peak['source']})")
    print(f"{'='*70}")

    # ─── Summary Table ────────────────────────────────────────────────────
    key_gens = [10, 20, 30, 50, 100, 200, 500]
    print(f"\n{'='*70}")
    print(f"SUMMARY TABLE: η² at key generations")
    print(f"{'='*70}")
    print(f"{'Gen':>5}  {'η²(fitness)':>14}  {'95% CI':>22}  {'η²(diversity)':>14}  {'95% CI':>22}  {'Source':>8}")
    print("-" * 95)

    for gen in key_gens:
        match = [r for r in all_results if r['gen'] == gen]
        if match:
            r = match[0]
            fit_ci = f"[{r['fit_lo']:.4f}, {r['fit_hi']:.4f}]"
            div_ci = f"[{r['div_lo']:.4f}, {r['div_hi']:.4f}]"
            print(f"{gen:>5}  {r['fit_eta']:>14.4f}  {fit_ci:>22}  {r['div_eta']:>14.4f}  {div_ci:>22}  {r['source']:>8}")
        else:
            print(f"{gen:>5}  {'N/A':>14}  {'':>22}  {'N/A':>14}  {'':>22}")

    # ─── Generate Figure ──────────────────────────────────────────────────
    print("\nGenerating publication figure...")

    # Separate data for the two regimes
    t_gens = np.array([r['gen'] for r in transient_results])
    t_fit = np.array([r['fit_eta'] for r in transient_results])
    t_fit_lo = np.array([r['fit_lo'] for r in transient_results])
    t_fit_hi = np.array([r['fit_hi'] for r in transient_results])
    t_div = np.array([r['div_eta'] for r in transient_results])
    t_div_lo = np.array([r['div_lo'] for r in transient_results])
    t_div_hi = np.array([r['div_hi'] for r in transient_results])

    p_gens = np.array([r['gen'] for r in pilot_results])
    p_fit = np.array([r['fit_eta'] for r in pilot_results])
    p_fit_lo = np.array([r['fit_lo'] for r in pilot_results])
    p_fit_hi = np.array([r['fit_hi'] for r in pilot_results])
    p_div = np.array([r['div_eta'] for r in pilot_results])
    p_div_lo = np.array([r['div_lo'] for r in pilot_results])
    p_div_hi = np.array([r['div_hi'] for r in pilot_results])

    # Color palette: professional, colorblind-friendly
    COLOR_FIT = '#2166ac'    # Blue
    COLOR_DIV = '#b2182b'    # Red
    COLOR_FIT_LIGHT = '#92c5de'
    COLOR_DIV_LIGHT = '#f4a582'

    # ─── Figure setup ─────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 4), dpi=300)

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.8)
    ax.spines['bottom'].set_linewidth(0.8)

    # ─── 10-run data (solid lines, filled CI) ────────────────────────────
    # Confidence intervals
    ax.fill_between(t_gens, t_fit_lo, t_fit_hi,
                    alpha=0.20, color=COLOR_FIT, linewidth=0)
    ax.fill_between(t_gens, t_div_lo, t_div_hi,
                    alpha=0.20, color=COLOR_DIV, linewidth=0)

    # Main lines
    ax.plot(t_gens, t_fit, color=COLOR_FIT, linewidth=1.8, marker='o',
            markersize=4, markeredgewidth=0, zorder=5)
    ax.plot(t_gens, t_div, color=COLOR_DIV, linewidth=1.8, marker='s',
            markersize=4, markeredgewidth=0, zorder=5)

    # ─── 3-run data (thin scatter + rolling mean, NO CI bands) ────────────
    # The 3-run CIs are too wide to be informative; show raw points + trend
    # Connect with last transient point for visual continuity
    bridge_fit_gens = np.concatenate([[t_gens[-1]], p_gens])
    bridge_fit_vals = np.concatenate([[t_fit[-1]], p_fit])
    bridge_div_gens = np.concatenate([[t_gens[-1]], p_gens])
    bridge_div_vals = np.concatenate([[t_div[-1]], p_div])

    # Raw points (small, semi-transparent)
    ax.scatter(p_gens, p_fit, color=COLOR_FIT, s=8, alpha=0.30, zorder=3,
               edgecolors='none')
    ax.scatter(p_gens, p_div, color=COLOR_DIV, s=8, alpha=0.30, zorder=3,
               edgecolors='none', marker='s')

    # Rolling mean (window=5 = 50 gens) for trend line
    def rolling_mean(vals, window=5):
        """Simple centered rolling mean with edge handling."""
        result = np.empty_like(vals)
        for i in range(len(vals)):
            lo = max(0, i - window // 2)
            hi = min(len(vals), i + window // 2 + 1)
            result[i] = np.mean(vals[lo:hi])
        return result

    p_fit_smooth = rolling_mean(p_fit)
    p_div_smooth = rolling_mean(p_div)

    # Smooth trend (dashed, connecting from last transient point)
    ax.plot(np.concatenate([[t_gens[-1]], p_gens]),
            np.concatenate([[t_fit[-1]], p_fit_smooth]),
            color=COLOR_FIT, linewidth=1.0, linestyle='--', alpha=0.55, zorder=4)
    ax.plot(np.concatenate([[t_gens[-1]], p_gens]),
            np.concatenate([[t_div[-1]], p_div_smooth]),
            color=COLOR_DIV, linewidth=1.0, linestyle='--', alpha=0.55, zorder=4)

    # ─── Vertical line at regime boundary ─────────────────────────────────
    ax.axvline(x=100, color='gray', linewidth=0.6, linestyle=':', alpha=0.5)
    ax.text(105, 0.94, r'$n\!=\!3$', fontsize=7, color='gray',
            ha='left', va='top')
    ax.text(95, 0.94, r'$n\!=\!10$', fontsize=7, color='gray',
            ha='right', va='top')

    # ─── Peak annotations ─────────────────────────────────────────────────
    # Fitness peak — annotate slightly above and to the right
    ax.annotate(
        r'$\eta^2 = $' + f'{fit_peak["fit_eta"]:.2f}',
        xy=(fit_peak['gen'], fit_peak['fit_eta']),
        xytext=(fit_peak['gen'] + 40, fit_peak['fit_eta'] + 0.06),
        fontsize=8,
        color=COLOR_FIT,
        fontweight='bold',
        arrowprops=dict(arrowstyle='->', color=COLOR_FIT, lw=0.8,
                        connectionstyle='arc3,rad=0.15'),
        ha='left', va='bottom',
    )

    # Diversity peak — annotate below to avoid overlap
    ax.annotate(
        r'$\eta^2 = $' + f'{div_peak["div_eta"]:.2f}',
        xy=(div_peak['gen'], div_peak['div_eta']),
        xytext=(div_peak['gen'] + 40, div_peak['div_eta'] - 0.15),
        fontsize=8,
        color=COLOR_DIV,
        fontweight='bold',
        arrowprops=dict(arrowstyle='->', color=COLOR_DIV, lw=0.8,
                        connectionstyle='arc3,rad=-0.15'),
        ha='left', va='top',
    )

    # ─── Axis labels and formatting ───────────────────────────────────────
    ax.set_xlabel('Generation', fontsize=10, labelpad=6)
    ax.set_ylabel(r'$\eta^2$ (effect size)', fontsize=10, labelpad=6)

    ax.set_xlim(-5, 510)
    ax.set_ylim(-0.02, 1.02)

    ax.tick_params(axis='both', which='major', labelsize=8.5, length=3, width=0.6)

    # X-axis ticks
    ax.set_xticks([0, 50, 100, 200, 300, 400, 500])

    # Y-axis ticks
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

    # ─── Legend ───────────────────────────────────────────────────────────
    legend_elements = [
        Line2D([0], [0], color=COLOR_FIT, linewidth=1.8, marker='o',
               markersize=4, markeredgewidth=0, label='Mean fitness'),
        Line2D([0], [0], color=COLOR_DIV, linewidth=1.8, marker='s',
               markersize=4, markeredgewidth=0, label='Diversity'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', frameon=False,
              fontsize=8.5, handlelength=2.5)

    # ─── Effect size reference lines ──────────────────────────────────────
    # Cohen's benchmarks for η²: small=0.01, medium=0.06, large=0.14
    for level, label in [(0.14, 'large'), (0.06, 'medium')]:
        ax.axhline(y=level, color='gray', linewidth=0.4, linestyle='-', alpha=0.25)
        ax.text(505, level, label, fontsize=6, color='gray', alpha=0.5,
                ha='left', va='center')

    plt.tight_layout()

    # ─── Save ─────────────────────────────────────────────────────────────
    out_dir = os.path.dirname(__file__)
    png_path = os.path.join(out_dir, "rsquared_over_time.png")
    pdf_path = os.path.join(out_dir, "rsquared_over_time.pdf")

    fig.savefig(png_path, dpi=300, bbox_inches='tight', pad_inches=0.1)
    fig.savefig(pdf_path, bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)

    print(f"\nFigure saved to:")
    print(f"  PNG: {png_path}")
    print(f"  PDF: {pdf_path}")

    # ─── Extended view analysis ───────────────────────────────────────────
    print(f"\n{'='*70}")
    print("EXTENDED VIEW ANALYSIS: Does the 500-gen window add information?")
    print(f"{'='*70}")

    # Check if there's meaningful η² beyond gen 100
    post100_fit = [r['fit_eta'] for r in pilot_results if not np.isnan(r['fit_eta'])]
    post100_div = [r['div_eta'] for r in pilot_results if not np.isnan(r['div_eta'])]

    if post100_fit:
        print(f"\nPost-100 η²(fitness):  mean={np.mean(post100_fit):.4f}, "
              f"max={np.max(post100_fit):.4f}, min={np.min(post100_fit):.4f}")
    if post100_div:
        print(f"Post-100 η²(diversity): mean={np.mean(post100_div):.4f}, "
              f"max={np.max(post100_div):.4f}, min={np.min(post100_div):.4f}")

    # Compare with peak
    print(f"\nPeak η²(fitness) = {fit_peak['fit_eta']:.4f} at gen {fit_peak['gen']}")
    if post100_fit:
        pct_of_peak = np.mean(post100_fit) / fit_peak['fit_eta'] * 100 if fit_peak['fit_eta'] > 0 else 0
        print(f"Post-100 mean is {pct_of_peak:.1f}% of peak")

    print(f"\nPeak η²(diversity) = {div_peak['div_eta']:.4f} at gen {div_peak['gen']}")
    if post100_div:
        pct_of_peak_d = np.mean(post100_div) / div_peak['div_eta'] * 100 if div_peak['div_eta'] > 0 else 0
        print(f"Post-100 mean is {pct_of_peak_d:.1f}% of peak")

    # ─── Write CSV of all results ─────────────────────────────────────────
    csv_path = os.path.join(out_dir, "rsquared_curve_data.csv")
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['generation', 'eta_sq_fitness', 'ci_lo_fitness', 'ci_hi_fitness',
                         'eta_sq_diversity', 'ci_lo_diversity', 'ci_hi_diversity', 'source'])
        for r in all_results:
            writer.writerow([
                r['gen'], f"{r['fit_eta']:.6f}", f"{r['fit_lo']:.6f}", f"{r['fit_hi']:.6f}",
                f"{r['div_eta']:.6f}", f"{r['div_lo']:.6f}", f"{r['div_hi']:.6f}",
                r['source'],
            ])
    print(f"\nData saved to: {csv_path}")

    print("\nDone.")


if __name__ == "__main__":
    main()
