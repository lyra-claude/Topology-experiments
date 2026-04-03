#!/usr/bin/env python3
"""
Analyze Foster census sweep: 13 Foster-graph topologies x 30 seeds.

Tests whether vertex count (n) and cycle rank (beta_1) predict GA performance
across a principled family of 3-regular graphs with increasing size.

Run when results/foster_sweep/ is fully populated.
"""

import os
import csv
import math
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from scipy import stats

RESULTS_DIR = "results/foster_sweep"
FIG_DIR = "results/foster_sweep/figures"

# Topology metadata: n = vertex count, beta1 = cycle rank = |E| - n + 1
# All Foster graphs are 3-regular, so |E| = 3n/2, giving beta1 = n/2 + 1
TOPOLOGY_META = {
    "k4":            {"n": 4,  "beta1": 3,  "density": 3.0 / (4  - 1), "girth": 3, "diameter": 1},
    "k33":           {"n": 6,  "beta1": 4,  "density": 3.0 / (6  - 1), "girth": 4, "diameter": 2},
    "cube":          {"n": 8,  "beta1": 5,  "density": 3.0 / (8  - 1), "girth": 4, "diameter": 3},
    "petersen":      {"n": 10, "beta1": 6,  "density": 3.0 / (10 - 1), "girth": 5, "diameter": 2},
    "heawood":       {"n": 14, "beta1": 8,  "density": 3.0 / (14 - 1), "girth": 6, "diameter": 3},
    "mobius-kantor": {"n": 16, "beta1": 9,  "density": 3.0 / (16 - 1), "girth": 6, "diameter": 4},
    "pappus":        {"n": 18, "beta1": 10, "density": 3.0 / (18 - 1), "girth": 6, "diameter": 4},
    "dodecahedron":  {"n": 20, "beta1": 11, "density": 3.0 / (20 - 1), "girth": 5, "diameter": 5},
    "desargues":     {"n": 20, "beta1": 11, "density": 3.0 / (20 - 1), "girth": 6, "diameter": 5},
    "nauru":         {"n": 24, "beta1": 13, "density": 3.0 / (24 - 1), "girth": 6, "diameter": 4},
    "f26a":          {"n": 26, "beta1": 14, "density": 3.0 / (26 - 1), "girth": 6, "diameter": 5},
    "coxeter":       {"n": 28, "beta1": 15, "density": 3.0 / (28 - 1), "girth": 7, "diameter": 4},
    "tutte-coxeter": {"n": 30, "beta1": 16, "density": 3.0 / (30 - 1), "girth": 8, "diameter": 4},
}

KEY_GENERATIONS = [0, 10, 30, 50, 100, 200, 500]


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_csv(path):
    """Load a run CSV, tolerating any comment lines before the header."""
    with open(path) as f:
        lines = [l for l in f if l.strip() and
                 (l.strip()[0].isdigit() or l.strip().startswith("generation"))]
    reader = csv.DictReader(lines)
    return list(reader)


def load_all_data():
    """Return dict: topo_name -> list of runs, each run is a list of row-dicts."""
    results = defaultdict(list)
    for fname in sorted(os.listdir(RESULTS_DIR)):
        if not fname.endswith(".csv"):
            continue
        topo = fname.rsplit("_seed", 1)[0]
        if topo not in TOPOLOGY_META:
            continue
        rows = load_csv(os.path.join(RESULTS_DIR, fname))
        if rows:
            results[topo].append(rows)
    return results


def get_all_generations(results):
    gens = set()
    for runs in results.values():
        for run in runs:
            for row in run:
                gens.add(int(row["generation"]))
    return sorted(gens)


# ---------------------------------------------------------------------------
# Statistics helpers (stdlib only, scipy used only in figures for p-values)
# ---------------------------------------------------------------------------

def _mean(xs):
    return sum(xs) / len(xs) if xs else float("nan")


def _std(xs):
    if len(xs) < 2:
        return 0.0
    m = _mean(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))


def _stderr(xs):
    return _std(xs) / math.sqrt(len(xs)) if len(xs) >= 2 else 0.0


def _pearson_r(xs, ys):
    n = len(xs)
    if n < 3:
        return float("nan")
    mx, my = _mean(xs), _mean(ys)
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if sx == 0 or sy == 0:
        return float("nan")
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / (sx * sy)


def _eta_squared(groups):
    """groups: list of lists of floats. Returns eta^2."""
    all_vals = [v for g in groups for v in g]
    if not all_vals:
        return float("nan")
    grand_mean = _mean(all_vals)
    ss_total = sum((x - grand_mean) ** 2 for x in all_vals)
    if ss_total == 0:
        return 0.0
    ss_between = sum(len(g) * (_mean(g) - grand_mean) ** 2 for g in groups if g)
    return ss_between / ss_total


def _one_way_f(groups):
    """Returns (F, p) for one-way ANOVA using scipy.stats.f_oneway."""
    non_empty = [g for g in groups if len(g) >= 2]
    if len(non_empty) < 2:
        return float("nan"), float("nan")
    return stats.f_oneway(*non_empty)


def extract_at_gen(results, topo, gen, metric):
    """Return list of per-seed values for a given topology/gen/metric."""
    vals = []
    for run in results[topo]:
        for row in run:
            if int(row["generation"]) == gen:
                try:
                    vals.append(float(row[metric]))
                except (KeyError, ValueError):
                    pass
                break
    return vals


def topos_sorted_by_n(results):
    """Return topology names present in data, sorted by vertex count then name."""
    present = [t for t in TOPOLOGY_META if t in results and results[t]]
    return sorted(present, key=lambda t: (TOPOLOGY_META[t]["n"], t))


# ---------------------------------------------------------------------------
# Trajectory builder (used by figures)
# ---------------------------------------------------------------------------

def get_metric_trajectories(results, metric="diversity"):
    """Return dict: topo -> (gens_array, means_array, se_array)."""
    all_gens = get_all_generations(results)
    trajectories = {}
    for topo in topos_sorted_by_n(results):
        gens_out, means_out, ses_out = [], [], []
        for gen in all_gens:
            vals = extract_at_gen(results, topo, gen, metric)
            if vals:
                gens_out.append(gen)
                means_out.append(_mean(vals))
                ses_out.append(_stderr(vals))
        trajectories[topo] = (
            np.array(gens_out),
            np.array(means_out),
            np.array(ses_out),
        )
    return trajectories


# ---------------------------------------------------------------------------
# Text analysis
# ---------------------------------------------------------------------------

def run_analysis(results):
    topos = topos_sorted_by_n(results)
    n_seeds = {t: len(results[t]) for t in topos}

    print("FOSTER CENSUS SWEEP — ANALYSIS")
    print(f"({len(topos)} topologies loaded, up to 30 seeds each)")
    print("=" * 100)
    print()

    # ------------------------------------------------------------------
    # 1. Summary table sorted by vertex count
    # ------------------------------------------------------------------
    print("SUMMARY TABLE (sorted by n)")
    header = (f"{'Topology':<16} {'n':>4} {'β₁':>4} {'Seeds':>6} "
              f"{'Div(gen30)':>12} {'±SE':>8} "
              f"{'Fit(gen100)':>12} {'±SE':>8} "
              f"{'Fit(gen500)':>12} {'±SE':>8}")
    print(header)
    print("-" * len(header))

    for topo in topos:
        meta = TOPOLOGY_META[topo]
        divs30  = extract_at_gen(results, topo, 30,  "diversity")
        fits100 = extract_at_gen(results, topo, 100, "meanFitness")
        fits500 = extract_at_gen(results, topo, 500, "meanFitness")

        def fmt(xs):
            if not xs:
                return f"{'N/A':>12}", f"{'':>8}"
            return f"{_mean(xs):>12.6f}", f"{_stderr(xs):>8.4f}"

        d30m, d30s   = fmt(divs30)
        f100m, f100s = fmt(fits100)
        f500m, f500s = fmt(fits500)

        print(f"{topo:<16} {meta['n']:>4} {meta['beta1']:>4} {n_seeds[topo]:>6} "
              f"{d30m}{d30s} {f100m}{f100s} {f500m}{f500s}")

    print()

    # ------------------------------------------------------------------
    # 2. Transient diversity table at key generations
    # ------------------------------------------------------------------
    print("DIVERSITY BY GENERATION (mean ± SE)")
    col_header = f"{'Gen':>5}"
    for topo in topos:
        col_header += f"  {topo[:12]:>14}"
    print(col_header)
    print("-" * len(col_header))

    for gen in KEY_GENERATIONS:
        row = f"{gen:>5}"
        for topo in topos:
            vals = extract_at_gen(results, topo, gen, "diversity")
            if vals:
                row += f"  {_mean(vals):6.4f}±{_stderr(vals):.4f}"
            else:
                row += f"  {'N/A':>14}"
        print(row)

    print()

    # ------------------------------------------------------------------
    # 3. Pearson correlations at key generations
    # ------------------------------------------------------------------
    print("PEARSON CORRELATIONS AT KEY GENERATIONS")
    print(f"{'Gen':>5} {'r(n,div)':>10} {'r(β₁,div)':>11} {'r(n,fit)':>10} {'r(β₁,fit)':>11}")
    print("-" * 55)

    for gen in KEY_GENERATIONS:
        ns, b1s, divs, fits = [], [], [], []
        for topo in topos:
            dv = extract_at_gen(results, topo, gen, "diversity")
            ft = extract_at_gen(results, topo, gen, "meanFitness")
            if dv and ft:
                ns.append(TOPOLOGY_META[topo]["n"])
                b1s.append(TOPOLOGY_META[topo]["beta1"])
                divs.append(_mean(dv))
                fits.append(_mean(ft))

        if len(ns) >= 3:
            r_n_div  = _pearson_r(ns, divs)
            r_b1_div = _pearson_r(b1s, divs)
            r_n_fit  = _pearson_r(ns, fits)
            r_b1_fit = _pearson_r(b1s, fits)

            def fmt_r(r):
                return f"{r:>10.4f}" if not math.isnan(r) else f"{'N/A':>10}"

            print(f"{gen:>5} {fmt_r(r_n_div)} {fmt_r(r_b1_div)} {fmt_r(r_n_fit)} {fmt_r(r_b1_fit)}")
        else:
            print(f"{gen:>5}  (insufficient data)")

    print()
    print("Note: n and β₁ are collinear for 3-regular graphs (β₁ = n/2 + 1).")
    print("Both columns shown for completeness but they measure the same underlying variable.")
    print()

    # ------------------------------------------------------------------
    # 4. ANOVA (one-way F test) at gen 30
    # ------------------------------------------------------------------
    print("ONE-WAY ANOVA AT GEN 30 (topology as factor)")
    for metric, label in [("diversity", "Diversity"), ("meanFitness", "Mean Fitness")]:
        groups = [extract_at_gen(results, t, 30, metric) for t in topos]
        F, p = _one_way_f(groups)
        eta2 = _eta_squared([g for g in groups if g])
        if not math.isnan(F):
            print(f"  {label:<15}: F = {F:.3f},  p = {p:.4f},  η² = {eta2:.4f}")
        else:
            print(f"  {label:<15}: insufficient data")

    print()

    # ------------------------------------------------------------------
    # 5. Effect size (eta^2) at key generations
    # ------------------------------------------------------------------
    print("EFFECT SIZE (η²) AT KEY GENERATIONS")
    print(f"{'Gen':>5} {'η²(div)':>10} {'η²(fit)':>10}")
    print("-" * 30)

    for gen in KEY_GENERATIONS:
        div_groups = [extract_at_gen(results, t, gen, "diversity")    for t in topos]
        fit_groups = [extract_at_gen(results, t, gen, "meanFitness")  for t in topos]
        e2d = _eta_squared([g for g in div_groups if g])
        e2f = _eta_squared([g for g in fit_groups if g])
        print(f"{gen:>5} {e2d:>10.4f} {e2f:>10.4f}")

    print()

    # ------------------------------------------------------------------
    # 6. Convergence speed: first generation reaching mean fitness >= 0.95
    # ------------------------------------------------------------------
    print("CONVERGENCE SPEED: first gen with mean fitness ≥ 0.95")
    print(f"{'Topology':<16} {'n':>4} {'β₁':>4} {'Conv.Gen':>10}")
    print("-" * 40)

    all_gens = get_all_generations(results)

    for topo in topos:
        meta = TOPOLOGY_META[topo]
        conv_gen = None
        for gen in all_gens:
            fits = extract_at_gen(results, topo, gen, "meanFitness")
            if fits and _mean(fits) >= 0.95:
                conv_gen = gen
                break
        if conv_gen is not None:
            print(f"{topo:<16} {meta['n']:>4} {meta['beta1']:>4} {conv_gen:>10}")
        else:
            print(f"{topo:<16} {meta['n']:>4} {meta['beta1']:>4} {'not reached':>10}")

    print()

    # ------------------------------------------------------------------
    # 7. Dodecahedron vs Desargues natural control
    # ------------------------------------------------------------------
    if "dodecahedron" in results and "desargues" in results:
        print("DODECAHEDRON vs DESARGUES — NATURAL CONTROL")
        print("Both: n=20, β₁=11, density=3/19 ≈ 0.158.  Differ in girth (5 vs 6) and structure.")
        print()

        for gen in KEY_GENERATIONS:
            dod_div = extract_at_gen(results, "dodecahedron", gen, "diversity")
            des_div = extract_at_gen(results, "desargues",    gen, "diversity")
            dod_fit = extract_at_gen(results, "dodecahedron", gen, "meanFitness")
            des_fit = extract_at_gen(results, "desargues",    gen, "meanFitness")

            row = f"  Gen {gen:>3}:  "
            if dod_div and des_div:
                row += f"Div: dod={_mean(dod_div):.4f}±{_stderr(dod_div):.4f}  des={_mean(des_div):.4f}±{_stderr(des_div):.4f}"
            if dod_fit and des_fit:
                row += f"  |  Fit: dod={_mean(dod_fit):.4f}±{_stderr(dod_fit):.4f}  des={_mean(des_fit):.4f}±{_stderr(des_fit):.4f}"
            print(row)

        # Mann-Whitney U at gen 30
        dod_div30 = extract_at_gen(results, "dodecahedron", 30, "diversity")
        des_div30 = extract_at_gen(results, "desargues",    30, "diversity")
        dod_fit30 = extract_at_gen(results, "dodecahedron", 30, "meanFitness")
        des_fit30 = extract_at_gen(results, "desargues",    30, "meanFitness")

        print()
        print("  Mann-Whitney U at gen 30:")
        if dod_div30 and des_div30 and len(dod_div30) >= 2 and len(des_div30) >= 2:
            U_div, p_div = stats.mannwhitneyu(dod_div30, des_div30, alternative='two-sided')
            print(f"    Diversity:  U = {U_div:.1f},  p = {p_div:.4f}  (n_dod={len(dod_div30)}, n_des={len(des_div30)})")
        else:
            print("    Diversity:  insufficient data")
        if dod_fit30 and des_fit30 and len(dod_fit30) >= 2 and len(des_fit30) >= 2:
            U_fit, p_fit = stats.mannwhitneyu(dod_fit30, des_fit30, alternative='two-sided')
            print(f"    Fitness:    U = {U_fit:.1f},  p = {p_fit:.4f}  (n_dod={len(dod_fit30)}, n_des={len(des_fit30)})")
        else:
            print("    Fitness:    insufficient data")

        print()
    else:
        print("DODECAHEDRON vs DESARGUES — skipped (data not available for both)")
        print()

    # ------------------------------------------------------------------
    # 8. Correlations with girth and diameter
    # ------------------------------------------------------------------
    print("CORRELATIONS WITH GIRTH AND DIAMETER AT KEY GENERATIONS")
    print(f"{'Gen':>5} {'r(girth,div)':>13} {'r(diam,div)':>13} {'r(girth,fit)':>13} {'r(diam,fit)':>13}")
    print("-" * 65)

    for gen in KEY_GENERATIONS:
        girths, diams, divs, fits = [], [], [], []
        for topo in topos:
            dv = extract_at_gen(results, topo, gen, "diversity")
            ft = extract_at_gen(results, topo, gen, "meanFitness")
            if dv and ft:
                girths.append(TOPOLOGY_META[topo]["girth"])
                diams.append(TOPOLOGY_META[topo]["diameter"])
                divs.append(_mean(dv))
                fits.append(_mean(ft))

        if len(girths) >= 3:
            r_g_div = _pearson_r(girths, divs)
            r_d_div = _pearson_r(diams, divs)
            r_g_fit = _pearson_r(girths, fits)
            r_d_fit = _pearson_r(diams, fits)

            def fmt_r(r):
                return f"{r:>13.4f}" if not math.isnan(r) else f"{'N/A':>13}"

            print(f"{gen:>5} {fmt_r(r_g_div)} {fmt_r(r_d_div)} {fmt_r(r_g_fit)} {fmt_r(r_d_fit)}")
        else:
            print(f"{gen:>5}  (insufficient data)")

    print()

    # ------------------------------------------------------------------
    # 9. Collinearity warning
    # ------------------------------------------------------------------
    print("NOTE: For 3-regular connected graphs, n, beta_1, and density are perfectly")
    print("collinear. Correlations with these variables cannot be distinguished.")
    print("The directed cycle experiment (constant n, m) provides the proper decorrelation.")
    print()

    print("=" * 100)


# ---------------------------------------------------------------------------
# Plotting style
# ---------------------------------------------------------------------------

def setup_style():
    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 11,
        'axes.linewidth': 0.8,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'figure.dpi': 150,
    })


# ---------------------------------------------------------------------------
# Figure 1: Diversity trajectories colored by n
# ---------------------------------------------------------------------------

def fig1_diversity_trajectories(results):
    trajectories = get_metric_trajectories(results, "diversity")
    topos = topos_sorted_by_n(results)

    ns = [TOPOLOGY_META[t]["n"] for t in topos]
    n_min, n_max = min(ns), max(ns)
    norm = Normalize(vmin=n_min, vmax=n_max)
    cmap = plt.cm.viridis

    fig, ax = plt.subplots(figsize=(9, 5.5))

    for topo in topos:
        gens, means, ses = trajectories[topo]
        n = TOPOLOGY_META[topo]["n"]
        color = cmap(norm(n))
        ax.plot(gens, means, color=color, linewidth=1.5, label=f"{topo} (n={n})")
        ax.fill_between(gens, means - ses, means + ses, color=color, alpha=0.15)

    ax.set_xlabel("Generation")
    ax.set_ylabel("Mean Population Diversity")
    ax.set_title("Diversity Trajectories — Foster Census Sweep (13 3-regular topologies)")
    ax.legend(fontsize=7.5, loc='upper right', frameon=False, ncol=2)

    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02, aspect=30)
    cbar.set_label("Vertex count (n)")

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Figure 2: Fitness trajectories colored by n
# ---------------------------------------------------------------------------

def fig2_fitness_trajectories(results):
    trajectories = get_metric_trajectories(results, "meanFitness")
    topos = topos_sorted_by_n(results)

    ns = [TOPOLOGY_META[t]["n"] for t in topos]
    n_min, n_max = min(ns), max(ns)
    norm = Normalize(vmin=n_min, vmax=n_max)
    cmap = plt.cm.viridis

    fig, ax = plt.subplots(figsize=(9, 5.5))

    for topo in topos:
        gens, means, ses = trajectories[topo]
        n = TOPOLOGY_META[topo]["n"]
        color = cmap(norm(n))
        ax.plot(gens, means, color=color, linewidth=1.5, label=f"{topo} (n={n})")
        ax.fill_between(gens, means - ses, means + ses, color=color, alpha=0.15)

    ax.set_xlabel("Generation")
    ax.set_ylabel("Mean Population Fitness")
    ax.set_title("Fitness Trajectories — Foster Census Sweep")
    ax.legend(fontsize=7.5, loc='lower right', frameon=False, ncol=2)

    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02, aspect=30)
    cbar.set_label("Vertex count (n)")

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Figure 3: Scatter — n vs mean diversity at gen 30 with labeled points
# ---------------------------------------------------------------------------

def fig3_scatter_n_vs_diversity(results):
    topos = topos_sorted_by_n(results)

    xs, ys, labels = [], [], []
    for topo in topos:
        vals = extract_at_gen(results, topo, 30, "diversity")
        if vals:
            xs.append(TOPOLOGY_META[topo]["n"])
            ys.append(_mean(vals))
            labels.append(topo)

    xs = np.array(xs)
    ys = np.array(ys)

    fig, ax = plt.subplots(figsize=(8, 5.5))

    ax.scatter(xs, ys, s=70, c='steelblue', edgecolors='black',
               linewidths=0.8, zorder=3)

    if len(xs) >= 3:
        slope, intercept, r, p, _ = stats.linregress(xs, ys)
        x_line = np.linspace(xs.min() - 1, xs.max() + 1, 200)
        ax.plot(x_line, slope * x_line + intercept, 'r--', linewidth=1.2,
                label=f"r = {r:.3f}, p = {p:.3f}")

    for xi, yi, lab in zip(xs, ys, labels):
        ax.annotate(lab, (xi, yi), textcoords="offset points",
                    xytext=(5, 6), fontsize=8, color='0.3')

    ax.set_xlabel("Vertex count (n)")
    ax.set_ylabel("Mean Diversity at Generation 30")
    ax.set_title("Population Size vs Diversity at Gen 30\n(all topologies 3-regular; β₁ = n/2 + 1)")
    if len(xs) >= 3:
        ax.legend(frameon=False, fontsize=10)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Figure 4: Convergence speed bar chart
# ---------------------------------------------------------------------------

def fig4_convergence_speed(results):
    topos = topos_sorted_by_n(results)
    all_gens = get_all_generations(results)

    conv_gens = []
    labels = []
    colors = []

    ns = [TOPOLOGY_META[t]["n"] for t in topos]
    n_min, n_max = (min(ns), max(ns)) if ns else (4, 30)
    norm = Normalize(vmin=n_min, vmax=n_max)
    cmap = plt.cm.viridis

    for topo in topos:
        n = TOPOLOGY_META[topo]["n"]
        conv_gen = None
        for gen in all_gens:
            fits = extract_at_gen(results, topo, gen, "meanFitness")
            if fits and _mean(fits) >= 0.95:
                conv_gen = gen
                break
        conv_gens.append(conv_gen if conv_gen is not None else all_gens[-1] if all_gens else 500)
        labels.append(f"{topo}\n(n={n})")
        colors.append(cmap(norm(n)))

    fig, ax = plt.subplots(figsize=(10, 5))
    x_pos = np.arange(len(labels))
    bars = ax.bar(x_pos, conv_gens, color=colors, edgecolor='black',
                  linewidth=0.6, width=0.7)

    # Mark topologies that did not converge
    all_gens_max = all_gens[-1] if all_gens else 500
    for i, (bar, topo) in enumerate(zip(bars, topos)):
        fits_final = extract_at_gen(results, topo, all_gens_max, "meanFitness")
        if fits_final and _mean(fits_final) < 0.95:
            ax.text(x_pos[i], conv_gens[i] + 5, "*", ha='center',
                    fontsize=12, color='red')

    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Generation of First Mean Fitness ≥ 0.95")
    ax.set_title("Convergence Speed by Topology\n(* = not reached by final generation)")

    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02, aspect=30)
    cbar.set_label("Vertex count (n)")

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Figure 5: Two-panel scatter — diversity vs n AND diversity vs density
# ---------------------------------------------------------------------------

def fig5_diversity_vs_n_and_density(results):
    topos = topos_sorted_by_n(results)

    ns, densities, divs, labels = [], [], [], []
    for topo in topos:
        vals = extract_at_gen(results, topo, 30, "diversity")
        if vals:
            ns.append(TOPOLOGY_META[topo]["n"])
            densities.append(TOPOLOGY_META[topo]["density"])
            divs.append(_mean(vals))
            labels.append(topo)

    ns = np.array(ns)
    densities = np.array(densities)
    divs = np.array(divs)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    # Left panel: diversity vs n
    ax1.scatter(ns, divs, s=70, c='steelblue', edgecolors='black',
                linewidths=0.8, zorder=3)
    if len(ns) >= 3:
        slope, intercept, r, p, _ = stats.linregress(ns, divs)
        x_line = np.linspace(ns.min() - 1, ns.max() + 1, 200)
        ax1.plot(x_line, slope * x_line + intercept, 'r--', linewidth=1.2,
                 label=f"r = {r:.3f}, p = {p:.3f}")
        ax1.legend(frameon=False, fontsize=10)
    for xi, yi, lab in zip(ns, divs, labels):
        ax1.annotate(lab, (xi, yi), textcoords="offset points",
                     xytext=(5, 6), fontsize=7.5, color='0.3')
    ax1.set_xlabel("Vertex count (n)")
    ax1.set_ylabel("Mean Diversity at Generation 30")
    ax1.set_title("Diversity increases with n (= β₁)")

    # Right panel: diversity vs density
    ax2.scatter(densities, divs, s=70, c='darkorange', edgecolors='black',
                linewidths=0.8, zorder=3)
    if len(densities) >= 3:
        slope, intercept, r, p, _ = stats.linregress(densities, divs)
        x_line = np.linspace(densities.min() - 0.02, densities.max() + 0.02, 200)
        ax2.plot(x_line, slope * x_line + intercept, 'r--', linewidth=1.2,
                 label=f"r = {r:.3f}, p = {p:.3f}")
        ax2.legend(frameon=False, fontsize=10)
    for xi, yi, lab in zip(densities, divs, labels):
        ax2.annotate(lab, (xi, yi), textcoords="offset points",
                     xytext=(5, 6), fontsize=7.5, color='0.3')
    ax2.set_xlabel("Density = 3/(n−1)")
    ax2.set_ylabel("Mean Diversity at Generation 30")
    ax2.set_title("Diversity decreases with density")

    fig.suptitle(
        "Same data, opposite narrative — the confound in 3-regular graphs",
        fontsize=12, fontweight='bold', y=1.02)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    setup_style()
    os.makedirs(FIG_DIR, exist_ok=True)

    print(f"Loading data from {RESULTS_DIR}/...")
    results = load_all_data()

    topos = topos_sorted_by_n(results)
    total_runs = sum(len(results[t]) for t in topos)
    print(f"  Loaded {total_runs} runs across {len(topos)} topologies")
    if len(topos) < len(TOPOLOGY_META):
        missing = sorted(set(TOPOLOGY_META) - set(topos))
        print(f"  WARNING: {len(missing)} topologies not yet in data: {missing}")
    print()

    # Text analysis
    run_analysis(results)

    # Figures
    figures = [
        ("fig1_diversity_trajectories", fig1_diversity_trajectories),
        ("fig2_fitness_trajectories",   fig2_fitness_trajectories),
        ("fig3_scatter_n_vs_diversity", fig3_scatter_n_vs_diversity),
        ("fig4_convergence_speed",      fig4_convergence_speed),
        ("fig5_diversity_vs_n_and_density", fig5_diversity_vs_n_and_density),
    ]

    for name, func in figures:
        print(f"Generating {name}...")
        fig = func(results)
        png_path = os.path.join(FIG_DIR, f"{name}.png")
        pdf_path = os.path.join(FIG_DIR, f"{name}.pdf")
        fig.savefig(png_path, dpi=300, bbox_inches='tight')
        fig.savefig(pdf_path, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved {png_path}")
        print(f"  Saved {pdf_path}")

    print(f"\nAll figures saved to {FIG_DIR}/")


if __name__ == "__main__":
    main()
