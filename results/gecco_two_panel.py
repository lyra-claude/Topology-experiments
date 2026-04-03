#!/usr/bin/env python3
"""
GECCO Two-Panel Figure: Cycle Rank Discovery + Effect Size Over Time

Panel A: Cycle rank (beta_1) vs peak mean fitness at gen 30 (OneMax)
Panel B: Eta-squared (effect size) over generations with 95% CI

Publication quality, 300 DPI, two-column paper format.
"""

import os
import csv
import io
import numpy as np
import networkx as nx
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from collections import defaultdict

# ─── Configuration ────────────────────────────────────────────────────────────

RESULTS_DIR = os.path.dirname(os.path.abspath(__file__))
PILOT_DIR = os.path.join(RESULTS_DIR, "pilot_onemax")
TRANSIENT_DIR = os.path.join(RESULTS_DIR, "transient_onemax_10runs")

TOPOLOGIES = [
    "disconnected", "ring", "star", "complete",
    "hypercube", "barbell", "watts-strogatz", "random-regular",
]

PILOT_SEEDS = [42, 137, 2718]
TRANSIENT_SEEDS = [42, 137, 2718, 314, 1618, 271, 577, 997, 1729, 4242]

N_BOOTSTRAP = 2000
RNG_SEED = 42

# ─── Graph Construction (faithful to src/Main.hs) ────────────────────────────

N_ISLANDS = 8

def build_graphs():
    """Build all 8 topologies matching the Haskell code exactly."""
    graphs = {}

    # Disconnected
    graphs['disconnected'] = nx.empty_graph(N_ISLANDS)

    # Ring
    graphs['ring'] = nx.cycle_graph(N_ISLANDS)

    # Star
    graphs['star'] = nx.star_graph(N_ISLANDS - 1)

    # Complete
    graphs['complete'] = nx.complete_graph(N_ISLANDS)

    # Hypercube (dim=3 -> 8 nodes)
    G = nx.hypercube_graph(3)
    mapping = {}
    for node in G.nodes():
        idx = 0
        for bit in node:
            idx = idx * 2 + bit
        mapping[node] = idx
    graphs['hypercube'] = nx.relabel_nodes(G, mapping)

    # Barbell
    half = N_ISLANDS // 2
    G = nx.Graph()
    G.add_nodes_from(range(N_ISLANDS))
    for i in range(half):
        for j in range(i + 1, half):
            G.add_edge(i, j)
    for i in range(half, N_ISLANDS):
        for j in range(i + 1, N_ISLANDS):
            G.add_edge(i, j)
    G.add_edge(half - 1, half)
    graphs['barbell'] = G

    # Watts-Strogatz (deterministic, matching Haskell)
    n, k, p, seed = N_ISLANDS, 4, 0.3, 42
    halfK = k // 2
    adj = [set() for _ in range(n)]
    for i in range(n):
        for d in range(1, halfK + 1):
            j_fwd = (i + d) % n
            j_bwd = (i - d) % n
            adj[i].add(j_fwd)
            adj[i].add(j_bwd)
            adj[j_fwd].add(i)
            adj[j_bwd].add(i)
    for i in range(n):
        for d in range(1, halfK + 1):
            j = (i + d) % n
            h = (i * 7919 + j * 6271 + seed * 1031) % 10000
            if (h / 10000.0) < p:
                newJ = ((i * 3571 + j * 2749 + seed * 947) % (n - 1))
                newJ_prime = newJ + 1 if newJ >= i else newJ
                if newJ_prime == j or newJ_prime in adj[i]:
                    continue
                adj[i].discard(j)
                adj[j].discard(i)
                adj[i].add(newJ_prime)
                adj[newJ_prime].add(i)
    G = nx.Graph()
    G.add_nodes_from(range(n))
    for i in range(n):
        for j in adj[i]:
            if j > i:
                G.add_edge(i, j)
    graphs['watts-strogatz'] = G

    # Random-regular (deterministic, matching Haskell)
    G = nx.cycle_graph(N_ISLANDS)
    d_reg, seed_rr = 3, 42
    for i in range(N_ISLANDS):
        if G.degree(i) >= d_reg:
            continue
        target = ((i * 5021 + seed_rr * 1733) % (N_ISLANDS - 2))
        target_prime = target + 1 if target >= i else target
        if G.has_edge(i, target_prime) or G.degree(target_prime) >= d_reg:
            continue
        G.add_edge(i, target_prime)
    graphs['random-regular'] = G

    return graphs


def compute_cycle_rank(G):
    """Cycle rank = |E| - |V| + c (first Betti number of the graph)."""
    return G.number_of_edges() - G.number_of_nodes() + nx.number_connected_components(G)


# ─── Data Loading ─────────────────────────────────────────────────────────────

def load_csv(filepath):
    """Load a CSV, handling potential cabal output prefix."""
    with open(filepath) as f:
        lines = f.readlines()
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
    """Load all CSVs -> dict[topology][seed_idx] = {gen: {metrics}}"""
    all_data = defaultdict(list)
    for topo in TOPOLOGIES:
        for seed in seeds:
            filepath = os.path.join(results_dir, f"{topo}_seed{seed}.csv")
            if not os.path.exists(filepath):
                print(f"WARNING: Missing {filepath}")
                continue
            all_data[topo].append(load_csv(filepath))
    return all_data


def get_gen30_peak_fitness(data):
    """Extract mean fitness at generation 30 (the peak effect window) for each topology."""
    stats = {}
    for topo in TOPOLOGIES:
        fitnesses = []
        for run in data[topo]:
            if 30 in run:
                fitnesses.append(run[30]['meanFitness'])
        stats[topo] = np.mean(fitnesses) if fitnesses else np.nan
    return stats


# ─── Eta-squared (reused from rsquared_curve.py) ─────────────────────────────

def compute_eta_squared(groups):
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


def bootstrap_eta_squared(data, gen, metric, n_bootstrap=N_BOOTSTRAP, rng=None):
    if rng is None:
        rng = np.random.default_rng(RNG_SEED)
    groups = {}
    for topo in TOPOLOGIES:
        vals = []
        for run in data[topo]:
            if gen in run:
                vals.append(run[gen][metric])
        groups[topo] = np.array(vals)
    if any(len(v) == 0 for v in groups.values()):
        return np.nan, np.nan, np.nan
    observed = compute_eta_squared([groups[t].tolist() for t in TOPOLOGIES])
    boot_vals = np.empty(n_bootstrap)
    for b in range(n_bootstrap):
        boot_groups = []
        for topo in TOPOLOGIES:
            g = groups[topo]
            idx = rng.integers(0, len(g), size=len(g))
            boot_groups.append(g[idx].tolist())
        boot_vals[b] = compute_eta_squared(boot_groups)
    return observed, np.nanpercentile(boot_vals, 2.5), np.nanpercentile(boot_vals, 97.5)


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    print("Building graphs...")
    graphs = build_graphs()

    print("Computing cycle ranks...")
    cycle_ranks = {}
    for name, G in graphs.items():
        cr = compute_cycle_rank(G)
        cycle_ranks[name] = cr
        print(f"  {name:20s}  |E|={G.number_of_edges():2d}  |V|={G.number_of_nodes()}  "
              f"c={nx.number_connected_components(G)}  beta_1={cr}")

    print("\nLoading pilot OneMax data...")
    pilot_data = load_dataset(PILOT_DIR, PILOT_SEEDS)

    print("Loading transient OneMax data...")
    transient_data = load_dataset(TRANSIENT_DIR, TRANSIENT_SEEDS)

    # ─── Panel A data: cycle rank vs peak fitness ─────────────────────────
    peak_fitness = get_gen30_peak_fitness(pilot_data)
    print("\nGen 30 peak fitness by topology:")
    for topo in TOPOLOGIES:
        print(f"  {topo:20s}  beta_1={cycle_ranks[topo]:2d}  fitness={peak_fitness[topo]:.6f}")

    # Correlation: all 8 topologies
    cr_arr_all = np.array([cycle_ranks[t] for t in TOPOLOGIES])
    fit_arr_all = np.array([peak_fitness[t] for t in TOPOLOGIES])
    rho_all, p_val_all = spearmanr(cr_arr_all, fit_arr_all)
    print(f"\nSpearman correlation (all 8): rho={rho_all:.3f}, p={p_val_all:.4f}")

    # Also compute for connected only (excluding disconnected)
    connected_topos = [t for t in TOPOLOGIES if t != 'disconnected']
    cr_arr_conn = np.array([cycle_ranks[t] for t in connected_topos])
    fit_arr_conn = np.array([peak_fitness[t] for t in connected_topos])
    rho_conn, p_val_conn = spearmanr(cr_arr_conn, fit_arr_conn)
    print(f"Spearman correlation (connected only): rho={rho_conn:.3f}, p={p_val_conn:.4f}")

    # Use all-8 correlation for the figure (matches the rho=0.893 from memory)
    rho, p_val = rho_all, p_val_all

    # ─── Panel B data: eta-squared over generations ───────────────────────
    print("\nComputing eta-squared over generations...")
    rng = np.random.default_rng(RNG_SEED)

    transient_gens = list(range(0, 101, 10))
    pilot_gens = list(range(110, 501, 10))

    transient_results = []
    for gen in transient_gens:
        fit_eta, fit_lo, fit_hi = bootstrap_eta_squared(transient_data, gen, 'meanFitness', rng=rng)
        div_eta, div_lo, div_hi = bootstrap_eta_squared(transient_data, gen, 'diversity', rng=rng)
        transient_results.append({
            'gen': gen,
            'fit_eta': fit_eta, 'fit_lo': fit_lo, 'fit_hi': fit_hi,
            'div_eta': div_eta, 'div_lo': div_lo, 'div_hi': div_hi,
        })

    pilot_results = []
    for gen in pilot_gens:
        fit_eta, fit_lo, fit_hi = bootstrap_eta_squared(pilot_data, gen, 'meanFitness', rng=rng)
        div_eta, div_lo, div_hi = bootstrap_eta_squared(pilot_data, gen, 'diversity', rng=rng)
        pilot_results.append({
            'gen': gen,
            'fit_eta': fit_eta, 'fit_lo': fit_lo, 'fit_hi': fit_hi,
            'div_eta': div_eta, 'div_lo': div_lo, 'div_hi': div_hi,
        })

    all_results = transient_results + pilot_results
    fit_peak = max(all_results, key=lambda r: r['fit_eta'] if not np.isnan(r['fit_eta']) else -1)

    # ═══════════════════════════════════════════════════════════════════════
    # FIGURE
    # ═══════════════════════════════════════════════════════════════════════

    print("\nGenerating two-panel GECCO figure...")

    # Color palette
    COLOR_FIT = '#2166ac'       # Blue
    COLOR_DIV = '#b2182b'       # Red
    COLOR_FIT_LIGHT = '#92c5de'
    COLOR_DIV_LIGHT = '#f4a582'
    COLOR_POINT = '#333333'
    COLOR_STAR = '#e66101'      # Orange for anomaly highlights
    COLOR_BARBELL = '#5e3c99'   # Purple for anomaly highlights
    COLOR_GOLDILOCKS = '#fdb863'  # Light orange for Goldilocks zone

    # Figure: two panels side by side, sized for two-column paper
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.2, 3.5), dpi=300)
    plt.rcParams.update({
        'font.size': 8,
        'axes.labelsize': 9,
        'axes.titlesize': 9.5,
        'xtick.labelsize': 7.5,
        'ytick.labelsize': 7.5,
    })

    # ─── PANEL A: Cycle Rank vs Peak Fitness ──────────────────────────────

    for spine in ['top', 'right']:
        ax_a.spines[spine].set_visible(False)
    ax_a.spines['left'].set_linewidth(0.8)
    ax_a.spines['bottom'].set_linewidth(0.8)

    # Display names and label positions (manually tuned to avoid overlap)
    display_names = {
        'disconnected': 'Disconnected',
        'ring': 'Ring',
        'star': 'Star',
        'complete': 'Complete',
        'hypercube': 'Hypercube',
        'barbell': 'Barbell',
        'watts-strogatz': 'Watts-Strogatz',
        'random-regular': 'Random-Regular',
    }

    # Plot all points
    for topo in TOPOLOGIES:
        cr = cycle_ranks[topo]
        fit = peak_fitness[topo]
        is_star = topo == 'star'
        is_barbell = topo == 'barbell'

        if is_star:
            color = COLOR_STAR
            marker = 'D'
            ms = 7
            zorder = 10
        elif is_barbell:
            color = COLOR_BARBELL
            marker = 'D'
            ms = 7
            zorder = 10
        else:
            color = COLOR_POINT
            marker = 'o'
            ms = 5
            zorder = 5

        ax_a.scatter(cr, fit, c=color, marker=marker, s=ms**2, zorder=zorder,
                     edgecolors='white', linewidths=0.5)

    # Label offset tuning to avoid overlaps
    label_offsets = {
        'disconnected': (0.3, -0.012),
        'ring':         (0.3, 0.003),
        'star':         (0.3, 0.003),
        'complete':     (-0.3, 0.005),
        'hypercube':    (0.3, 0.003),
        'barbell':      (0.3, 0.003),
        'watts-strogatz': (-0.3, -0.008),
        'random-regular': (0.3, -0.008),
    }
    label_ha = {
        'complete': 'right',
        'watts-strogatz': 'right',
    }

    for topo in TOPOLOGIES:
        cr = cycle_ranks[topo]
        fit = peak_fitness[topo]
        dx, dy = label_offsets.get(topo, (0.3, 0.0))
        ha = label_ha.get(topo, 'left')
        is_anomaly = topo in ('star', 'barbell')
        color = COLOR_STAR if topo == 'star' else (COLOR_BARBELL if topo == 'barbell' else '#555555')
        weight = 'bold' if is_anomaly else 'normal'
        fontsize = 7 if is_anomaly else 6.5

        ax_a.annotate(
            display_names[topo],
            xy=(cr, fit),
            xytext=(cr + dx, fit + dy),
            fontsize=fontsize,
            color=color,
            fontweight=weight,
            ha=ha, va='center',
        )

    # Correlation annotation
    ax_a.text(0.05, 0.05,
              f'Spearman $\\rho$ = {rho:.3f}\n$p$ = {p_val:.3f}',
              transform=ax_a.transAxes,
              fontsize=7,
              verticalalignment='bottom',
              bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                        edgecolor='#cccccc', alpha=0.9))

    # Resolved anomaly annotations
    # Star: beta_1=0, lowest among connected -> correctly predicts bottleneck
    ax_a.annotate(
        r'$\beta_1\!=\!0$ (tree)',
        xy=(cycle_ranks['star'], peak_fitness['star']),
        xytext=(cycle_ranks['star'] + 1.5, peak_fitness['star'] - 0.025),
        fontsize=6, color=COLOR_STAR, fontstyle='italic',
        arrowprops=dict(arrowstyle='->', color=COLOR_STAR, lw=0.7,
                        connectionstyle='arc3,rad=0.2'),
        ha='left', va='top',
    )

    # Barbell: beta_1=6, highest -> correctly predicts rich cycles
    ax_a.annotate(
        r'$\beta_1\!=\!6$ (rich cycles)',
        xy=(cycle_ranks['barbell'], peak_fitness['barbell']),
        xytext=(cycle_ranks['barbell'] - 1.5, peak_fitness['barbell'] + 0.015),
        fontsize=6, color=COLOR_BARBELL, fontstyle='italic',
        arrowprops=dict(arrowstyle='->', color=COLOR_BARBELL, lw=0.7,
                        connectionstyle='arc3,rad=-0.2'),
        ha='right', va='bottom',
    )

    ax_a.set_xlabel(r'Cycle rank ($\beta_1 = |E| - |V| + c$)', fontsize=9, labelpad=4)
    ax_a.set_ylabel('Peak mean fitness (gen 30)', fontsize=9, labelpad=4)
    ax_a.set_title('(a) Cycle Rank vs GA Performance', fontsize=9.5, fontweight='bold', pad=8)

    # X-axis: range from -0.5 to max+1
    max_cr = max(cycle_ranks.values())
    ax_a.set_xlim(-0.8, max_cr + 1.5)
    ax_a.set_xticks(range(0, max_cr + 1))

    ax_a.tick_params(axis='both', which='major', length=3, width=0.6)

    # ─── PANEL B: Eta-squared Over Generations ────────────────────────────

    for spine in ['top', 'right']:
        ax_b.spines[spine].set_visible(False)
    ax_b.spines['left'].set_linewidth(0.8)
    ax_b.spines['bottom'].set_linewidth(0.8)

    # Transient data (10-run, solid with CI)
    t_gens = np.array([r['gen'] for r in transient_results])
    t_fit = np.array([r['fit_eta'] for r in transient_results])
    t_fit_lo = np.array([r['fit_lo'] for r in transient_results])
    t_fit_hi = np.array([r['fit_hi'] for r in transient_results])
    t_div = np.array([r['div_eta'] for r in transient_results])
    t_div_lo = np.array([r['div_lo'] for r in transient_results])
    t_div_hi = np.array([r['div_hi'] for r in transient_results])

    # Pilot data (3-run, dashed trend)
    p_gens = np.array([r['gen'] for r in pilot_results])
    p_fit = np.array([r['fit_eta'] for r in pilot_results])
    p_div = np.array([r['div_eta'] for r in pilot_results])

    # Goldilocks zone highlight (gen 20-40 where eta-squared peaks)
    ax_b.axvspan(20, 40, alpha=0.12, color=COLOR_GOLDILOCKS, zorder=0)
    ax_b.text(30, 0.96, 'Goldilocks\nzone', fontsize=6, color='#b35900',
              ha='center', va='top', fontstyle='italic', alpha=0.8)

    # CI bands (10-run only)
    ax_b.fill_between(t_gens, t_fit_lo, t_fit_hi,
                      alpha=0.20, color=COLOR_FIT, linewidth=0)
    ax_b.fill_between(t_gens, t_div_lo, t_div_hi,
                      alpha=0.20, color=COLOR_DIV, linewidth=0)

    # Main lines (10-run)
    ax_b.plot(t_gens, t_fit, color=COLOR_FIT, linewidth=1.8, marker='o',
              markersize=3.5, markeredgewidth=0, zorder=5)
    ax_b.plot(t_gens, t_div, color=COLOR_DIV, linewidth=1.8, marker='s',
              markersize=3.5, markeredgewidth=0, zorder=5)

    # 3-run data: small scatter + rolling mean trend
    ax_b.scatter(p_gens, p_fit, color=COLOR_FIT, s=5, alpha=0.25, zorder=3,
                 edgecolors='none')
    ax_b.scatter(p_gens, p_div, color=COLOR_DIV, s=5, alpha=0.25, zorder=3,
                 edgecolors='none', marker='s')

    def rolling_mean(vals, window=5):
        result = np.empty_like(vals)
        for i in range(len(vals)):
            lo = max(0, i - window // 2)
            hi = min(len(vals), i + window // 2 + 1)
            result[i] = np.mean(vals[lo:hi])
        return result

    p_fit_smooth = rolling_mean(p_fit)
    p_div_smooth = rolling_mean(p_div)

    ax_b.plot(np.concatenate([[t_gens[-1]], p_gens]),
              np.concatenate([[t_fit[-1]], p_fit_smooth]),
              color=COLOR_FIT, linewidth=1.0, linestyle='--', alpha=0.5, zorder=4)
    ax_b.plot(np.concatenate([[t_gens[-1]], p_gens]),
              np.concatenate([[t_div[-1]], p_div_smooth]),
              color=COLOR_DIV, linewidth=1.0, linestyle='--', alpha=0.5, zorder=4)

    # Regime boundary
    ax_b.axvline(x=100, color='gray', linewidth=0.6, linestyle=':', alpha=0.5)
    ax_b.text(105, 0.96, '$n$=3', fontsize=6, color='gray', ha='left', va='top')
    ax_b.text(95, 0.96, '$n$=10', fontsize=6, color='gray', ha='right', va='top')

    # Peak annotation
    ax_b.annotate(
        r'$\eta^2 = $' + f'{fit_peak["fit_eta"]:.2f}',
        xy=(fit_peak['gen'], fit_peak['fit_eta']),
        xytext=(fit_peak['gen'] + 50, fit_peak['fit_eta'] + 0.06),
        fontsize=7,
        color=COLOR_FIT,
        fontweight='bold',
        arrowprops=dict(arrowstyle='->', color=COLOR_FIT, lw=0.7,
                        connectionstyle='arc3,rad=0.15'),
        ha='left', va='bottom',
    )

    # Cohen's benchmarks
    for level, label in [(0.14, 'large'), (0.06, 'medium')]:
        ax_b.axhline(y=level, color='gray', linewidth=0.4, linestyle='-', alpha=0.25)
        ax_b.text(505, level, label, fontsize=5.5, color='gray', alpha=0.5,
                  ha='left', va='center')

    # Legend
    legend_elements = [
        Line2D([0], [0], color=COLOR_FIT, linewidth=1.5, marker='o',
               markersize=3.5, markeredgewidth=0, label='Mean fitness'),
        Line2D([0], [0], color=COLOR_DIV, linewidth=1.5, marker='s',
               markersize=3.5, markeredgewidth=0, label='Diversity'),
    ]
    ax_b.legend(handles=legend_elements, loc='upper right', frameon=False,
                fontsize=7, handlelength=2)

    ax_b.set_xlabel('Generation', fontsize=9, labelpad=4)
    ax_b.set_ylabel(r'$\eta^2$ (effect size)', fontsize=9, labelpad=4)
    ax_b.set_title(r'(b) Topology Effect Size Over Time', fontsize=9.5, fontweight='bold', pad=8)

    ax_b.set_xlim(-5, 520)
    ax_b.set_ylim(-0.02, 1.02)
    ax_b.set_xticks([0, 50, 100, 200, 300, 400, 500])
    ax_b.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax_b.tick_params(axis='both', which='major', length=3, width=0.6)

    # ─── Save ─────────────────────────────────────────────────────────────

    fig.subplots_adjust(left=0.08, right=0.95, bottom=0.14, top=0.90, wspace=0.32)

    png_path = os.path.join(RESULTS_DIR, "gecco_two_panel.png")
    pdf_path = os.path.join(RESULTS_DIR, "gecco_two_panel.pdf")

    fig.savefig(png_path, dpi=300, bbox_inches='tight', pad_inches=0.15)
    fig.savefig(pdf_path, bbox_inches='tight', pad_inches=0.15)
    plt.close(fig)

    print(f"\nFigure saved:")
    print(f"  PNG: {png_path}")
    print(f"  PDF: {pdf_path}")
    print("\nDone.")


if __name__ == "__main__":
    main()
