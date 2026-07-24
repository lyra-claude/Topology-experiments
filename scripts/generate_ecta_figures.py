#!/usr/bin/env python3
"""
Generate all 6 ECTA 2026 figures from REAL experimental data.

Reads CSV data extracted from git branches:
  - Directed cycle experiment (240 runs): /tmp/ecta_data/directed/
  - Foster sweep (390 runs): /tmp/ecta_data/foster/
  - NK pilot (90 runs): /tmp/ecta_data/nk/

Outputs PDF + PNG to: /home/lyra/projects/Topology-experiments/paper/ecta2026/figures/

Author: Lyra (code agent)
Date: 2026-04-07
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# ============================================================
# Publication settings for LNCS
# ============================================================
plt.rcParams.update({
    'font.size': 10,
    'font.family': 'serif',
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'legend.fontsize': 8.5,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'lines.linewidth': 1.5,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
    'axes.grid': False,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# ============================================================
# Configuration
# ============================================================
DIRECTED_DIR = '/tmp/ecta_data/directed'
FOSTER_DIR = '/tmp/ecta_data/foster'
NK_DIR = '/tmp/ecta_data/nk'
OUTDIR = '/home/lyra/projects/Topology-experiments/paper/ecta2026/figures'
os.makedirs(OUTDIR, exist_ok=True)

# Topology families — from Table 2 of the paper
# (name, kappa, color, linestyle)
FAMILIES = [
    ('dag-layer',         0, '#1f77b4', '-'),
    ('dag-wide',          0, '#2ca02c', '--'),
    ('lowcyc-1',          3, '#ff7f0e', '-'),
    ('bidir-ring',       10, '#d62728', '-'),
    ('two-cliques',      14, '#9467bd', '--'),
    ('mesh-cyclic',      20, '#8c564b', '-'),
    ('dense-triangles',  29, '#e377c2', '--'),
    ('ring-skip2',       47, '#7f7f7f', '-'),
]

NAMES = [f[0] for f in FAMILIES]
KAPPAS = np.array([f[1] for f in FAMILIES], dtype=float)
COLORS = {f[0]: f[2] for f in FAMILIES}
STYLES = {f[0]: f[3] for f in FAMILIES}

# Foster census graphs: name -> n (number of vertices)
FOSTER_GRAPHS = {
    'k4': 4,
    'k33': 6,
    'cube': 8,
    'petersen': 10,
    'heawood': 14,
    'mobius-kantor': 16,
    'pappus': 18,
    'dodecahedron': 20,
    'desargues': 20,
    'f26a': 26,    # also known as F026A
    'nauru': 24,
    'coxeter': 28,
    'tutte-coxeter': 30,
}

# NK pilot topologies and K values
NK_TOPOS = ['dag-layer', 'bidir-ring', 'ring-skip2']
NK_K_VALUES = [0, 4, 6]


# ============================================================
# Data loading
# ============================================================
def load_csvs(directory, pattern):
    """Load all CSVs matching pattern from directory.
    Returns dict: key -> list of DataFrames (one per seed).
    Key is extracted from filename before _seed.
    Handles files with 'Up to date' prefix or build log contamination.
    """
    data = {}
    files = sorted(glob.glob(os.path.join(directory, pattern)))
    for f in files:
        basename = os.path.basename(f)
        # Parse: topology_seedNNN.csv or nkK_topology_NNN.csv
        parts = basename.replace('.csv', '')
        key = parts.rsplit('_seed', 1)[0] if '_seed' in parts else parts.rsplit('_', 1)[0]
        try:
            # Read raw lines to find the CSV header
            with open(f, 'r') as fh:
                lines = fh.readlines()
            # Find the line that starts with 'generation,'
            header_idx = None
            for i, line in enumerate(lines):
                if line.startswith('generation,'):
                    header_idx = i
                    break
            if header_idx is None:
                print(f"  Warning: no CSV header in {basename}, skipping")
                continue
            # Parse from the header line onward
            from io import StringIO
            csv_text = ''.join(lines[header_idx:])
            df = pd.read_csv(StringIO(csv_text))
            if 'generation' not in df.columns:
                print(f"  Warning: no generation column in {basename}, skipping")
                continue
            if key not in data:
                data[key] = []
            data[key].append(df)
        except Exception as e:
            print(f"  Warning: failed to read {basename}: {e}")
    return data


def load_directed_data():
    """Load directed cycle experiment data."""
    print("Loading directed cycle data...")
    data = load_csvs(DIRECTED_DIR, '*.csv')
    result = {}
    for name in NAMES:
        if name in data:
            dfs = data[name]
            # Stack into arrays: (n_seeds, n_gens) for each metric
            gens = dfs[0]['generation'].values
            diversity = np.array([df['diversity'].values for df in dfs])
            mean_fit = np.array([df['meanFitness'].values for df in dfs])
            best_fit = np.array([df['bestFitness'].values for df in dfs])
            result[name] = {
                'generations': gens,
                'diversity': diversity,      # shape (n_seeds, n_gens)
                'meanFitness': mean_fit,
                'bestFitness': best_fit,
                'n_seeds': len(dfs),
            }
            print(f"  {name}: {len(dfs)} seeds, {len(gens)} generations")
        else:
            print(f"  WARNING: no data for {name}")
    return result


def load_foster_data():
    """Load Foster sweep data."""
    print("Loading Foster sweep data...")
    data = load_csvs(FOSTER_DIR, '*.csv')
    result = {}
    for name in FOSTER_GRAPHS:
        if name in data:
            dfs = data[name]
            gens = dfs[0]['generation'].values
            diversity = np.array([df['diversity'].values for df in dfs])
            mean_fit = np.array([df['meanFitness'].values for df in dfs])
            result[name] = {
                'generations': gens,
                'diversity': diversity,
                'meanFitness': mean_fit,
                'n_seeds': len(dfs),
                'n': FOSTER_GRAPHS[name],
            }
            print(f"  {name} (n={FOSTER_GRAPHS[name]}): {len(dfs)} seeds")
        else:
            print(f"  WARNING: no data for {name}")
    return result


def load_nk_data():
    """Load NK pilot data."""
    print("Loading NK pilot data...")
    data = load_csvs(NK_DIR, '*.csv')
    result = {}
    for k in NK_K_VALUES:
        for topo in NK_TOPOS:
            key = f'nk{k}_{topo}'
            if key in data:
                dfs = data[key]
                gens = dfs[0]['generation'].values
                diversity = np.array([df['diversity'].values for df in dfs])
                mean_fit = np.array([df['meanFitness'].values for df in dfs])
                best_fit = np.array([df['bestFitness'].values for df in dfs])
                result[key] = {
                    'generations': gens,
                    'diversity': diversity,
                    'meanFitness': mean_fit,
                    'bestFitness': best_fit,
                    'n_seeds': len(dfs),
                    'K': k,
                    'topo': topo,
                }
                print(f"  K={k}, {topo}: {len(dfs)} seeds")
            else:
                print(f"  WARNING: no data for {key}")
    return result


# ============================================================
# Statistical helpers
# ============================================================
def compute_eta2_at_gen(data_dict, gen_idx, names):
    """One-way ANOVA eta-squared at a given generation index."""
    groups = []
    for n in names:
        if n in data_dict:
            groups.append(data_dict[n]['diversity'][:, gen_idx])
    if len(groups) < 2:
        return 0.0
    all_vals = np.concatenate(groups)
    grand_mean = all_vals.mean()
    ss_between = sum(len(g) * (g.mean() - grand_mean)**2 for g in groups)
    ss_total = np.sum((all_vals - grand_mean)**2)
    return ss_between / ss_total if ss_total > 1e-15 else 0.0


def compute_eta2_at_gen_metric(data_dict, gen_idx, names, metric='diversity'):
    """One-way ANOVA eta-squared at a given generation index for any metric."""
    groups = []
    for n in names:
        if n in data_dict:
            groups.append(data_dict[n][metric][:, gen_idx])
    if len(groups) < 2:
        return 0.0
    all_vals = np.concatenate(groups)
    grand_mean = all_vals.mean()
    ss_between = sum(len(g) * (g.mean() - grand_mean)**2 for g in groups)
    ss_total = np.sum((all_vals - grand_mean)**2)
    return ss_between / ss_total if ss_total > 1e-15 else 0.0


def compute_r_at_gen(data_dict, gen_idx, names, kappas):
    """Pearson r between kappa and mean diversity at a given generation."""
    means = []
    valid_kappas = []
    for i, n in enumerate(names):
        if n in data_dict:
            means.append(data_dict[n]['diversity'][:, gen_idx].mean())
            valid_kappas.append(kappas[i])
    means = np.array(means)
    valid_kappas = np.array(valid_kappas)
    if means.std() < 1e-12:
        return 0.0
    return np.corrcoef(valid_kappas, means)[0, 1]


def find_gen_index(generations, target_gen):
    """Find the index of the closest generation to target_gen."""
    return np.argmin(np.abs(generations - target_gen))


def save_fig(fig, name):
    """Save figure as both PDF and PNG."""
    fig.savefig(os.path.join(OUTDIR, f'{name}.pdf'))
    fig.savefig(os.path.join(OUTDIR, f'{name}.png'))
    plt.close(fig)
    print(f"  Saved {name}.pdf and {name}.png")


# ============================================================
# Load all data
# ============================================================
directed = load_directed_data()
foster = load_foster_data()
nk = load_nk_data()

# Get generation arrays
gens = directed[NAMES[0]]['generations']
n_seeds = directed[NAMES[0]]['n_seeds']

print(f"\nData summary:")
print(f"  Directed: {sum(d['n_seeds'] for d in directed.values())} total runs")
print(f"  Foster: {sum(d['n_seeds'] for d in foster.values())} total runs")
print(f"  NK: {sum(d['n_seeds'] for d in nk.values())} total runs")

# ============================================================
# Compute key statistics
# ============================================================
gen30_idx = find_gen_index(gens, 30)
gen50_idx = find_gen_index(gens, 50)
gen100_idx = find_gen_index(gens, 100)

div_means_30 = np.array([directed[n]['diversity'][:, gen30_idx].mean() for n in NAMES])
r_30 = np.corrcoef(KAPPAS, div_means_30)[0, 1]
eta2_30 = compute_eta2_at_gen(directed, gen30_idx, NAMES)

print(f"\nKey statistics at gen 30 (idx={gen30_idx}, actual gen={gens[gen30_idx]}):")
print(f"  r(kappa, diversity) = {r_30:.3f}")
print(f"  eta2(diversity) = {eta2_30:.3f}")
for n in NAMES:
    k = KAPPAS[NAMES.index(n)]
    m = directed[n]['diversity'][:, gen30_idx].mean()
    print(f"    {n} (kappa={int(k)}): mean_div={m:.4f}")

# ============================================================
# FIGURE 1: Diversity trajectories for all 8 topologies
# ============================================================
print("\nGenerating Fig 1: Diversity trajectories...")
fig, ax = plt.subplots(figsize=(6.5, 3.8))

for name, kappa, color, style in FAMILIES:
    d = directed[name]
    mean = d['diversity'].mean(axis=0)
    se = d['diversity'].std(axis=0) / np.sqrt(d['n_seeds'])
    label = f'{name} ($\\kappa$={kappa})'
    ax.plot(gens, mean, color=color, linestyle=style, label=label, zorder=3)
    ax.fill_between(gens, mean - se, mean + se, alpha=0.12, color=color, zorder=2)

ax.set_xlabel('Generation')
ax.set_ylabel('Mean pairwise Hamming diversity')
ax.set_xlim(0, 200)
# Auto y-limits based on data
ymax = max(directed[n]['diversity'].mean(axis=0).max() for n in NAMES) * 1.05
ax.set_ylim(0, min(ymax, 0.55))
ax.legend(loc='upper right', ncol=2, framealpha=0.9, edgecolor='none', fontsize=7.5)
ax.axvline(x=30, color='gray', linestyle=':', linewidth=0.8, alpha=0.5)
ax.text(32, ax.get_ylim()[1] * 0.96, 'gen 30', fontsize=7, color='gray')

fig.tight_layout()
save_fig(fig, 'fig1_diversity_trajectories')

# ============================================================
# FIGURE 2: Correlation r(kappa, diversity) over time
# ============================================================
print("\nGenerating Fig 2: Correlation over time...")

# Compute r at each generation
gen_range = gens[gens >= 5]
gen_indices = [find_gen_index(gens, g) for g in gen_range]
r_over_time = np.array([compute_r_at_gen(directed, idx, NAMES, KAPPAS) for idx in gen_indices])

fig, ax = plt.subplots(figsize=(5.0, 3.2))

ax.plot(gen_range, r_over_time, color='#1f77b4', linewidth=2)
ax.axhline(y=0, color='gray', linewidth=0.5, linestyle='-')

# Reference markers
ax.axvline(x=30, color='gray', linestyle=':', linewidth=0.8, alpha=0.5)
ax.text(32, 0.25, 'gen 30', fontsize=7, color='gray')

# Annotate key values
r_at_30 = compute_r_at_gen(directed, gen30_idx, NAMES, KAPPAS)
ax.plot(30, r_at_30, 'o', color='#d62728', markersize=5, zorder=5)
ax.annotate(f'$r = {r_at_30:.3f}$', xy=(30, r_at_30),
            xytext=(50, r_at_30 - 0.08), fontsize=8, color='#d62728',
            arrowprops=dict(arrowstyle='->', color='#d62728', lw=0.8))

ax.set_xlabel('Generation')
ax.set_ylabel('Pearson $r(\\kappa, \\mathrm{diversity})$')
ax.set_xlim(5, 200)
ax.set_ylim(-1.0, 0.5)

fig.tight_layout()
save_fig(fig, 'fig2_correlation_over_time')

# ============================================================
# FIGURE 3: Scatter kappa vs diversity at gen 30
# ============================================================
print("\nGenerating Fig 3: Scatter kappa vs diversity at gen 30...")
fig, ax = plt.subplots(figsize=(4.5, 3.5))

se_at_30 = np.array([directed[n]['diversity'][:, gen30_idx].std() / np.sqrt(directed[n]['n_seeds']) for n in NAMES])

for idx, (name, kappa, color, style) in enumerate(FAMILIES):
    ax.errorbar(kappa, div_means_30[idx], yerr=se_at_30[idx],
                fmt='o', color=color, markersize=6, capsize=3, capthick=1,
                markeredgecolor='white', markeredgewidth=0.5, zorder=4)

# Regression line
slope, intercept = np.polyfit(KAPPAS, div_means_30, 1)
x_line = np.linspace(-2, 52, 100)
ax.plot(x_line, slope * x_line + intercept, color='gray', linewidth=1,
        linestyle='--', zorder=2)

# Labels with manual offsets to avoid overlap
label_offsets = {
    'dag-layer': (4, 6),
    'dag-wide': (-15, -12),
    'lowcyc-1': (4, 5),
    'bidir-ring': (4, -10),
    'two-cliques': (-28, 6),
    'mesh-cyclic': (4, 5),
    'dense-triangles': (-22, -12),
    'ring-skip2': (4, 5),
}
for idx, (name, kappa, _, _) in enumerate(FAMILIES):
    ax.annotate(name, (kappa, div_means_30[idx]),
                xytext=label_offsets[name], textcoords='offset points',
                fontsize=6.5, color='#444444')

# Compute p-value for the correlation
r_val = np.corrcoef(KAPPAS, div_means_30)[0, 1]
n_points = len(KAPPAS)
t_stat = r_val * np.sqrt((n_points - 2) / (1 - r_val**2))
p_val = 2 * stats.t.sf(abs(t_stat), df=n_points - 2)

ax.text(0.97, 0.95, f'$r = {r_val:.3f}$\n$p = {p_val:.3f}$\n$n = {n_points}$',
        transform=ax.transAxes, ha='right', va='top', fontsize=9,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                  edgecolor='gray', alpha=0.8))

ax.set_xlabel('Directed simple cycle count $\\kappa(D)$')
ax.set_ylabel('Mean diversity at gen. 30')
ax.set_xlim(-3, 52)

fig.tight_layout()
save_fig(fig, 'fig3_scatter_gen30')
print(f"  r = {r_val:.3f}, p = {p_val:.3f}")

# ============================================================
# FIGURE 4: Effect size (eta-squared) over time
# ============================================================
print("\nGenerating Fig 4: Effect size over time...")

eta2_div_over_time = np.array([compute_eta2_at_gen_metric(directed, idx, NAMES, 'diversity') for idx in gen_indices])
eta2_fit_over_time = np.array([compute_eta2_at_gen_metric(directed, idx, NAMES, 'meanFitness') for idx in gen_indices])

fig, ax = plt.subplots(figsize=(5.0, 3.2))

ax.plot(gen_range, eta2_div_over_time, color='#1f77b4', linewidth=2,
        label='Diversity $\\eta^2$')
ax.plot(gen_range, eta2_fit_over_time, color='#d62728', linewidth=2,
        label='Fitness $\\eta^2$', linestyle='--')

# Cohen's thresholds
for thresh, lbl in [(0.01, 'small'), (0.06, 'medium'), (0.14, 'large')]:
    ax.axhline(y=thresh, color='gray', linewidth=0.5, linestyle=':', alpha=0.4)
    ax.text(198, thresh + 0.004, lbl, fontsize=6.5, color='gray', ha='right', alpha=0.6)

ax.axvline(x=30, color='gray', linestyle=':', linewidth=0.8, alpha=0.5)

ax.set_xlabel('Generation')
ax.set_ylabel('$\\eta^2$ (effect size)')
ax.set_xlim(5, 200)
ymax = max(eta2_div_over_time.max(), eta2_fit_over_time.max()) * 1.15
ax.set_ylim(0, min(ymax, 0.5))
ax.legend(loc='upper right', framealpha=0.9, edgecolor='none')

fig.tight_layout()
save_fig(fig, 'fig4_effect_size')

# Print key temporal stats
print(f"\n  Temporal stats:")
print(f"  Gen 30: eta2_div={eta2_div_over_time[gen_range == 30][0] if 30 in gen_range else 'N/A':.3f}, eta2_fit={eta2_fit_over_time[gen_range == 30][0] if 30 in gen_range else 'N/A':.3f}")

# ============================================================
# FIGURE 5: Foster census — same data, opposite narrative
# ============================================================
print("\nGenerating Fig 5: Foster census confound illustration...")

# Use generation 100 diversity means for each Foster graph
# Pick closest gen index to 100
foster_gens = list(foster.values())[0]['generations']
foster_gen100_idx = find_gen_index(foster_gens, 100)
actual_gen = foster_gens[foster_gen100_idx]
print(f"  Using generation {actual_gen} (closest to 100)")

# Collect data: for each graph, get mean diversity and SE at gen 100
foster_names_sorted = sorted(FOSTER_GRAPHS.keys(), key=lambda x: FOSTER_GRAPHS[x])
foster_ns = np.array([FOSTER_GRAPHS[n] for n in foster_names_sorted])
foster_div_means = np.array([foster[n]['diversity'][:, foster_gen100_idx].mean() for n in foster_names_sorted])
foster_div_se = np.array([foster[n]['diversity'][:, foster_gen100_idx].std() / np.sqrt(foster[n]['n_seeds']) for n in foster_names_sorted])

# Density for 3-regular graphs: m = 3n/2, density = m / C(n,2) = (3n/2) / (n(n-1)/2) = 3/(n-1)
foster_density = 3.0 / (foster_ns - 1)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.5, 3.2), sharey=True)

# Left: diversity vs n (should show increasing trend)
ax1.errorbar(foster_ns, foster_div_means, yerr=foster_div_se,
             fmt='o', color='#1f77b4', markersize=5, capsize=2.5, capthick=0.8,
             markeredgecolor='white', markeredgewidth=0.4)
s1, i1 = np.polyfit(foster_ns, foster_div_means, 1)
xl = np.linspace(foster_ns.min() - 1, foster_ns.max() + 1, 100)
ax1.plot(xl, s1 * xl + i1, color='#1f77b4', linewidth=1, linestyle='--', alpha=0.7)
r1 = np.corrcoef(foster_ns, foster_div_means)[0, 1]
ax1.text(0.05, 0.95, f'$r = {r1:.2f}$', transform=ax1.transAxes, ha='left', va='top',
         fontsize=9, bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                               edgecolor='gray', alpha=0.8))
ax1.set_xlabel('Number of islands $n$')
ax1.set_ylabel(f'Mean diversity (gen. {actual_gen})')
ax1.set_title('"$\\beta_1$ drives diversity"', fontsize=10, style='italic')

# Label the n=20 pair (dodecahedron and desargues)
dodec_idx = foster_names_sorted.index('dodecahedron')
desar_idx = foster_names_sorted.index('desargues')
ax1.annotate('Dodec.', (foster_ns[dodec_idx], foster_div_means[dodec_idx]),
             xytext=(-28, 7), textcoords='offset points', fontsize=6, color='#666')
ax1.annotate('Desar.', (foster_ns[desar_idx], foster_div_means[desar_idx]),
             xytext=(4, -10), textcoords='offset points', fontsize=6, color='#666')

# Mann-Whitney U test for dodecahedron vs desargues
dodec_div = foster['dodecahedron']['diversity'][:, foster_gen100_idx]
desar_div = foster['desargues']['diversity'][:, foster_gen100_idx]
u_stat, p_mann = stats.mannwhitneyu(dodec_div, desar_div, alternative='two-sided')
print(f"  Dodecahedron vs Desargues Mann-Whitney p = {p_mann:.3f}")

# Right: diversity vs density (should show decreasing trend — SAME data!)
ax2.errorbar(foster_density, foster_div_means, yerr=foster_div_se,
             fmt='s', color='#d62728', markersize=5, capsize=2.5, capthick=0.8,
             markeredgecolor='white', markeredgewidth=0.4)
s2, i2 = np.polyfit(foster_density, foster_div_means, 1)
xl2 = np.linspace(foster_density.min() - 0.02, foster_density.max() + 0.02, 100)
ax2.plot(xl2, s2 * xl2 + i2, color='#d62728', linewidth=1, linestyle='--', alpha=0.7)
r2 = np.corrcoef(foster_density, foster_div_means)[0, 1]
ax2.text(0.95, 0.95, f'$r = {r2:.2f}$', transform=ax2.transAxes, ha='right', va='top',
         fontsize=9, bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                               edgecolor='gray', alpha=0.8))
ax2.set_xlabel('Density $m / \\binom{n}{2}$')
ax2.set_title('"Density drives diversity"', fontsize=10, style='italic')

fig.tight_layout(w_pad=2.0)
save_fig(fig, 'fig5_foster_same_data')
print(f"  r(n, div) = {r1:.3f}, r(density, div) = {r2:.3f}")

# Compute eta2 within the 3-regular family (should be low — confound!)
foster_eta2 = compute_eta2_at_gen(
    {n: {'diversity': foster[n]['diversity']} for n in foster_names_sorted},
    foster_gen100_idx,
    foster_names_sorted
)
print(f"  Foster overall eta2 at gen {actual_gen}: {foster_eta2:.3f}")

# ============================================================
# FIGURE 6: NK pilot — eta2 by K value
# ============================================================
print("\nGenerating Fig 6: NK landscape effect size...")

# Find generation 200 index in NK data
nk_key = list(nk.keys())[0]
nk_gens = nk[nk_key]['generations']
nk_gen200_idx = find_gen_index(nk_gens, 200)
actual_nk_gen = nk_gens[nk_gen200_idx]
print(f"  Using generation {actual_nk_gen} (closest to 200)")

# Compute eta2 for each K value across the 3 topologies
nk_eta2_div = []
nk_eta2_fit = []
nk_eta2_best = []

for k in NK_K_VALUES:
    topo_keys = [f'nk{k}_{t}' for t in NK_TOPOS]

    # Diversity eta2
    groups_div = [nk[key]['diversity'][:, nk_gen200_idx] for key in topo_keys if key in nk]
    all_div = np.concatenate(groups_div)
    gm = all_div.mean()
    ssb = sum(len(g) * (g.mean() - gm)**2 for g in groups_div)
    sst = np.sum((all_div - gm)**2)
    nk_eta2_div.append(ssb / sst if sst > 1e-15 else 0)

    # Mean fitness eta2
    groups_fit = [nk[key]['meanFitness'][:, nk_gen200_idx] for key in topo_keys if key in nk]
    all_fit = np.concatenate(groups_fit)
    gm = all_fit.mean()
    ssb = sum(len(g) * (g.mean() - gm)**2 for g in groups_fit)
    sst = np.sum((all_fit - gm)**2)
    nk_eta2_fit.append(ssb / sst if sst > 1e-15 else 0)

    # Best fitness eta2
    groups_best = [nk[key]['bestFitness'][:, nk_gen200_idx] for key in topo_keys if key in nk]
    all_best = np.concatenate(groups_best)
    gm = all_best.mean()
    ssb = sum(len(g) * (g.mean() - gm)**2 for g in groups_best)
    sst = np.sum((all_best - gm)**2)
    nk_eta2_best.append(ssb / sst if sst > 1e-15 else 0)

    # ANOVA p-value for diversity
    if len(groups_div) >= 2:
        f_stat, p_val = stats.f_oneway(*groups_div)
        print(f"  K={k}: eta2_div={nk_eta2_div[-1]:.3f}, eta2_fit={nk_eta2_fit[-1]:.3f}, eta2_best={nk_eta2_best[-1]:.3f}, ANOVA p={p_val:.4f}")

fig, ax = plt.subplots(figsize=(4.5, 3.5))

x = np.arange(len(NK_K_VALUES))
width = 0.25

ax.bar(x - width, nk_eta2_div, width, label='Diversity', color='#1f77b4',
       edgecolor='white', linewidth=0.5)
ax.bar(x, nk_eta2_fit, width, label='Mean fitness', color='#d62728',
       edgecolor='white', linewidth=0.5)
ax.bar(x + width, nk_eta2_best, width, label='Best fitness', color='#2ca02c',
       edgecolor='white', linewidth=0.5)

# Cohen's thresholds
for thresh, lbl in [(0.01, 'small'), (0.06, 'medium'), (0.14, 'large')]:
    ax.axhline(y=thresh, color='gray', linewidth=0.5, linestyle=':', alpha=0.4)
    ax.text(2.42, thresh + 0.008, lbl, fontsize=6.5, color='gray', ha='right', alpha=0.6)

# Significance stars based on ANOVA p-values
for i, k in enumerate(NK_K_VALUES):
    topo_keys = [f'nk{k}_{t}' for t in NK_TOPOS]
    groups_div = [nk[key]['diversity'][:, nk_gen200_idx] for key in topo_keys if key in nk]
    groups_fit = [nk[key]['meanFitness'][:, nk_gen200_idx] for key in topo_keys if key in nk]

    if len(groups_div) >= 2:
        _, p_div = stats.f_oneway(*groups_div)
        _, p_fit = stats.f_oneway(*groups_fit)

        def sig_stars(p):
            if p < 0.001: return '***'
            if p < 0.01: return '**'
            if p < 0.05: return '*'
            return ''

        s_div = sig_stars(p_div)
        s_fit = sig_stars(p_fit)
        if s_div:
            ax.text(i - width, nk_eta2_div[i] + 0.02, s_div, ha='center', va='bottom',
                    fontsize=8, fontweight='bold')
        if s_fit:
            ax.text(i, nk_eta2_fit[i] + 0.02, s_fit, ha='center', va='bottom',
                    fontsize=8, fontweight='bold')

# Amplification annotation
if nk_eta2_div[0] > 0.001:
    amplification = nk_eta2_div[-1] / nk_eta2_div[0]
    ax.annotate(f'${amplification:.0f}\\times$',
                xy=(len(NK_K_VALUES) - 1 - width, nk_eta2_div[-1]),
                xytext=(len(NK_K_VALUES) - 0.7, nk_eta2_div[-1] * 0.85),
                fontsize=9, ha='center',
                arrowprops=dict(arrowstyle='->', color='gray', lw=0.8),
                color='#444444')

ax.set_xlabel('Epistasis $K$')
ax.set_ylabel(f'$\\eta^2$ (effect size) at gen. {actual_nk_gen}')
ax.set_xticks(x)
ax.set_xticklabels([f'$K = {k}$' for k in NK_K_VALUES])
ymax_nk = max(max(nk_eta2_div), max(nk_eta2_fit), max(nk_eta2_best)) * 1.2
ax.set_ylim(0, min(ymax_nk, 1.0))
ax.legend(loc='upper left', framealpha=0.9, edgecolor='none')

fig.tight_layout()
save_fig(fig, 'fig6_nk_effect_size')

# ============================================================
# Final summary
# ============================================================
print("\n" + "=" * 60)
print("ALL 6 FIGURES GENERATED FROM REAL DATA")
print("=" * 60)
print(f"Output directory: {OUTDIR}/")
for i in range(1, 7):
    name = ['fig1_diversity_trajectories', 'fig2_correlation_over_time',
            'fig3_scatter_gen30', 'fig4_effect_size',
            'fig5_foster_same_data', 'fig6_nk_effect_size'][i-1]
    print(f"  {name}.pdf / .png")

print(f"\nKey statistics (from REAL data):")
print(f"  Directed cycle experiment:")
print(f"    r(kappa, diversity) at gen 30: {r_30:.3f}")
print(f"    eta2(diversity) at gen 30: {eta2_30:.3f}")
print(f"    n_seeds per topology: {n_seeds}")
print(f"  Foster census:")
print(f"    r(n, diversity): {r1:.3f}")
print(f"    r(density, diversity): {r2:.3f}")
print(f"    Dodecahedron vs Desargues p = {p_mann:.3f}")
print(f"  NK pilot:")
for i, k in enumerate(NK_K_VALUES):
    print(f"    K={k}: eta2_div={nk_eta2_div[i]:.3f}, eta2_fit={nk_eta2_fit[i]:.3f}, eta2_best={nk_eta2_best[i]:.3f}")
