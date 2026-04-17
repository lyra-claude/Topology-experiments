#!/usr/bin/env python3
"""
Generate publication-quality figures for star-vs-cycle edge-controlled experiments.

360 runs: star_8 vs cycle_7 (7 edges), star_12 vs cycle_11 (11 edges),
NK K=0, K=2, K=4, 30 seeds each.
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# ── Configuration ──────────────────────────────────────────────────────
RESULTS_DIR = '/home/lyra/projects/Topology-experiments/results/star_vs_cycle'
FIG_DIR = os.path.join(RESULTS_DIR, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

# Publication style
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'axes.labelsize': 13,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# Colors
STAR_COLOR = '#2166ac'   # blue
CYCLE_COLOR = '#b2182b'  # red

NK_LEVELS = [0, 2, 4]
# Edge-controlled pairs: star_8/cycle_7 (7 edges), star_12/cycle_11 (11 edges)
PAIRS = [('star8', 'cycle7'), ('star12', 'cycle11')]
PAIR_LABELS = ['7-edge (star-8 vs cycle-7)', '11-edge (star-12 vs cycle-11)']


def load_condition(prefix, nk):
    """Load all seeds for a given condition, return DataFrame with seed column."""
    pattern = os.path.join(RESULTS_DIR, f'{prefix}_nk{nk}_seed*.csv')
    files = sorted(glob.glob(pattern))
    dfs = []
    for f in files:
        seed = int(os.path.basename(f).split('seed')[1].replace('.csv', ''))
        df = pd.read_csv(f)
        df['seed'] = seed
        dfs.append(df)
    if not dfs:
        raise ValueError(f"No files for {prefix}_nk{nk}")
    return pd.concat(dfs, ignore_index=True)


def aggregate(df):
    """Aggregate across seeds: mean and SE of diversity."""
    agg = df.groupby('generation')['diversity'].agg(['mean', 'std', 'count']).reset_index()
    agg['se'] = agg['std'] / np.sqrt(agg['count'])
    return agg


def compute_eta_squared(star_df, cycle_df):
    """Compute eta-squared (proportion of variance explained by topology) per generation."""
    gens = sorted(star_df['generation'].unique())
    results = []
    for g in gens:
        s = star_df[star_df['generation'] == g]['diversity'].values
        c = cycle_df[cycle_df['generation'] == g]['diversity'].values
        combined = np.concatenate([s, c])
        grand_mean = combined.mean()
        ss_between = len(s) * (s.mean() - grand_mean)**2 + len(c) * (c.mean() - grand_mean)**2
        ss_total = np.sum((combined - grand_mean)**2)
        eta2 = ss_between / ss_total if ss_total > 0 else 0
        results.append({'generation': g, 'eta_squared': eta2})
    return pd.DataFrame(results)


def compute_cohens_d(star_df, cycle_df, gen_range=None):
    """Compute Cohen's d for diversity (cycle - star) over a generation range."""
    if gen_range is not None:
        star_vals = star_df[(star_df['generation'] >= gen_range[0]) &
                           (star_df['generation'] <= gen_range[1])].groupby('seed')['diversity'].mean().values
        cycle_vals = cycle_df[(cycle_df['generation'] >= gen_range[0]) &
                             (cycle_df['generation'] <= gen_range[1])].groupby('seed')['diversity'].mean().values
    else:
        star_vals = star_df.groupby('seed')['diversity'].mean().values
        cycle_vals = cycle_df.groupby('seed')['diversity'].mean().values

    n1, n2 = len(star_vals), len(cycle_vals)
    pooled_std = np.sqrt(((n1 - 1) * star_vals.std(ddof=1)**2 +
                          (n2 - 1) * cycle_vals.std(ddof=1)**2) / (n1 + n2 - 2))
    if pooled_std == 0:
        return 0
    return (cycle_vals.mean() - star_vals.mean()) / pooled_std


# ── Load all data ──────────────────────────────────────────────────────
data = {}
for pair_idx, (star_prefix, cycle_prefix) in enumerate(PAIRS):
    for nk in NK_LEVELS:
        data[(star_prefix, nk)] = load_condition(star_prefix, nk)
        data[(cycle_prefix, nk)] = load_condition(cycle_prefix, nk)

print("Data loaded. Conditions:", len(data))
for key, df in data.items():
    seeds = df['seed'].nunique()
    print(f"  {key}: {seeds} seeds, {len(df)} rows")


# ── Figure 1: Diversity over time (3-panel, pooled across edge counts) ──
fig, axes = plt.subplots(1, 3, figsize=(12, 3.5), sharey=True)
nk_titles = ['NK $K$=0 (smooth)', 'NK $K$=2 (moderate)', 'NK $K$=4 (rugged)']

for i, nk in enumerate(NK_LEVELS):
    ax = axes[i]

    # Pool both edge-count pairs
    star_dfs = [data[(sp, nk)] for sp, _ in PAIRS]
    cycle_dfs = [data[(cp, nk)] for _, cp in PAIRS]
    star_all = pd.concat(star_dfs, ignore_index=True)
    cycle_all = pd.concat(cycle_dfs, ignore_index=True)

    star_agg = aggregate(star_all)
    cycle_agg = aggregate(cycle_all)

    ax.plot(star_agg['generation'], star_agg['mean'], color=STAR_COLOR,
            linewidth=1.5, label='Star')
    ax.fill_between(star_agg['generation'],
                    star_agg['mean'] - star_agg['se'],
                    star_agg['mean'] + star_agg['se'],
                    color=STAR_COLOR, alpha=0.2)

    ax.plot(cycle_agg['generation'], cycle_agg['mean'], color=CYCLE_COLOR,
            linewidth=1.5, label='Cycle')
    ax.fill_between(cycle_agg['generation'],
                    cycle_agg['mean'] - cycle_agg['se'],
                    cycle_agg['mean'] + cycle_agg['se'],
                    color=CYCLE_COLOR, alpha=0.2)

    ax.set_title(nk_titles[i])
    ax.set_xlabel('Generation')
    if i == 0:
        ax.set_ylabel('Population diversity')
    ax.legend(frameon=False, loc='upper right')

plt.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'diversity_over_time.{ext}'))
print("Saved: diversity_over_time.png/pdf")
plt.close()


# ── Figure 2: η² over time (3 NK levels on one plot) ──────────────────
fig, ax = plt.subplots(figsize=(6, 3.5))
nk_colors = {0: '#4daf4a', 2: '#ff7f00', 4: '#984ea3'}
nk_styles = {0: '--', 2: '-', 4: '-.'}

for nk in NK_LEVELS:
    star_all = pd.concat([data[(sp, nk)] for sp, _ in PAIRS], ignore_index=True)
    cycle_all = pd.concat([data[(cp, nk)] for _, cp in PAIRS], ignore_index=True)
    eta2_df = compute_eta_squared(star_all, cycle_all)

    ax.plot(eta2_df['generation'], eta2_df['eta_squared'],
            color=nk_colors[nk], linestyle=nk_styles[nk], linewidth=1.8,
            label=f'$K$={nk}')

ax.set_xlabel('Generation')
ax.set_ylabel('$\\eta^2$ (topology effect size)')
ax.set_title('Effect of topology on diversity over time')
ax.legend(frameon=False)
ax.set_ylim(bottom=0)

plt.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'eta_squared_over_time.{ext}'))
print("Saved: eta_squared_over_time.png/pdf")
plt.close()


# ── Figure 3: Cohen's d bar chart ─────────────────────────────────────
fig, ax = plt.subplots(figsize=(5, 3.5))

# Compute Cohen's d for full run and for early/late phases
d_full = []
d_early = []  # gen 0-100
d_late = []   # gen 200-500
for nk in NK_LEVELS:
    star_all = pd.concat([data[(sp, nk)] for sp, _ in PAIRS], ignore_index=True)
    cycle_all = pd.concat([data[(cp, nk)] for _, cp in PAIRS], ignore_index=True)
    d_full.append(compute_cohens_d(star_all, cycle_all))
    d_early.append(compute_cohens_d(star_all, cycle_all, gen_range=(0, 100)))
    d_late.append(compute_cohens_d(star_all, cycle_all, gen_range=(200, 500)))

x = np.arange(len(NK_LEVELS))
width = 0.25

bars_early = ax.bar(x - width, d_early, width, label='Early (gen 0-100)',
                    color='#66c2a5', edgecolor='black', linewidth=0.5)
bars_full = ax.bar(x, d_full, width, label='Full run',
                   color='#fc8d62', edgecolor='black', linewidth=0.5)
bars_late = ax.bar(x + width, d_late, width, label='Late (gen 200-500)',
                   color='#8da0cb', edgecolor='black', linewidth=0.5)

ax.axhline(y=0, color='black', linewidth=0.5)
ax.set_xticks(x)
ax.set_xticklabels([f'$K$={k}' for k in NK_LEVELS])
ax.set_ylabel("Cohen's $d$ (cycle $-$ star)")
ax.set_title('Topology effect by landscape ruggedness')
ax.legend(frameon=False, fontsize=9)

# Annotate direction
ax.text(0, min(d_full[0], d_early[0], d_late[0]) - 0.08, 'star wins',
        ha='center', fontsize=8, fontstyle='italic', color='gray')
ax.text(1, max(d_full[1], d_early[1], d_late[1]) + 0.05, 'cycle wins',
        ha='center', fontsize=8, fontstyle='italic', color='gray')

plt.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'cohens_d_by_nk.{ext}'))
print("Saved: cohens_d_by_nk.png/pdf")
plt.close()


# ── Figure 4: Diversity by edge-count pair (supplementary) ─────────────
fig, axes = plt.subplots(2, 3, figsize=(12, 6), sharey='row')

for pair_idx, (star_prefix, cycle_prefix) in enumerate(PAIRS):
    for col, nk in enumerate(NK_LEVELS):
        ax = axes[pair_idx, col]

        star_agg = aggregate(data[(star_prefix, nk)])
        cycle_agg = aggregate(data[(cycle_prefix, nk)])

        ax.plot(star_agg['generation'], star_agg['mean'], color=STAR_COLOR,
                linewidth=1.5, label='Star')
        ax.fill_between(star_agg['generation'],
                        star_agg['mean'] - star_agg['se'],
                        star_agg['mean'] + star_agg['se'],
                        color=STAR_COLOR, alpha=0.2)

        ax.plot(cycle_agg['generation'], cycle_agg['mean'], color=CYCLE_COLOR,
                linewidth=1.5, label='Cycle')
        ax.fill_between(cycle_agg['generation'],
                        cycle_agg['mean'] - cycle_agg['se'],
                        cycle_agg['mean'] + cycle_agg['se'],
                        color=CYCLE_COLOR, alpha=0.2)

        if pair_idx == 0:
            ax.set_title(nk_titles[col])
        if col == 0:
            ax.set_ylabel(PAIR_LABELS[pair_idx] + '\nDiversity')
        ax.set_xlabel('Generation')
        if pair_idx == 0 and col == 2:
            ax.legend(frameon=False, loc='upper right')

plt.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'diversity_by_pair.{ext}'))
print("Saved: diversity_by_pair.png/pdf")
plt.close()


# ── Print summary statistics ───────────────────────────────────────────
print("\n=== Summary Statistics ===")
for nk in NK_LEVELS:
    star_all = pd.concat([data[(sp, nk)] for sp, _ in PAIRS], ignore_index=True)
    cycle_all = pd.concat([data[(cp, nk)] for _, cp in PAIRS], ignore_index=True)
    d = compute_cohens_d(star_all, cycle_all)
    eta2_df = compute_eta_squared(star_all, cycle_all)
    peak_eta2 = eta2_df['eta_squared'].max()
    peak_gen = eta2_df.loc[eta2_df['eta_squared'].idxmax(), 'generation']
    final_eta2 = eta2_df[eta2_df['generation'] == 500]['eta_squared'].values[0]
    print(f"NK K={nk}: Cohen's d={d:.3f}, peak η²={peak_eta2:.3f} (gen {peak_gen}), final η²={final_eta2:.3f}")

print(f"\nFigures saved to: {FIG_DIR}")
