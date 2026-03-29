#!/usr/bin/env python3
"""
Deep transient analysis - per-run correlations, effect sizes, and anomaly investigation.
"""

import csv
import os
import math
from collections import defaultdict

RESULTS_DIR = "results/transient_onemax_10runs"
TOPOLOGIES = ["disconnected", "ring", "star", "complete", "hypercube", "barbell", "watts-strogatz", "random-regular"]
SEEDS = [42, 137, 2718, 314, 1618, 271, 577, 997, 1729, 4242]

LAMBDA2 = {
    "disconnected": 0.0, "barbell": 0.07, "random-regular": 0.27,
    "ring": 0.59, "star": 1.0, "watts-strogatz": 1.5,
    "hypercube": 2.0, "complete": 8.0,
}

def mean(xs): return sum(xs) / len(xs) if xs else 0
def std(xs):
    m = mean(xs)
    return math.sqrt(sum((x-m)**2 for x in xs) / (len(xs)-1)) if len(xs)>1 else 0
def se(xs): return std(xs) / math.sqrt(len(xs)) if xs else 0
def ci95(xs):
    t_val = 2.262  # t_{0.025, 9}
    return t_val * se(xs)

def normal_cdf(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))

def spearman_rank(xs, ys):
    n = len(xs)
    if n < 3: return 0.0, 1.0
    def rank_data(vals):
        indexed = sorted(enumerate(vals), key=lambda x: x[1])
        ranks = [0.0]*n
        i = 0
        while i < n:
            j = i
            while j < n-1 and indexed[j+1][1] == indexed[j][1]: j += 1
            avg_rank = (i+j)/2.0 + 1.0
            for k in range(i, j+1): ranks[indexed[k][0]] = avg_rank
            i = j+1
        return ranks
    rx, ry = rank_data(xs), rank_data(ys)
    mx, my = mean(rx), mean(ry)
    num = sum((rx[i]-mx)*(ry[i]-my) for i in range(n))
    den_x = math.sqrt(sum((rx[i]-mx)**2 for i in range(n)))
    den_y = math.sqrt(sum((ry[i]-my)**2 for i in range(n)))
    if den_x == 0 or den_y == 0: return 0.0, 1.0
    rho = num / (den_x * den_y)
    if abs(rho) >= 1.0: return rho, 0.0
    t_stat = rho * math.sqrt((n-2)/(1-rho**2))
    p_val = 2 * (1 - normal_cdf(abs(t_stat)))
    return rho, p_val

def cohens_d(xs, ys):
    """Compute Cohen's d effect size."""
    nx, ny = len(xs), len(ys)
    mx, my = mean(xs), mean(ys)
    vx = sum((x-mx)**2 for x in xs) / (nx-1) if nx>1 else 0
    vy = sum((y-my)**2 for y in ys) / (ny-1) if ny>1 else 0
    pooled_sd = math.sqrt(((nx-1)*vx + (ny-1)*vy) / (nx+ny-2))
    return (mx - my) / pooled_sd if pooled_sd > 0 else 0

def load_all():
    all_data = defaultdict(list)
    for topo in TOPOLOGIES:
        for seed in SEEDS:
            filepath = os.path.join(RESULTS_DIR, f"{topo}_seed{seed}.csv")
            data = {}
            with open(filepath) as f:
                for row in csv.DictReader(f):
                    gen = int(row['generation'])
                    data[gen] = {
                        'meanFitness': float(row['meanFitness']),
                        'bestFitness': float(row['bestFitness']),
                        'diversity': float(row['diversity']),
                    }
            all_data[topo].append(data)
    return all_data

def main():
    data = load_all()

    print("=" * 80)
    print("  DEEP TRANSIENT ANALYSIS — ANOMALIES & EFFECT SIZES")
    print("=" * 80)

    # ─── 1. Per-run Spearman (using all 80 data points) ─────────────────────
    print("\n" + "=" * 80)
    print("  PER-RUN SPEARMAN: λ₂ vs metrics (80 data points: 8 topos × 10 seeds)")
    print("=" * 80)

    for gen in [20, 30, 40, 50]:
        # Each data point is one run from one topology
        l2_all = []
        fit_all = []
        div_all = []
        for topo in TOPOLOGIES:
            for run in data[topo]:
                l2_all.append(LAMBDA2[topo])
                fit_all.append(run[gen]['meanFitness'])
                div_all.append(run[gen]['diversity'])

        rho_fit, p_fit = spearman_rank(l2_all, fit_all)
        rho_div, p_div = spearman_rank(l2_all, div_all)
        sig_f = "***" if p_fit < 0.001 else "**" if p_fit < 0.01 else "*" if p_fit < 0.05 else "n.s."
        sig_d = "***" if p_div < 0.001 else "**" if p_div < 0.01 else "*" if p_div < 0.05 else "n.s."

        print(f"\n  Gen {gen} (n=80):")
        print(f"    λ₂ vs Fitness:   ρ = {rho_fit:.4f}, p = {p_fit:.8f}  {sig_f}")
        print(f"    λ₂ vs Diversity: ρ = {rho_div:.4f}, p = {p_div:.8f}  {sig_d}")

    # ─── 2. Cohen's d effect sizes ───────────────────────────────────────────
    print("\n" + "=" * 80)
    print("  COHEN'S d EFFECT SIZES: Complete vs Disconnected")
    print("=" * 80)

    for gen in [20, 30, 40, 50]:
        comp_fit = [r[gen]['meanFitness'] for r in data['complete']]
        disc_fit = [r[gen]['meanFitness'] for r in data['disconnected']]
        comp_div = [r[gen]['diversity'] for r in data['complete']]
        disc_div = [r[gen]['diversity'] for r in data['disconnected']]

        d_fit = cohens_d(comp_fit, disc_fit)
        d_div = cohens_d(comp_div, disc_div)

        mag_f = "huge" if abs(d_fit)>1.2 else "very large" if abs(d_fit)>1.0 else "large" if abs(d_fit)>0.8 else "medium" if abs(d_fit)>0.5 else "small"
        mag_d = "huge" if abs(d_div)>1.2 else "very large" if abs(d_div)>1.0 else "large" if abs(d_div)>0.8 else "medium" if abs(d_div)>0.5 else "small"

        print(f"\n  Gen {gen}:")
        print(f"    Fitness:   d = {d_fit:+.3f} ({mag_f})")
        print(f"    Diversity: d = {d_div:+.3f} ({mag_d})")

    # ─── 3. Star anomaly deep dive ──────────────────────────────────────────
    print("\n" + "=" * 80)
    print("  STAR ANOMALY: Detailed Investigation")
    print("=" * 80)

    print("\n  Star (λ₂=1.0) ranks BELOW ring (λ₂=0.59) and random-regular (λ₂=0.27)")
    print("  This is the main source of rank disruption.\n")

    for gen in [20, 30]:
        print(f"  --- Gen {gen} ---")
        star_vals = [r[gen]['meanFitness'] for r in data['star']]
        ring_vals = [r[gen]['meanFitness'] for r in data['ring']]
        rr_vals = [r[gen]['meanFitness'] for r in data['random-regular']]

        # Star vs Ring
        d_sr = cohens_d(star_vals, ring_vals)
        from_stat = "Star < Ring" if mean(star_vals) < mean(ring_vals) else "Star > Ring"
        print(f"    Star vs Ring: Star={mean(star_vals):.4f}, Ring={mean(ring_vals):.4f}, "
              f"d={d_sr:+.3f} ({from_stat})")

        # Star vs Random-Regular
        d_srr = cohens_d(star_vals, rr_vals)
        from_stat2 = "Star < RR" if mean(star_vals) < mean(rr_vals) else "Star > RR"
        print(f"    Star vs RandReg: Star={mean(star_vals):.4f}, RR={mean(rr_vals):.4f}, "
              f"d={d_srr:+.3f} ({from_stat2})")

        # Star diversity: is it preserving too much?
        star_div = [r[gen]['diversity'] for r in data['star']]
        ring_div = [r[gen]['diversity'] for r in data['ring']]
        print(f"    Star diversity={mean(star_div):.4f}, Ring diversity={mean(ring_div):.4f}")

    print("\n  Star topology interpretation:")
    print("  Hub-and-spoke creates an asymmetric migration pattern:")
    print("  - Hub island receives from ALL others (over-mixing at center)")
    print("  - Spoke islands receive ONLY from hub (under-mixing at periphery)")
    print("  - Net effect: diversity preserved in spokes, but convergence slowed")
    print("  - This is a STRUCTURAL anomaly: λ₂ measures algebraic connectivity")
    print("    but doesn't capture the ASYMMETRY of information flow.")

    # ─── 4. Barbell overperformance ──────────────────────────────────────────
    print("\n" + "=" * 80)
    print("  BARBELL OVERPERFORMANCE: Detailed Investigation")
    print("=" * 80)

    print("\n  Barbell (λ₂=0.07) ranks 3rd in fitness at gen 20!")
    print("  Expected rank by λ₂: 7th (near disconnected)\n")

    for gen in [20, 30]:
        barbell_fit = [r[gen]['meanFitness'] for r in data['barbell']]
        disc_fit = [r[gen]['meanFitness'] for r in data['disconnected']]
        hyper_fit = [r[gen]['meanFitness'] for r in data['hypercube']]

        d_bd = cohens_d(barbell_fit, disc_fit)
        d_bh = cohens_d(barbell_fit, hyper_fit)

        print(f"  Gen {gen}:")
        print(f"    Barbell={mean(barbell_fit):.4f} vs Disconnected={mean(disc_fit):.4f}, d={d_bd:+.3f}")
        print(f"    Barbell={mean(barbell_fit):.4f} vs Hypercube={mean(hyper_fit):.4f}, d={d_bh:+.3f}")

        barbell_div = [r[gen]['diversity'] for r in data['barbell']]
        disc_div = [r[gen]['diversity'] for r in data['disconnected']]
        print(f"    Barbell diversity={mean(barbell_div):.4f} vs Disconnected diversity={mean(disc_div):.4f}")

    print("\n  Barbell interpretation:")
    print("  Two 4-island cliques with a single bridge.")
    print("  Within-clique: COMPLETE connectivity (λ₂=4 for K4)")
    print("  Between-clique: single bridge (bottleneck)")
    print("  Net effect: fast mixing within each clique + selective transfer between cliques")
    print("  This is the GOLDILOCKS topology — λ₂ underestimates its effective connectivity")
    print("  because the INTRA-clique connectivity is what drives early convergence,")
    print("  while the bridge preserves inter-clique diversity.")

    # ─── 5. Goldilocks Zone Analysis ─────────────────────────────────────────
    print("\n" + "=" * 80)
    print("  GOLDILOCKS ZONE: Middle-connectivity vs Extremes")
    print("=" * 80)

    # Group topologies
    low = ["disconnected"]
    mid = ["ring", "random-regular", "barbell", "star", "watts-strogatz", "hypercube"]
    high = ["complete"]

    for gen in [20, 30]:
        low_fit = [r[gen]['meanFitness'] for t in low for r in data[t]]
        mid_fit = [r[gen]['meanFitness'] for t in mid for r in data[t]]
        high_fit = [r[gen]['meanFitness'] for t in high for r in data[t]]

        low_div = [r[gen]['diversity'] for t in low for r in data[t]]
        mid_div = [r[gen]['diversity'] for t in mid for r in data[t]]
        high_div = [r[gen]['diversity'] for t in high for r in data[t]]

        print(f"\n  Gen {gen}:")
        print(f"    Low (disconnected):  fit={mean(low_fit):.4f}, div={mean(low_div):.4f}")
        print(f"    Mid (6 topologies):  fit={mean(mid_fit):.4f}, div={mean(mid_div):.4f}")
        print(f"    High (complete):     fit={mean(high_fit):.4f}, div={mean(high_div):.4f}")
        print(f"    Fitness ordering: {'Low < Mid < High' if mean(low_fit) < mean(mid_fit) < mean(high_fit) else 'UNEXPECTED'}")

    # ─── 6. Removing Star: Does Spearman improve? ───────────────────────────
    print("\n" + "=" * 80)
    print("  SENSITIVITY: Spearman WITHOUT Star")
    print("=" * 80)

    topos_no_star = [t for t in TOPOLOGIES if t != "star"]
    for gen in [20, 30]:
        l2_ns = [LAMBDA2[t] for t in topos_no_star]
        fit_ns = [mean([r[gen]['meanFitness'] for r in data[t]]) for t in topos_no_star]
        div_ns = [mean([r[gen]['diversity'] for r in data[t]]) for t in topos_no_star]

        rho_f, p_f = spearman_rank(l2_ns, fit_ns)
        rho_d, p_d = spearman_rank(l2_ns, div_ns)
        sig_f = "***" if p_f < 0.001 else "**" if p_f < 0.01 else "*" if p_f < 0.05 else "n.s."
        sig_d = "***" if p_d < 0.001 else "**" if p_d < 0.01 else "*" if p_d < 0.05 else "n.s."

        print(f"\n  Gen {gen} (7 topologies, excluding Star):")
        print(f"    λ₂ vs Fitness:   ρ = {rho_f:.4f}, p = {p_f:.6f}  {sig_f}")
        print(f"    λ₂ vs Diversity: ρ = {rho_d:.4f}, p = {p_d:.6f}  {sig_d}")

    # ─── 7. Per-run Spearman without Star ────────────────────────────────────
    print("\n" + "=" * 80)
    print("  PER-RUN SPEARMAN WITHOUT STAR (n=70)")
    print("=" * 80)

    for gen in [20, 30]:
        l2_all = []
        fit_all = []
        div_all = []
        for topo in topos_no_star:
            for run in data[topo]:
                l2_all.append(LAMBDA2[topo])
                fit_all.append(run[gen]['meanFitness'])
                div_all.append(run[gen]['diversity'])

        rho_fit, p_fit = spearman_rank(l2_all, fit_all)
        rho_div, p_div = spearman_rank(l2_all, div_all)
        sig_f = "***" if p_fit < 0.001 else "**" if p_fit < 0.01 else "*" if p_fit < 0.05 else "n.s."
        sig_d = "***" if p_div < 0.001 else "**" if p_div < 0.01 else "*" if p_div < 0.05 else "n.s."

        print(f"\n  Gen {gen} (n=70, excluding Star):")
        print(f"    λ₂ vs Fitness:   ρ = {rho_fit:.4f}, p = {p_fit:.8f}  {sig_f}")
        print(f"    λ₂ vs Diversity: ρ = {rho_div:.4f}, p = {p_div:.8f}  {sig_d}")

    # ─── 8. Variance decomposition: between-topology vs within-topology ─────
    print("\n" + "=" * 80)
    print("  VARIANCE DECOMPOSITION: Between vs Within Topology")
    print("=" * 80)

    for gen in [20, 30]:
        # Overall mean
        all_fit = [r[gen]['meanFitness'] for t in TOPOLOGIES for r in data[t]]
        grand_mean = mean(all_fit)

        # SS_between
        ss_between = sum(len(data[t]) * (mean([r[gen]['meanFitness'] for r in data[t]]) - grand_mean)**2
                        for t in TOPOLOGIES)
        # SS_within
        ss_within = sum((r[gen]['meanFitness'] - mean([r2[gen]['meanFitness'] for r2 in data[t]]))**2
                       for t in TOPOLOGIES for r in data[t])
        # SS_total
        ss_total = sum((v - grand_mean)**2 for v in all_fit)

        eta_sq = ss_between / ss_total if ss_total > 0 else 0
        # F-statistic
        k = len(TOPOLOGIES)
        N = len(all_fit)
        ms_between = ss_between / (k - 1)
        ms_within = ss_within / (N - k)
        F = ms_between / ms_within if ms_within > 0 else 0

        print(f"\n  Gen {gen} (Mean Fitness):")
        print(f"    SS_between = {ss_between:.8f}")
        print(f"    SS_within  = {ss_within:.8f}")
        print(f"    SS_total   = {ss_total:.8f}")
        print(f"    η² = {eta_sq:.4f} ({eta_sq*100:.1f}% of variance explained by topology)")
        print(f"    F({k-1},{N-k}) = {F:.3f}")

        # Same for diversity
        all_div = [r[gen]['diversity'] for t in TOPOLOGIES for r in data[t]]
        grand_mean_d = mean(all_div)
        ss_between_d = sum(len(data[t]) * (mean([r[gen]['diversity'] for r in data[t]]) - grand_mean_d)**2
                          for t in TOPOLOGIES)
        ss_within_d = sum((r[gen]['diversity'] - mean([r2[gen]['diversity'] for r2 in data[t]]))**2
                         for t in TOPOLOGIES for r in data[t])
        ss_total_d = sum((v - grand_mean_d)**2 for v in all_div)
        eta_sq_d = ss_between_d / ss_total_d if ss_total_d > 0 else 0
        ms_between_d = ss_between_d / (k - 1)
        ms_within_d = ss_within_d / (N - k)
        F_d = ms_between_d / ms_within_d if ms_within_d > 0 else 0

        print(f"\n  Gen {gen} (Diversity):")
        print(f"    SS_between = {ss_between_d:.8f}")
        print(f"    SS_within  = {ss_within_d:.8f}")
        print(f"    SS_total   = {ss_total_d:.8f}")
        print(f"    η² = {eta_sq_d:.4f} ({eta_sq_d*100:.1f}% of variance explained by topology)")
        print(f"    F({k-1},{N-k}) = {F_d:.3f}")

    # ─── 9. Window of maximum differentiation ───────────────────────────────
    print("\n" + "=" * 80)
    print("  WINDOW OF MAXIMUM DIFFERENTIATION")
    print("=" * 80)

    print("\n  Range of mean fitness across topologies at each generation:")
    for gen in [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
        fits = [mean([r[gen]['meanFitness'] for r in data[t]]) for t in TOPOLOGIES]
        divs = [mean([r[gen]['diversity'] for r in data[t]]) for t in TOPOLOGIES]
        spread_f = max(fits) - min(fits)
        spread_d = max(divs) - min(divs)
        print(f"    Gen {gen:>3}: Fitness range = {spread_f:.6f}, Diversity range = {spread_d:.6f}")

    print()

if __name__ == "__main__":
    main()
