#!/usr/bin/env python3
"""
Analyze interference experiment across NK roughness levels.

Key question from Clio: does the 10:1 ratio (primary beta_1 effect vs
arrangement effect) hold across NK roughness, or does arrangement effect
GROW with landscape ruggedness?

Prediction from cactus group theory: non-abelian interference should
amplify on rugged landscapes (higher NK K), because the fitness landscape
makes path-dependent exploration more consequential.
"""

import os
import csv
import math
import json
from collections import defaultdict

TOPOLOGIES = [
    "interference-adjacent",
    "interference-separated",
    "interference-nested",
]

TOPO_LABELS = {
    "interference-adjacent": "Adjacent",
    "interference-separated": "Separated",
    "interference-nested": "Nested",
}

# Map NK level to results directory
NK_LEVELS = {
    "NK0": "results/interference-nk0",
    "NK2": "results/interference-nk2",
    "NK4": "results/interference",       # original runs
    "NK6": "results/interference-nk6",
}

CHECKPOINTS = list(range(0, 510, 10))


def load_csv(path):
    with open(path) as f:
        lines = [l for l in f if l.strip() and
                 (l.strip()[0].isdigit() or l.strip().startswith("generation"))]
    reader = csv.DictReader(lines)
    return list(reader)


def mean(xs):
    return sum(xs) / len(xs) if xs else 0.0


def std(xs):
    if len(xs) < 2:
        return 0.0
    m = mean(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))


def stderr(xs):
    return std(xs) / math.sqrt(len(xs)) if xs else 0.0


def eta_squared(groups):
    all_vals = [x for g in groups for x in g]
    if not all_vals:
        return 0.0
    grand_mean = mean(all_vals)
    ss_total = sum((x - grand_mean) ** 2 for x in all_vals)
    if ss_total == 0:
        return 0.0
    ss_between = sum(len(g) * (mean(g) - grand_mean) ** 2 for g in groups if g)
    return ss_between / ss_total


def f_statistic(groups):
    k = len(groups)
    n_total = sum(len(g) for g in groups)
    if k < 2 or n_total <= k:
        return 0.0
    grand_mean = mean([x for g in groups for x in g])
    ss_between = sum(len(g) * (mean(g) - grand_mean) ** 2 for g in groups if g)
    ss_within = sum(sum((x - mean(g)) ** 2 for x in g) for g in groups if g)
    df_between = k - 1
    df_within = n_total - k
    if ss_within == 0 or df_within == 0:
        return float('inf') if ss_between > 0 else 0.0
    return (ss_between / df_between) / (ss_within / df_within)


def cohens_d(xs, ys):
    if len(xs) < 2 or len(ys) < 2:
        return 0.0
    mx, my = mean(xs), mean(ys)
    sx, sy = std(xs), std(ys)
    pooled = math.sqrt(((len(xs) - 1) * sx**2 + (len(ys) - 1) * sy**2) /
                       (len(xs) + len(ys) - 2))
    if pooled == 0:
        return 0.0
    return (mx - my) / pooled


def load_nk_data(results_dir):
    """Load all interference data from a results directory."""
    data = defaultdict(lambda: defaultdict(list))
    fitness_data = defaultdict(lambda: defaultdict(list))
    file_counts = defaultdict(int)

    if not os.path.isdir(results_dir):
        return data, fitness_data, file_counts

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
                fitness_data[topo][gen].append(float(row["bestFitness"]))

    return data, fitness_data, file_counts


def main():
    print("=" * 75)
    print("INTERFERENCE NK SWEEP — Arrangement Effect vs Landscape Ruggedness")
    print("=" * 75)
    print()

    all_results = {}
    summary = {}

    for nk_label, results_dir in NK_LEVELS.items():
        data, fitness_data, file_counts = load_nk_data(results_dir)

        total_files = sum(file_counts.values())
        if total_files == 0:
            print(f"  {nk_label}: NO DATA (directory: {results_dir})")
            continue

        print(f"\n{'─' * 75}")
        print(f"  {nk_label}: {total_files} runs ({', '.join(f'{t}: {file_counts[t]}' for t in TOPOLOGIES)})")
        print(f"{'─' * 75}")

        # Compute eta-squared trajectory
        eta2_trajectory = {}
        for gen in CHECKPOINTS:
            groups = [data[topo].get(gen, []) for topo in TOPOLOGIES]
            if all(len(g) > 0 for g in groups):
                eta2_trajectory[gen] = eta_squared(groups)

        if not eta2_trajectory:
            print("  No complete generation data found.")
            continue

        # Key statistics
        peak_eta2 = max(eta2_trajectory.values())
        peak_gen = max(eta2_trajectory, key=eta2_trajectory.get)
        final_eta2 = eta2_trajectory.get(500, eta2_trajectory.get(max(eta2_trajectory.keys()), 0))
        mean_eta2 = mean(list(eta2_trajectory.values()))

        # Diversity at key generations
        key_gens = [50, 100, 200, 300, 500]
        print(f"\n  {'Gen':>5}  ", end="")
        for topo in TOPOLOGIES:
            print(f"  {TOPO_LABELS[topo]:>14}", end="")
        print(f"  {'eta^2':>8}  {'F':>8}")
        print(f"  {'─' * 65}")

        for gen in key_gens:
            if gen not in eta2_trajectory:
                continue
            groups = [data[topo].get(gen, []) for topo in TOPOLOGIES]
            f_val = f_statistic(groups)
            print(f"  {gen:>5}  ", end="")
            for topo in TOPOLOGIES:
                vals = data[topo].get(gen, [])
                if vals:
                    print(f"  {mean(vals):>10.4f}({stderr(vals):.3f})", end="")
                else:
                    print(f"  {'N/A':>14}", end="")
            print(f"  {eta2_trajectory[gen]:>8.4f}  {f_val:>8.2f}")

        # Pairwise at peak
        print(f"\n  Peak eta^2 = {peak_eta2:.4f} at gen {peak_gen}")
        for i, t1 in enumerate(TOPOLOGIES):
            for t2 in TOPOLOGIES[i + 1:]:
                d = cohens_d(data[t1].get(peak_gen, []), data[t2].get(peak_gen, []))
                print(f"    {TOPO_LABELS[t1]} vs {TOPO_LABELS[t2]}: d = {d:+.3f}")

        summary[nk_label] = {
            "peak_eta2": peak_eta2,
            "peak_gen": peak_gen,
            "final_eta2": final_eta2,
            "mean_eta2": mean_eta2,
            "n_runs": total_files,
        }
        all_results[nk_label] = eta2_trajectory

    # --- CROSS-NK COMPARISON ---
    print(f"\n\n{'=' * 75}")
    print("CROSS-NK COMPARISON: Does arrangement effect scale with roughness?")
    print(f"{'=' * 75}\n")

    if len(summary) < 2:
        print("Need at least 2 NK levels with data for comparison.")
        return

    print(f"  {'NK Level':>10}  {'Peak eta^2':>12}  {'Peak Gen':>10}  {'Final eta^2':>12}  {'Mean eta^2':>12}  {'Runs':>6}")
    print(f"  {'─' * 70}")
    for nk_label in ["NK0", "NK2", "NK4", "NK6"]:
        if nk_label in summary:
            s = summary[nk_label]
            print(f"  {nk_label:>10}  {s['peak_eta2']:>12.4f}  {s['peak_gen']:>10}  {s['final_eta2']:>12.4f}  {s['mean_eta2']:>12.4f}  {s['n_runs']:>6}")

    # Compute ratio: arrangement effect at each NK level
    if "NK4" in summary and len(summary) > 1:
        nk4_peak = summary["NK4"]["peak_eta2"]
        print(f"\n  Ratio relative to NK4 (peak eta^2 = {nk4_peak:.4f}):")
        for nk_label in ["NK0", "NK2", "NK4", "NK6"]:
            if nk_label in summary:
                ratio = summary[nk_label]["peak_eta2"] / nk4_peak if nk4_peak > 0 else 0
                print(f"    {nk_label}: {ratio:.2f}x")

    # Interpretation
    print(f"\n{'─' * 75}")
    print("INTERPRETATION")
    print(f"{'─' * 75}")

    nk_keys = sorted(summary.keys(), key=lambda x: int(x[2:]))
    peak_vals = [summary[k]["peak_eta2"] for k in nk_keys]

    if len(peak_vals) >= 2:
        # Check if arrangement effect increases with K
        increasing = all(peak_vals[i] <= peak_vals[i + 1] for i in range(len(peak_vals) - 1))
        if increasing:
            print("MONOTONIC INCREASE: Arrangement effect grows with landscape ruggedness.")
            print("Consistent with Clio's prediction: non-abelian interference amplifies")
            print("on rugged landscapes where path-dependent exploration matters more.")
        else:
            # Check for general trend
            from_first_to_last = peak_vals[-1] - peak_vals[0]
            if from_first_to_last > 0.01:
                print("GENERAL INCREASE (non-monotonic): Arrangement effect generally grows")
                print("with roughness, but not strictly monotonically.")
            elif abs(from_first_to_last) < 0.01:
                print("FLAT: Arrangement effect is roughly constant across NK levels.")
                print("The 10:1 ratio appears to be landscape-independent.")
            else:
                print("DECREASE: Arrangement effect shrinks with roughness.")
                print("Counter to Clio's prediction. Needs investigation.")

    if len(peak_vals) >= 2:
        ratio_range = max(peak_vals) / min(peak_vals) if min(peak_vals) > 0 else float('inf')
        print(f"\nRange: {min(peak_vals):.4f} to {max(peak_vals):.4f} (ratio: {ratio_range:.1f}x)")

    # Save summary JSON
    output_json = "results/interference_nk_sweep_summary.json"
    with open(output_json, 'w') as f:
        json.dump({
            "summary": summary,
            "trajectories": {k: {str(g): v for g, v in traj.items()} for k, traj in all_results.items()},
        }, f, indent=2)
    print(f"\nSummary saved to: {output_json}")

    # --- Generate plot ---
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        nk_colors = {"NK0": "#95a5a6", "NK2": "#3498db", "NK4": "#e74c3c", "NK6": "#8e44ad"}

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Left: eta^2 trajectories per NK level
        ax1 = axes[0]
        for nk_label in ["NK0", "NK2", "NK4", "NK6"]:
            if nk_label in all_results:
                traj = all_results[nk_label]
                gens = sorted(traj.keys())
                vals = [traj[g] for g in gens]
                ax1.plot(gens, vals, color=nk_colors[nk_label], linewidth=2, label=nk_label)

        ax1.set_xlabel("Generation", fontsize=12)
        ax1.set_ylabel("$\\eta^2$ (arrangement effect on diversity)", fontsize=12)
        ax1.set_title("Interference Effect Across NK Roughness", fontsize=13)
        ax1.legend(fontsize=11)
        ax1.grid(True, alpha=0.3)
        ax1.axhline(y=0.06, color='gray', linestyle=':', alpha=0.5)
        ax1.set_ylim(bottom=0)

        # Right: peak eta^2 vs K
        ax2 = axes[1]
        k_vals = [int(k[2:]) for k in nk_keys]
        peaks = [summary[k]["peak_eta2"] for k in nk_keys]
        finals = [summary[k]["final_eta2"] for k in nk_keys]

        ax2.bar([k - 0.15 for k in k_vals], peaks, width=0.3,
                color=[nk_colors[f"NK{k}"] for k in k_vals], alpha=0.8, label="Peak $\\eta^2$")
        ax2.bar([k + 0.15 for k in k_vals], finals, width=0.3,
                color=[nk_colors[f"NK{k}"] for k in k_vals], alpha=0.4, label="Final $\\eta^2$ (gen 500)")

        ax2.set_xlabel("NK K value (landscape ruggedness)", fontsize=12)
        ax2.set_ylabel("$\\eta^2$ (arrangement effect)", fontsize=12)
        ax2.set_title("Arrangement Effect Scales with Ruggedness?", fontsize=13)
        ax2.set_xticks(k_vals)
        ax2.legend(fontsize=10)
        ax2.grid(True, alpha=0.3, axis='y')

        plt.tight_layout()
        plot_path = "results/interference_nk_sweep.png"
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        print(f"Plot saved to: {plot_path}")

    except ImportError:
        print("matplotlib not available; skipping plot.")


if __name__ == "__main__":
    os.chdir("/home/lyra/projects/Topology-experiments")
    main()
