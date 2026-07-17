# Results: Volume-vs-Spectrum Pilot

**Status:** Complete  
**Date:** 2026-07-17  
**Pre-registration:** `experiments/PREREGISTRATION-volume-vs-spectrum.md`  
**Data:** `experiments/results/volume_vs_spectrum/pilot.csv` (original), `pilot_with_fitness.csv` (fitness-instrumented re-run)  
**Analysis scripts:** `experiments/analyze_volume_vs_spectrum.py`, `experiments/analysis_partial_lyra_20260717.py`, `experiments/refute_lyra_20260717.py`

---

## 1. Setup

**Question:** Does migration-graph spectrum (algebraic connectivity λ₂) order final within-island genetic diversity once total migration volume is held exactly constant?

**Design:** 10 connected cubic (3-regular) graphs on 12 vertices, with λ₂ spanning [0.17, 1.47]. Because every graph is 3-regular, every island has exactly 3 neighbours, and total migration volume per generation is:

    12 islands × 3 neighbours × 5 migrants = 180  (identical for every graph by construction)

Volume is not a free variable — it cannot explain any diversity ordering across this graph set. Any remaining signal must come from structure, not quantity.

**Graph set (from `cubic12_graphset.json`):**

| Graph ID        | λ₂       |
|-----------------|----------|
| rand3_12_s148   | 0.167742 |
| rand3_12_s114   | 0.290725 |
| rand3_12_s13    | 0.447119 |
| rand3_12_s8     | 0.606910 |
| rand3_12_s43    | 0.732886 |
| rand3_12_s36    | 0.885092 |
| rand3_12_s89    | 1.032717 |
| rand3_12_s17    | 1.096788 |
| franklin        | 1.267949 |
| rand3_12_s168   | 1.467911 |

**GA parameters (matching the paper):** 12 islands, population 50/island, migration interval 10 generations, 5 migrants/edge/event, 500 generations.  
**Domains:** OneMax (negative control) and NK4 (N=100, K=4, rugged landscape).  
**Seeds:** 20 per (graph, domain) combination = 400 total runs.

**Primary metric:** Mean final within-island diversity (mean pairwise Hamming distance at generation 500), averaged over 12 islands then over 20 seeds, yielding one value per graph per domain.

---

## 2. Pre-registered hypotheses and decision rule

From the pre-registration (unchanged):

- **H_spectrum:** Higher λ₂ ⇒ faster mixing ⇒ less standing diversity. Prediction: strong negative Spearman ρ (ρ ≤ −0.7, p < 0.05) in NK4.
- **H_volume-was-everything:** Once volume is fixed, spectrum does nothing. Prediction: ρ ≈ 0, flat across the graph set.

**Decision rule:**
- Support H_spectrum if ρ ≤ −0.7 with rank test p < 0.05, in at least NK4.
- Refute the spectrum mechanism if |ρ| < 0.3 OR sign inconsistent across seeds/domains OR test fails.
- Ambiguous if 0.3 ≤ |ρ| < 0.7.

---

## 3. Results

### 3.1 Primary: Spearman ρ (λ₂ vs mean-seed diversity per graph)

| Domain  | Graph-level n | ρ (graph-level) | p-value | Run-level n | ρ (run-level) | Kendall W |
|---------|--------------|-----------------|---------|-------------|---------------|-----------|
| NK4     | 10           | **−0.842**      | 0.002   | 200         | −0.40         | 0.22      |
| OneMax  | 10           | +0.006          | —       | 200         | ≈0            | —         |

**NK4** clears the pre-registered threshold: ρ = −0.842 ≤ −0.7, p = 0.002 < 0.05. H_spectrum is supported.  
**OneMax** is flat (ρ ≈ 0) — the negative control passes cleanly.

The run-level ρ = −0.40 (n=200) is weaker and noisier, as expected: individual runs have high within-graph variance, and the spectral signal operates at the graph level. Kendall's W = 0.22 reflects modest but non-trivial agreement across seeds.

### 3.2 Robustness checks

**Leave-one-out (LOO):** Graph-level ρ was recomputed 10 times, each time dropping one graph. The LOO range is [−0.95, −0.78]. No single graph drives the result; the ordering is stable.

**Mean vs median aggregation:** Aggregating per-graph diversity with the median rather than the mean yields ρ = −0.879. The sign and magnitude are unchanged.

---

## 4. Fitness confound test

### 4.1 Why it was necessary

The original pilot runner (`run_pilot.py`) logged diversity but not fitness. This left open an alternative explanation: perhaps high-λ₂ graphs simply converge faster (reach higher fitness), and lower genetic diversity is a side-effect of being closer to the optimum — not of spectral mixing per se. Without fitness data, this confound was unmeasurable.

### 4.2 Instrumentation and determinism check

Script `run_pilot_with_fitness.py` was added to log `final_best_fitness` and `final_mean_fitness` at generation 500. The experiment was re-run with the exact same 400 seeds. Comparing the new `pilot_with_fitness.csv` against the original `pilot.csv`:

    Max |diversity_new − diversity_original| = 0.00  (across all 400 runs)

Diversity reproduced bit-for-bit. The GA is deterministic given seed; re-running with fitness logging does not alter the diversity trajectories. This also serves as an implicit integrity check on the codebase.

### 4.3 Partial correlation results

With fitness now available, partial Spearman correlation was computed:

| Level       | ρ(λ₂, diversity) | ρ(λ₂, diversity | final_mean_fitness) | Attenuation |
|-------------|------------------|--------------------------------------|-------------|
| Graph-level | −0.84            | −0.82                                | <3%         |
| Run-level   | −0.40            | −0.41                                | 0% (slight increase) |

Controlling for fitness leaves the diversity-spectrum association essentially unchanged.

**λ₂ vs fitness:** graph-level ρ = +0.39 (not significant), run-level ρ = +0.05. Spectral connectivity does not strongly predict the level of fitness achieved — the graphs are producing the same optimization outcomes at different standing diversities.

**Interpretation:** Higher-λ₂ graphs produce less within-island diversity at the same final fitness level. The mechanism is deme homogenization by spectral structure, not faster convergence to the optimum.

---

## 5. Verdict

**H_spectrum is supported** (NK4: graph-level ρ = −0.842, p = 0.002; LOO range [−0.95, −0.78]; partial ρ = −0.82 after fitness control; negative control clean).

The paper's mechanism claim survives the volume-control test. λ₂ orders diversity in NK4 even when total migration volume is exactly constant across graphs. The spectral bridge is not simply a proxy for migration volume.

---

## 6. Caveats and scope

1. **n = 10 graphs is a pilot.** The graph-level Spearman ρ is based on 10 data points. The result is statistically significant and robust to LOO, but a larger graph set would further harden the conclusion.

2. **Run-level signal is modest.** ρ = −0.40 at the run level reflects genuine within-graph variance. Seeds with identical spectral structure can produce different diversity outcomes; the spectral ordering is a tendency, not a deterministic rule.

3. **Uniform stationary distribution (π uniform) by construction.** All 10 graphs are 3-regular, so the random walk on each has uniform stationary distribution. This design successfully controls for volume but it cannot say anything about the role of **edge direction** or **non-uniform π**. Directed graphs with asymmetric degree (e.g., asymmetric star topologies) can break the λ₂ ↔ diversity link through direction alone — this is a qualitatively different mechanism.

4. **Experiment #2 (planned): asymmetric star, non-uniform π.** The next experiment will vary edge direction at fixed volume, testing whether the directed Laplacian's spectral structure adds explanatory power beyond what λ₂ of the symmetrised graph captures. The current pilot is silent on this question by design.

5. **Domain scope.** OneMax is a clean negative control but not a challenging landscape. NK4 is the relevant test case. Results should not be extrapolated to domains with qualitatively different fitness landscapes without further experiments.

---

## 7. Files

| File | Description |
|------|-------------|
| `experiments/PREREGISTRATION-volume-vs-spectrum.md` | Pre-registration (unchanged) |
| `experiments/results/volume_vs_spectrum/pilot.csv` | Original 400-run results (diversity only) |
| `experiments/results/volume_vs_spectrum/pilot_with_fitness.csv` | Fitness-instrumented re-run (400 rows, diversity reproduced bit-for-bit) |
| `experiments/run_pilot_with_fitness.py` | Runner that adds fitness logging |
| `experiments/analyze_volume_vs_spectrum.py` | Primary analysis: graph-level ρ, LOO, Kendall W |
| `experiments/analysis_partial_lyra_20260717.py` | Partial Spearman controlling for fitness |
| `experiments/refute_lyra_20260717.py` | Adversarial refutation check |
