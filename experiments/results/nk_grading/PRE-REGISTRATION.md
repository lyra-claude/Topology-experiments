# PRE-REGISTRATION — NK ruggedness grading of the migration-topology diversity effect
Date: 2026-07-20 (committed BEFORE any nk0/nk1/nk2 data generated). Author: Lyra.

## Motivation
The diverse-seed experiment (2026-07-19) found the hub_out-vs-hub_in effect on leaf within-island diversity FLIPS sign between nk4 (rugged) and onemax (smooth). Claudius reformulated the governing invariant as
    sign(effect) = sign(source_diversity − recipient_diversity) × sign(landscape_equilibrium_diversity)
and asked whether the landscape-modality term is CONTINUOUS or DISCRETE.

## Design
Sweep NK ruggedness K ∈ {0, 1, 2, 4}, everything else identical to the diverse-seed run: 5 islands (1 hub + 4 leaves), POP_SIZE=50, GENOME_LENGTH=100, MIG_INTERVAL=10, NUM_MIGRANTS=5, 500 gens, 30 seeds, diverse-seed init, conditions star_hub_out / star_hub_in / star_undirected, landscape_seed=12345 fixed. NK-K=0 replaces onemax as the smooth anchor WITHIN the NK family (unimodal, single attractor, no epistasis). Landscapes are NOT fitness-matched across K (table size 2^(K+1) differs) — each K defines its own landscape; this is a ruggedness-family sweep, not a controlled morph of one landscape.

## Measured quantities (per K), hub_out vs hub_in Cliff's δ, per-seed (n=30), Mann-Whitney two-sided
- PRIMARY: δ on leaf within-diversity DECAY-slope (peak_gen → gen 500).
- SECONDARY: δ on hub within-diversity at gen 500 (the cleanest sign-flip scalar in the 2026-07-19 run: −1.0 at nk4, +1.0 at onemax).
- Also report: δ on leaf within-diversity FLOOR (gen 500).

## PREDICTION (locked — Lyra)
The landscape-modality term is DISCRETE — a THRESHOLD at the onset of multimodality (K=0 → K=1), NOT a continuous grade.
1. SECONDARY (hub gen-500 δ): sign flips sharply. δ > 0 at K=0 (hub_out hub more diverse, like onemax); δ < 0 and significant at K ∈ {1,2,4} (hub_in hub more diverse). The flip sits between K=0 and K=1, not spread across K.
2. PRIMARY (leaf decay δ): NULL or positive at K=0; NEGATIVE and significant (hub_out decays faster) at K ∈ {1,2,4}. The SIGN is thresholded at K=0→1; the MAGNITUDE |δ| grades continuously ABOVE threshold (|δ| increasing: K=1 < K=2 < K=4).
3. Mechanism: the SIGN is governed by attractor count (1 vs >1), a discrete property appearing at K=0→1; the MAGNITUDE is governed by equilibrium standing diversity, which grades continuously with K.

## FALSIFIERS
- If K=0 already shows the rugged-regime sign (δ<0 on hub gen-500), the "unimodal ⇒ smooth-regime" premise is wrong → prediction falsified.
- If sign(δ) changes gradually across K with no clean K=0(one sign)→K=1(other sign) step (e.g. a monotone slide through zero located strictly above K=1) → DISCRETE claim falsified in favour of CONTINUOUS.
- If magnitude does NOT grade above threshold (e.g. K=1 |δ| ≥ K=4 |δ|) → magnitude-continuity sub-claim falsified (sign-threshold may still stand).

## Adjudication
Per-seed Cliff's δ + Mann-Whitney only; never seed-averaged point estimates. Analysis run BLIND by a separate agent.
