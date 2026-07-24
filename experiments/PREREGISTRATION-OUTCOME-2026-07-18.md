# Outcome vs pre-registration — within/between decomposition

**Adjudicated 2026-07-18 against `PREREGISTRATION-within-between-2026-07-18.md` (locked before the run).**
Numbers from `results/direction/within_between_decomposition.csv` / `DECOMPOSITION.md`.
Reconstruction gate PASSED (within+between → pooled to 5.6e-17; seeds reproduced bit-for-bit; fitness matched, nk4 KW p=0.152).

## nk4 (rugged) — the decisive landscape

| Component | hub_in | hub_out | undirected | ordering |
|---|---|---|---|---|
| between-island | 0.345 | 0.203 | 0.064 | hub_in > hub_out > undirected, all sig (δ up to 1.00) |
| within-island  | 0.075 | 0.062 | 0.040 | hub_in > hub_out > undirected, all sig (δ 0.54–0.88) |

between spread 0.06→0.35 ≈ 4× the within spread 0.04→0.08.

## Verdict
- **Shared prediction (between carries hub_in ≫ hub_out ≥ undirected): CONFIRMED.** Between is strongly ordered and dominates in magnitude.
- **Lyra's "within is FLAT" (whole effect between-island): REFUTED.** Within-island is significantly ordered, all three pairs.
- **Claudius's "within is INVERTED" (hub_out/undirected higher within): REFUTED.** Within runs the SAME direction as between, not opposed.
- **Truth = a third, un-pre-registered option:** within-island tracks between-island in *direction* (hub_in > hub_out > undirected) but is dominated by between in *magnitude*. Both components pull the pooled ordering the same way.

## Mechanism (post-hoc, flagged as post-hoc)
"Who receives a homogenizing consensus signal" drives both components together:
- **hub_in** (leaves→hub): the 4 leaves are pure sources, never homogenized → high between; the hub receives 4 diverse migrant streams → high within. Highest of both.
- **hub_out** (hub→leaves): leaves all receive the *same* hub genotypes → converge toward hub → mid.
- **undirected**: everyone blends to consensus → lowest of both.

The pooled ordering is between+within pulling the SAME way, NOT (Lyra) between-only, NOT (Claudius) between-minus-within.

## What survives of each side
- Lyra: the *magnitude* claim ("most of the effect lives between islands") survives — between spread ≈4× within. The *categorical* claim ("within is flat") does not.
- Claudius: the intuition that migration perturbs within-island dispersion survives as "within is affected," but the *sign* (inverted) is wrong.

## onemax (smooth) — near-null
Same direction in both components, ~10× smaller (within 0.024–0.030, between 0.024–0.031). Ruggedness amplifies, does not create — as expected.

## Caveat
between-island here = mean cross-island pairwise Hamming (AMOVA partition), which relates to Wright F_ST but is not identical to variance-of-island-means. The AMOVA identity is exact and standard; the mechanism paragraph above is POST-HOC, not pre-registered.
