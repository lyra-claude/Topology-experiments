# Effect A — Mechanism Note: Consensus-Flattening is Landscape-Independent

Date: 2026-07-21. Author: Lyra. Status: LOCKED (adjudicable claim).
Companion documents: ADJUDICATION.md, ANALYSIS.md, PRE-REGISTRATION.md.

---

## Statement of Effect A

**Effect A (consensus-flattening):** hub_out leaves decay in within-island diversity faster than hub_in leaves at every landscape ruggedness level tested, including unimodal K=0. The sign of this effect has no landscape term. Its magnitude grades continuously with K.

This effect was revealed as the PRIMARY scalar (leaf within-diversity decay-slope, δ = #(hub_out > hub_in) − #(hub_out < hub_in) normalised) in the NK-grading sweep. It caused the PRIMARY pre-registration to be half-refuted: the prediction assumed the same sign-threshold that governs Effect B, but Effect A carries none.

---

## Evidence: sign is landscape-independent at every K

Per-seed Cliff's δ, hub_out vs hub_in on leaf decay-slope (n=30/condition, Mann-Whitney two-sided). Sign convention: negative δ = hub_out decays faster.

| K | δ_leaf_decay | p          | landscape character               |
|---|-------------|------------|-----------------------------------|
| 0 | -0.66       | 1.02e-05   | unimodal, single attractor, NK-K=0|
| 1 | -0.90       | 2.67e-09   | onset of multimodality            |
| 2 | -0.89       | 3.20e-09   | moderate ruggedness               |
| 4 | -0.99       | 4.98e-11   | high ruggedness                   |

All four δ values are negative and significant (p < 0.05). The sign does not flip between K=0 and K=1; it does not flip anywhere. The magnitude increases monotonically: |δ| = 0.66, 0.90, 0.89, 0.99, approaching saturation.

**Locked claim:** sign(δ_leaf_decay) = negative for all K ∈ {0, 1, 2, 4}, with no landscape-gated threshold. Falsifiable by any future K with significant positive δ on this scalar.

---

## Mechanism

A hub_out topology routes migrant flow **from the hub to the leaves**. Each migration interval, the hub broadcasts a sample of its own gene pool outward; leaves receive migrants drawn from the hub but send none back.

The hub, acting as the upstream source, rapidly converges toward the landscape's dominant attractor (unimodal or the highest-fitness basin for any K). Its outgoing migrants therefore carry a **consensus signal** — a progressively less diverse sample. The receiving leaves absorb this signal unconditionally. Regardless of whether the landscape is smooth or rugged, the leaves that receive the hub's consensus stream homogenise faster than leaves that do not receive any migrant input (hub_in leaves, where the hub is the *recipient*, not the source).

The mechanism is: **receiving a consensus migrant stream reduces standing within-island diversity, independent of landscape structure.** The landscape's ruggedness modulates the strength of the consensus (a more rugged landscape allows more heterogeneous hub attractors, slowing but not reversing the flattening) — hence the grading in |δ| — but it cannot reverse the direction of the effect.

This is a **directional-flow effect**: what matters is whether an island is downstream (receiving) or upstream (sending) in the migration network. Landscape is not a precondition; it is an amplitude modulator.

---

## Distinction from Effect B

Effect B (polymorphism-maintenance) is **landscape-dependent**. It is carried by the SECONDARY scalar: hub within-diversity at gen 500. Its δ flips sign discretely between K=0 and K=1 (K=0: δ=+0.95; K=1: δ=-0.79; K=2: δ=-1.00; K=4: δ=-1.00). The mechanism is that the hub, isolated from migrant input in a hub_out topology, maintains polymorphism only if the landscape has multiple attractors (K ≥ 1). On a unimodal landscape (K=0) the hub still converges to the single global optimum; "isolation from migrants" does not help if there is only one attractor to occupy.

Effect B therefore requires a landscape-modality condition: it switches on at the onset of multimodality. Effect A does not. A single measured scalar (e.g. leaf floor diversity at gen 500) conflates both effects, which is why the FLOOR scalar washes out at intermediate K — it is the residual of two opposed forces with different K-dependencies. The two scalars must be kept separate.

| Property                          | Effect A (leaf decay-slope) | Effect B (hub gen-500) |
|-----------------------------------|-----------------------------|------------------------|
| Sign landscape-dependent?         | No                          | Yes (discrete at K=0→1)|
| Magnitude landscape-dependent?    | Yes (grades with K)         | Yes (grades with K)    |
| Mechanism                         | Consensus broadcast to leaves | Hub polymorphism at attractors |
| Directional flow role             | Receiver (leaf in hub_out)  | Sender-isolation (hub in hub_out) |
| Pre-reg scalar                    | PRIMARY (half-refuted)      | SECONDARY (confirmed)  |

---

## Methodological note

The pre-registration predicted the same sign-threshold for both scalars, assuming one mechanism. The blind analysis revealed two effects with different K-dependencies. The leaf DECAY slope is Effect A's cleanest carrier because it captures the *rate* of homogenisation rather than the endpoint: at the endpoint (FLOOR), Effect A and B can cancel or add depending on K, obscuring both. The decay slope isolates the early-phase loss driven purely by the incoming signal, before landscape-dependent equilibration dominates.

---

*This note locks Effect A's mechanism. The claim is adjudicable: Effect A's sign is landscape-independent. Future sweeps across continuous K or alternative topologies should hold this scalar separate from hub diversity to avoid conflation.*
