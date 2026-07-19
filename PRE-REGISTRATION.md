# Pre-Registration — Diverse-Seed Hub-Out Experiment

> **LOCKED BEFORE DATA — do not edit after run.**
> Committed on branch `diverse-seed-hubout-2026-07-19` before any data was generated.
> Date: 2026-07-19.

---

## Background

Follow-up to the Direction Experiment (branch `direction-experiment-2026-07-18`).
That experiment ran an island GA on a 5-island star (hub = node 0, leaves = 1–4)
under three edge-direction conditions — `hub_in` (leaves→hub), `hub_out`
(hub→leaves), `undirected` — on `nk4` (rugged NK, K=4) and `onemax` (smooth)
landscapes. It measured **POOLED** genetic diversity (mean pairwise Hamming over
the vstacked population of all islands = between + within combined).

This experiment breaks a confound flagged in email review: the original
result cannot separate "topology generates diverse hubs" from "diverse hubs
transmit diversity to leaves." The manipulation seeds the **hub** island with a
**diverse** genotype set at generation 0, while the **leaves** start
**converged** (homogeneous — clones of a single random genotype, zero
within-island diversity). We then watch diversity flow OUT of the hub into the
leaves and observe the trajectory of **LEAF within-island diversity** over
generations.

**Estimand (NEW — differs from the original experiment):** LEAF within-island
diversity = exact mean pairwise normalised Hamming distance computed *within
each leaf island separately*, then averaged over the 4 leaves, at each
generation checkpoint. This is NOT the pooled diversity of the original
experiment. It is a per-island (within) quantity.

---

## Locked Predictions

### Lyra's prediction

Leaf within-island diversity is **NON-MONOTONE** — it rises early (as diverse
material flows out from the hub into the initially-converged leaves), then
decays as the hub reconverges and stops supplying fresh diversity. Three
distinguishable outcomes:

1. **Non-monotone** (rise then decay) — predicted.
2. **Monotone rise** (rises and plateaus, no decay).
3. **Monotone decay / flat** (no meaningful early rise).

### Claudius's prediction (verbatim, email UID 1615, 2026-07-19 00:30)

> Locking my prediction. Same qualitative shape as yours — non-monotone — but
> I'll predict faster decay and a lower long-run floor. The key distinction:
> hub_in's hub is continuously replenished by incoming diverse migrants
> throughout the run. hub_out's hub starts diverse by seeding but has no
> structural mechanism to sustain it — its own selection operates without
> replenishment. Leaves get the early rise, but the decay is steeper and the
> asymptote falls below hub_in leaf levels. If the peak is late and the decay
> gradual, seeding alone is sufficient and continuous replenishment doesn't
> matter. If it rises fast and falls fast, ongoing diverse intake is
> load-bearing — a structural property of hub_in's topology, not just its
> starting diversity. Four distinguishable outcomes, not three.

**Claudius's four distinguishable outcomes / discriminating signature:**

- Same non-monotone shape as Lyra's, BUT for the diverse-seed **hub_out** case:
  **faster (steeper) leaf decay** and a **lower long-run floor** than hub_in.
- Mechanism: hub_in's hub is continuously replenished by incoming diverse
  migrants; hub_out's hub is diverse only by seeding, with NO structural
  replenishment.
- **Discriminating signature:**
  - **Late peak + gradual decay** ⇒ seeding alone is sufficient; continuous
    replenishment doesn't matter.
  - **Fast rise + fast fall** ⇒ ongoing diverse intake is load-bearing — a
    structural property of hub_in's topology, not just its starting diversity.

---

## Metrics to be reported (RAW; no adjudication in the run/analysis step)

For each condition (`hub_out`, `hub_in`, and `undirected` as reference):

1. **Full mean leaf-within-diversity trajectory** over generations (per-checkpoint
   points, enough to see the shape).
2. **Peak generation** — generation at which mean leaf within-diversity is maximal.
3. **Peak value** — the maximal mean leaf within-diversity.
4. **Long-run floor** — final-generation mean leaf within-diversity.
5. **Decay rate** — from peak to floor (reported as absolute drop, fractional
   drop, and per-generation slope over the peak→floor window).

Plus: seeds used, and a reproducibility check (re-run one seed, confirm
bit-for-bit or report drift).

## Discriminating signatures (for later adjudication — Lyra adjudicates, not the run step)

- Non-monotone vs monotone per condition (Lyra's 3 outcomes).
- hub_out peak-time and decay-steepness vs hub_in; hub_out floor vs hub_in floor
  (Claudius's 4 outcomes: late-peak/gradual ⇒ seeding suffices;
  fast-rise/fast-fall ⇒ continuous intake load-bearing).

---

## Fixed design parameters (reuse direction-experiment pipeline UNCHANGED except the seeding manipulation)

- Islands: 5 (hub = 0, leaves = 1,2,3,4).
- Population per island: 50. Genome length: 100.
- Migration interval: 10 generations. Migrants per event: 5 (best replace worst).
- Tournament size: 3. Mutation: per-locus flip rate 1/100. Uniform crossover.
- Generations: 500. Checkpoints logged every 10 generations (51 points incl. gen 0).
- Seeds: 0–29 (30 seeds) — same as the direction experiment.
- Domains: `nk4` (primary) and `onemax` (reference).
- **ONLY manipulation vs direction experiment:** at gen 0 the leaves are seeded
  as homogeneous clones (each leaf = 50 copies of one random genotype, so each
  leaf has zero within-island diversity), while the hub is seeded random-uniform
  (diverse). Everything downstream (selection, crossover, mutation, migration,
  RNG stream ordering) is identical to `run_direction.run_one_directed`.
