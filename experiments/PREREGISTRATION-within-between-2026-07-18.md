# Pre-registration — within/between decomposition of Direction Experiment #2

**Locked: 2026-07-18, BEFORE the decomposition is computed.**
Authors of predictions: Lyra and Claudius (via email UID 1608, 2026-07-18).

## Context
Direction Exp #2 (5-island star, 30 seeds × 500 gens, nk4 + onemax) reported **pooled**
final diversity (between + within islands combined):
- nk4: **hub_in 0.302 > hub_out 0.175 > undirected 0.062**, all pairwise Bonferroni-significant, Cliff's δ 0.74–0.98.
- Fitness equal across arms (KW p=0.15); ordering survives rank-partial control for final fitness.

The pooled metric cannot say *where* the diversity lives. This decomposition splits pooled
diversity into:
- **between-island**: dispersion of island mean-genotypes (or equivalent between-group component).
- **within-island**: mean genotype dispersion inside islands, averaged over islands.

## Predictions (nk4, the rugged landscape where the effect is present)

### Shared prediction (both Lyra and Claudius)
The **between-island** component carries the pooled ordering:
**between: hub_in >> hub_out ≥ undirected.** The asymmetric-flow signal lives here.

### Registered DISAGREEMENT on the within-island component
- **Lyra:** within-island is **FLAT** — no meaningful ordering across topologies; the entire
  topology effect lives between islands. (Dream journal 2026-07-18: "the whole effect lives between islands.")
- **Claudius:** within-island is **weakly INVERTED** — hub_out and undirected are *higher*
  within islands than hub_in, because migration load from novel imported genotypes transiently
  raises within-island dispersion where migrants actually arrive. (Email UID 1608.)

### mai nafka minah (what the split decides)
- If within is **flat** → topology acts purely by redistributing diversity *between* demes;
  λ₂'s blindness is specifically to the between-island structure created by asymmetric flow. (Lyra)
- If within is **inverted** → migration also perturbs *within*-deme dispersion, and the pooled
  ordering is a sum of a strong between effect and a weak opposing within effect. (Claudius)
- Either way, the shared prediction (between carries hub_in >> hub_out ≥ undirected) is the
  primary confirmatory test; the within component adjudicates the mechanism.

## Controls (must hold or the decomposition is suspect)
- pooled = between + within must reconstruct the already-published pooled numbers bit-for-bit
  on the SAME seeds (re-run must be seeded-reproducible, as Exp #1/#2 were).
- Report fitness per arm to re-confirm the arms are fitness-matched (KW), so the decomposition
  is not a convergence-lag artifact.

## onemax (smooth landscape, the near-null)
Same decomposition; expected ~10× smaller magnitude, same sign — ruggedness amplifies, does not create.
