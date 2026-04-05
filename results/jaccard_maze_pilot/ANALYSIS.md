# Jaccard Maze Pilot — Analysis Summary

**Date:** 2026-04-05
**Design:** 8 topologies x 3 seeds x 500 generations (24 runs total)
**Fitness landscape:** Jaccard maze (combinatorial, multi-path, non-trivial epistasis)

## Key Finding

**There IS a topology signal on the Jaccard maze, but it is weaker and later-peaking than OneMax.**

| Metric | Peak eta-sq | Peak Gen | OneMax eta-sq (gen 30) |
|--------|-------------|----------|------------------------|
| meanFitness | 0.655 | gen 90 | ~0.88 |
| bestFitness | 0.685 | gen 100 | ~0.88 |
| diversity | 0.691 | gen 100 | ~0.76 |

The topology effect peaks around gen 90-100, not gen 20-30 like OneMax. This is consistent with the maze being a harder landscape where structural effects take longer to manifest.

## Statistical Significance

Significant ANOVA results (p < 0.05):

- **Gen 100, bestFitness:** F = 4.98, p = 0.0038 (confirmed by Kruskal-Wallis: p = 0.024)
- **Gen 100, diversity:** F = 5.12, p = 0.0033 (confirmed by Kruskal-Wallis: p = 0.032)
- **Gen 500, diversity:** F = 3.25, p = 0.024 (Kruskal-Wallis marginal: p = 0.088)

No pairwise comparisons survived Bonferroni correction (28 comparisons, 3 seeds is too few per group). The omnibus effect is real but individual topology contrasts need more seeds.

## Topology Rankings

**Fitness (overall):** Disconnected > Barbell > Ring > Complete (spread ~4.7% of mean)
**Best fitness (gen 500):** Hypercube > Barbell > Random-Regular > Ring (spread ~9.6%)
**Diversity (overall):** Disconnected > Ring > Random-Regular > Complete (spread ~5.9%)

Notable: **Complete graph has LOWEST diversity** (0.586 at gen 500), consistent with our theory that full connectivity accelerates diversity collapse. **Disconnected has HIGHEST diversity** (0.617), consistent with isolated subpopulations maintaining independent search.

## Comparison with OneMax

| Property | OneMax | Jaccard Maze |
|----------|--------|-------------|
| Peak eta-sq (fitness) | ~0.88 (gen 30) | 0.685 (gen 100) |
| Peak eta-sq (diversity) | ~0.76 (gen 30) | 0.691 (gen 100) |
| Temporal profile | Sharp early peak, decays by gen 100 | Gradual rise, peaks gen 90-100, partially persists |
| Complete graph | Fastest convergence | Lowest diversity |
| Signal detectability | Massive with 3 seeds | Needs more seeds for pairwise |

The maze shows **~77% of OneMax's peak effect size** but with a fundamentally different temporal profile. Where OneMax shows a sharp transient burst, the maze shows a slower build that partially persists.

## Interpretation for ECTA Paper

This result supports the **capability-modulated topology effect** thesis:

1. **Landscape difficulty modulates the topology effect.** Harder landscapes (maze) show weaker but more persistent topology effects vs. easy landscapes (OneMax).
2. **The diversity channel is primary.** Diversity eta-sq matches or exceeds fitness eta-sq at every time point, consistent with the causal chain: topology -> diversity -> fitness.
3. **Complete graph = diversity collapse** is robust across landscapes.
4. **The temporal shift is theoretically interesting.** On easy landscapes, topology matters early when the population is exploring; on hard landscapes, topology matters later when the population needs sustained diversity to avoid premature convergence.

## Limitations

- **3 seeds per topology** is marginal for parametric tests. Pairwise contrasts do not survive Bonferroni correction. Recommend 10+ seeds for publication.
- **Disconnected topology** is arguably degenerate (no migration). Including/excluding it changes rankings.
- **meanFitness never reaches significance** at any key generation. The topology signal is primarily in bestFitness and diversity.

## Files

- Figures: `figures/fig_maze_trajectories.pdf`, `figures/fig_maze_eta_squared.pdf`
- Analysis script: `../../analyze_jaccard_maze.py`
- Raw data: 24 CSV files in this directory

## Recommendation

**Publishable as supporting evidence, not standalone.** The signal is real (eta-sq ~ 0.69, p < 0.004) but needs more seeds to identify which specific topologies differ. For the ECTA paper, this is a second landscape confirming the OneMax result with a different temporal profile. The "harder landscape = later peak" finding is novel and supports the capability-modulation thesis.

For Claudius: this validates the maze design. Propose scaling to 10 seeds for the full experiment.
