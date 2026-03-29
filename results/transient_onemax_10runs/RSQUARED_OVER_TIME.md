# Temporal Profile of Topology Effects on OneMax

> Eta-squared and Spearman rho across 8 topologies, 10 runs each, 100 generations (OneMax, pop=200)

## Key Findings

- **Eta-squared (fitness) peaks at generation 30** with eta^2 = 0.8799
- **Eta-squared (diversity) peaks at generation 30** with eta^2 = 0.7567
- **Spearman rho (fitness) peaks at generation 30** with rho = 0.7619 (p = 0.0280)
- **Spearman rho (diversity) peaks at generation 0** with rho = 0.8729 (p = 0.0047)

## Signal Windows

| Metric | Window (eta^2 > 0.10) | Peak Gen | Peak Value |
|--------|----------------------|----------|------------|
| eta^2 fitness | gen 20--50 | 30 | 0.8799 |
| eta^2 diversity | gen 20--30 | 30 | 0.7567 |

Both eta-squared metrics drop below 0.10 at generation **10**.


## Full Table

| Gen | eta^2 fitness | eta^2 diversity | rho fitness | p(fit) | rho diversity | p(div) |
|-----|--------------|----------------|-------------|--------|--------------|--------|
| 0 | 0.0000 | 0.0000 | -0.3273 | 0.4287 | +0.8729 | 0.0047 |
| 10 | 0.0000 | 0.0000 | +0.0825 | 0.8461 | +0.3944 | 0.3336 |
| 20 | 0.7713 | 0.5764 | +0.6190 | 0.1017 | -0.5238 | 0.1827 |
| 30 | 0.8799 | 0.7567 | +0.7619 | 0.0280 | -0.5476 | 0.1600 |
| 40 | 0.5511 | 0.0558 | +0.5000 | 0.2070 | -0.1667 | 0.6932 |
| 50 | 0.1452 | 0.0199 | +0.5000 | 0.2070 | -0.2619 | 0.5309 |
| 60 | 0.0037 | 0.0021 | +0.6905 | 0.0580 | +0.4910 | 0.2166 |
| 70 | 0.0280 | 0.0070 | +0.7381 | 0.0366 | -0.2619 | 0.5309 |
| 80 | 0.0364 | 0.0116 | +0.7381 | 0.0366 | -0.8623 | 0.0059 |
| 90 | 0.0046 | 0.0018 | +0.0476 | 0.9108 | +0.0599 | 0.8880 |
| 100 | 0.0120 | 0.0071 | +0.0238 | 0.9554 | -0.4431 | 0.2715 |

## Curve Shape

### Fitness
- Peak: gen 30 (eta^2 = 0.8799)
- Half-peak on rise: gen 20
- Half-peak on fall: gen 40
- Rise width: 10 gens, Fall width: 10 gens
- Asymmetry (fall/rise): 1.00x

### Diversity
- Peak: gen 30 (eta^2 = 0.7567)
- Half-peak on rise: gen 20
- Half-peak on fall: gen N/A

## Per-Topology Means at Key Generations

### Generation 0

| Topology | lambda_2 | Mean Fitness (sd) | Diversity (sd) |
|----------|---------|-------------------|----------------|
| barbell | 0.07 | 0.4996 (0.0029) | 0.5010 (0.0108) |
| complete | 8.0 | 0.4996 (0.0029) | 0.5010 (0.0108) |
| disconnected | 0.0 | 0.4996 (0.0029) | 0.5010 (0.0108) |
| hypercube | 2.0 | 0.4996 (0.0029) | 0.5010 (0.0108) |
| random-regular | 0.27 | 0.4996 (0.0029) | 0.5010 (0.0108) |
| ring | 0.59 | 0.4996 (0.0029) | 0.5010 (0.0108) |
| star | 1.0 | 0.4996 (0.0029) | 0.5010 (0.0108) |
| watts-strogatz | 1.5 | 0.4996 (0.0029) | 0.5010 (0.0108) |

### Generation 10

| Topology | lambda_2 | Mean Fitness (sd) | Diversity (sd) |
|----------|---------|-------------------|----------------|
| barbell | 0.07 | 0.7778 (0.0056) | 0.3232 (0.0160) |
| complete | 8.0 | 0.7778 (0.0056) | 0.3232 (0.0160) |
| disconnected | 0.0 | 0.7778 (0.0056) | 0.3232 (0.0160) |
| hypercube | 2.0 | 0.7778 (0.0056) | 0.3232 (0.0160) |
| random-regular | 0.27 | 0.7778 (0.0056) | 0.3232 (0.0160) |
| ring | 0.59 | 0.7778 (0.0056) | 0.3232 (0.0160) |
| star | 1.0 | 0.7778 (0.0056) | 0.3232 (0.0160) |
| watts-strogatz | 1.5 | 0.7778 (0.0056) | 0.3232 (0.0160) |

### Generation 20

| Topology | lambda_2 | Mean Fitness (sd) | Diversity (sd) |
|----------|---------|-------------------|----------------|
| barbell | 0.07 | 0.9351 (0.0034) | 0.1204 (0.0099) |
| complete | 8.0 | 0.9417 (0.0030) | 0.1102 (0.0078) |
| disconnected | 0.0 | 0.9161 (0.0037) | 0.1499 (0.0059) |
| hypercube | 2.0 | 0.9342 (0.0045) | 0.1210 (0.0096) |
| random-regular | 0.27 | 0.9314 (0.0046) | 0.1248 (0.0091) |
| ring | 0.59 | 0.9292 (0.0040) | 0.1292 (0.0113) |
| star | 1.0 | 0.9251 (0.0050) | 0.1360 (0.0131) |
| watts-strogatz | 1.5 | 0.9356 (0.0031) | 0.1206 (0.0088) |

### Generation 30

| Topology | lambda_2 | Mean Fitness (sd) | Diversity (sd) |
|----------|---------|-------------------|----------------|
| barbell | 0.07 | 0.9786 (0.0013) | 0.0437 (0.0035) |
| complete | 8.0 | 0.9808 (0.0008) | 0.0394 (0.0035) |
| disconnected | 0.0 | 0.9665 (0.0020) | 0.0643 (0.0051) |
| hypercube | 2.0 | 0.9789 (0.0014) | 0.0439 (0.0031) |
| random-regular | 0.27 | 0.9778 (0.0015) | 0.0435 (0.0034) |
| ring | 0.59 | 0.9771 (0.0015) | 0.0465 (0.0037) |
| star | 1.0 | 0.9761 (0.0021) | 0.0477 (0.0058) |
| watts-strogatz | 1.5 | 0.9789 (0.0008) | 0.0426 (0.0035) |

### Generation 50

| Topology | lambda_2 | Mean Fitness (sd) | Diversity (sd) |
|----------|---------|-------------------|----------------|
| barbell | 0.07 | 0.9822 (0.0006) | 0.0328 (0.0041) |
| complete | 8.0 | 0.9826 (0.0009) | 0.0325 (0.0040) |
| disconnected | 0.0 | 0.9815 (0.0010) | 0.0347 (0.0043) |
| hypercube | 2.0 | 0.9823 (0.0009) | 0.0331 (0.0048) |
| random-regular | 0.27 | 0.9826 (0.0007) | 0.0330 (0.0046) |
| ring | 0.59 | 0.9823 (0.0005) | 0.0329 (0.0045) |
| star | 1.0 | 0.9824 (0.0008) | 0.0332 (0.0042) |
| watts-strogatz | 1.5 | 0.9825 (0.0008) | 0.0336 (0.0041) |

### Generation 100

| Topology | lambda_2 | Mean Fitness (sd) | Diversity (sd) |
|----------|---------|-------------------|----------------|
| barbell | 0.07 | 0.9817 (0.0013) | 0.0355 (0.0041) |
| complete | 8.0 | 0.9818 (0.0011) | 0.0348 (0.0036) |
| disconnected | 0.0 | 0.9815 (0.0014) | 0.0354 (0.0037) |
| hypercube | 2.0 | 0.9816 (0.0012) | 0.0351 (0.0039) |
| random-regular | 0.27 | 0.9819 (0.0009) | 0.0349 (0.0044) |
| ring | 0.59 | 0.9819 (0.0013) | 0.0349 (0.0037) |
| star | 1.0 | 0.9817 (0.0014) | 0.0358 (0.0038) |
| watts-strogatz | 1.5 | 0.9817 (0.0012) | 0.0351 (0.0043) |

## Interpretation

*Topology effects are transient.* The signal appears, peaks, and decays as the population converges.
The window of maximum topology influence is the critical regime for practitioners to understand.
