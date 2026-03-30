# Timing Summary

Per-topology wall-clock time across seeds (from .time files, microsecond precision).

| Topology | Domain | Seeds | Mean (s) | Std Dev (s) | Min (s) | Max (s) |
|----------|--------|------:|--------:|-----------:|-------:|-------:|
| barbell | 15x15 maze | 3 | 2.17 | 0.01 | 2.17 | 2.18 |
| complete | 15x15 maze | 3 | 2.18 | 0.00 | 2.17 | 2.18 |
| disconnected | 15x15 maze | 3 | 2.17 | 0.01 | 2.16 | 2.18 |
| hypercube | 15x15 maze | 3 | 2.18 | 0.01 | 2.16 | 2.18 |
| random-regular | 15x15 maze | 3 | 2.18 | 0.00 | 2.18 | 2.18 |
| ring | 15x15 maze | 3 | 2.18 | 0.00 | 2.17 | 2.18 |
| star | 15x15 maze | 3 | 2.17 | 0.01 | 2.16 | 2.18 |
| watts-strogatz | 15x15 maze | 3 | 2.18 | 0.01 | 2.17 | 2.19 |
| barbell | 8x8 maze | 3 | 0.11 | 0.00 | 0.11 | 0.11 |
| complete | 8x8 maze | 3 | 0.11 | 0.00 | 0.11 | 0.11 |
| disconnected | 8x8 maze | 3 | 0.11 | 0.00 | 0.11 | 0.11 |
| hypercube | 8x8 maze | 3 | 0.11 | 0.00 | 0.11 | 0.11 |
| random-regular | 8x8 maze | 3 | 0.11 | 0.00 | 0.11 | 0.11 |
| ring | 8x8 maze | 3 | 0.11 | 0.00 | 0.11 | 0.11 |
| star | 8x8 maze | 3 | 0.11 | 0.00 | 0.11 | 0.11 |
| watts-strogatz | 8x8 maze | 3 | 0.11 | 0.00 | 0.11 | 0.11 |

**Note:** No .time files found for OneMax (pilot_onemax/ or transient_onemax_10runs/). OneMax timing data was not captured.

## Cross-Domain Comparison (mean seconds)

| Topology | 15x15 maze | 8x8 maze | Speedup |
|----------|----------:|--------:|--------:|
| barbell | 2.17 | 0.11 | 19.8x |
| complete | 2.18 | 0.11 | 19.9x |
| disconnected | 2.17 | 0.11 | 19.6x |
| hypercube | 2.18 | 0.11 | 19.6x |
| random-regular | 2.18 | 0.11 | 19.8x |
| ring | 2.18 | 0.11 | 19.5x |
| star | 2.17 | 0.11 | 19.6x |
| watts-strogatz | 2.18 | 0.11 | 19.5x |

