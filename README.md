# Topology-Experiments: How Communication Topology Affects Multi-Agent Performance

Empirical experiments measuring how the communication graph topology of island-model genetic algorithms affects population diversity and fitness. We use two topological invariants -- the first Betti number (beta_1, cycle rank) and algebraic connectivity (lambda_2) -- to predict and explain GA performance across multiple fitness domains.

## Key Findings

- **beta_1 outpredicts model choice.** The cycle rank of the migration graph explains more variance in diversity (eta-squared up to 0.69) than population size, migration rate, or other hyperparameters.
- **NK landscape ruggedness amplifies topology effects.** Smooth landscapes (K=0) show near-zero topology effect; rugged landscapes (K=6) show eta-squared = 0.69. Topology IS landscape-dependent.
- **Two-invariant temporal separation.** beta_1 drives transient diversity (peaks early, fades by generation 500), while lambda_2 drives persistent diversity (emerges later, persists). Two independent channels confirmed by iso-spectral bridge experiments (320 runs).
- **Directed cycles decorrelate from density.** Undirected beta_1 is confounded with edge density (beta_1 = |E| - n + 1). Directed cycle count escapes this confound -- 8 digraphs at constant density show 0 to 47 directed cycles. Directed cycle count predicts diversity at constant density (r = -0.68, eta-squared = 0.17).
- **Non-abelian interference.** Different cycle arrangements at constant beta_1 produce different eta-squared (0.045 difference). Shared paths create correlated information flow. Arrangement is ~10% of the primary beta_1 effect.
- **Goldilocks zone.** Too few cycles = fragile convergence; too many = wandering without selection pressure. Intermediate beta_1 is optimal, and the optimum depends on landscape ruggedness.

## Architecture

The GA engine is written in **Haskell** (`src/`). Analysis and plotting scripts are in **Python**.

### Source (Haskell)

| File | Purpose |
|------|---------|
| `src/Main.hs` | CLI entry point, experiment runner |
| `src/IslandGA.hs` | Island-model GA with configurable topology |
| `src/Domain.hs` | Typeclass for fitness domains |
| `src/NKLandscape.hs` | NK landscape fitness domain (parameterized ruggedness) |
| `src/Maze.hs` / `src/Maze8.hs` / `src/MazeCore.hs` | Maze generation and solving domain |
| `src/OneMax.hs` | OneMax (bit-counting) baseline domain |
| `src/SudokuSolver.hs` | Sudoku fitness domain |

### Analysis (Python)

| File | Purpose |
|------|---------|
| `compute_lambda2.py` | Compute lambda_2 for all topologies |
| `topologies.py` | Adjacency matrices for candidate topologies |
| `analyze_pilot.py` | Statistical analysis of pilot experiments |
| `analyze_directed_full.py` | Directed cycle experiment analysis |
| `analyze_deep_transient.py` | Transient dynamics and temporal separation |
| `plot_star_vs_cycle.py` | Star vs cycle comparison figures |

### Experiment Scripts

Shell scripts in the root directory (`run_pilot.sh`, `run_directed_full.sh`, etc.) orchestrate batch experiment runs.

## Experiments (by branch)

| Branch | Experiment | Runs |
|--------|-----------|------|
| `feat/pilot-run-batch1` | Initial topology sweep (8 topologies x NK0/2/4/6) | ~160 |
| `feat/bridge-experiment` | Iso-spectral bridge: beta_1 vs lambda_2 separation | 320 |
| `feat/directed-cycle-experiment` | Directed cycles at constant density | 240 |
| `feat/directed-cycles-multidomain` | Directed cycles across NK0/2/4/6 | 960 |
| `feat/interference-experiment` | Cycle arrangement at constant beta_1 | 180+ |
| `feat/star-vs-cycle` | Star vs cycle comparison | 270 |
| `feat/multi-domain-experiments` | Maze and Sudoku cross-domain validation | 480+ |
| `paper/ecta2026` | ECTA 2026 paper draft and figures | -- |

Results (CSV data and figures) are stored in `results/`.

## Running Experiments

### Prerequisites

- GHC 9.6+ and Cabal (Haskell toolchain)
- Python 3 with `numpy`, `scipy`, `matplotlib`, `pandas`, `seaborn`

### Build and run

```bash
# Build the GA engine
cabal build

# Run a single experiment (example: NK landscape, K=2, ring topology, 8 islands)
cabal run topology-experiments -- --domain nk --nk-k 2 --topology ring --islands 8 --generations 500 --runs 10

# Run a batch experiment
bash run_pilot.sh
```

### Analysis

```bash
# Compute spectral properties
python3 compute_lambda2.py

# Analyze pilot results
python3 analyze_pilot.py
```

## Related Publications

- [The One Number That Predicts Whether Your AI Agent Team Will Work](https://medium.com/@lyraclaude20/the-one-number-that-predicts-whether-your-ai-agent-team-will-work-4554a2e92a11) -- introducing beta_1 for multi-agent systems
- [Why 85% x 85% x 85% Is the Wrong Math for Your AI Agent Team](https://medium.com/@lyraclaude20/why-85-x-85-x-85-is-the-wrong-math-for-your-ai-agent-team-d4c6a6a0808a) -- error propagation and topology
- [The Topology of Agent Attacks](https://medium.com/@lyraclaude20/the-topology-of-agent-attacks-ac2af55e3304) -- security implications of communication topology

## Collaboration

This project is a collaboration between **Lyra** (GA simulations, empirical analysis) and **Claudius** (spectral theory, topology selection), with **Robin** providing the mathematical framework and research direction.

## License

Research code. Not yet formally licensed.
