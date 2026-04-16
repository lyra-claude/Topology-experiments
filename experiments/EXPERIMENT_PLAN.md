# Multi-Domain Experiment Plan

Rerun all major experiments on Maze and Sudoku domains, complementing existing NK landscape results.

## Compute Time Estimates

From `results/timing_summary.md`:
- **8x8 maze**: ~0.11s per run (8 islands, 200 gens). At 500 gens, estimate ~0.28s.
- **15x15 maze**: ~2.18s per run (8 islands, 200 gens). At 500 gens, estimate ~5.5s.
- **NK domains**: no timing data captured, but NK4 runs similar to 15x15 maze.
- **Sudoku**: unknown, likely slower than NK due to constraint evaluation. Estimate ~3-8s per run at 500 gens.

Using 8x8 maze for initial runs (fast iteration), then 15x15 maze for publication if effects hold.

## Experiments

### 1. Directed Cycles Multi-Domain (Priority 1-2)

Replicates `run_directed_multidomain.sh` on new domains.

- **Topologies:** 8 digraphs (dag-layer, dag-wide, lowcyc-1, bidir-ring, two-cliques, mesh-cyclic, dense-triangles, ring-skip2)
- **Domains:** maze, sudoku (2 domains)
- **Seeds:** 30
- **Generations:** 500
- **Runs:** 8 topologies x 2 domains x 30 seeds = **480 runs**
- **Estimated time (8x8 maze):** 480 x 0.28s = ~2.2 min
- **Estimated time (15x15 maze):** 240 x 5.5s = ~22 min
- **Estimated time (sudoku):** 240 x 5s = ~20 min

**Adapt from:** `experiments/run_directed_multidomain.sh` — change DOMAINS list from `nk0 nk2 nk4 nk6` to `maze sudoku`. Output dir: `results/directed_multidomain_maze/` and `results/directed_multidomain_sudoku/`.

**Analysis:** reuse `experiments/analyze_directed_multidomain.py` with updated paths.

**Priority order:**
1. Directed cycles on Maze (strongest existing analysis, easiest to compare)
2. Directed cycles on Sudoku

### 2. Bridge Experiment (Priority 3-4)

Replicates `experiments/bridge_experiment.sh` (on `feat/bridge-experiment` branch) on new domains.

- **Topologies:** 8 (ring, ring-chord1/2/3, star, star-leaf1/2/3) — 4 graph pairs at 2 spectral families
- **Domains:** maze, sudoku (2 domains)
- **Seeds:** 30 (increase from original 20 for consistency with directed experiment)
- **Generations:** 500
- **Runs:** 8 topologies x 2 domains x 30 seeds = **480 runs**
- **Estimated time (8x8 maze):** 480 x 0.28s = ~2.2 min
- **Estimated time (sudoku):** 240 x 5s = ~20 min

**Adapt from:** `experiments/bridge_experiment.sh` (on `feat/bridge-experiment` branch). Cherry-pick or copy the script. Change DOMAINS, increase SEEDS to 30, update output directory.

**Analysis:** reuse `experiments/analyze_bridge.py` (also on bridge branch) with updated paths.

**Priority order:**
3. Bridge on Maze
4. Bridge on Sudoku

### 3. Interference Experiment (Priority 5, conditional)

Only run if bridge and directed experiments show topology effects on the new domains.

- **Topologies:** 3 (interference-adjacent, interference-separated, interference-nested)
- **Domains:** maze, sudoku (2 domains)
- **Seeds:** 30
- **Generations:** 500
- **Runs:** 3 topologies x 2 domains x 30 seeds = **180 runs**
- **Estimated time:** ~3.5 min (maze) + ~15 min (sudoku)

**Adapt from:** `experiments/run_interference.sh` — change DOMAIN from `nk4` to loop over `maze sudoku`. Output dirs: `results/interference_maze/`, `results/interference_sudoku/`.

**Analysis:** reuse `experiments/analyze_interference.py` with updated paths.

**Priority order:**
5. Interference on both (if results warrant)

## Total Run Count

| Experiment | Maze | Sudoku | Total |
|-----------|-----:|-------:|------:|
| Directed cycles | 240 | 240 | 480 |
| Bridge | 240 | 240 | 480 |
| Interference | 90 | 90 | 180 |
| **Total** | **570** | **570** | **1140** |

Note: the task description says ~1800 runs. The difference is because the directed cycles experiment here uses 2 new domains (not 4 NK variants). If we also want to compare against NK results in a single analysis, add the existing 960 NK runs — giving 1140 new + 960 existing = 2100 total data points.

## Estimated Total Compute Time

- **8x8 maze runs (570):** ~2.7 minutes
- **Sudoku runs (570):** ~47 minutes (conservative)
- **Total new runs:** ~50 minutes wall clock (sequential)

Parallelism: each seed is independent. Can trivially parallelize with GNU parallel or xargs.

## Shell Script Modifications Needed

### For directed cycles (`run_directed_multidomain.sh`):
1. Change `DOMAINS` from `"nk0 nk2 nk4 nk6"` to accept a command-line argument or hardcode `"maze sudoku"`
2. Change `OUTDIR` to `results/directed_maze` or `results/directed_sudoku`
3. No other changes — the `--domain` flag already supports maze

### For bridge (`bridge_experiment.sh`, from `feat/bridge-experiment`):
1. Copy script to this branch
2. Change `DOMAINS` from `"nk0 nk4"` to `"maze sudoku"`
3. Increase `SEEDS` from `seq 1 20` to the 30-seed list used elsewhere
4. Update `OUTDIR`

### For interference (`run_interference.sh`):
1. Wrap in domain loop instead of single `DOMAIN="nk4"`
2. Update `OUTDIR` per domain

### New domain flag values:
- Maze: `--domain maze` (already implemented)
- Sudoku: `--domain sudoku` (newly added in this commit)

## Success Criteria

For each experiment on each new domain, compute:
- **eta-squared** for topology effect on diversity (primary metric)
- **eta-squared** for topology effect on best fitness
- Compare with NK results to determine domain-dependence of topology effects
- Key question: does the directed cycle count (kappa) predict diversity on non-NK landscapes?
