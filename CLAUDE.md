# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

SigProfilerSimulator generates realistic null hypothesis models of cancer mutational landscapes. It randomly shuffles observed somatic mutations while preserving mutational signatures at a user-selected resolution (SBS-6 through SBS-24576, plus indel classifications ID-83/ID-415). Uses a Monte Carlo approach with PCG random number generator (via `fastrand`). Works in unison with SigProfilerMatrixGenerator for mutation classification.

**Reference:** Bergstrom et al. (2020) BMC Bioinformatics 21:438. Detailed summary in `dev_notes/bergstrom_2020_summary.md`.

## Build and test commands

```bash
# Install the package (editable)
pip install -e .

# Install reference genome (required for tests and usage)
SigProfilerMatrixGenerator install GRCh37 --local_genome /path/to/genome/

# Run smoke test (just imports the package)
python test.py

# Run full test suite (requires GRCh37 genome installed)
pytest -s -rw tests

# Run a single test
pytest -s -rw tests/SPS_test.py::test_maf_file_behavior[seed1234-chromT-True]
```

Tests are parameterized across 4 configs varying seed presence and `chrom_based` flag. A session-scoped fixture runs the simulator 8 times (4 configs x 2 copies) before any test. Tests verify that seeded runs produce identical output and unseeded runs produce different output.

**Important:** Reproducibility requires the same number of CPU cores, because seeds are distributed one-per-worker via `np.random.SeedSequence.spawn()`.

## Architecture

The package has only two source modules:

- **`SigProfilerSimulator/SigProfilerSimulator.py`** (~1000 lines) — Main entry point `SigProfilerSimulator()` function. Handles seed setup (`setup_seeds()`), input validation, reference genome checks, context distribution file generation, mutation preparation, and multiprocessing pool orchestration. Delegates per-chromosome simulation to the simulator module.

- **`SigProfilerSimulator/mutational_simulator.py`** (~5900 lines) — Core simulation engine. Key functions:
  - `simulator()` — Main Monte Carlo loop. Runs per-chromosome, placing mutations at random positions using `fastrand.pcg32()`. Supports updating mode (mutations modify reference), TSB (transcriptional strand bias), BED regions, probability masks.
  - `mutation_preparation()` / `mutation_preparation_chromosomes()` / `mutation_preparation_region()` — Read input VCF/MAF files, parse mutation counts per context per sample.
  - `mut_tracker()` — Allocate mutations proportionally to chromosomes based on nucleotide context distributions.
  - `combine_simulation_files()` — Merge per-chromosome output files into final MAF/VCF.

### Simulation flow

1. `setup_seeds()` — Generate or read master seed, spawn child seeds (one per CPU core)
2. Validate reference genome files and generate context distribution files if missing (via SigProfilerMatrixGenerator)
3. `mutation_preparation*()` — Read input mutations, classify by context
4. `mut_tracker()` — Distribute mutations across chromosomes proportionally
5. Multiprocessing pool calls `simulator()` per chromosome — place mutations at random positions using PCG32 RNG
6. `combine_simulation_files()` — Merge per-chromosome results into per-sample MAF/VCF output

### Key dependencies

- `SigProfilerMatrixGenerator` — Mutation classification, reference genome management, context distributions
- `sigProfilerPlotting` — Visualization
- `fastrand` — PCG32 random number generator for deterministic, fast mutation placement

### Output structure

```
project_path/
├── input/              # User-provided VCF/MAF files (copied here)
├── output/
│   ├── SBS/, DBS/, ID/ # Mutation count matrices
│   ├── simulations/    # Simulated MAF/VCF files (1.maf, 2.maf, ...)
│   └── Simulator_seeds.txt
└── logs/               # .err and .out log files
```

## Supported mutation contexts

SBS: 6, 24, 96, 384, 1536, 6144, 24576 (increasing sequence context and optional transcriptional strand bias). Indels: ID-83, ID-415. DBS: 78, 186.
