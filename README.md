# Bosonic DQMC On A Triangular Lattice

This repository contains a finite-temperature grand-canonical determinant quantum Monte Carlo code for a two-flavor Bose-Hubbard model on a frustrated triangular lattice.

The primary executable is built from Fortran sources under `src/`. Runtime input and output use fixed filenames in the current working directory, so run the executable from a directory containing `paramC_sets.txt`, `confin.txt`, and `seeds.txt`.

## Build

```bash
make build
```

The executable is written to `build/bosonDQMC.out`. The Makefile detects local and server library roots:

- `/home/yyk/Lib_90_new`
- `/home/zxli_1/Lib_90_new`

Override paths when needed:

```bash
make build FC=mpiifort LIB_ROOT=/path/to/Lib_90_new
```

## Local Run

Run the default example:

```bash
make run-example
```

This delegates to `scripts/run_local.sh`, which builds the executable, changes into
`runs/examples/triangle_3x2`, and runs:

```bash
mpirun -np 1 <repo>/build/bosonDQMC.out
```

Run another directory:

```bash
bash scripts/run_local.sh runs/examples/triangle_3x2 1
```

## Server Submission

Submit the template with an `yyk_` job name. The available queues are `node6348` and `fat6348`.

```bash
sbatch -p node6348 --job-name yyk_triangle scripts/sbatch_triangle.sh runs/examples/triangle_3x2
sbatch -p fat6348 --job-name yyk_triangle_fat scripts/sbatch_triangle.sh runs/examples/triangle_3x2
```

## Benchmark

Run the strict live DQMC benchmark suite:

```bash
make benchmark
```

`make benchmark` is an alias for the live suite. `make benchmark-dqmc` is kept as the explicit DQMC target name.

This runs `mpirun -np 1` through the local helper in fresh temporary directories. The default live suite includes:

- one free-boson `U1=U2=0` analytic case
- the four ED reference cases recorded in `temp/benchmark.txt`

All live suite cases compare `total_NE`, `total_kinetic`, `doubleOcc`, `squareOcc`, `numsquare_up`, and `numsquare_do`. All live suite cases use `dtau = beta / Ltrot = 0.01`; the interacting cases use `Nbin = 100000`. On the current WSL workstation with `MPI_NP=1`, the full suite was observed at `real 604.50` seconds, about 10 minutes 5 seconds; allow 15 minutes or more under load. Do not reduce `Ltrot` or `Nbin` just to make a substantive algorithm benchmark faster.

For a fast live benchmark, run only the no-interaction `U1=U2=0` DQMC case against the analytic reference:

```bash
make benchmark-fast
```

For a fast check of JSON parsing, fixture normalization, and comparison-script behavior without running DQMC, use `make check-fixtures`. This is not physics validation and must not replace the live suite for algorithm changes.

Optional ED recomputation:

```bash
make benchmark-ed
```

This runs the dense fixed-particle-number ED script in `benchmarks/ed/EDtriangle_symm_NEblock.py` for the checked-in 3x2 ED parameter file.

## Documentation

- `docs/physics.md` maps physical symbols and observables to code variables and output files.
- `docs/development.md` documents repository layout, Makefile targets, run directories, and benchmark extension.
- `AGENTS.md` gives instructions for coding agents working in this repository.
