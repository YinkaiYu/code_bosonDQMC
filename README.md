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

Run the fast comparison against checked-in fixture data:

```bash
make benchmark
```

The default benchmark does not recompute ED. It validates the comparison semantics for total particle number and total kinetic energy.

Optional ED recomputation:

```bash
make benchmark-ed
```

This requires the Python dependencies used by `benchmarks/ed/EDtriangle_symm_NEblock.py`, including QuSpin and Numba.

## Documentation

- `docs/physics.md` maps physical symbols and observables to code variables and output files.
- `docs/development.md` documents repository layout, Makefile targets, run directories, and benchmark extension.
- `AGENTS.md` gives instructions for coding agents working in this repository.
