# Agent Instructions

This repository contains a finite-temperature grand-canonical determinant quantum Monte Carlo code for a two-flavor Bose-Hubbard model on a triangular lattice.

Start by reading these files before changing code:

- `README.md`
- `docs/physics.md`
- `docs/development.md`
- `benchmarks/README.md`

## Repository Map

- `src/` contains the Fortran source files.
- `Makefile` is the primary build and workflow entrypoint.
- `runs/examples/triangle_3x2/` contains the committed example run inputs.
- `scripts/run_local.sh` runs a local MPI job from a run directory.
- `scripts/sbatch_triangle.sh` is the Slurm submission template.
- `benchmarks/compare.py` performs the fast fixture/reference comparison.
- `benchmarks/ed/EDtriangle_symm_NEblock.py` is optional and can be slow because it recomputes ED.

## Runtime Contract

The Fortran executable uses fixed filenames in the current working directory. Always run it from an explicit run directory containing:

- `paramC_sets.txt`
- `confin.txt`
- `seeds.txt`

The executable writes outputs into that same current working directory. Keep generated outputs inside run directories; do not move fixed filename I/O casually unless the user explicitly asks for an I/O redesign.

Useful commands:

```bash
make print-config
make build
make run-example
make benchmark
```

Local runs use:

```bash
bash scripts/run_local.sh runs/examples/triangle_3x2 1
```

Server submissions should use an `yyk_` job name and one of the available queues:

```bash
sbatch -p node6348 --job-name yyk_triangle scripts/sbatch_triangle.sh runs/examples/triangle_3x2
sbatch -p fat6348 --job-name yyk_triangle_fat scripts/sbatch_triangle.sh runs/examples/triangle_3x2
```

## Physics And Observable Rules

Do not change observable normalization as a mechanical refactor. Treat changes to these outputs as physics changes that require benchmark and documentation updates:

- `density_up`
- `density_do`
- `num_up`
- `num_do`
- `kinetic`
- `doubleOcc`
- `squareOcc`
- `numsquare_up`
- `numsquare_do`

The program samples the `b` flavor explicitly and reconstructs the `c` flavor by complex conjugation. Preserve that convention unless the requested physics change explicitly says otherwise.

Important benchmark conversions:

```text
total_NE_DQMC = last(num_up) + last(num_do)
total_kinetic_DQMC = last(kinetic) * Lq
```

The ED reference script reports total particle number and total kinetic expectation. The DQMC `kinetic` output is normalized by `Lq` in `src/obser_equal.f90`.

## Build And Benchmark Notes

The active executable is built from the `SOURCES` list in `Makefile`. `globalK.f90` and `global_update.f90` are kept in `src/`, but they are not part of the active executable. Do not re-enable global update unless explicitly requested.

Default benchmark:

```bash
make benchmark
```

This is a fast fixture comparison and must not recompute ED. Do not run:

```bash
make benchmark-ed
```

unless the user explicitly asks for ED recomputation and accepts the Python dependencies and runtime cost.

## Verification Before Completion

Before claiming a code change is complete, run the relevant verification commands. For most changes:

```bash
make print-config
make build
python3 -m py_compile benchmarks/compare.py benchmarks/ed/EDtriangle_symm_NEblock.py
make benchmark
git diff --check
```

Run `make run-example` when the change affects build, runtime scripts, input layout, or Fortran execution. If sandboxed MPI fails with socket permission errors, report that exact failure and rerun with approved elevated permissions if needed.

If compiler, MPI, external libraries, or Python ED dependencies are unavailable, report the exact command and failure. Do not claim success from partial verification.
