# Development Guide

## Layout

- `src/` contains Fortran source files.
- `build/` contains generated objects, modules, and `bosonDQMC.out`.
- `runs/examples/` contains committed run input directories.
- `scripts/` contains local and server run helpers.
- `benchmarks/` contains comparison scripts, fixture outputs, ED references, and optional ED tooling.
- `docs/` contains physics and development documentation.

## Build Targets

```bash
make print-config
make build
make run-example
make benchmark
make benchmark-ed
make clean
```

`make build` compiles the active object list in `Makefile`. `globalK.f90` and `global_update.f90` are kept in `src/` but are not part of the active executable.

## Run Directory Contract

The executable reads these files from the current working directory:

- `paramC_sets.txt`
- `confin.txt`
- `seeds.txt`

The executable writes scalar observables and logs into the current working directory. Run from a dedicated directory to avoid mixing outputs from unrelated parameter sets.

## Adding A Run Case

Create a directory under `runs/examples/` or an untracked working run directory with:

- `paramC_sets.txt`
- `confin.txt`
- `seeds.txt`

Run it locally:

```bash
bash scripts/run_local.sh path/to/run_dir 1
```

## Adding A Benchmark Case

Create:

- `benchmarks/references/<case>.json`
- optional fixture outputs under `benchmarks/fixtures/<case>_mc_outputs/`
- optional ED params under `benchmarks/ed/params_<case>.txt`

The reference JSON must document whether each expected value is total, per flavor, per site, or divided by `Lq`.

## Generated Files

Generated build files are under `build/` and removed by:

```bash
make clean
```

Generated run outputs inside `runs/**` are ignored by git. Fixture outputs under `benchmarks/fixtures/` are committed when they are intentionally used by `make benchmark`.
