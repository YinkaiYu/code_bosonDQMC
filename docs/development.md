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
make benchmark-fast
make benchmark-dqmc
make check-fixtures
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

Output filenames are fixed by the Fortran executable. Keep generated outputs inside run directories unless the runtime I/O contract is being deliberately changed.

## Discrete-HS Input Contract

This branch uses discrete local auxiliary-field flips. Every stored field value
in `confin.txt` or `confout.txt` must be one of the legal labels:

```text
-2, -1, 1, 2
```

Old continuous-HS `confout` or `confin` files that contain Gaussian field
values are not valid restarts for this branch unless they are converted to
legal discrete labels before reuse. The restart reader loads raw values from
`confin.txt`; discrete label values are strictly checked when they are used by
the discrete HS eta/gamma, log-weight, or local-update routines. Invalid
restart labels will stop the run during those operations, not necessarily at
file-read time, so debug restart failures by checking every stored field value
against the legal label set above.

The fixed `paramC_sets.txt` format still includes `shiftLoc` and
`shiftWarm(1:2)` for compatibility with existing input files and scripts. The
discrete local proposal ignores those values.

## Adding A Run Case

Create a directory under `runs/examples/` or an untracked working run directory with:

- `paramC_sets.txt`
- `confin.txt`
- `seeds.txt`

Run it locally:

```bash
bash scripts/run_local.sh path/to/run_dir 1
```

The default example target writes scalar observables, logs, and pole diagnostic files into `runs/examples/triangle_3x2`:

```bash
make run-example
```

Validate the pole diagnostic file structure and cross-file consistency for any completed run with:

```bash
python3 benchmarks/check_pole_diagnostics.py <run_dir>
```

`benchmarks/check_pole_diagnostics.py` remains valid for this discrete-HS branch
because the diagnostic file names and row semantics are unchanged.

## Adding A Benchmark Case

Create:

- `benchmarks/references/<case>.json`
- fixture outputs under `benchmarks/fixtures/<case>_mc_outputs/`
- optional ED params under `benchmarks/ed/params_<case>.txt`

The reference JSON must include `dqmc_fixture` and document whether each expected value is total, per flavor, per site, or divided by `Lq`. After adding a case, run:

```bash
python3 -m unittest benchmarks/test_compare.py -v
make check-fixtures
```

For changes that substantively alter the DQMC algorithm, run the live benchmark:

```bash
make benchmark
```

`make benchmark` is an alias for `make benchmark-dqmc`. It uses fresh temporary run directories, so generated output from older runs cannot contaminate the comparison. The default suite is listed in `benchmarks/dqmc_suite.json` and includes one free analytic case plus all four ED reference cases from `temp/benchmark.txt`.

All live suite cases use `dtau = beta / Ltrot = 0.01`. Interacting live cases use `Nbin = 100000` and are compared to ED with block-estimated standard errors of the Monte Carlo mean. On the current WSL workstation with `MPI_NP=1`, the full suite was observed at `real 604.50` seconds, about 10 minutes 5 seconds; budget at least 15 minutes and do not shorten `Ltrot` or `Nbin` unless the benchmark definition is being deliberately changed.

## Generated Files

Generated build files are under `build/` and removed by:

```bash
make clean
```

Generated run outputs inside `runs/**` are ignored by git. Fixture outputs under `benchmarks/fixtures/` are committed when they are intentionally used by `make check-fixtures`.
