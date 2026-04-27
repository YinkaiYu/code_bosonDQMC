# Benchmarks

The default benchmark is the strict live DQMC-vs-ED suite:

```bash
make benchmark
```

`make benchmark` is an alias for `make benchmark-dqmc`. It runs the manifest in `benchmarks/dqmc_suite.json`, which includes the free-boson analytic case plus all four ED cases recorded in `/mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/temp/benchmark.txt`. All cases must pass.

All live benchmark inputs use `dtau = beta / Ltrot = 0.01`. The interacting cases use `Nbin = 100000`, `Nsweep = 1`, `shiftLoc = 1.5`, and warm-up enabled. On the current WSL workstation with `MPI_NP=1`, the full live suite was observed at `real 604.50` seconds, about 10 minutes 5 seconds; budget at least 15 minutes, and longer if the machine is busy.

The fast benchmark is a real DQMC run for the no-interaction `U1=U2=0` case:

```bash
make benchmark-fast
```

It compares against the analytic reference in `benchmarks/dqmc_references/triangle_3x2_free_beta3_mu-2.5.json`. This is useful for quick executable, Green-function, particle-number, and kinetic-normalization checks, but it does not exercise the interacting Monte Carlo path.

`make check-fixtures` is a fast comparison-script check over checked-in fixture files. It discovers every `*.json` file under `benchmarks/references/`. Each reference points to a fixture DQMC-style output directory under `benchmarks/fixtures/`. This validates parsing and documented DQMC-to-ED normalization semantics, but it is not physics validation and must not replace `make benchmark` for algorithm changes.

For one-off debugging of a single case, override both `DQMC_BENCHMARK_INPUT_DIR` and `DQMC_BENCHMARK_REFERENCE`. Do not use a single-case override as the final benchmark evidence for substantive algorithm changes.

## Live Pass/Fail Semantics

For an interacting live DQMC observable, `benchmarks/compare.py` first converts each bin to the ED scale. It then groups the 100000 samples into 10 blocks of 10000 samples. The reported `stderr` is the standard error of the Monte Carlo mean:

```text
stderr = std(block_means) / sqrt(number_of_blocks)
z = (DQMC_block_mean - ED_value) / stderr
```

The field `stderr_tolerance` is the allowed number of standard errors of the mean, not the standard deviation of individual Monte Carlo samples. The current references use `stderr_tolerance = 3.0`; a case passes when:

```text
abs(DQMC_block_mean - ED_value) <= max(atol, rtol * abs(ED_value), stderr_tolerance * stderr)
```

The free-boson case is deterministic and uses a tight absolute tolerance instead of block statistics.

## Observable Conversions

For the current DQMC output:

```text
total_NE = last(num_up) + last(num_do)
total_kinetic = last(kinetic) * Lq
```

The ED script accumulates `NE = NE_b + NE_c` and its kinetic operator includes both b and c hopping layers. The DQMC `kinetic` file includes both flavors and is divided by `Lq` in `src/obser_equal.f90`.

Interacting real DQMC benchmarks use sample means instead of the last bin:

```text
total_NE = mean(num_up) + mean(num_do)
total_kinetic = mean(kinetic) * Lq
```

## Adding A Case

Add one JSON file under `benchmarks/references/` and one fixture directory under `benchmarks/fixtures/`. The reference JSON must include:

- `case`: stable case identifier, matching the JSON filename stem
- `dqmc_fixture`: repository-relative path to the fixture output directory
- `parameters`: `Lx`, `Ly`, `t`, `U1`, `U2`, `beta`, `mu`
- `observables`: at least the particle and kinetic observables needed by the case, with explicit operation names such as `last`, `mean`, `sum_mean`, or `mean_times_lq`

The fixture directory must contain scalar output files with DQMC names, currently `num_up`, `num_do`, and `kinetic`.

Run the suite after adding a case:

```bash
python3 -m unittest benchmarks/test_compare.py -v
make check-fixtures
```

If the case is meant to exercise the live DQMC executable, add a JSON under `benchmarks/dqmc_references/`, add or reuse an input directory under `runs/benchmarks/`, wire it through `benchmarks/dqmc_suite.json`, and run `make benchmark`.

## Optional ED Recompute

Run:

```bash
make benchmark-ed
```

This runs `benchmarks/ed/EDtriangle_symm_NEblock.py` from inside `benchmarks/ed/`. It requires QuSpin and Numba and can take much longer than the default benchmark.
