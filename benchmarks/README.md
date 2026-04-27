# Benchmarks

The default benchmark is a fast comparison path:

```bash
make benchmark
```

It reads `benchmarks/references/triangle_3x2.json` and fixture DQMC-style output files under `benchmarks/fixtures/triangle_3x2_mc_outputs/`.

The fixture validates comparison semantics. It is not a replacement for a production Monte Carlo run.

## Observable Conversions

For the current DQMC output:

```text
total_NE = last(num_up) + last(num_do)
total_kinetic = last(kinetic) * Lq
```

The ED script reports total particle number and total kinetic expectation. The DQMC `kinetic` file is normalized by `Lq` in `src/obser_equal.f90`.

## Optional ED Recompute

Run:

```bash
make benchmark-ed
```

This runs `benchmarks/ed/EDtriangle_symm_NEblock.py` from inside `benchmarks/ed/`. It requires QuSpin and Numba and can take much longer than the default benchmark.
