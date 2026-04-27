# Agent Instructions

Start by reading:

- `README.md`
- `docs/physics.md`
- `docs/development.md`

The Fortran executable uses fixed filenames in the current working directory. Keep runs inside explicit run directories containing:

- `paramC_sets.txt`
- `confin.txt`
- `seeds.txt`

Do not change observable normalization as a mechanical refactor. Treat changes to `density_up`, `density_do`, `num_up`, `num_do`, `kinetic`, `doubleOcc`, `squareOcc`, `numsquare_up`, or `numsquare_do` as physics changes that require benchmark updates.

The program currently samples the `b` flavor explicitly and reconstructs the `c` flavor by complex conjugation. Preserve that convention unless the requested physics change explicitly says otherwise.

Do not re-enable `global_update.f90` or `globalK.f90` unless explicitly requested. The active executable is built from the object list in `Makefile`.

Before claiming a code change is complete, run the relevant verification commands:

```bash
make print-config
make build
python3 -m py_compile benchmarks/compare.py
make benchmark
```

If compiler, MPI, or external library availability prevents a command from running locally, report the exact command and failure.
