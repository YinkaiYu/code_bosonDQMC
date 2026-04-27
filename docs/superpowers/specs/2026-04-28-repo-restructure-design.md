# Repository Restructure Design

## Goal

Restructure the current bosonic DQMC repository into a Makefile-led, collaboration-friendly layout without changing the core Fortran algorithm. The new layout should make build, local execution, server submission, physics documentation, and benchmark comparison clear enough for both human collaborators and coding agents.

## Scope

This is a medium engineering refactor. It may break the old root-level and `test/` entrypoints, but it should preserve the numerical behavior of the Fortran code unless a change is explicitly required for relocation.

In scope:

- Move Fortran source files under `src/`.
- Replace the root build flow with a clear Makefile target set.
- Move runnable input examples under `runs/examples/`.
- Add local run and server submission scripts.
- Add default fast benchmark comparison against checked-in reference values.
- Add optional ED recomputation support.
- Add physics and contributor documentation.
- Preserve current fixed filename I/O inside each run directory.

Out of scope:

- Rewriting the DQMC update algorithm.
- Refactoring global runtime state in `CalcBasic`.
- Changing observable definitions or normalizations.
- Making ED recomputation mandatory for default tests.
- Migrating to `fpm`.

## Current Project Findings

The current repository is a flat Fortran project:

- `main.f90` is the executable entrypoint.
- `CalcBasic` reads `paramC_sets.txt` from the current working directory and stores global runtime parameters.
- `Fields_mod` reads `confin.txt` and `seeds.txt` from the current working directory.
- `FourierTrans_mod` writes observables such as `density_up`, `density_do`, `kinetic`, `num_up`, and `num_do` into the current working directory.
- `global_update.f90` and `globalK.f90` exist, but the main program currently leaves global update disabled.
- `Makefile` delegates the build to `Compile`, which hard-codes the object order.
- `test/` currently mixes example inputs, a local run helper, and a Slurm submission script.

The existing fixed filename I/O is part of the runtime contract. This refactor should document it and make run directories explicit instead of changing those filenames.

## Target Layout

```text
.
├── AGENTS.md
├── README.md
├── Makefile
├── src/
│   ├── main.f90
│   ├── calc_basic.f90
│   ├── lattice.f90
│   ├── fields.f90
│   └── ...
├── build/
│   ├── obj/
│   ├── mod/
│   └── bosonDQMC.out
├── runs/
│   └── examples/
│       └── triangle_3x2/
│           ├── paramC_sets.txt
│           ├── confin.txt
│           └── seeds.txt
├── scripts/
│   ├── run_local.sh
│   └── sbatch_triangle.sh
├── benchmarks/
│   ├── README.md
│   ├── compare.py
│   ├── fixtures/
│   │   └── triangle_3x2_mc_outputs/
│   ├── references/
│   │   └── triangle_3x2.json
│   └── ed/
│       └── EDtriangle_symm_NEblock.py
└── docs/
    ├── physics.md
    └── development.md
```

`build/` remains generated output and should be ignored by git. The committed example run directory is intentionally small and contains input files only.

## Build Design

The root `Makefile` remains the single primary build interface. It should:

- Detect the Fortran compiler and external library root using the current local/server conventions.
- Compile all source files from `src/`.
- Write object files to `build/obj/`.
- Write module files to `build/mod/`.
- Link the executable to `build/bosonDQMC.out`.
- Preserve the existing external library assumptions: `Modules/modules_90.a`, `MyEis/libeis.a`, `MyNag/libnag.a`, `MyLin/liblin.a`, and `Ran/libran.a`.

Required targets:

- `make build`: compile `build/bosonDQMC.out`.
- `make clean`: remove generated build artifacts.
- `make print-config`: print compiler, library root, flags, and executable path.
- `make run-example`: build and run the default example locally.
- `make benchmark`: run the fast benchmark comparison against reference values.
- `make benchmark-ed`: explicitly attempt ED recomputation.

The old `Compile` file can be removed once its object ordering is reproduced in the new Makefile.

## Runtime Design

The Fortran executable should still be launched from a run directory containing:

- `paramC_sets.txt`
- `confin.txt`
- `seeds.txt`

Local execution should use:

```bash
mpirun -np 1 /absolute/or/relative/path/to/build/bosonDQMC.out
```

`scripts/run_local.sh` should accept a run directory argument, defaulting to `runs/examples/triangle_3x2`, build if needed, then execute from inside that run directory so the current Fortran I/O behavior remains unchanged.

The server submission template should live at `scripts/sbatch_triangle.sh` and should:

- Use an `yyk_` job name.
- Support `node6348` and `fat6348` partitions.
- Run from a chosen run directory.
- Launch the same executable with `mpirun`.

## Benchmark Design

Default benchmark behavior should be fast and deterministic. It should compare DQMC output files against checked-in reference values derived from prior ED runs, without recomputing ED.

`benchmarks/references/triangle_3x2.json` should include:

- Lattice and Hamiltonian parameters: `Lx`, `Ly`, `t`, `U1`, `U2`, `beta`, `mu`.
- Reference observables from ED.
- Explicit normalization metadata for each observable.
- Tolerances for stochastic DQMC comparison.

`benchmarks/compare.py` should:

- Read DQMC output files from a run directory.
- Read the JSON reference.
- Compare the late-run or averaged DQMC values to the reference according to documented semantics.
- Print a concise pass/fail report.
- Avoid recomputing ED by default.

The default `make benchmark` target should validate the comparison path without requiring a long Monte Carlo run. It may compare the default run directory after `make run-example` has produced outputs, or it may compare a committed lightweight fixture under `benchmarks/fixtures/`. Fixture values should be labeled as script-validation data if they are not production MC statistics.

`benchmarks/ed/EDtriangle_symm_NEblock.py` should be included as an optional reference generator. `make benchmark-ed` or an explicit script flag should be required to run it. Missing Python dependencies such as `quspin` or `numba` should produce a clear skip/error message for the ED-specific command only, not for the default benchmark.

## Physics Documentation Design

`docs/physics.md` should be the canonical map between the physical model and code. It should include:

- The two-flavor triangular-lattice Bose-Hubbard Hamiltonian.
- The sign convention used by the code: `t = 1 > 0`, `U1 < 0`, `U2 > 0`.
- Continuous Hubbard-Stratonovich decoupling for the `U1` and `U2` terms.
- The fact that the code samples the `b` flavor explicitly and reconstructs the `c` flavor by complex conjugation.
- The equal-time Green function convention `G_ij = <b_i b_j^+>`.
- A table mapping physics symbols to code variables and files.
- A table mapping observables to Fortran fields and output files.
- Benchmark normalization notes, especially total versus per-flavor particle number and total kinetic energy versus kinetic energy density.

Required mapping examples:

```text
Physical symbol        Code name / location
Lx, Ly                 Nlx, Nly in CalcBasic, read from paramC_sets.txt
beta                   Beta in CalcBasic
Delta tau              Dtau = Beta / Ltrot
t                      RT in CalcBasic, currently set to 1.d0 in Params_set
U1, U2                 RU1, RU2 in CalcBasic, read from paramC_sets.txt
mu                     mu in CalcBasic, read from paramC_sets.txt
auxiliary fields       Conf%phi_list(ns, ii, nt) in Fields_mod
lattice bonds          Latt%L_bonds in MyLattice
b flavor Green matrix  Prop%Gr
c flavor Green matrix  dconjg(Prop%Gr)
```

Observable mapping should include:

```text
DQMC output file       Fortran field             Meaning
density_up            Obs%density_up            per-site density for b flavor
density_do            Obs%density_do            per-site density for c flavor
num_up                Obs%num_up                total b-flavor particle number
num_do                Obs%num_do                total c-flavor particle number
kinetic               Obs%kinetic               current code divides by Lq in Obs_equal_calc
doubleOcc             Obs%doubleOcc             per-site cross-flavor density product
squareOcc             Obs%squareOcc             per-site same-flavor square contribution
```

The benchmark documentation should state that the ED script reports total `NE = NE_b + NE_c` and total kinetic expectation, while DQMC may expose both per-flavor particle numbers and a kinetic observable normalized by `Lq`. The comparison script must make these conversions explicit.

## Contributor And Agent Documentation

`README.md` should explain:

- What the code simulates.
- How to build locally.
- How to run the default example locally.
- How to run the fast benchmark.
- How to optionally rerun ED.
- How to submit the server script.

`AGENTS.md` should explain:

- Start with `README.md`, `docs/physics.md`, and `docs/development.md`.
- Keep the executable running from run directories because of fixed filename I/O.
- Treat observable normalization changes as physics changes requiring benchmark updates.
- Do not re-enable global update unless explicitly requested.
- Prefer small, reviewable changes and run `make build` plus relevant benchmark commands before claiming success.

`docs/development.md` should document:

- Source/module layout.
- Makefile targets.
- Run directory contract.
- Generated files and cleanup expectations.
- How to add a new benchmark case.

## Error Handling

The build should fail clearly when no Fortran compiler is found. The Makefile should continue to support explicit `FC=...` and `LIB_ROOT=...` overrides.

Run scripts should fail with a clear message if:

- The requested run directory is missing.
- Required input files are absent.
- `mpirun` is not available.
- The executable cannot be built.

Benchmark scripts should fail clearly if required DQMC output files are missing. Optional ED recomputation should distinguish dependency failure from numerical comparison failure.

## Testing And Verification

Minimum verification after implementation:

- `make print-config`
- `make build`
- `make run-example` if local MPI and external libraries are available
- `make benchmark` after a run directory has generated output

If local compilation is impossible because compiler or external libraries are unavailable, the refactor should still be validated by:

- Checking that all Makefile source paths match files under `src/`.
- Running syntax checks for Python scripts with `python -m py_compile`.
- Running benchmark comparison against a small fixture if included.

## Migration Notes

The old root-level source paths and `test/` run scripts are allowed to break or disappear. The new documented entrypoints are the supported interface.

The Fortran code should continue to use fixed runtime filenames in the current working directory during this refactor. Any future change to configurable input/output paths should be handled as a separate design because it touches runtime behavior across several modules.
