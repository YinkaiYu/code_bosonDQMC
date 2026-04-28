# Pole Diagnostics And Discrete HS Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add once-per-bin bosonic pole diagnostics to the continuous-HS baseline and build a discrete-HS branch that reuses the same diagnostic output contract.

**Architecture:** Add a shared pole-diagnostics path around the bin-level `Prop%Gr`, not inside the per-tau equal-time observable loop. Continuous HS computes `log_P_HS` from the Gaussian field sum; the discrete branch stores HS labels and maps labels to `eta(l)` and `gamma(l)` for propagation, acceptance, and log-weight diagnostics. MPI ranks are treated as independent Monte Carlo samples for pole diagnostics: diagnostics are not averaged across ranks, and rank 0 writes all rank-local samples in deterministic rank order once per bin.

**Tech Stack:** Fortran 90, MPI, BLAS/LAPACK through the existing Intel/MKL link, Python benchmark helper scripts and unittest, Makefile workflow.

---

## File Structure

Main continuous-HS changes:

- Modify: `src/fields.f90`
  Add `AuxConf%log_weight()` for continuous Gaussian HS.
- Create: `src/pole_diagnostics.f90`
  Compute `z_a`, `d_pole`, `X_pole`, `rho_G`, `smax_G`, and `log_weight` from `Prop%Gr`.
- Modify: `src/main.f90`
  Call pole diagnostics once per bin after `Sweep_local%sweep`.
- Modify: `Makefile`
  Compile `pole_diagnostics.f90` before `main.f90`.
- Create: `benchmarks/check_pole_diagnostics.py`
  Validate pole diagnostic output files.
- Modify: `benchmarks/test_compare.py`
  Unit-test the pole diagnostic checker.
- Modify: `.gitignore`
  Ignore generated pole diagnostic files under `runs/**`.
- Modify: `docs/physics.md`, `docs/development.md`, `benchmarks/README.md`
  Document continuous pole diagnostics.

Discrete branch changes in `.worktrees/discrete-hs-pole-diagnostics`:

- Modify: `src/operator_Hubbard.f90`
  Use discrete HS `alpha_disc`, map stored labels to `eta(l)` and `gamma(l)`.
- Modify: `src/fields.f90`
  Initialize and output discrete labels; compute discrete `log_P_HS`.
- Modify: `src/localU.f90`
  Replace continuous shifts with symmetric discrete label proposals.
- Modify: `docs/physics.md`, `docs/development.md`, `benchmarks/README.md`
  Document discrete HS constants, labels, and log-weight semantics.

---

### Task 1: Add Pole Diagnostic Checker Tests

**Files:**
- Modify: `benchmarks/test_compare.py`
- Create: `benchmarks/check_pole_diagnostics.py`

- [ ] **Step 1: Write the failing tests**

Add imports near the top of `benchmarks/test_compare.py`:

```python
import math
```

Add this constant after `COMPARE`:

```python
POLE_CHECK = REPO_ROOT / "benchmarks" / "check_pole_diagnostics.py"
```

Add these tests before `if __name__ == "__main__":`:

```python
    def test_pole_diagnostic_checker_accepts_consistent_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            (run_dir / "paramC_sets.txt").write_text(
                "-1.0 1.0 -2.5\n"
                "3 2 60 6.0\n"
                "3 2 60\n"
                "10 2 1 0.3\n"
                ".false. 0\n"
                ".false. 500 1.0 1.0\n"
                "2 0.1 0.0 0.0\n",
                encoding="utf-8",
            )
            z_row = " ".join(["1.0 0.0"] * 6)
            (run_dir / "info.txt").write_text("# Cores                                        : 2\n", encoding="utf-8")
            (run_dir / "pole_z").write_text(f"{z_row}\n{z_row}\n{z_row}\n{z_row}\n", encoding="utf-8")
            (run_dir / "pole_distance").write_text("1.0\n0.5\n0.1\n0.2\n", encoding="utf-8")
            (run_dir / "pole_x").write_text(
                "0.0\n0.3010299956639812\n1.0\n0.6989700043360187\n",
                encoding="utf-8",
            )
            (run_dir / "green_spectral_radius").write_text("1.0\n2.0\n10.0\n5.0\n", encoding="utf-8")
            (run_dir / "green_smax").write_text("1.0\n2.1\n11.0\n5.5\n", encoding="utf-8")
            (run_dir / "log_weight").write_text("-3.0\n-2.5\n-2.0\n-1.5\n", encoding="utf-8")

            result = subprocess.run(
                [sys.executable, str(POLE_CHECK), str(run_dir)],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=False,
            )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("Pole diagnostic files passed", result.stdout)

    def test_pole_diagnostic_checker_rejects_inconsistent_radius(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            (run_dir / "paramC_sets.txt").write_text(
                "-1.0 1.0 -2.5\n"
                "3 2 60 6.0\n"
                "3 2 60\n"
                "10 1 1 0.3\n"
                ".false. 0\n"
                ".false. 500 1.0 1.0\n"
                "2 0.1 0.0 0.0\n",
                encoding="utf-8",
            )
            (run_dir / "pole_z").write_text(" ".join(["1.0 0.0"] * 6) + "\n", encoding="utf-8")
            (run_dir / "pole_distance").write_text("0.1\n", encoding="utf-8")
            (run_dir / "pole_x").write_text("1.0\n", encoding="utf-8")
            (run_dir / "green_spectral_radius").write_text("9.0\n", encoding="utf-8")
            (run_dir / "green_smax").write_text("9.0\n", encoding="utf-8")
            (run_dir / "log_weight").write_text("-2.0\n", encoding="utf-8")

            result = subprocess.run(
                [sys.executable, str(POLE_CHECK), str(run_dir)],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=False,
            )

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("green_spectral_radius", result.stdout + result.stderr)
```

- [ ] **Step 2: Run tests and verify RED**

Run:

```bash
python3 -m unittest benchmarks/test_compare.py -v
```

Expected: fails because `benchmarks/check_pole_diagnostics.py` does not exist.

- [ ] **Step 3: Implement the checker**

Create `benchmarks/check_pole_diagnostics.py`:

```python
#!/usr/bin/env python3
"""Validate once-per-bin pole diagnostic files produced by the DQMC executable."""

from __future__ import annotations

import argparse
import math
from pathlib import Path


SCALAR_FILES = (
    "pole_distance",
    "pole_x",
    "green_spectral_radius",
    "green_smax",
    "log_weight",
)


def parse_params(path: Path) -> tuple[int, int]:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped[0].isalpha():
            continue
        rows.append(stripped.split())
        if len(rows) >= 4:
            break
    if len(rows) < 4:
        raise ValueError(f"could not read dimensions and Nbin from {path}")
    nlx = int(rows[1][0])
    nly = int(rows[1][1])
    nbin = int(rows[3][1])
    return nlx * nly, nbin


def parse_rank_count(run_dir: Path, explicit_ranks: int | None) -> int:
    if explicit_ranks is not None:
        if explicit_ranks <= 0:
            raise ValueError("--ranks must be positive")
        return explicit_ranks
    info = run_dir / "info.txt"
    if not info.exists():
        return 1
    for line in info.read_text(encoding="utf-8").splitlines():
        if "# Cores" in line:
            return int(line.split(":")[-1])
    return 1


def read_scalar(path: Path) -> list[float]:
    if not path.is_file():
        raise FileNotFoundError(f"missing diagnostic file: {path}")
    values = [float(line.split()[0]) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not values:
        raise ValueError(f"no values in diagnostic file: {path}")
    for value in values:
        if not math.isfinite(value):
            raise ValueError(f"non-finite value in {path}: {value}")
    return values


def read_pole_z(path: Path, ndim: int) -> list[list[float]]:
    if not path.is_file():
        raise FileNotFoundError(f"missing diagnostic file: {path}")
    rows: list[list[float]] = []
    expected = 2 * ndim
    for line_no, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        values = [float(item) for item in line.split()]
        if len(values) != expected:
            raise ValueError(f"{path}:{line_no} has {len(values)} columns, expected {expected}")
        if not all(math.isfinite(value) for value in values):
            raise ValueError(f"{path}:{line_no} contains a non-finite value")
        rows.append(values)
    if not rows:
        raise ValueError(f"no rows in diagnostic file: {path}")
    return rows


def validate(run_dir: Path, rtol: float, atol: float, explicit_ranks: int | None) -> None:
    ndim, nbin = parse_params(run_dir / "paramC_sets.txt")
    ranks = parse_rank_count(run_dir, explicit_ranks)
    expected_samples = nbin * ranks
    pole_z = read_pole_z(run_dir / "pole_z", ndim)
    scalars = {name: read_scalar(run_dir / name) for name in SCALAR_FILES}

    counts = {"pole_z": len(pole_z), **{name: len(values) for name, values in scalars.items()}}
    if len(set(counts.values())) != 1:
        rendered = ", ".join(f"{name}={count}" for name, count in sorted(counts.items()))
        raise ValueError(
            f"diagnostic sample counts must match across files: {rendered}"
        )
    sample_count = len(pole_z)
    if sample_count < expected_samples or sample_count % expected_samples != 0:
        rendered = ", ".join(f"{name}={count}" for name, count in sorted(counts.items()))
        raise ValueError(
            f"diagnostic sample counts must be a positive multiple of Nbin*ranks={expected_samples}: {rendered}"
        )

    for index, distance in enumerate(scalars["pole_distance"]):
        if distance <= 0.0:
            raise ValueError(f"pole_distance[{index}] must be positive, got {distance}")
        expected_x = -math.log10(distance)
        actual_x = scalars["pole_x"][index]
        if not math.isclose(actual_x, expected_x, rel_tol=rtol, abs_tol=atol):
            raise ValueError(f"pole_x[{index}]={actual_x} is inconsistent with pole_distance={distance}")
        expected_radius = 1.0 / distance
        actual_radius = scalars["green_spectral_radius"][index]
        if not math.isclose(actual_radius, expected_radius, rel_tol=rtol, abs_tol=atol):
            raise ValueError(
                f"green_spectral_radius[{index}]={actual_radius} is inconsistent with pole_distance={distance}"
            )
        if scalars["green_smax"][index] + atol < actual_radius:
            raise ValueError(f"green_smax[{index}] is smaller than green_spectral_radius[{index}]")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--ranks", type=int, help="Expected MPI rank count; defaults to info.txt or 1")
    parser.add_argument("--rtol", type=float, default=1e-8)
    parser.add_argument("--atol", type=float, default=1e-10)
    args = parser.parse_args()
    try:
        validate(args.run_dir, args.rtol, args.atol, args.ranks)
    except Exception as exc:
        print(f"ERROR: {exc}")
        return 1
    print(f"Pole diagnostic files passed: {args.run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run tests and verify GREEN**

Run:

```bash
python3 -m unittest benchmarks/test_compare.py -v
python3 -m py_compile benchmarks/check_pole_diagnostics.py
```

Expected: both commands pass.

- [ ] **Step 5: Commit**

Run:

```bash
git add benchmarks/test_compare.py benchmarks/check_pole_diagnostics.py
git commit -m "添加pole诊断输出检查脚本"
```

---

### Task 2: Implement Continuous-HS Pole Diagnostics On Main

**Files:**
- Modify: `src/fields.f90`
- Create: `src/pole_diagnostics.f90`
- Modify: `src/main.f90`
- Modify: `Makefile`

- [ ] **Step 1: Verify RED against the executable**

Run:

```bash
make run-example
python3 benchmarks/check_pole_diagnostics.py runs/examples/triangle_3x2
```

Expected: `make run-example` may pass, but the checker fails with `missing diagnostic file: .../pole_z`.

- [ ] **Step 2: Add continuous HS log weight helper**

In `src/fields.f90`, extend `type :: AuxConf`:

```fortran
    contains
        procedure   :: make       => AuxConf_make
        procedure   :: log_weight => AuxConf_log_weight
        final       ::              AuxConf_clear
```

Add this procedure before `conf_in`:

```fortran
    real(kind=8) function AuxConf_log_weight(this) result(log_weight)
        class(AuxConf), intent(in) :: this
        log_weight = -0.5d0 * sum(this%phi_list * this%phi_list)
        return
    end function AuxConf_log_weight
```

- [ ] **Step 3: Add the pole diagnostics module**

Create `src/pole_diagnostics.f90`:

```fortran
module PoleDiagnostics_mod
    use ProcessMatrix
    use DQMC_Model_mod
    implicit none
    private
    public :: PoleDiagnostics

    type :: PoleDiagnostics
    contains
        procedure :: write => PoleDiagnostics_write
    end type PoleDiagnostics

contains

    subroutine compute_eigenvalues(Gr, evals)
        complex(kind=8), dimension(Ndim, Ndim), intent(in) :: Gr
        complex(kind=8), dimension(Ndim), intent(out) :: evals
        complex(kind=8), dimension(Ndim, Ndim) :: amat
        complex(kind=8), dimension(1, 1) :: vl, vr
        complex(kind=8), dimension(1) :: work_query
        complex(kind=8), dimension(:), allocatable :: work
        real(kind=8), dimension(max(1, 2*Ndim)) :: rwork
        integer :: info, lwork

        amat = Gr
        call ZGEEV('N', 'N', Ndim, amat, Ndim, evals, vl, 1, vr, 1, &
                   work_query, -1, rwork, info)
        if (info .ne. 0) then
            write(6,*) 'ZGEEV workspace query failed in pole diagnostics, info =', info
            stop
        endif
        lwork = max(1, int(real(work_query(1))))
        allocate(work(lwork))
        amat = Gr
        call ZGEEV('N', 'N', Ndim, amat, Ndim, evals, vl, 1, vr, 1, &
                   work, lwork, rwork, info)
        deallocate(work)
        if (info .ne. 0) then
            write(6,*) 'ZGEEV failed in pole diagnostics, info =', info
            stop
        endif
        return
    end subroutine compute_eigenvalues

    subroutine compute_singular_values(Gr, svals)
        complex(kind=8), dimension(Ndim, Ndim), intent(in) :: Gr
        real(kind=8), dimension(Ndim), intent(out) :: svals
        complex(kind=8), dimension(Ndim, Ndim) :: amat
        complex(kind=8), dimension(1, 1) :: u, vt
        complex(kind=8), dimension(1) :: work_query
        complex(kind=8), dimension(:), allocatable :: work
        real(kind=8), dimension(max(1, 5*Ndim)) :: rwork
        integer :: info, lwork

        amat = Gr
        call ZGESVD('N', 'N', Ndim, Ndim, amat, Ndim, svals, u, 1, vt, 1, &
                    work_query, -1, rwork, info)
        if (info .ne. 0) then
            write(6,*) 'ZGESVD workspace query failed in pole diagnostics, info =', info
            stop
        endif
        lwork = max(1, int(real(work_query(1))))
        allocate(work(lwork))
        amat = Gr
        call ZGESVD('N', 'N', Ndim, Ndim, amat, Ndim, svals, u, 1, vt, 1, &
                    work, lwork, rwork, info)
        deallocate(work)
        if (info .ne. 0) then
            write(6,*) 'ZGESVD failed in pole diagnostics, info =', info
            stop
        endif
        return
    end subroutine compute_singular_values

    subroutine calc_local_diagnostics(Prop, zvals, scalar_local)
        class(Propagator), intent(in) :: Prop
        complex(kind=8), dimension(Ndim), intent(out) :: zvals
        real(kind=8), dimension(5), intent(out) :: scalar_local
        complex(kind=8), dimension(Ndim) :: evals
        real(kind=8), dimension(Ndim) :: svals
        real(kind=8) :: pole_distance, pole_x, rho_g, smax_g, log_weight
        integer :: ia

        call compute_eigenvalues(Prop%Gr, evals)
        call compute_singular_values(Prop%Gr, svals)

        do ia = 1, Ndim
            if (abs(evals(ia)) .gt. tiny(1.d0)) then
                zvals(ia) = dcmplx(1.d0, 0.d0) / evals(ia)
            else
                zvals(ia) = dcmplx(huge(1.d0), 0.d0)
            endif
        enddo

        pole_distance = minval(abs(zvals))
        if (pole_distance .gt. tiny(1.d0)) then
            pole_x = -log10(pole_distance)
        else
            pole_x = -log10(tiny(1.d0))
        endif
        rho_g = maxval(abs(evals))
        smax_g = maxval(svals)
        log_weight = Conf%log_weight() + 2.d0 * sum(log(max(svals, tiny(1.d0))))

        scalar_local(1) = pole_distance
        scalar_local(2) = pole_x
        scalar_local(3) = rho_g
        scalar_local(4) = smax_g
        scalar_local(5) = log_weight
        return
    end subroutine calc_local_diagnostics

    subroutine append_pole_z(z_collect)
        complex(kind=8), dimension(Ndim, ISIZE), intent(in) :: z_collect
        integer :: ia, irank_file

        open(unit=81, file='pole_z', status='unknown', action='write', position='append')
        do irank_file = 1, ISIZE
            do ia = 1, Ndim
                write(81, '(1X,ES24.16,1X,ES24.16)', advance='no') &
                    real(z_collect(ia, irank_file)), aimag(z_collect(ia, irank_file))
            enddo
            write(81,*)
        enddo
        close(81)
        return
    end subroutine append_pole_z

    subroutine append_scalar_column(filename, scalar_collect, column)
        character(len=*), intent(in) :: filename
        real(kind=8), dimension(5, ISIZE), intent(in) :: scalar_collect
        integer, intent(in) :: column
        integer :: irank_file

        open(unit=82, file=filename, status='unknown', action='write', position='append')
        do irank_file = 1, ISIZE
            write(82,*) scalar_collect(column, irank_file)
        enddo
        close(82)
        return
    end subroutine append_scalar_column

    subroutine PoleDiagnostics_write(this, Prop)
        include 'mpif.h'
        class(PoleDiagnostics), intent(inout) :: this
        class(Propagator), intent(in) :: Prop
        complex(kind=8), dimension(Ndim) :: zvals
        complex(kind=8), dimension(:,:), allocatable :: z_collect
        real(kind=8), dimension(5) :: scalar_local
        real(kind=8), dimension(:,:), allocatable :: scalar_collect

        call calc_local_diagnostics(Prop, zvals, scalar_local)

        allocate(z_collect(Ndim, ISIZE), scalar_collect(5, ISIZE))
        call MPI_GATHER(zvals, Ndim, MPI_complex16, z_collect, Ndim, &
                        MPI_complex16, 0, MPI_COMM_WORLD, IERR)
        call MPI_GATHER(scalar_local, 5, MPI_Real8, scalar_collect, 5, &
                        MPI_Real8, 0, MPI_COMM_WORLD, IERR)

        if (IRANK == 0) then
            call append_pole_z(z_collect)
            call append_scalar_column('pole_distance', scalar_collect, 1)
            call append_scalar_column('pole_x', scalar_collect, 2)
            call append_scalar_column('green_spectral_radius', scalar_collect, 3)
            call append_scalar_column('green_smax', scalar_collect, 4)
            call append_scalar_column('log_weight', scalar_collect, 5)
        endif
        deallocate(z_collect, scalar_collect)
        return
    end subroutine PoleDiagnostics_write
end module PoleDiagnostics_mod
```

Rank 0 appends `ISIZE` rows per bin to each diagnostic file, in rank order. Existing physics observables continue to use their current MPI averages; only pole diagnostics use all-rank sample pooling.

- [ ] **Step 4: Wire diagnostics into the build and main loop**

In `Makefile`, add the new source before `main.f90`:

```make
  fourier_trans.f90 \
  pole_diagnostics.f90 \
  main.f90
```

In `src/main.f90`, add:

```fortran
    use PoleDiagnostics_mod
```

Declare:

```fortran
    type(PoleDiagnostics) :: PoleDiag
```

After `call Fourier%preq(Obs_equal)`, add:

```fortran
        call PoleDiag%write(Prop)
```

This placement emits one pole diagnostic sample per bin.

- [ ] **Step 5: Verify GREEN**

Run:

```bash
make build
make run-example
python3 benchmarks/check_pole_diagnostics.py runs/examples/triangle_3x2
```

Expected: build and run pass; checker reports `Pole diagnostic files passed`.

- [ ] **Step 6: Commit**

Run:

```bash
git add src/fields.f90 src/pole_diagnostics.f90 src/main.f90 Makefile
git commit -m "添加连续场pole诊断输出"
```

---

### Task 3: Document Continuous Pole Diagnostics

**Files:**
- Modify: `.gitignore`
- Modify: `docs/physics.md`
- Modify: `docs/development.md`
- Modify: `benchmarks/README.md`

- [ ] **Step 1: Write the failing documentation/output check**

Run:

```bash
python3 benchmarks/check_pole_diagnostics.py runs/examples/triangle_3x2
git status --short
```

Expected before `.gitignore` update: the checker passes, and `git status --short` shows generated pole output files under `runs/examples/triangle_3x2`.

- [ ] **Step 2: Ignore generated pole output files**

Add to `.gitignore` after existing generated observable files:

```gitignore
/runs/**/pole_z
/runs/**/pole_distance
/runs/**/pole_x
/runs/**/green_spectral_radius
/runs/**/green_smax
/runs/**/log_weight
```

- [ ] **Step 3: Update physics documentation**

In `docs/physics.md`, add a section covering:

```text
Continuous HS alpha:
alpha_cont(U1) = sqrt(-2 * U1 * Dtau)
alpha_cont(U2) = i * sqrt(2 * U2 * Dtau)

Pole outputs:
pole_z: Re/Im pairs for z_a = 1 / mu_a(G)
pole_distance: min_a |z_a|
pole_x: -log10(pole_distance)
green_spectral_radius: max_a |mu_a(G)|
green_smax: largest singular value of G
log_weight: log_P_HS + 2 * sum log(s_a(G))
```

Also state that pole diagnostics do not average across MPI ranks. For `MPI_NP > 1`, each rank contributes one independent configuration sample per bin, and rank 0 writes `ISIZE` rows per bin in rank order. Use `MPI_NP=1` when a single Markov-chain time series is needed for tracking individual pole spikes.

- [ ] **Step 4: Update development and benchmark docs**

In `docs/development.md`, document that `make run-example` writes pole diagnostics into the run directory and that `benchmarks/check_pole_diagnostics.py <run_dir>` validates them.

In `benchmarks/README.md`, document that live benchmark physics pass/fail remains based on `num_up`, `num_do`, and `kinetic`; pole diagnostics are structural and analysis outputs.

- [ ] **Step 5: Verify docs and ignore behavior**

Run:

```bash
python3 benchmarks/check_pole_diagnostics.py runs/examples/triangle_3x2
git diff --check
git status --short
```

Expected: checker passes; no whitespace errors; generated pole files do not appear in git status.

- [ ] **Step 6: Commit**

Run:

```bash
git add .gitignore docs/physics.md docs/development.md benchmarks/README.md
git commit -m "记录连续场pole诊断输出"
```

---

### Task 4: Verify Main Continuous-HS Baseline

**Files:**
- Verification task only.

- [ ] **Step 1: Run standard verification**

Run:

```bash
make print-config
make build
python3 -m py_compile benchmarks/compare.py benchmarks/ed/EDtriangle_symm_NEblock.py benchmarks/run_dqmc_suite.py benchmarks/check_pole_diagnostics.py
python3 -m unittest benchmarks/test_compare.py -v
make check-fixtures
make run-example
python3 benchmarks/check_pole_diagnostics.py runs/examples/triangle_3x2
git diff --check
```

Expected: all commands pass. The Intel `-mkl` deprecation warning is acceptable if it remains the only compiler warning.

- [ ] **Step 2: Handle failures**

Expected: no fix commit is needed. If any command fails, stop at the first failing command, inspect the failure, and return to the task that introduced that behavior before continuing.

---

### Task 5: Bring Main Pole Diagnostics Into The Discrete Worktree

**Files:**
- Worktree: `.worktrees/discrete-hs-pole-diagnostics`

- [ ] **Step 1: Update the worktree from main**

Run from the main repository root:

```bash
git status --short
git -C .worktrees/discrete-hs-pole-diagnostics status --short
git -C .worktrees/discrete-hs-pole-diagnostics merge main
```

Expected: both status checks are clean before merge; merge succeeds with the continuous pole diagnostic commits.

- [ ] **Step 2: Verify inherited diagnostics in the worktree**

Run:

```bash
make -C .worktrees/discrete-hs-pole-diagnostics build
make -C .worktrees/discrete-hs-pole-diagnostics run-example
python3 .worktrees/discrete-hs-pole-diagnostics/benchmarks/check_pole_diagnostics.py .worktrees/discrete-hs-pole-diagnostics/runs/examples/triangle_3x2
```

Expected: inherited continuous diagnostics build and validate before discrete changes begin.

---

### Task 6: Convert The Worktree To Discrete HS

**Files:**
- Modify in `.worktrees/discrete-hs-pole-diagnostics`: `src/operator_Hubbard.f90`
- Modify in `.worktrees/discrete-hs-pole-diagnostics`: `src/fields.f90`
- Modify in `.worktrees/discrete-hs-pole-diagnostics`: `src/localU.f90`

- [ ] **Step 1: Write RED check for discrete mode marker**

Run in `.worktrees/discrete-hs-pole-diagnostics`:

```bash
make run-example
grep -q "Auxiliary field scheme.*discrete" runs/examples/triangle_3x2/info.txt
```

Expected: grep fails because the branch still reports no discrete scheme marker.

- [ ] **Step 2: Add discrete HS constants and label helpers**

In `src/operator_Hubbard.f90`, add these helpers in the module `contains` section before `opU_set`:

```fortran
    real(kind=8) function hs_eta(label) result(eta)
        real(kind=8), intent(in) :: label
        select case (nint(label))
        case (-1)
            eta = -sqrt(2.d0 * (3.d0 - sqrt(6.d0)))
        case (1)
            eta = sqrt(2.d0 * (3.d0 - sqrt(6.d0)))
        case (-2)
            eta = -sqrt(2.d0 * (3.d0 + sqrt(6.d0)))
        case (2)
            eta = sqrt(2.d0 * (3.d0 + sqrt(6.d0)))
        case default
            write(6,*) 'illegal discrete HS label for eta:', label, 'rank', IRANK
            stop
        end select
        return
    end function hs_eta

    real(kind=8) function hs_gamma(label) result(gamma)
        real(kind=8), intent(in) :: label
        select case (abs(nint(label)))
        case (1)
            gamma = 1.d0 + sqrt(6.d0) / 3.d0
        case (2)
            gamma = 1.d0 - sqrt(6.d0) / 3.d0
        case default
            write(6,*) 'illegal discrete HS label for gamma:', label, 'rank', IRANK
            stop
        end select
        return
    end function hs_gamma
```

Replace `opU_set` with the discrete coupling constants:

```fortran
    subroutine opU_set(this, RU)
        class(OperatorHubbard), intent(inout) :: this
        real(kind=8), intent(in) :: RU
        this%alpha = dcmplx( 0.d0, 0.d0 )
        if ( RU < -Zero ) this%alpha = dcmplx( sqrt(-RU * Dtau), 0.d0 )
        if ( RU >  Zero ) this%alpha = dcmplx( 0.d0, sqrt( RU * Dtau) )
        return
    end subroutine opU_set
```

Replace `opU_get_exp` with the discrete label mapping:

```fortran
    subroutine opU_get_exp(this, phi, nflag)
        class(OperatorHubbard), intent(inout) :: this
        integer, intent(in) :: nflag
        real(kind=8), intent(in) :: phi
        this%gaussian = dcmplx(hs_gamma(phi), 0.d0)
        this%expalpha = exp( this%alpha * hs_eta(phi) * nflag )
        return
    end subroutine opU_get_exp
```

The existing `opU_get_delta` then yields `gamma(new)/gamma(old)` and the correct rank-1 `Delta`.

- [ ] **Step 3: Store and initialize discrete labels**

In `src/fields.f90`, keep `phi_list` as `real(kind=8)` storage for compatibility but store only labels `-2.d0`, `-1.d0`, `1.d0`, `2.d0`.

Add this helper before `AuxConf_log_weight`:

```fortran
    real(kind=8) function log_gamma_label(label) result(log_gamma)
        real(kind=8), intent(in) :: label
        select case (abs(nint(label)))
        case (1)
            log_gamma = log(1.d0 + sqrt(6.d0) / 3.d0)
        case (2)
            log_gamma = log(1.d0 - sqrt(6.d0) / 3.d0)
        case default
            write(6,*) 'illegal discrete HS label for log gamma:', label, 'rank', IRANK
            stop
        end select
        return
    end function log_gamma_label
```

Replace `AuxConf_log_weight` with:

```fortran
    real(kind=8) function AuxConf_log_weight(this) result(log_weight)
        class(AuxConf), intent(in) :: this
        integer :: ns, ii, nt

        log_weight = 0.d0
        do nt = 1, Ltrot
            do ii = 1, Ndim
                do ns = 1, Naux
                    log_weight = log_weight + log_gamma_label(this%phi_list(ns, ii, nt))
                enddo
            enddo
        enddo
        return
    end function AuxConf_log_weight
```

Replace the body of `conf_set` with a symmetric random label draw:

```fortran
        do nt = 1, LtrotTherm
            do ii = 1, NdimTherm
                do ns = 1, Naux
                    select case (nranf(itmp, 4))
                    case (1)
                        phi_list(ns, ii, nt) = -2.d0
                    case (2)
                        phi_list(ns, ii, nt) = -1.d0
                    case (3)
                        phi_list(ns, ii, nt) = 1.d0
                    case (4)
                        phi_list(ns, ii, nt) = 2.d0
                    end select
                enddo
            enddo
        enddo
        return
```

The discrete branch ignores `iniType`, `iniAmpl`, and `iniBias` for generated initial configurations because labels are drawn from the four-point HS support.

- [ ] **Step 4: Replace local continuous shifts with discrete proposals**

In `src/localU.f90`, add this helper before `LocalU_metro`:

```fortran
    real(kind=8) function propose_discrete_label(old_label, iseed) result(new_label)
        real(kind=8), intent(in) :: old_label
        integer, intent(inout) :: iseed
        real(kind=8), dimension(4), parameter :: labels = (/ -2.d0, -1.d0, 1.d0, 2.d0 /)
        integer :: old_index, draw_index, proposed_index, ilabel

        old_index = 0
        do ilabel = 1, 4
            if (nint(old_label) == nint(labels(ilabel))) old_index = ilabel
        enddo
        if (old_index == 0) then
            write(6,*) 'illegal old discrete HS label:', old_label, 'rank', IRANK
            stop
        endif

        draw_index = nranf(iseed, 3)
        proposed_index = draw_index
        if (proposed_index >= old_index) proposed_index = proposed_index + 1
        new_label = labels(proposed_index)
        return
    end function propose_discrete_label
```

In both `LocalU_metro` and `LocalU_metro_therm`, replace:

```fortran
xflip = ranf(iseed)
Xdif = dble((xflip - 0.5) * abs(shiftLoc))
phi_new = phi_old + Xdif
```

with:

```fortran
phi_new = propose_discrete_label(phi_old, iseed)
```

Remove unused local variables `xflip` and `Xdif` from those two subroutines. Because proposals are symmetric over the three alternate labels, no proposal-probability ratio is added.

- [ ] **Step 5: Add the info marker**

In `src/calc_basic.f90`, write:

```fortran
write(50,*) 'Auxiliary field scheme                         : discrete'
```

near the other model/update metadata.

- [ ] **Step 6: Verify GREEN**

Run in `.worktrees/discrete-hs-pole-diagnostics`:

```bash
make build
make run-example
grep -q "Auxiliary field scheme.*discrete" runs/examples/triangle_3x2/info.txt
python3 benchmarks/check_pole_diagnostics.py runs/examples/triangle_3x2
```

Expected: all commands pass.

- [ ] **Step 7: Commit**

Run in `.worktrees/discrete-hs-pole-diagnostics`:

```bash
git add src/operator_Hubbard.f90 src/fields.f90 src/localU.f90 src/calc_basic.f90
git commit -m "实现离散辅助场局域更新"
```

---

### Task 7: Document Discrete HS Branch

**Files:**
- Modify in `.worktrees/discrete-hs-pole-diagnostics`: `docs/physics.md`
- Modify in `.worktrees/discrete-hs-pole-diagnostics`: `docs/development.md`
- Modify in `.worktrees/discrete-hs-pole-diagnostics`: `benchmarks/README.md`

- [ ] **Step 1: Update physics docs**

Document:

```text
alpha_disc(U1) = sqrt(-U1 * Dtau)
alpha_disc(U2) = i * sqrt(U2 * Dtau)
stored labels l in {-2, -1, 1, 2}
eta(l) mapping
gamma(l) mapping
log_P_HS = sum log(gamma(l))
```

State that `shiftLoc` and `shiftWarm` remain in the input file for format compatibility but are not used by the discrete proposal.

- [ ] **Step 2: Update development and benchmark docs**

Document that this branch uses discrete local flips and the same pole diagnostic output files as continuous main. State that full live benchmark is required for algorithm validation.

- [ ] **Step 3: Verify and commit**

Run in `.worktrees/discrete-hs-pole-diagnostics`:

```bash
git diff --check
git add docs/physics.md docs/development.md benchmarks/README.md
git commit -m "记录离散辅助场数据约定"
```

---

### Task 8: Verify Discrete HS Branch

**Files:**
- Verification task only.

- [ ] **Step 1: Run standard fast verification**

Run in `.worktrees/discrete-hs-pole-diagnostics`:

```bash
make print-config
make build
python3 -m py_compile benchmarks/compare.py benchmarks/ed/EDtriangle_symm_NEblock.py benchmarks/run_dqmc_suite.py benchmarks/check_pole_diagnostics.py
python3 -m unittest benchmarks/test_compare.py -v
make check-fixtures
make run-example
python3 benchmarks/check_pole_diagnostics.py runs/examples/triangle_3x2
git diff --check
```

Expected: all commands pass.

- [ ] **Step 2: Run the live benchmark**

Run in `.worktrees/discrete-hs-pole-diagnostics`:

```bash
make benchmark
```

Expected: the free analytic case and all four ED reference cases pass. Budget at least 15 minutes. If MPI fails with socket permission errors, rerun the same command with approved elevated permissions and report the exact first failure if it still fails.

- [ ] **Step 3: Handle failures**

Expected: no fix commit is needed. If any command fails, stop at the first failing command, inspect the failure, and return to the task that introduced that behavior before continuing.

---

### Task 9: Final Review

**Files:**
- Main repository root
- `.worktrees/discrete-hs-pole-diagnostics`

- [ ] **Step 1: Check branch states**

Run:

```bash
git status --short
git log --oneline -5
git -C .worktrees/discrete-hs-pole-diagnostics status --short
git -C .worktrees/discrete-hs-pole-diagnostics log --oneline -8
```

Expected: both working trees are clean except ignored generated files; logs show the continuous pole diagnostics on main and discrete HS commits on the branch.

- [ ] **Step 2: Summarize verification evidence**

Record:

```text
main verification commands and outcomes
discrete branch verification commands and outcomes
whether make benchmark passed for the discrete branch
location of worktree
new output file names
```

No code change is needed in this step.
