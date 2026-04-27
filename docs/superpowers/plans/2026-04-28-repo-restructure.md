# Repository Restructure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restructure the bosonic DQMC repository into a Makefile-led, collaboration-friendly layout with documented physics mappings, explicit run directories, and fast benchmark comparison.

**Architecture:** Keep the Fortran numerical code behavior unchanged while moving source files into `src/` and generated build products into `build/`. Preserve fixed filename runtime I/O by running the executable from a run directory. Add Python benchmark tooling that compares DQMC output semantics to ED reference values without recomputing ED by default.

**Tech Stack:** Fortran 90, MPI, Intel/MPI Fortran-compatible Makefile, shell scripts, Python 3 standard library, optional QuSpin/Numba for ED recomputation.

---

## File Structure Map

Create or modify these files:

- Modify: `Makefile` - replace `Compile` delegation with direct Makefile build, run, and benchmark targets.
- Modify: `.gitignore` - ignore `build/` and generated run outputs.
- Create: `AGENTS.md` - agent-facing repository instructions.
- Create: `README.md` - human-facing build, run, benchmark, and server usage guide.
- Create: `docs/physics.md` - physical model to code variable mapping.
- Create: `docs/development.md` - engineering layout and workflow guide.
- Create: `scripts/run_local.sh` - local `mpirun -np 1` wrapper.
- Create: `scripts/sbatch_triangle.sh` - Slurm submission template for `node6348` and `fat6348`.
- Create: `benchmarks/README.md` - benchmark semantics and ED notes.
- Create: `benchmarks/compare.py` - fast DQMC-vs-reference comparison.
- Create: `benchmarks/references/triangle_3x2.json` - ED reference metadata.
- Create: `benchmarks/fixtures/triangle_3x2_mc_outputs/num_up` - fixture DQMC b-flavor particle number.
- Create: `benchmarks/fixtures/triangle_3x2_mc_outputs/num_do` - fixture DQMC c-flavor particle number.
- Create: `benchmarks/fixtures/triangle_3x2_mc_outputs/kinetic` - fixture DQMC kinetic density.
- Create: `benchmarks/ed/params_triangle_3x2.txt` - optional ED parameter file.
- Create: `benchmarks/ed/EDtriangle_symm_NEblock.py` - copy from `/mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/temp/EDtriangle_symm_NEblock.py`.
- Move: `calc_basic.f90`, `lattice.f90`, `fields.f90`, `non_interact.f90`, `operator_Hubbard.f90`, `model.f90`, `process_matrix.f90`, `multiply.f90`, `obser_equal.f90`, `obser_tau.f90`, `stabilization.f90`, `localU.f90`, `dynamics.f90`, `local_sweep.f90`, `fourier_trans.f90`, `main.f90`, `globalK.f90`, `global_update.f90` into `src/`.
- Move: `test/paramC_sets.txt`, `test/confin.txt`, and `test/seeds.txt` into `runs/examples/triangle_3x2/`.
- Delete: `Compile`, `auto.sh`, `test/auto.sh`, and `test/dqmc`.

---

### Task 1: Move Source And Example Inputs

**Files:**
- Create: `src/`
- Create: `runs/examples/triangle_3x2/`
- Move: root `*.f90` files to `src/`
- Move: `test/paramC_sets.txt`, `test/confin.txt`, `test/seeds.txt` to `runs/examples/triangle_3x2/`
- Delete: `Compile`
- Delete: `auto.sh`
- Delete: `test/auto.sh`
- Delete: `test/dqmc`

- [ ] **Step 1: Create target directories**

Run:

```bash
mkdir -p src runs/examples/triangle_3x2 scripts benchmarks/fixtures/triangle_3x2_mc_outputs benchmarks/references benchmarks/ed docs
```

Expected: command exits with status 0.

- [ ] **Step 2: Move Fortran sources into `src/`**

Run:

```bash
git mv calc_basic.f90 dynamics.f90 fields.f90 fourier_trans.f90 globalK.f90 global_update.f90 lattice.f90 localU.f90 local_sweep.f90 main.f90 model.f90 multiply.f90 non_interact.f90 obser_equal.f90 obser_tau.f90 operator_Hubbard.f90 process_matrix.f90 stabilization.f90 src/
```

Expected: command exits with status 0.

- [ ] **Step 3: Move runnable example inputs**

Run:

```bash
git mv test/paramC_sets.txt runs/examples/triangle_3x2/paramC_sets.txt
git mv test/confin.txt runs/examples/triangle_3x2/confin.txt
git mv test/seeds.txt runs/examples/triangle_3x2/seeds.txt
```

Expected: each command exits with status 0.

- [ ] **Step 4: Remove legacy entrypoints**

Run:

```bash
git rm Compile auto.sh test/auto.sh test/dqmc
```

Expected: command exits with status 0.

- [ ] **Step 5: Verify the moved file inventory**

Run:

```bash
find src -maxdepth 1 -type f -name '*.f90' | sort
find runs/examples/triangle_3x2 -maxdepth 1 -type f | sort
```

Expected output includes all 18 Fortran files under `src/` and the three required run input files under `runs/examples/triangle_3x2/`.

- [ ] **Step 6: Commit the move**

Run:

```bash
git add src runs
git commit -m "重排源码和示例输入目录"
```

Expected: git creates a commit containing only source/input moves and legacy script deletions.

---

### Task 2: Replace The Build System With A Direct Makefile

**Files:**
- Modify: `Makefile`
- Modify: `.gitignore`

- [ ] **Step 1: Replace `Makefile` with the direct build implementation**

Write `Makefile` exactly as follows:

```makefile
SERVER_LIB_ROOT := /home/zxli_1/Lib_90_new
LOCAL_LIB_ROOT  := /home/yyk/Lib_90_new

DETECTED_LIB_ROOT := $(firstword $(wildcard $(SERVER_LIB_ROOT) $(LOCAL_LIB_ROOT)))
LIB_ROOT ?= $(if $(DETECTED_LIB_ROOT),$(DETECTED_LIB_ROOT),$(SERVER_LIB_ROOT))

MPIIFORT := $(shell command -v mpiifort 2>/dev/null)
MPIIFX   := $(shell command -v mpiifx   2>/dev/null)
MPIFORT  := $(shell command -v mpifort  2>/dev/null)
MPIF90   := $(shell command -v mpif90   2>/dev/null)
GFORTRAN := $(shell command -v gfortran 2>/dev/null)
IFORT    := $(shell command -v ifort    2>/dev/null)
IFX      := $(shell command -v ifx      2>/dev/null)

ifneq ($(filter command% environment,$(origin FC)),)
  FC_CMD := $(FC)
else ifneq ($(MPIIFORT),)
  ifneq ($(IFORT),)
    FC_CMD := $(MPIIFORT) -fc=ifort
  else ifneq ($(IFX),)
    FC_CMD := $(MPIIFORT) -fc=ifx
  else
    FC_CMD := $(MPIIFORT)
  endif
else ifneq ($(MPIIFX),)
  ifneq ($(IFX),)
    FC_CMD := $(MPIIFX) -fc=ifx
  else
    FC_CMD := $(MPIIFX)
  endif
else ifneq ($(MPIFORT),)
  FC_CMD := $(MPIFORT)
else ifneq ($(MPIF90),)
  FC_CMD := $(MPIF90)
else ifneq ($(GFORTRAN),)
  FC_CMD := $(GFORTRAN)
endif

ifneq ($(filter clean,$(MAKECMDGOALS)),clean)
  ifneq ($(filter help,$(MAKECMDGOALS)),help)
    ifeq ($(strip $(FC_CMD)),)
      $(error No Fortran compiler found. Set FC explicitly, for example: make FC=mpiifort)
    endif
  endif
endif

SRC_DIR := src
BUILD_DIR := build
OBJ_DIR := $(BUILD_DIR)/obj
MOD_DIR := $(BUILD_DIR)/mod
TARGET := $(BUILD_DIR)/bosonDQMC.out

FFLAGS ?= -O3 -I$(LIB_ROOT)/Modules -I$(MOD_DIR)
SUFFIX ?=
LDFLAGS ?= -mkl
PYTHON ?= python3
MPI_NP ?= 1
RUN_DIR ?= runs/examples/triangle_3x2
BENCHMARK_RUN_DIR ?= benchmarks/fixtures/triangle_3x2_mc_outputs
BENCHMARK_REFERENCE ?= benchmarks/references/triangle_3x2.json

ifneq (,$(findstring gfortran,$(FC_CMD)))
  MOD_FLAG ?= -J$(MOD_DIR)
else
  MOD_FLAG ?= -module $(MOD_DIR)
endif

LDLIBS := $(LIB_ROOT)/Modules/modules_90.a \
          $(LIB_ROOT)/MyEis/libeis.a \
          $(LIB_ROOT)/MyNag/libnag.a \
          $(LIB_ROOT)/MyLin/liblin.a \
          $(LIB_ROOT)/Ran/libran.a

SOURCES := \
  calc_basic.f90 \
  lattice.f90 \
  fields.f90 \
  non_interact.f90 \
  operator_Hubbard.f90 \
  model.f90 \
  process_matrix.f90 \
  multiply.f90 \
  obser_equal.f90 \
  obser_tau.f90 \
  stabilization.f90 \
  localU.f90 \
  dynamics.f90 \
  local_sweep.f90 \
  fourier_trans.f90 \
  main.f90

OBJECTS := $(addprefix $(OBJ_DIR)/,$(SOURCES:.f90=.o))

.PHONY: all build clean print-config run-example benchmark benchmark-ed help

all: build

build: $(TARGET)

$(TARGET): $(OBJECTS) | $(BUILD_DIR)
	$(FC_CMD) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC_CMD) $(SUFFIX) -c $(FFLAGS) $(MOD_FLAG) -o $@ $<

$(BUILD_DIR) $(OBJ_DIR) $(MOD_DIR):
	mkdir -p $@

print-config:
	@echo "FC=$(FC_CMD)"
	@echo "LIB_ROOT=$(LIB_ROOT)"
	@echo "FFLAGS=$(FFLAGS)"
	@echo "LDFLAGS=$(LDFLAGS)"
	@echo "MOD_FLAG=$(MOD_FLAG)"
	@echo "TARGET=$(TARGET)"
	@echo "RUN_DIR=$(RUN_DIR)"
	@echo "BENCHMARK_RUN_DIR=$(BENCHMARK_RUN_DIR)"
	@echo "BENCHMARK_REFERENCE=$(BENCHMARK_REFERENCE)"

run-example: build
	bash scripts/run_local.sh $(RUN_DIR) $(MPI_NP)

benchmark:
	$(PYTHON) benchmarks/compare.py --reference $(BENCHMARK_REFERENCE) --run-dir $(BENCHMARK_RUN_DIR)

benchmark-ed:
	cd benchmarks/ed && cp params_triangle_3x2.txt params.txt && $(PYTHON) EDtriangle_symm_NEblock.py

clean:
	rm -rf $(BUILD_DIR)
	rm -f *.mod *.o *.lst *.opt-report bosonDQMC.out

help:
	@echo "Targets:"
	@echo "  make build          Build build/bosonDQMC.out"
	@echo "  make run-example    Run the default example with mpirun"
	@echo "  make benchmark      Compare fixture or run output to ED reference"
	@echo "  make benchmark-ed   Optionally recompute ED reference in benchmarks/ed"
	@echo "  make clean          Remove generated build artifacts"
	@echo "  make print-config   Print compiler and path configuration"
```

- [ ] **Step 2: Update `.gitignore`**

Write `.gitignore` exactly as follows:

```gitignore
.codex
build/
*.o
*.mod
*.lst
*.opt-report
bosonDQMC.out

/runs/**/info.txt
/runs/**/confout.txt
/runs/**/density_up
/runs/**/density_do
/runs/**/kinetic
/runs/**/doubleOcc
/runs/**/squareOcc
/runs/**/num_up
/runs/**/num_do
/runs/**/numsquare_up
/runs/**/numsquare_do
/runs/**/den_upup_sub*
/runs/**/den_dodo_sub*
/runs/**/den_updo
/benchmarks/ed/params.txt
/benchmarks/ed/results.txt
/benchmarks/ed/calc.log
```

- [ ] **Step 3: Verify Makefile path references**

Run:

```bash
make print-config
```

Expected: output includes `TARGET=build/bosonDQMC.out`, `RUN_DIR=runs/examples/triangle_3x2`, and `BENCHMARK_REFERENCE=benchmarks/references/triangle_3x2.json`.

- [ ] **Step 4: Verify the build graph reaches files under `src/`**

Run:

```bash
make -n build
```

Expected: dry-run compile commands reference `src/calc_basic.f90`, `src/lattice.f90`, and `src/main.f90`; link command writes `build/bosonDQMC.out`.

- [ ] **Step 5: Commit the build system change**

Run:

```bash
git add Makefile .gitignore
git commit -m "改用Makefile直接构建"
```

Expected: git creates a commit containing `Makefile` and `.gitignore` changes.

---

### Task 3: Add Local And Server Run Scripts

**Files:**
- Create: `scripts/run_local.sh`
- Create: `scripts/sbatch_triangle.sh`

- [ ] **Step 1: Create `scripts/run_local.sh`**

Write `scripts/run_local.sh` exactly as follows:

```bash
#!/usr/bin/env bash
set -euo pipefail

run_dir="${1:-runs/examples/triangle_3x2}"
np="${2:-${MPI_NP:-1}}"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exe="${repo_root}/build/bosonDQMC.out"

case "${run_dir}" in
  /*) ;;
  *) run_dir="${repo_root}/${run_dir}" ;;
esac

if [[ ! -d "${run_dir}" ]]; then
  echo "Run directory does not exist: ${run_dir}" >&2
  exit 2
fi

for required in paramC_sets.txt confin.txt seeds.txt; do
  if [[ ! -f "${run_dir}/${required}" ]]; then
    echo "Missing required input file: ${run_dir}/${required}" >&2
    exit 2
  fi
done

if ! command -v mpirun >/dev/null 2>&1; then
  echo "mpirun is not available on PATH" >&2
  exit 2
fi

make -C "${repo_root}" build

cd "${run_dir}"
exec mpirun -np "${np}" "${exe}"
```

- [ ] **Step 2: Create `scripts/sbatch_triangle.sh`**

Write `scripts/sbatch_triangle.sh` exactly as follows:

```bash
#!/usr/bin/env bash
#SBATCH -J yyk_triangle
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p node6348

set -euo pipefail

run_dir="${1:-${RUN_DIR:-runs/examples/triangle_3x2}}"
np="${SLURM_NTASKS:-1}"

repo_root="${SLURM_SUBMIT_DIR:-$(pwd)}"
exe="${repo_root}/build/bosonDQMC.out"

case "${run_dir}" in
  /*) ;;
  *) run_dir="${repo_root}/${run_dir}" ;;
esac

if [[ ! "${SLURM_JOB_NAME:-yyk_triangle}" == yyk_* ]]; then
  echo "Slurm job name should start with yyk_: ${SLURM_JOB_NAME}" >&2
  exit 2
fi

if [[ ! -d "${run_dir}" ]]; then
  echo "Run directory does not exist: ${run_dir}" >&2
  exit 2
fi

for required in paramC_sets.txt confin.txt seeds.txt; do
  if [[ ! -f "${run_dir}/${required}" ]]; then
    echo "Missing required input file: ${run_dir}/${required}" >&2
    exit 2
  fi
done

make -C "${repo_root}" build

cd "${run_dir}"
exec mpirun -np "${np}" "${exe}"
```

- [ ] **Step 3: Mark scripts executable**

Run:

```bash
chmod +x scripts/run_local.sh scripts/sbatch_triangle.sh
```

Expected: command exits with status 0.

- [ ] **Step 4: Validate script syntax**

Run:

```bash
bash -n scripts/run_local.sh
bash -n scripts/sbatch_triangle.sh
```

Expected: both commands exit with status 0.

- [ ] **Step 5: Commit run scripts**

Run:

```bash
git add scripts/run_local.sh scripts/sbatch_triangle.sh
git commit -m "添加本地和服务器运行脚本"
```

Expected: git creates a commit containing both scripts with executable mode.

---

### Task 4: Add Fast Benchmark Comparison And Optional ED Inputs

**Files:**
- Create: `benchmarks/compare.py`
- Create: `benchmarks/references/triangle_3x2.json`
- Create: `benchmarks/fixtures/triangle_3x2_mc_outputs/num_up`
- Create: `benchmarks/fixtures/triangle_3x2_mc_outputs/num_do`
- Create: `benchmarks/fixtures/triangle_3x2_mc_outputs/kinetic`
- Create: `benchmarks/ed/params_triangle_3x2.txt`
- Create: `benchmarks/ed/EDtriangle_symm_NEblock.py`

- [ ] **Step 1: Create `benchmarks/compare.py`**

Write `benchmarks/compare.py` exactly as follows:

```python
#!/usr/bin/env python3
"""Compare DQMC scalar output files against documented ED reference values."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable


def read_series(path: Path) -> list[float]:
    if not path.exists():
        raise FileNotFoundError(f"missing DQMC output file: {path}")
    values: list[float] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        values.append(float(stripped.split()[0]))
    if not values:
        raise ValueError(f"no numeric values found in {path}")
    return values


def last_values(run_dir: Path, files: Iterable[str]) -> list[float]:
    return [read_series(run_dir / name)[-1] for name in files]


def compute_actual(run_dir: Path, operation: str, files: list[str], lq: int) -> float:
    values = last_values(run_dir, files)
    if operation == "last":
        if len(values) != 1:
            raise ValueError("operation 'last' requires exactly one file")
        return values[0]
    if operation == "sum_last":
        return sum(values)
    if operation == "last_times_lq":
        if len(values) != 1:
            raise ValueError("operation 'last_times_lq' requires exactly one file")
        return values[0] * float(lq)
    raise ValueError(f"unsupported comparison operation: {operation}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path, help="JSON reference file")
    parser.add_argument("--run-dir", required=True, type=Path, help="Directory containing DQMC scalar output files")
    args = parser.parse_args()

    reference = json.loads(args.reference.read_text(encoding="utf-8"))
    params = reference["parameters"]
    lq = int(params["Lx"]) * int(params["Ly"])

    failures: list[str] = []
    print(f"Reference: {args.reference}")
    print(f"DQMC run directory: {args.run_dir}")
    print("")

    for name, spec in reference["observables"].items():
        dqmc = spec["dqmc"]
        actual = compute_actual(args.run_dir, dqmc["operation"], dqmc["files"], lq)
        expected = float(spec["value"])
        atol = float(spec.get("atol", 0.0))
        rtol = float(spec.get("rtol", 0.0))
        passed = math.isclose(actual, expected, abs_tol=atol, rel_tol=rtol)
        status = "PASS" if passed else "FAIL"
        print(
            f"{status} {name}: actual={actual:.16g} expected={expected:.16g} "
            f"abs_diff={abs(actual - expected):.3g} atol={atol:.3g} rtol={rtol:.3g}"
        )
        if not passed:
            failures.append(name)

    if failures:
        print("")
        print("Failed observables: " + ", ".join(failures))
        return 1

    print("")
    print("All benchmark comparisons passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: Create the ED reference JSON**

Write `benchmarks/references/triangle_3x2.json` exactly as follows:

```json
{
  "case": "triangle_3x2_beta6_mu-2.5_u1_0_u2_1",
  "source": "Reference values copied from temp/benchmark.txt; default benchmark uses fixture data to validate comparison semantics.",
  "parameters": {
    "Lx": 3,
    "Ly": 2,
    "t": 1.0,
    "U1": 0.0,
    "U2": 1.0,
    "beta": 6.0,
    "mu": -2.5
  },
  "normalization": {
    "ed_total_NE": "ED reports total particle number NE_b + NE_c.",
    "dqmc_total_NE": "DQMC total particle number is num_up + num_do.",
    "ed_total_kinetic": "ED reports total kinetic expectation.",
    "dqmc_kinetic": "DQMC file kinetic is accumulated with division by Lq in Obs_equal_calc, so compare kinetic * Lq to ED total kinetic."
  },
  "observables": {
    "total_NE": {
      "value": 0.0008330248173669302,
      "atol": 1e-12,
      "rtol": 0.0,
      "dqmc": {
        "operation": "sum_last",
        "files": ["num_up", "num_do"]
      }
    },
    "total_kinetic": {
      "value": -0.001659465125418983,
      "atol": 1e-12,
      "rtol": 0.0,
      "dqmc": {
        "operation": "last_times_lq",
        "files": ["kinetic"]
      }
    }
  }
}
```

- [ ] **Step 3: Create fixture DQMC output files**

Write `benchmarks/fixtures/triangle_3x2_mc_outputs/num_up` exactly as follows:

```text
0.0004165124086834651
```

Write `benchmarks/fixtures/triangle_3x2_mc_outputs/num_do` exactly as follows:

```text
0.0004165124086834651
```

Write `benchmarks/fixtures/triangle_3x2_mc_outputs/kinetic` exactly as follows:

```text
-0.00027657752090316385
```

- [ ] **Step 4: Add the optional ED parameter file**

Write `benchmarks/ed/params_triangle_3x2.txt` exactly as follows:

```text
Lx = 3
Ly = 2
t = 1.0
U1 = 0.0
U2 = 1.0
beta = 6.0
mu = -2.5
```

- [ ] **Step 5: Copy the optional ED script into the repository**

Run:

```bash
cp /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/temp/EDtriangle_symm_NEblock.py benchmarks/ed/EDtriangle_symm_NEblock.py
```

Expected: command exits with status 0 and creates `benchmarks/ed/EDtriangle_symm_NEblock.py`.

- [ ] **Step 6: Normalize the ED script shebang**

Change only the first line of `benchmarks/ed/EDtriangle_symm_NEblock.py` to:

```python
#!/usr/bin/env python3
```

Expected: the rest of the ED script remains copied from the external reference.

- [ ] **Step 7: Validate Python syntax and benchmark behavior**

Run:

```bash
python3 -m py_compile benchmarks/compare.py benchmarks/ed/EDtriangle_symm_NEblock.py
make benchmark
```

Expected: Python compilation exits with status 0. `make benchmark` prints `PASS total_NE`, `PASS total_kinetic`, and `All benchmark comparisons passed.`

- [ ] **Step 8: Commit benchmark tooling**

Run:

```bash
git add benchmarks
git commit -m "添加快速benchmark比较工具"
```

Expected: git creates a commit containing benchmark scripts, fixture outputs, reference values, and optional ED inputs.

---

### Task 5: Add Physics, Development, And Agent Documentation

**Files:**
- Create: `README.md`
- Create: `AGENTS.md`
- Create: `docs/physics.md`
- Create: `docs/development.md`
- Create: `benchmarks/README.md`

- [ ] **Step 1: Create `README.md`**

Write `README.md` exactly as follows:

```markdown
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

This calls:

```bash
mpirun -np 1 build/bosonDQMC.out
```

from inside `runs/examples/triangle_3x2`.

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
```

- [ ] **Step 2: Create `AGENTS.md`**

Write `AGENTS.md` exactly as follows:

```markdown
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
```

- [ ] **Step 3: Create `docs/physics.md`**

Write `docs/physics.md` exactly as follows:

```markdown
# Physics And Code Map

## Model

The code simulates the two-flavor Bose-Hubbard model on a triangular lattice:

```text
H = t sum_<ij> (b_i^dagger b_j + c_i^dagger c_j + h.c.)
  + U1 sum_i (n_{b,i} + n_{c,i})^2
  + U2 sum_i (n_{b,i} - n_{c,i})^2
```

The current convention is:

- `t = 1 > 0`
- `U1 < 0`
- `U2 > 0`
- flavors `b` and `c` are also called `up` and `down`

## Hubbard-Stratonovich Fields

The Trotter decomposition uses two continuous auxiliary fields:

- `phi_1` couples to `n_b + n_c` with real coefficient `sqrt(-2 * Dtau * U1)`.
- `phi_2` couples to `n_b - n_c` with imaginary coefficient `i * sqrt(2 * Dtau * U2)`.

After decoupling, the two flavor Hamiltonians are complex conjugates. The code samples the `b` flavor explicitly. The `c` flavor Green matrix is reconstructed with complex conjugation.

## Green Function Convention

For a fixed auxiliary-field configuration, the equal-time Green matrix is:

```text
G_ij(tau) = < b_i b_j^dagger >
```

In code, `Prop%Gr` stores the `b` flavor matrix. `ObserEqual_mod` constructs:

```text
Grup  = Prop%Gr
Grupc = transpose(Grup) - ZKRON
Grdo  = dconjg(Prop%Gr)
Grdoc = dconjg(transpose(Grdo)) - ZKRON
```

## Physical Symbols To Code Variables

| Physical symbol | Code variable or location |
| --- | --- |
| `Lx`, `Ly` | `Nlx`, `Nly` in `src/calc_basic.f90`, read from `paramC_sets.txt` |
| `Lq = Lx * Ly` | `Lq` in `CalcBasic` |
| `beta` | `Beta` in `CalcBasic` |
| `Delta tau` | `Dtau = Beta / Ltrot` in `Params_set` |
| `Ltrot` | `Ltrot` in `CalcBasic`, read from `paramC_sets.txt` |
| `t` | `RT` in `CalcBasic`, currently set to `1.d0` in `Params_set` |
| `U1`, `U2` | `RU1`, `RU2` in `CalcBasic`, read from `paramC_sets.txt` |
| `mu` | `mu` in `CalcBasic`, read from `paramC_sets.txt` |
| auxiliary field flavor index | `ns = 1` for `U1`, `ns = 2` for `U2` |
| auxiliary fields | `Conf%phi_list(ns, ii, nt)` in `src/fields.f90` |
| triangular lattice nearest-neighbor bonds | `Latt%L_bonds(ii, nb)` in `src/lattice.f90` |
| space-time bonds | `Latt%LT_bonds(iit, nb)` in `src/lattice.f90` |
| `b` flavor Green matrix | `Prop%Gr` |
| `c` flavor Green matrix | `dconjg(Prop%Gr)` |
| local update shift | `shiftLoc` in `CalcBasic`, read from `paramC_sets.txt` |
| warm-up shift | `shiftWarm(1:2)` in `CalcBasic`, read from `paramC_sets.txt` |

## Observable Map

| DQMC output file | Fortran field | Meaning |
| --- | --- | --- |
| `density_up` | `Obs%density_up` | per-site density for the `b` flavor |
| `density_do` | `Obs%density_do` | per-site density for the `c` flavor |
| `num_up` | `Obs%num_up` | total `b`-flavor particle number |
| `num_do` | `Obs%num_do` | total `c`-flavor particle number |
| `kinetic` | `Obs%kinetic` | kinetic observable accumulated with division by `Lq` in `Obs_equal_calc` |
| `doubleOcc` | `Obs%doubleOcc` | per-site cross-flavor density product |
| `squareOcc` | `Obs%squareOcc` | per-site same-flavor square contribution |
| `numsquare_up` | `Obs%numsquare_up` | total `b`-flavor number-square estimator |
| `numsquare_do` | `Obs%numsquare_do` | total `c`-flavor number-square estimator |
| `den_upup_sub11` | `Obs%den_corr_up` after Fourier transform | `b-b` density correlation for the single orbital case |
| `den_dodo_sub11` | `Obs%den_corr_do` after Fourier transform | `c-c` density correlation for the single orbital case |
| `den_updo` | `Obs%den_corr_updo` after Fourier transform | cross-flavor density correlation |

## Benchmark Normalization

The ED reference script reports:

- total particle number `NE = NE_b + NE_c`
- total kinetic expectation value

The DQMC code writes:

- `num_up` for the `b` flavor
- `num_do` for the `c` flavor
- `kinetic` divided by `Lq`

Therefore the benchmark comparison should use:

```text
total_NE_DQMC = last(num_up) + last(num_do)
total_kinetic_DQMC = last(kinetic) * Lq
```

Changing this conversion is a physics-level change and must be reflected in `benchmarks/references/*.json` and `benchmarks/README.md`.
```

- [ ] **Step 4: Create `docs/development.md`**

Write `docs/development.md` exactly as follows:

```markdown
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
```

- [ ] **Step 5: Create `benchmarks/README.md`**

Write `benchmarks/README.md` exactly as follows:

```markdown
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
```

- [ ] **Step 6: Verify documentation links and forbidden marker strings**

Run:

```bash
rg -n "T[B]D|T[O]DO|待[定]" README.md AGENTS.md docs benchmarks/README.md
```

Expected: command exits with status 1 because no matches are found.

- [ ] **Step 7: Commit documentation**

Run:

```bash
git add README.md AGENTS.md docs/physics.md docs/development.md benchmarks/README.md
git commit -m "补充物理和开发文档"
```

Expected: git creates a commit containing documentation only.

---

### Task 6: Full Verification And Cleanup

**Files:**
- Modify only if verification finds a path, syntax, or generated-file ignore problem.

- [ ] **Step 1: Check repository status before verification**

Run:

```bash
git status --short
```

Expected: no output after the previous task commits.

- [ ] **Step 2: Verify generated Makefile configuration**

Run:

```bash
make print-config
```

Expected: output includes non-empty `FC=`, `LIB_ROOT=`, `TARGET=build/bosonDQMC.out`, and `RUN_DIR=runs/examples/triangle_3x2`.

- [ ] **Step 3: Verify Python benchmark tooling**

Run:

```bash
python3 -m py_compile benchmarks/compare.py benchmarks/ed/EDtriangle_symm_NEblock.py
make benchmark
```

Expected: Python compilation exits with status 0. `make benchmark` exits with status 0 and prints both benchmark observables as `PASS`.

- [ ] **Step 4: Attempt the Fortran build**

Run:

```bash
make build
```

Expected if compiler and libraries are available: command exits with status 0 and creates `build/bosonDQMC.out`.

Expected if local compiler or external libraries are unavailable: command exits nonzero with a clear compiler or library error. Record the exact error in the final implementation summary; do not claim build success.

- [ ] **Step 5: Run the local example if the build succeeded**

Run:

```bash
make run-example
```

Expected if build and MPI are available: command runs from `runs/examples/triangle_3x2` and writes `info.txt`, scalar observables, and `confout.txt` in that run directory.

Expected if MPI runtime is unavailable: command exits nonzero with a clear `mpirun` or runtime error. Record the exact error in the final implementation summary.

- [ ] **Step 6: Check ignored generated files**

Run:

```bash
git status --short
```

Expected: generated `build/` and generated run outputs do not appear. Any modified tracked file should be intentional.

- [ ] **Step 7: Run whitespace checks**

Run:

```bash
git diff --check
```

Expected: no output.

- [ ] **Step 8: Final commit for verification fixes if needed**

If Steps 1-7 required fixes, run:

```bash
git add Makefile .gitignore scripts benchmarks README.md AGENTS.md docs
git commit -m "修复重构后的验证问题"
```

Expected: git creates a commit only if there were fixes after prior commits. If there were no fixes, skip this step.

---

## Spec Coverage Self-Review

- Move Fortran source files under `src/`: covered by Task 1.
- Replace root build flow with Makefile target set: covered by Task 2.
- Move runnable examples under `runs/examples/`: covered by Task 1.
- Add local run and server submission scripts: covered by Task 3.
- Add default fast benchmark comparison: covered by Task 4.
- Add optional ED recomputation support: covered by Task 4.
- Add physics and contributor documentation: covered by Task 5.
- Preserve fixed filename I/O inside run directories: covered by Tasks 2, 3, and 5.
- Avoid algorithm and observable normalization changes: stated in Tasks 1, 4, and 5; no Fortran behavior edits are planned except file relocation.

## Execution Notes

Implement tasks in order. Do not combine Task 1 with Task 2 because the file move should remain reviewable separately from Makefile logic. Do not run `make benchmark-ed` during routine verification unless the user explicitly asks for ED recomputation.
