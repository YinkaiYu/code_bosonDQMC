#!/usr/bin/env bash
#SBATCH -J yyk_triangle
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p node6348

set -euo pipefail

run_dir="${1:-${RUN_DIR:-runs/examples/triangle_3x2}}"
np="${SLURM_NTASKS:-1}"

script_repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ -n "${REPO_ROOT:-}" ]]; then
  repo_root="${REPO_ROOT}"
elif [[ -f "${script_repo}/Makefile" && -d "${script_repo}/src" ]]; then
  repo_root="${script_repo}"
else
  repo_root="${SLURM_SUBMIT_DIR:-$(pwd)}"
fi
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

if ! command -v mpirun >/dev/null 2>&1; then
  echo "mpirun is not available on PATH" >&2
  exit 2
fi

make -C "${repo_root}" build

cd "${run_dir}"
exec mpirun -np "${np}" "${exe}"
