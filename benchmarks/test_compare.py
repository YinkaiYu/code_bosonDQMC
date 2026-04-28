#!/usr/bin/env python3
"""Tests for the benchmark comparison suite."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
COMPARE = REPO_ROOT / "benchmarks" / "compare.py"
POLE_CHECK = REPO_ROOT / "benchmarks" / "check_pole_diagnostics.py"
REFERENCES = REPO_ROOT / "benchmarks" / "references"
DQMC_REFERENCES = REPO_ROOT / "benchmarks" / "dqmc_references"
DQMC_SUITE = REPO_ROOT / "benchmarks" / "dqmc_suite.json"

EXPECTED_CASES = {
    "triangle_3x2_beta3_mu-2.5_u1_0_u2_1": {
        "beta": 3.0,
        "mu": -2.5,
        "U1": 0.0,
        "U2": 1.0,
        "total_NE": 0.08999902580923487,
        "total_kinetic": -0.17786109516236884,
    },
    "triangle_3x2_beta6_mu-2.5_u1_0_u2_1": {
        "beta": 6.0,
        "mu": -2.5,
        "U1": 0.0,
        "U2": 1.0,
        "total_NE": 0.0008330248173669302,
        "total_kinetic": -0.001659465125418983,
    },
    "triangle_3x2_beta1_mu-5_u1_-0.1_u2_1": {
        "beta": 1.0,
        "mu": -5.0,
        "U1": -0.1,
        "U2": 1.0,
        "total_NE": 0.13997061039816858,
        "total_kinetic": -0.2558942387788746,
    },
    "triangle_3x2_beta1.4_mu-5_u1_-0.1_u2_1": {
        "beta": 1.4,
        "mu": -5.0,
        "U1": -0.1,
        "U2": 1.0,
        "total_NE": 0.027122698028982088,
        "total_kinetic": -0.052062742470074974,
    },
}


class BenchmarkSuiteTests(unittest.TestCase):
    def test_reference_directory_runs_all_committed_cases(self) -> None:
        result = subprocess.run(
            [
                sys.executable,
                str(COMPARE),
                "--reference-dir",
                str(REFERENCES),
            ],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        for case_name in EXPECTED_CASES:
            self.assertIn(f"Case: {case_name}", result.stdout)
        self.assertIn("All benchmark cases passed.", result.stdout)

    def test_committed_references_cover_benchmark_note_cases(self) -> None:
        references = {
            path.stem: json.loads(path.read_text(encoding="utf-8"))
            for path in REFERENCES.glob("*.json")
        }

        self.assertLessEqual(set(EXPECTED_CASES), set(references))
        for case_name, expected in EXPECTED_CASES.items():
            reference = references[case_name]
            params = reference["parameters"]
            observables = reference["observables"]

            self.assertEqual(reference["case"], case_name)
            self.assertEqual(params["Lx"], 3)
            self.assertEqual(params["Ly"], 2)
            self.assertEqual(params["t"], 1.0)
            self.assertEqual(params["beta"], expected["beta"])
            self.assertEqual(params["mu"], expected["mu"])
            self.assertEqual(params["U1"], expected["U1"])
            self.assertEqual(params["U2"], expected["U2"])
            self.assertEqual(observables["total_NE"]["value"], expected["total_NE"])
            self.assertEqual(
                observables["total_kinetic"]["value"], expected["total_kinetic"]
            )
            self.assertEqual(
                observables["total_NE"]["dqmc"]["operation"], "sum_last"
            )
            self.assertEqual(
                observables["total_NE"]["dqmc"]["files"], ["num_up", "num_do"]
            )
            self.assertEqual(
                observables["total_kinetic"]["dqmc"]["operation"], "last_times_lq"
            )
            self.assertIn("dqmc_fixture", reference)

    def test_makefile_exposes_real_dqmc_benchmark_target(self) -> None:
        result = subprocess.run(
            ["make", "-n", "benchmark-dqmc"],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("benchmarks/run_dqmc_suite.py", result.stdout)
        self.assertIn("benchmarks/dqmc_suite.json", result.stdout)

    def test_makefile_benchmark_is_strict_live_suite(self) -> None:
        result = subprocess.run(
            ["make", "-n", "benchmark"],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("benchmarks/run_dqmc_suite.py", result.stdout)
        self.assertIn("benchmarks/dqmc_suite.json", result.stdout)

    def test_makefile_exposes_fast_noninteracting_benchmark(self) -> None:
        result = subprocess.run(
            ["make", "-n", "benchmark-fast"],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("scripts/run_dqmc_benchmark.sh", result.stdout)
        self.assertIn("triangle_3x2_free_beta3_mu-2.5", result.stdout)
        self.assertIn("benchmarks/dqmc_references/triangle_3x2_free_beta3_mu-2.5.json", result.stdout)

    def test_makefile_keeps_fixture_check_outside_benchmark_targets(self) -> None:
        result = subprocess.run(
            ["make", "-n", "check-fixtures"],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("benchmarks/compare.py", result.stdout)
        self.assertIn("benchmarks/references", result.stdout)
        self.assertNotIn("benchmarks/run_dqmc_suite.py", result.stdout)

        removed_alias = subprocess.run(
            ["make", "-n", "benchmark-fixture"],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertNotEqual(removed_alias.returncode, 0)

    def test_dqmc_suite_covers_all_ed_benchmark_cases(self) -> None:
        suite = json.loads(DQMC_SUITE.read_text(encoding="utf-8"))
        cases = {case["name"]: case for case in suite["cases"]}

        self.assertIn("triangle_3x2_free_beta3_mu-2.5", cases)
        for case_name in EXPECTED_CASES:
            self.assertIn(case_name, cases)

        for case in cases.values():
            reference = REPO_ROOT / case["reference"]
            input_dir = REPO_ROOT / case["input_dir"]
            self.assertTrue(reference.is_file(), reference)
            self.assertTrue(input_dir.is_dir(), input_dir)
            live_reference = json.loads(reference.read_text(encoding="utf-8"))
            params = live_reference["parameters"]
            self.assertEqual(params["Ltrot"], int(round(params["beta"] * 100)))

        for case_name in EXPECTED_CASES:
            reference_path = REPO_ROOT / cases[case_name]["reference"]
            reference = json.loads(reference_path.read_text(encoding="utf-8"))
            params = reference["parameters"]
            observables = reference["observables"]
            self.assertGreaterEqual(params["Nbin"], 100000)
            self.assertEqual(observables["total_NE"]["dqmc"]["operation"], "sum_mean")
            self.assertEqual(
                observables["total_kinetic"]["dqmc"]["operation"], "mean_times_lq"
            )
            for observable in observables.values():
                self.assertIn("statistics", observable["dqmc"])
                self.assertIn("stderr_tolerance", observable["dqmc"]["statistics"])
                self.assertEqual(
                    observable["dqmc"]["statistics"]["stderr_tolerance"], 3.0
                )
                self.assertGreaterEqual(
                    observable["dqmc"]["statistics"]["block_size"], 10000
                )

    def test_compare_supports_mean_observables_for_real_dqmc_runs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            run_dir = tmp_path / "run"
            run_dir.mkdir()
            (run_dir / "num_up").write_text("0.08\n0.09\n0.10\n", encoding="utf-8")
            (run_dir / "kinetic").write_text("-0.05\n-0.06\n-0.07\n", encoding="utf-8")

            reference = tmp_path / "reference.json"
            reference.write_text(
                json.dumps(
                    {
                        "case": "mean_dqmc_case",
                        "parameters": {"Lx": 3, "Ly": 2},
                        "observables": {
                            "total_NE": {
                                "value": 0.09,
                                "atol": 1e-12,
                                "dqmc": {
                                    "operation": "mean",
                                    "files": ["num_up"],
                                    "min_samples": 3,
                                },
                            },
                            "total_kinetic": {
                                "value": -0.36,
                                "atol": 1e-12,
                                "dqmc": {
                                    "operation": "mean_times_lq",
                                    "files": ["kinetic"],
                                    "min_samples": 3,
                                },
                            },
                        },
                    }
                ),
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(COMPARE),
                    "--reference",
                    str(reference),
                    "--run-dir",
                    str(run_dir),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=False,
            )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("PASS total_NE", result.stdout)
        self.assertIn("PASS total_kinetic", result.stdout)

    def test_compare_supports_block_statistical_error_for_real_dqmc_runs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            run_dir = tmp_path / "run"
            run_dir.mkdir()
            (run_dir / "num_up").write_text("1\n3\n5\n7\n", encoding="utf-8")
            (run_dir / "num_do").write_text("2\n4\n6\n8\n", encoding="utf-8")

            reference = tmp_path / "reference.json"
            reference.write_text(
                json.dumps(
                    {
                        "case": "block_stats_case",
                        "parameters": {"Lx": 1, "Ly": 1},
                        "observables": {
                            "total_NE": {
                                "value": 9.1,
                                "atol": 0.0,
                                "dqmc": {
                                    "operation": "sum_mean",
                                    "files": ["num_up", "num_do"],
                                    "min_samples": 4,
                                    "statistics": {
                                        "block_size": 2,
                                        "stderr_tolerance": 1.0,
                                    },
                                },
                            }
                        },
                    }
                ),
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(COMPARE),
                    "--reference",
                    str(reference),
                    "--run-dir",
                    str(run_dir),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=False,
            )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("PASS total_NE", result.stdout)
        self.assertIn("stderr=", result.stdout)
        self.assertIn("z=", result.stdout)

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

    def test_pole_diagnostic_checker_rejects_malformed_required_param_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            (run_dir / "paramC_sets.txt").write_text(
                "-1.0 1.0 -2.5\n"
                "BROKEN LATTICE ROW\n"
                "3 2 60\n"
                "10 2 1 0.3\n"
                ".false. 0\n"
                ".false. 500 1.0 1.0\n"
                "2 2 0.0 0.0\n",
                encoding="utf-8",
            )
            z_row = " ".join(["1.0 0.0"] * 6)
            (run_dir / "pole_z").write_text(f"{z_row}\n{z_row}\n", encoding="utf-8")
            (run_dir / "pole_distance").write_text("1.0\n0.5\n", encoding="utf-8")
            (run_dir / "pole_x").write_text("0.0\n0.3010299956639812\n", encoding="utf-8")
            (run_dir / "green_spectral_radius").write_text("1.0\n2.0\n", encoding="utf-8")
            (run_dir / "green_smax").write_text("1.0\n2.1\n", encoding="utf-8")
            (run_dir / "log_weight").write_text("-3.0\n-2.5\n", encoding="utf-8")

            result = subprocess.run(
                [sys.executable, str(POLE_CHECK), str(run_dir)],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=False,
            )

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("ERROR:", result.stdout + result.stderr)
        self.assertIn("paramC_sets.txt", result.stdout + result.stderr)

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


if __name__ == "__main__":
    unittest.main()
