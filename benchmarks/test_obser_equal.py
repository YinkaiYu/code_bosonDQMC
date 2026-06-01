#!/usr/bin/env python3
"""Tests for equal-time observable Green-function conventions."""

from __future__ import annotations

import json
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
OBSER_EQUAL = REPO_ROOT / "src" / "obser_equal.f90"
PHYSICS_DOC = REPO_ROOT / "docs" / "physics.md"


def _compact(text: str) -> str:
    return re.sub(r"\s+", "", text.lower())


class ObserEqualGreenConventionTests(unittest.TestCase):
    def test_down_density_matrix_uses_down_green_without_extra_conjugation(self) -> None:
        source = _compact(OBSER_EQUAL.read_text(encoding="utf-8"))

        self.assertIn("grdo=dconjg(prop%gr)", source)
        self.assertIn("grdoc=transpose(grdo)-zkron", source)
        self.assertNotIn("grdoc=dconjg(transpose(grdo))-zkron", source)

    def test_physics_doc_matches_down_green_convention(self) -> None:
        doc = _compact(PHYSICS_DOC.read_text(encoding="utf-8"))

        self.assertIn("grdo=dconjg(prop%gr)", doc)
        self.assertIn("grdoc=transpose(grdo)-zkron", doc)
        self.assertNotIn("grdoc=dconjg(transpose(grdo))-zkron", doc)

    def test_ed_script_emits_equal_time_scalar_observables(self) -> None:
        ed_script = REPO_ROOT / "benchmarks" / "ed" / "EDtriangle_symm_NEblock.py"
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            (run_dir / "params.txt").write_text(
                "Lx = 3\n"
                "Ly = 2\n"
                "t = 1.0\n"
                "U1 = 0.0\n"
                "U2 = 1.0\n"
                "beta = 6.0\n"
                "mu = -2.5\n",
                encoding="utf-8",
            )
            result = subprocess.run(
                [sys.executable, str(ed_script)],
                cwd=run_dir,
                text=True,
                capture_output=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            data = json.loads((run_dir / "results.json").read_text(encoding="utf-8"))

        observables = data["observables"]
        self.assertAlmostEqual(observables["doubleOcc"], 5.0104309962475375e-06)
        self.assertAlmostEqual(observables["squareOcc"], 2.3674394261440157e-09)
        self.assertAlmostEqual(observables["numsquare_up"], 0.00041666206338322737)
        self.assertAlmostEqual(observables["numsquare_do"], 0.00041666206338322737)


if __name__ == "__main__":
    unittest.main()
