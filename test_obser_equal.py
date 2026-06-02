#!/usr/bin/env python3
"""Tests for equal-time observable Green-function conventions."""

from __future__ import annotations

import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent
OBSER_EQUAL = REPO_ROOT / "obser_equal.f90"


def _compact(text: str) -> str:
    return re.sub(r"\s+", "", text.lower())


class ObserEqualGreenConventionTests(unittest.TestCase):
    def test_down_density_matrix_uses_down_green_without_extra_conjugation(self) -> None:
        source = _compact(OBSER_EQUAL.read_text(encoding="utf-8"))

        self.assertIn("grdo=dconjg(prop%gr)", source)
        self.assertIn("grdoc=transpose(grdo)-zkron", source)
        self.assertNotIn("grdoc=dconjg(transpose(grdo))-zkron", source)


if __name__ == "__main__":
    unittest.main()
