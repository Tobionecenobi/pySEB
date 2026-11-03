"""Optional scientific checks, kept outside the straight-through example."""

from copy import deepcopy
import importlib.util
from pathlib import Path
import runpy
import unittest

import numpy as np


@unittest.skipUnless(importlib.util.find_spec("pyPRISM"), "optional pyPRISM is not installed")
class TestPyprismPrototype(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).resolve().parents[3] / "Examples" / "pyprism_seb.py"
        cls.data = runpy.run_path(str(path))

    def test_chain_normalization_and_independent_bond_sum(self):
        d = self.data
        q = np.array([0.0, 1e-4, 0.1, 0.5, 1.0, 4.0, 100.0, 1e6 / d["bond_length"]])
        n = d["n_sites"]
        e = np.sinc(q * d["bond_length"] / np.pi)
        # N-s site pairs are separated by s independently oriented bonds.
        expected = 1 + sum(2 * (n - s) / n * e**s for s in range(1, n))
        actual = n * np.asarray(d["world"].EvaluateFormFactor("chain", d["parameters"], q))
        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
        self.assertAlmostEqual(actual[0], n)
        self.assertAlmostEqual(actual[-1], 1, delta=3e-6)

    def test_matches_builtin_input_and_solution(self):
        d = self.data
        system = deepcopy(d["system"])
        reference = d["pyPRISM"].omega.FreelyJointedChain(
            length=d["n_sites"], l=d["bond_length"]
        )
        np.testing.assert_allclose(d["omega"], reference.calculate(d["q"]), rtol=1e-9)
        system.omega["site", "site"] = reference
        prism = system.createPRISM()
        result = prism.solve(options={"disp": False, "fatol": 1e-9, "maxiter": 200})
        self.assertTrue(result.success, result.message)
        self.assertLess(np.max(np.abs(prism.cost(result.x))), 1e-7)
        structure = d["pyPRISM"].calculate.structure_factor(prism)["site", "site"]
        np.testing.assert_allclose(d["structure"], structure, rtol=2e-6, atol=1e-8)
        self.assertTrue(np.all(np.isfinite(d["structure"]) & (d["structure"] > 0)))
        self.assertTrue(np.all(np.isfinite(d["pair"])) & np.all(d["pair"] >= -1e-5))
        np.testing.assert_allclose(d["pair"][system.domain.r < 1], 0, atol=1e-6)
        self.assertGreater(np.max(np.abs(d["structure"] - d["omega"])), 0.1)

    def test_dilute_solution_approaches_isolated_chain(self):
        d = self.data
        system = deepcopy(d["system"])
        system.density["site"] = 1e-8
        prism = system.createPRISM()
        result = prism.solve(options={"disp": False, "fatol": 1e-9, "maxiter": 200})
        self.assertTrue(result.success, result.message)
        self.assertLess(np.max(np.abs(prism.cost(result.x))), 1e-7)
        structure = d["pyPRISM"].calculate.structure_factor(prism)["site", "site"]
        np.testing.assert_allclose(structure, d["omega"], rtol=1e-5, atol=1e-7)


if __name__ == "__main__":
    unittest.main()
