import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pyseb
import numpy as np
import sympy


class TestPySEBSmoke(unittest.TestCase):
    def test_import_and_world_creation(self):
        world = pyseb.World()
        self.assertIsNotNone(world)

    def test_build_diblock_structure_and_form_factor(self):
        world = pyseb.World()
        graph_id = world.Add("GaussianPolymer", "poly1")
        world.Link("GaussianPolymer", "poly2.end1", "poly1.end2")
        world.Add(graph_id, "diblockcopolymer")

        form_factor = world.FormFactor("diblockcopolymer")
        self.assertIsNotNone(form_factor)

    def test_form_factor_string_is_generated(self):
        world = pyseb.World()
        graph_id = world.Add("GaussianPolymer", "poly1")
        world.Link("GaussianPolymer", "poly2.end1", "poly1.end2")
        world.Add(graph_id, "diblockcopolymer")

        form_factor = world.FormFactor("diblockcopolymer")
        expr_string = str(form_factor)
        self.assertTrue(len(expr_string) > 0)
        self.assertIn("beta_", expr_string)

    def test_print_form_factor_as_sympy_equation(self):
        world = pyseb.World()
        graph_id = world.Add("GaussianPolymer", "poly1")
        world.Link("GaussianPolymer", "poly2.end1", "poly1.end2")
        world.Add(graph_id, "diblockcopolymer")

        form_factor = world.FormFactor("diblockcopolymer")
        sympy_expr = pyseb.SymPyExpression(str(form_factor))
        simplified_expr = sympy.simplify(sympy_expr.expr)
        latex_expr = sympy.latex(simplified_expr)
        print(f"SymPy equation: {simplified_expr}")
        print(f"LaTeX equation: {latex_expr}")

        self.assertIsNotNone(sympy_expr.expr)

    def test_sympy_evaluate_expression(self):
        world = pyseb.World()
        graph_id = world.Add("GaussianPolymer", "poly1")
        world.Link("GaussianPolymer", "poly2.end1", "poly1.end2")
        world.Add(graph_id, "diblockcopolymer")

        form_factor = world.FormFactor("diblockcopolymer")
        params = {
            "beta_poly1": 1.0,
            "beta_poly2": 1.0,
            "Rg_poly1": 5.0,
            "Rg_poly2": 7.0,
        }

        value_at_scalar_q = pyseb.evaluate_expression(world, form_factor, params, 0.1)
        values_at_q_list = pyseb.evaluate_expression(world, form_factor, params, [0.05, 0.1, 0.2])

        self.assertTrue(isinstance(value_at_scalar_q, float))
        self.assertTrue(math.isfinite(value_at_scalar_q))
        self.assertEqual(len(values_at_q_list), 3)
        self.assertIsInstance(values_at_q_list, np.ndarray)
        for value in values_at_q_list:
            self.assertTrue(math.isfinite(float(value)))

    def test_add_type_and_name_overload_matches_examples(self):
        world = pyseb.World()
        graph_id = world.Add("GaussianPolymer", "poly1")

        world.Link("GaussianPolymer", "poly2.end1", "poly1.end2")
        world.Add(graph_id, "diblockcopolymer")

        form_factor = world.FormFactor("diblockcopolymer")
        self.assertIn("beta_poly1", str(form_factor))
        self.assertIn("beta_poly2", str(form_factor))

    def test_small_micelle_example_builds_and_evaluates(self):
        world = pyseb.World()
        graph_id = world.Add("SolidSphere", "sphere")

        for index in range(3):
            world.Link(
                "GaussianPolymer",
                f"poly{index}.end1",
                f"sphere.surface#r{index}",
                "poly",
            )

        world.Add(graph_id, "micelle")
        form_factor = world.FormFactor("micelle")

        params = {
            "beta_sphere": 1.0,
            "beta_poly": 0.5,
            "R_sphere": 50.0,
            "Rg_poly": 20.0,
        }
        value = pyseb.evaluate_expression(world, form_factor, params, 0.1)

        self.assertTrue(math.isfinite(value))

    def test_sympy_and_helper_evaluation_match_for_diblock(self):
        world = pyseb.World()
        graph_id = world.Add("GaussianPolymer", "poly1")
        world.Link("GaussianPolymer", "poly2.end1", "poly1.end2")
        world.Add(graph_id, "diblockcopolymer")

        form_factor = world.FormFactor("diblockcopolymer")
        sympy_form_factor = pyseb.to_sympy(form_factor)
        simplified_form_factor = sympy.simplify(sympy_form_factor)
        params = {
            "beta_poly1": 1.0,
            "beta_poly2": 0.8,
            "Rg_poly1": 5.0,
            "Rg_poly2": 8.0,
        }

        self.assertEqual(sympy.simplify(sympy_form_factor - simplified_form_factor), 0)

        q_values = np.array([0.01, 0.05, 0.1, 0.2])
        q = sympy.Symbol("q")
        sympy_fn = sympy.lambdify(q, sympy_form_factor.subs(params), modules=["numpy"])
        sympy_values = np.asarray(sympy_fn(q_values), dtype=float)
        helper_values = pyseb.evaluate_expression(world, form_factor, params, q_values)

        np.testing.assert_allclose(sympy_values, helper_values, rtol=1e-12, atol=1e-12)

    def test_sympy_and_helper_evaluation_match_for_micelle(self):
        world = pyseb.World()
        graph_id = world.Add("SolidSphere", "core")

        for index in range(6):
            world.Link(
                "GaussianPolymer",
                f"chain{index}.end1",
                f"core.surface#anchor{index}",
                "chain",
            )

        world.Add(graph_id, "micelle")
        form_factor = world.FormFactor("micelle")
        sympy_form_factor = pyseb.to_sympy(form_factor)
        simplified_form_factor = sympy.simplify(sympy_form_factor)
        params = {
            "beta_core": 1.0,
            "beta_chain": 0.25,
            "R_core": 30.0,
            "Rg_chain": 12.0,
        }

        self.assertEqual(sympy.simplify(sympy_form_factor - simplified_form_factor), 0)

        q = sympy.Symbol("q")
        numeric_expr = sympy_form_factor.subs(params)

        for q_value in (0.02, 0.05, 0.1, 0.2):
            sympy_value = float(numeric_expr.subs(q, q_value).evalf())
            helper_value = pyseb.evaluate_expression(world, form_factor, params, q_value)
            self.assertAlmostEqual(sympy_value, helper_value, places=12)

    def test_sympy_example_runs_when_launched_by_absolute_path(self):
        repo_root = Path(__file__).resolve().parents[3]
        script = repo_root / "pyseb" / "examples" / "python" / "sympy_example.py"

        with tempfile.TemporaryDirectory() as tmpdir:
            result = subprocess.run(
                [sys.executable, str(script)],
                cwd=tmpdir,
                text=True,
                capture_output=True,
                timeout=30,
            )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Form Factor Expression:", result.stdout)

    def test_python_examples_run_when_launched_by_absolute_path(self):
        repo_root = Path(__file__).resolve().parents[3]
        examples_dir = repo_root / "pyseb" / "examples" / "python"

        expected_output = {
            "backend_simplification_comparison.py": "Backend simplification comparison",
            "diblock_evaluation.py": "Diblock copolymer form factor",
            "figure13_star_chain.py": "Figure 13 star-chain form factor",
            "micelle_evaluation.py": "Micelle form factor",
            "nested_structures.py": "Nested diblock-star chain form factor",
            "symbolic_integral.py": "Solid cylinder symbolic integral",
            "symbolic_backend.py": "Available backends:",
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            for filename, expected in expected_output.items():
                result = subprocess.run(
                    [sys.executable, str(examples_dir / filename)],
                    cwd=tmpdir,
                    text=True,
                    capture_output=True,
                    timeout=30,
                )

                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn(expected, result.stdout)

    def test_python_symbolic_construction_converts_to_sympy(self):
        pyseb.set_backend("portable")
        x = pyseb.symbol("x")
        expr = pyseb.sin(x * pyseb.constant(2.0) + 3.0)

        self.assertIsInstance(expr, pyseb.Expression)
        self.assertIn("sin", expr.to_python())

        sympy_expr = pyseb.to_sympy(expr)
        self.assertIsInstance(sympy_expr, sympy.Expr)
        self.assertEqual(float(sympy_expr.subs("x", 4.0).evalf()), math.sin(11.0))
        with self.assertRaises(RuntimeError):
            expr.subs("x", 4.0)

    def test_symbolic_backend_exports_to_sympy_for_numeric_evaluation(self):
        results = {}

        pyseb.set_backend("portable")
        x = pyseb.symbol("x")
        expr = pyseb.sin(x * 2.0 + 3.0)
        results["sympy_from_portable"] = float(pyseb.to_sympy(expr).subs("x", 4.0).evalf())

        if "ginac" in pyseb.available_backends():
            pyseb.set_backend("ginac")
            x = pyseb.symbol("x")
            expr = pyseb.sin(x * 2.0 + 3.0)
            results["ginac"] = expr.subs("x", 4.0).eval()

        for value in results.values():
            self.assertAlmostEqual(value, math.sin(11.0), places=12)

        if "ginac" in results:
            self.assertAlmostEqual(results["sympy_from_portable"], results["ginac"], places=12)

    def test_symbolic_backend_common_function_parity(self):
        def build_expression():
            x = pyseb.symbol("x")
            y = pyseb.symbol("y")
            return (
                pyseb.sin(x * 2.0 + y)
                + pyseb.exp(x / 3.0)
                + pyseb.log(y + 4.0)
                + pyseb.sqrt(x + pyseb.e())
                + pyseb.sin(pyseb.pi())
            )

        results = {}
        pyseb.set_backend("portable")
        expr = build_expression()
        sympy_expr = pyseb.to_sympy(expr)
        results["sympy_from_portable"] = float(sympy_expr.subs({"x": 1.25, "y": 2.5}).evalf())

        if "ginac" in pyseb.available_backends():
            pyseb.set_backend("ginac")
            expr = build_expression()
            results["ginac"] = expr.subs("x", 1.25).subs("y", 2.5).eval()

        for value in results.values():
            self.assertTrue(math.isfinite(value))
            self.assertAlmostEqual(value, results["sympy_from_portable"], places=12)

    def test_symbolic_backend_special_function_parity(self):
        def build_expression():
            x = pyseb.symbol("x")
            return (
                (x / 2.0).erf()
                + (x / 5.0).erfc()
                + x.bessel_j0()
                + x.bessel_j1()
                + (x / 3.0).dawson()
            )

        pyseb.set_backend("portable")
        expr = build_expression()

        sympy_expr = pyseb.to_sympy(expr)
        value = float(sympy_expr.subs("x", 1.25).evalf())
        expected = (
            math.erf(1.25 / 2.0)
            + math.erfc(1.25 / 5.0)
            + float(sympy.besselj(0, 1.25).evalf())
            + float(sympy.besselj(1, 1.25).evalf())
            + float(
                (
                    sympy.exp(-(1.25 / 3.0) ** 2)
                    * sympy.sqrt(sympy.pi)
                    / 2
                    * sympy.erfi(1.25 / 3.0)
                ).evalf()
            )
        )

        self.assertAlmostEqual(value, expected, places=12)
        results = {"sympy_from_portable": value}

        if "ginac" in pyseb.available_backends():
            pyseb.set_backend("ginac")
            expr = build_expression()
            results["ginac"] = expr.subs("x", 1.25).eval()

        for result in results.values():
            self.assertTrue(math.isfinite(result))
            self.assertAlmostEqual(result, results["sympy_from_portable"], places=12)

    def test_portable_scattering_special_functions_export_to_sympy(self):
        pyseb.set_backend("portable")

        rod_world = pyseb.World()
        rod_world.Add("ThinRod", "rod")
        rod_form_factor = pyseb.to_sympy(rod_world.FormFactor("rod"))

        self.assertIn("Si", str(rod_form_factor))
        self.assertNotIn("sin(q*L_rod)/(q*L_rod)", str(rod_form_factor))

        cylinder_world = pyseb.World()
        cylinder_world.Add("SolidCylinder", "cylinder")
        cylinder_form_factor = pyseb.to_sympy(cylinder_world.FormFactor("cylinder"))

        self.assertIn("Integral", str(cylinder_form_factor))

    def test_sympy_dawson_method_matches_parser_mapping(self):
        x = sympy.Symbol("x")
        from_method = pyseb.SymPyExpression(x).dawson().expr
        from_string = pyseb.SymPyExpression("dawson(x)").expr

        self.assertEqual(sympy.simplify(from_method - from_string), 0)


if __name__ == "__main__":
    unittest.main()
