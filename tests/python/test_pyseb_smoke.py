import math
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
        graph_id = world.Add("GaussianPolymer", "poly1", "")
        world.Link("GaussianPolymer", "poly2.end1", "poly1.end2")
        world.Add(graph_id, "diblockcopolymer")

        form_factor = world.FormFactor("diblockcopolymer")
        self.assertIsNotNone(form_factor)

    def test_form_factor_string_is_generated(self):
        world = pyseb.World()
        graph_id = world.Add("GaussianPolymer", "poly1", "")
        world.Link("GaussianPolymer", "poly2.end1", "poly1.end2")
        world.Add(graph_id, "diblockcopolymer")

        form_factor = world.FormFactor("diblockcopolymer")
        expr_string = str(form_factor)
        self.assertTrue(len(expr_string) > 0)
        self.assertIn("beta_", expr_string)

    def test_print_form_factor_as_sympy_equation(self):
        world = pyseb.World()
        graph_id = world.Add("GaussianPolymer", "poly1", "")
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
        graph_id = world.Add("GaussianPolymer", "poly1", "")
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

    def test_python_symbolic_construction_converts_to_sympy(self):
        pyseb.set_backend("portable")
        x = pyseb.symbol("x")
        expr = pyseb.sin(x * pyseb.constant(2.0) + 3.0)

        self.assertIsInstance(expr, pyseb.Expression)
        self.assertIn("sin", expr.to_python())

        sympy_expr = pyseb.to_sympy(expr)
        self.assertIsInstance(sympy_expr, sympy.Expr)
        self.assertEqual(float(sympy_expr.subs("x", 4.0).evalf()), expr.subs("x", 4.0).eval())

    def test_symbolic_backend_numeric_parity(self):
        results = {}
        for backend in ("portable", "ginac"):
            if backend not in pyseb.available_backends():
                continue

            pyseb.set_backend(backend)
            x = pyseb.symbol("x")
            expr = pyseb.sin(x * 2.0 + 3.0)
            results[backend] = expr.subs("x", 4.0).eval()

        self.assertIn("portable", results)

        pyseb.set_backend("portable")
        x = pyseb.symbol("x")
        expr = pyseb.sin(x * 2.0 + 3.0)
        sympy_value = float(pyseb.to_sympy(expr).subs("x", 4.0).evalf())
        results["sympy"] = sympy_value

        for value in results.values():
            self.assertAlmostEqual(value, math.sin(11.0), places=12)

        if "ginac" in results:
            self.assertAlmostEqual(results["portable"], results["ginac"], places=12)

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
        for backend in ("portable", "ginac"):
            if backend not in pyseb.available_backends():
                continue

            pyseb.set_backend(backend)
            expr = build_expression()
            results[backend] = expr.subs("x", 1.25).subs("y", 2.5).eval()

        self.assertIn("portable", results)

        pyseb.set_backend("portable")
        expr = build_expression()
        sympy_expr = pyseb.to_sympy(expr)
        sympy_value = float(sympy_expr.subs({"x": 1.25, "y": 2.5}).evalf())
        results["sympy"] = sympy_value

        for value in results.values():
            self.assertTrue(math.isfinite(value))
            self.assertAlmostEqual(value, results["portable"], places=12)

    def test_portable_special_functions_convert_to_sympy(self):
        pyseb.set_backend("portable")
        x = pyseb.symbol("x")
        expr = (x / 2.0).erf() + (x / 5.0).erfc() + x.bessel_j0() + x.bessel_j1()

        sympy_expr = pyseb.to_sympy(expr)
        value = float(sympy_expr.subs("x", 1.25).evalf())
        expected = (
            math.erf(1.25 / 2.0)
            + math.erfc(1.25 / 5.0)
            + float(sympy.besselj(0, 1.25).evalf())
            + float(sympy.besselj(1, 1.25).evalf())
        )

        self.assertAlmostEqual(value, expected, places=12)


if __name__ == "__main__":
    unittest.main()
