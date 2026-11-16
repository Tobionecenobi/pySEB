import math
import unittest

import numpy as np
import pyseb
import sympy


class TestSymbolic(unittest.TestCase):
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

    def test_add_type_and_name_overload_matches_examples(self):
        world = pyseb.World()
        graph_id = world.Add("GaussianPolymer", "poly1")

        world.Link("GaussianPolymer", "poly2.end1", "poly1.end2")
        world.Add(graph_id, "diblockcopolymer")

        form_factor = world.FormFactor("diblockcopolymer")
        self.assertIn("beta_poly1", str(form_factor))
        self.assertIn("beta_poly2", str(form_factor))

    def test_python_symbolic_construction_converts_to_sympy(self):
        pyseb.set_backend("portable")
        x = pyseb.symbol("x")
        expr = pyseb.sin(x * pyseb.constant(2.0) + 3.0)

        self.assertIsInstance(expr, pyseb.Expression)
        self.assertIn("sin", expr.to_python())

        sympy_expr = pyseb.to_sympy(expr)
        self.assertIsInstance(sympy_expr, sympy.Expr)
        self.assertEqual(float(sympy_expr.subs("x", 4.0).evalf()), math.sin(11.0))
        self.assertAlmostEqual(expr.subs("x", 4.0).eval(), math.sin(11.0))

    def test_portable_symbolic_output_is_stable(self):
        pyseb.set_backend("portable")
        x = pyseb.symbol("portable_x")
        y = pyseb.symbol("portable_y")

        self.assertEqual(str(x + y), "(portable_x + portable_y)")
        self.assertEqual(str(-x), "(-portable_x)")
        self.assertEqual(str(x**pyseb.constant(2.0)), "(portable_x ** 2)")
        self.assertEqual(x.abs().to_python(), "Abs(portable_x)")
        self.assertEqual(x.bessel_j0().to_python(), "besselj(0, portable_x)")
        self.assertEqual(x.bessel_j1().to_python(), "besselj(1, portable_x)")
        self.assertEqual(x.dawson().to_python(), "dawson(portable_x)")

        combined = x.sin() + y.cos()
        expected = "(sin(portable_x) + cos(portable_y))"
        self.assertEqual(str(combined), expected)
        self.assertEqual(combined.to_latex(), expected)
        self.assertEqual(combined.to_cform(), expected)

    def test_portable_symbolic_output_preserves_parentheses(self):
        pyseb.set_backend("portable")
        x = pyseb.symbol("precedence_x")
        y = pyseb.symbol("precedence_y")
        z = pyseb.symbol("precedence_z")

        self.assertEqual(
            ((x + y) * z).to_python(),
            "((precedence_x + precedence_y) * precedence_z)",
        )
        self.assertEqual(
            (x + y * z).to_python(),
            "(precedence_x + (precedence_y * precedence_z))",
        )
        self.assertEqual(
            (x - (y - z)).to_python(),
            "(precedence_x - (precedence_y - precedence_z))",
        )
        self.assertEqual(
            pyseb.pow(x + y, pyseb.constant(2.0)).to_python(),
            "((precedence_x + precedence_y) ** 2)",
        )

    def test_portable_symbolic_output_round_trips_through_sympy(self):
        pyseb.set_backend("portable")
        x = pyseb.symbol("roundtrip_x")
        y = pyseb.symbol("roundtrip_y")
        expr = (x.sin() + y.sqrt()) / (x + 2.0)

        parsed = sympy.sympify(expr.to_python())
        expected = (
            sympy.sin(sympy.Symbol("roundtrip_x"))
            + sympy.sqrt(sympy.Symbol("roundtrip_y"))
        ) / (sympy.Symbol("roundtrip_x") + 2)

        self.assertEqual(sympy.simplify(parsed - expected), 0)

    def test_portable_special_function_output_round_trips_through_sympy(self):
        pyseb.set_backend("portable")
        x = pyseb.symbol("special_x")
        expr = (
            x.abs()
            + x.bessel_j0()
            + x.bessel_j1()
            + x.dawson()
            + x.erf()
            + x.erfc()
        )

        parsed = pyseb.SymPyExpression(expr.to_python()).expr
        self.assertIsInstance(parsed, sympy.Expr)
        self.assertIn("special_x", str(parsed))

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
