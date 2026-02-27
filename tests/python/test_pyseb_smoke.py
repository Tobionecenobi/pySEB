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


if __name__ == "__main__":
    unittest.main()
