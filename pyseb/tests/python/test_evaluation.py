import math
import unittest

import numpy as np
import pyseb
import sympy


class TestEvaluation(unittest.TestCase):
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

    def test_numerical_subunit_callbacks_and_symbolic_placeholder(self):
        unit = pyseb.NumericalSubunit(pyseb.NormalizationMode.Normalized)
        unit.addReferencePoint("center")
        unit.setTotalBeta(2.0)
        unit.setFormFactorFunction(lambda q, params: math.exp(-q * q))
        unit.setFormFactorAmplitudeFunction(
            "center", lambda q, params: math.exp(-0.5 * q * q)
        )
        unit.setRadiusOfGyration2(3.0)
        unit.setSigmaMSDRef2Scat("center", 3.0)
        unit.ValidateNumerically()

        world = pyseb.World()
        world.Add(unit, "numeric")

        self.assertAlmostEqual(
            world.EvaluateFormFactor("numeric", {}, 0.2),
            math.exp(-0.04),
        )
        self.assertAlmostEqual(
            world.EvaluateFormFactorAmplitude("numeric.center", {}, 0.2),
            math.exp(-0.02),
        )
        self.assertIn("F_numeric", str(world.FormFactor("numeric")))

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

    def test_numpy_helpers_evaluate_integrated_subunits(self):
        world = pyseb.World()
        world.Add("SolidCylinder", "cylinder")
        params = {
            "beta_cylinder": 2.0,
            "R_cylinder": 1.0,
            "L_cylinder": 1.5,
        }
        q_values = np.array([[0.0, 0.1], [0.3, 0.8]])

        form_factor = pyseb.evaluate_form_factor(
            world, "cylinder", params, q_values
        )
        amplitude = pyseb.evaluate_form_factor_amplitude(
            world, "cylinder.center", params, q_values
        )
        phase = pyseb.evaluate_phase_factor(
            world,
            "cylinder.hull",
            "cylinder.ends",
            params,
            q_values,
        )

        self.assertIsInstance(form_factor, np.ndarray)
        self.assertEqual(form_factor.shape, q_values.shape)
        self.assertEqual(amplitude.shape, q_values.shape)
        self.assertEqual(phase.shape, q_values.shape)

        expected = np.array(
            world.EvaluateFormFactor(
                "cylinder", params, q_values.reshape(-1).tolist()
            )
        ).reshape(q_values.shape)
        np.testing.assert_allclose(form_factor, expected, rtol=1e-12, atol=1e-12)
        self.assertEqual(
            pyseb.evaluate_form_factor(world, "cylinder", params, 0.0),
            1.0,
        )

        unnormalized = pyseb.evaluate_form_factor(
            world,
            "cylinder",
            params,
            q_values,
            normalized=False,
        )
        np.testing.assert_allclose(
            unnormalized,
            4.0 * form_factor,
            rtol=1e-12,
            atol=1e-12,
        )

    def test_numpy_helpers_preserve_empty_shape(self):
        world = pyseb.World()
        world.Add("ThinDisk", "disk")
        params = {"beta_disk": 1.0, "R_disk": 1.0}
        q_values = np.empty((2, 0, 3))

        result = pyseb.evaluate_form_factor(
            world, "disk", params, q_values
        )
        self.assertIsInstance(result, np.ndarray)
        self.assertEqual(result.shape, q_values.shape)

    def test_numpy_helpers_cover_scalar_empty_negative_and_nonfinite_q(self):
        world = pyseb.World()
        world.Add("SolidCylinder", "cylinder")
        parameters = {
            "beta_cylinder": 2.0,
            "R_cylinder": 1.0,
            "L_cylinder": 1.5,
        }

        evaluators = (
            lambda q: pyseb.evaluate_form_factor(
                world, "cylinder", parameters, q
            ),
            lambda q: pyseb.evaluate_form_factor_amplitude(
                world, "cylinder.center", parameters, q
            ),
            lambda q: pyseb.evaluate_phase_factor(
                world,
                "cylinder.hull",
                "cylinder.ends",
                parameters,
                q,
            ),
        )

        empty_q = np.empty((2, 0, 3))
        for evaluate in evaluators:
            with self.subTest(case="empty", evaluator=evaluate):
                result = evaluate(empty_q)
                self.assertIsInstance(result, np.ndarray)
                self.assertEqual(result.shape, empty_q.shape)

            with self.subTest(case="scalar", evaluator=evaluate):
                result = evaluate(np.asarray(0.2))
                self.assertIsInstance(result, float)
                self.assertTrue(math.isfinite(result))

            with self.subTest(case="negative", evaluator=evaluate):
                symmetric_q = np.array([-0.4, 0.4])
                original_q = symmetric_q.copy()
                values = evaluate(symmetric_q)
                np.testing.assert_array_equal(symmetric_q, original_q)
                np.testing.assert_allclose(
                    values[0], values[1], rtol=1e-12, atol=1e-12
                )

            with self.subTest(case="nonfinite", evaluator=evaluate):
                with self.assertRaises(RuntimeError):
                    evaluate(np.array([0.1, math.inf]))


if __name__ == "__main__":
    unittest.main()
