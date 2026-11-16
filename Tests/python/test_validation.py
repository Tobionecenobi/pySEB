import math
import unittest

import pyseb


class TestValidation(unittest.TestCase):
    def test_evaluators_reject_invalid_inputs_and_parameters(self):
        world = pyseb.World()
        world.Add("SolidCylinder", "cylinder")
        parameters = {
            "beta_cylinder": 2.0,
            "R_cylinder": 1.0,
            "L_cylinder": 1.5,
        }

        for invalid_q in (math.nan, math.inf, -math.inf):
            with self.subTest(q=invalid_q):
                with self.assertRaises(RuntimeError):
                    world.EvaluateFormFactor(
                        "cylinder", parameters, invalid_q
                    )

        with self.assertRaises(RuntimeError):
            world.EvaluateFormFactor(
                "cylinder", parameters, [0.1, math.nan]
            )

        missing_radius = dict(parameters)
        del missing_radius["R_cylinder"]
        with self.assertRaises(RuntimeError):
            world.EvaluateFormFactor(
                "cylinder", missing_radius, 0.1
            )

        for name, value in (
            ("R_cylinder", math.nan),
            ("beta_cylinder", math.inf),
        ):
            invalid = dict(parameters)
            invalid[name] = value
            with self.subTest(parameter=name, value=value):
                with self.assertRaises(RuntimeError):
                    world.EvaluateFormFactor(
                        "cylinder", invalid, 0.1
                    )

        with self.assertRaises(RuntimeError):
            world.EvaluateFormFactor("unknown", parameters, 0.1)
        with self.assertRaises(RuntimeError):
            world.EvaluateFormFactorAmplitude(
                "cylinder.unknown", parameters, 0.1
            )

    def test_numerical_callbacks_fail_fast(self):
        missing = pyseb.NumericalSubunit(
            pyseb.NormalizationMode.Normalized
        )
        missing.setTotalBeta(1.0)
        missing_world = pyseb.World()
        missing_world.Add(missing, "missing")
        with self.assertRaises(RuntimeError):
            missing_world.EvaluateFormFactor("missing", {}, 0.1)

        non_finite = pyseb.NumericalSubunit(
            pyseb.NormalizationMode.Normalized
        )
        non_finite.setTotalBeta(1.0)
        non_finite.setFormFactorFunction(
            lambda q, parameters: math.nan
        )
        non_finite_world = pyseb.World()
        non_finite_world.Add(non_finite, "nonfinite")
        with self.assertRaises(RuntimeError):
            non_finite_world.EvaluateFormFactor("nonfinite", {}, 0.1)

        callback_error = pyseb.NumericalSubunit(
            pyseb.NormalizationMode.Normalized
        )
        callback_error.setTotalBeta(1.0)

        def raise_from_callback(q, parameters):
            raise ValueError("callback failure")

        callback_error.setFormFactorFunction(raise_from_callback)
        callback_error_world = pyseb.World()
        callback_error_world.Add(callback_error, "callbackerror")
        with self.assertRaises(ValueError):
            callback_error_world.EvaluateFormFactor(
                "callbackerror", {}, 0.1
            )

        incorrectly_normalized = pyseb.NumericalSubunit(
            pyseb.NormalizationMode.Normalized
        )
        incorrectly_normalized.setTotalBeta(2.0)
        incorrectly_normalized.setFormFactorFunction(
            lambda q, parameters: 0.5
        )
        incorrectly_normalized.setRadiusOfGyration2(1.0)
        with self.assertRaises(RuntimeError):
            incorrectly_normalized.ValidateNumerically()
        with self.assertRaises(RuntimeError):
            incorrectly_normalized.ValidateNumerically(
                parameters={}, tolerance=0.0
            )

    def test_mixed_near_zero_contrast_only_allows_raw_scattering(self):
        def build_world(residual_beta):
            cloud = pyseb.DebyeSphereCloud(
                [
                    pyseb.SphereScatterer(
                        0.0, 0.0, 0.0, 0.0, 1.0
                    ),
                    pyseb.SphereScatterer(
                        1.0, 0.0, 0.0, 0.0, 1.0
                    ),
                ]
            )
            cloud.addReferencePoint("left", 0.0, 0.0, 0.0)
            cloud.addReferencePoint("join", 1.0, 0.0, 0.0)

            world = pyseb.World()
            graph = world.Add(cloud, "cloud")
            world.Link(
                "SolidSphere",
                "sphere.center",
                "cloud.join",
            )
            world.Add(graph, "mixed")
            parameters = {
                "beta_sphere": -2.0 + residual_beta,
                "R_sphere": 1.0,
            }
            return world, parameters

        for residual_beta in (0.0, 5e-15, -5e-15):
            world, parameters = build_world(residual_beta)
            with self.subTest(residual_beta=residual_beta):
                self.assertTrue(
                    math.isfinite(
                        world.EvaluateFormFactorUnnormalized(
                            "mixed", parameters, 0.2
                        )
                    )
                )
                self.assertTrue(
                    math.isfinite(
                        world.EvaluateFormFactorAmplitudeUnnormalized(
                            "mixed:cloud.left", parameters, 0.2
                        )
                    )
                )
                with self.assertRaises(RuntimeError):
                    world.EvaluateFormFactor(
                        "mixed", parameters, 0.2
                    )
                with self.assertRaises(RuntimeError):
                    world.EvaluateFormFactorAmplitude(
                        "mixed:cloud.left", parameters, 0.2
                    )
                with self.assertRaises(RuntimeError):
                    world.EvaluateRadiusOfGyration2(
                        "mixed", parameters
                    )


if __name__ == "__main__":
    unittest.main()
