"""Exact point-scattering identities and invisible-junction compatibility."""

import unittest

import numpy as np
import sympy
import pyseb


class TestPointSubunits(unittest.TestCase):
    def test_point_self_scattering_amplitude_and_size(self):
        world = pyseb.World()
        world.Add("Point", "site")
        for beta in (-2.0, 1.0, 3.0):
            parameters = {"beta_site": beta}
            for q in (0.0, 0.1, 100.0):
                self.assertEqual(world.EvaluateFormFactor("site", parameters, q), 1.0)
                self.assertEqual(world.EvaluateFormFactorUnnormalized("site", parameters, q), beta**2)
                self.assertEqual(world.EvaluateFormFactorAmplitude("site.point", parameters, q), 1.0)
                self.assertEqual(world.EvaluateFormFactorAmplitudeUnnormalized("site.point", parameters, q), beta)
                self.assertEqual(world.EvaluatePhaseFactor("site.point", "site.point", parameters, q), 1.0)
            self.assertEqual(world.EvaluateRadiusOfGyration2("site", parameters), 0.0)
            self.assertEqual(world.EvaluateSMSDRef2Scat("site.point", parameters), 0.0)
        self.assertEqual(world.EvaluateFormFactorUnnormalized("site", {"beta_site": 0.0}, 1.0), 0.0)

    def test_two_points_match_exact_interference_and_cloud(self):
        world = pyseb.World()
        graph = world.Add("Point", "left")
        world.Link("ThinRod", "bond.end1", "left.point")
        world.Link("Point", "right.point", "bond.end2")
        world.Add(graph, "pair")
        q = np.array([0.0, 1e-4, 0.1, 1.0, 5.0, 100.0])
        for left, right in ((1.0, 1.0), (-2.0, 3.0), (-2.0, 2.0)):
            parameters = {"beta_left": left, "beta_right": right, "beta_bond": 0.0, "L_bond": 1.5}
            expected = left**2 + right**2 + 2 * left * right * np.sinc(q * 1.5 / np.pi)
            actual = world.EvaluateFormFactorUnnormalized("pair", parameters, q)
            np.testing.assert_allclose(actual, expected, atol=1e-12)
            cloud = pyseb.DebyeSphereCloud([
                pyseb.SphereScatterer(0, 0, 0, 0, left),
                pyseb.SphereScatterer(1.5, 0, 0, 0, right),
            ])
            reference = pyseb.World()
            reference.Add(cloud, "cloud")
            np.testing.assert_allclose(
                actual, reference.EvaluateFormFactorUnnormalized("cloud", {}, q), atol=1e-12
            )
            expression = pyseb.to_sympy(world.FormFactor_Unnormalized("pair")).subs(parameters)
            function = sympy.lambdify(sympy.Symbol("q"), expression, modules="numpy")
            np.testing.assert_allclose(function(q[1:]), expected[1:], atol=1e-12)

    def test_invisible_junction_adds_no_scattering_or_parameters(self):
        world = pyseb.World()
        graph = world.Add("InvisiblePoint", "core")
        world.Link("GaussianPolymer", "arm.end1", "core.point")
        world.Add(graph, "star")
        parameters = {"beta_arm": 2.0, "Rg_arm": 3.0}
        q = [0.0, 0.1, 1.0]
        np.testing.assert_allclose(
            world.EvaluateFormFactorUnnormalized("star", parameters, q),
            world.EvaluateFormFactorUnnormalized("arm", parameters, q),
        )
        self.assertEqual(world.EvaluateFormFactorUnnormalized("core", {}, 1.0), 0.0)
        self.assertEqual(world.EvaluatePhaseFactor("core.point", "core.point", {}, 1.0), 1.0)
        reloaded = pyseb.world_from_yaml(world.to_yaml())
        self.assertIn("pyseb/InvisiblePoint", reloaded.to_yaml())
        np.testing.assert_allclose(
            world.EvaluateFormFactorUnnormalized("star", parameters, q),
            reloaded.EvaluateFormFactorUnnormalized("star", parameters, q),
        )


if __name__ == "__main__":
    unittest.main()
