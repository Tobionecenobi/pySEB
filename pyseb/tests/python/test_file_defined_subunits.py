import math
from pathlib import Path
import tempfile
import unittest

import pyseb

from pyseb.cli import subunit_main


ROOT = Path(__file__).resolve().parents[3]
MODELS = ROOT / "models"


class TestFileDefinedSubunits(unittest.TestCase):
    def test_spheroid_surface_moments_and_center_phase(self):
        world = pyseb.World("spheroid-area")
        world.Add("Spheroid", "shape")
        # Independent 60-digit surface-area integrals, including the sphere limit.
        for a, c, surface_moment in (
            (2.0, 3.0, 5.507228480390443),
            (3.0, 1.0, 5.5240378791693595),
            (2.0, 2.0, 4.0),
        ):
            with self.subTest(a=a, c=c):
                parameters = {"a_shape": a, "c_shape": c, "beta_shape": 1.0}
                volume_moment = (2 * a*a + c*c) / 5
                for first, second, expected in (
                    ("center", "surface", surface_moment),
                    ("pole", "surface", c*c + surface_moment),
                    ("surface", "surface", 2 * surface_moment),
                ):
                    self.assertAlmostEqual(
                        world.EvaluateSMSDRef2Ref(
                            "shape." + first, "shape." + second, parameters
                        ), expected, places=9,
                    )
                self.assertAlmostEqual(
                    world.EvaluateSMSDRef2Scat("shape.surface", parameters),
                    volume_moment + surface_moment, places=9,
                )
                q = 0.04 / max(a, c)
                for method, refs, moment in (
                    (world.EvaluateFormFactorAmplitude, ("shape.surface",),
                     volume_moment + surface_moment),
                    (world.EvaluatePhaseFactor, ("shape.center", "shape.surface"),
                     surface_moment),
                ):
                    slopes = [6 * (1 - method(*refs, parameters, step)) / step**2
                              for step in (q, q / 2)]
                    self.assertAlmostEqual(
                        (4 * slopes[1] - slopes[0]) / 3, moment, delta=2e-5,
                    )
        parameters = {"a_shape": 2.0, "c_shape": 2.0, "beta_shape": 1.0}
        self.assertAlmostEqual(
            world.EvaluatePhaseFactor("shape.center", "shape.surface", parameters, 0.5),
            math.sin(1.0), places=10,
        )

    def test_surface_linked_spheroid_pair_size_matches_guinier(self):
        world = pyseb.World("spheroid-pair")
        graph = world.Add("Spheroid", "first", "shape")
        world.Link("Spheroid", "second.surface#join", "first.surface#join", "shape")
        world.Add(graph, "pair")
        parameters = {"a_shape": 2.0, "c_shape": 3.0, "beta_shape": 1.0}
        # Equal-volume bodies: Rg_pair^2 = Rg_body^2 + <center separation^2>/4.
        # Independent surface draws give <center separation^2> = 2*Ms.
        expected = 17/5 + 5.507228480390443/2
        reported = world.EvaluateRadiusOfGyration2("pair", parameters)
        self.assertAlmostEqual(reported, expected, places=9)
        slopes = [3 * (1 - world.EvaluateFormFactor("pair", parameters, q)) / q**2
                  for q in (0.01, 0.005)]
        self.assertAlmostEqual((4 * slopes[1] - slopes[0]) / 3, reported, delta=2e-5)

    def test_size_integrals_require_q_independence(self):
        source = (MODELS / "Point.pyseb.yaml").read_text()
        source = source.split("\nvalidation:\n", 1)[0]
        source = source.replace(
            'radius_of_gyration_squared: "0"',
            'radius_of_gyration_squared: "moment"',
        )
        source += ('\ndefinitions: {moment: "geometry / 2"}\n'
                   'integrals:\n'
                   '  geometry: {variable: t, lower: "0", upper: "2", integrand: "t^2"}\n')
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "integral-size.pyseb.yaml"
            path.write_text(source)
            world = pyseb.World("integral-size")
            world.add_subunit_file(str(path), "unit")
            self.assertAlmostEqual(
                world.EvaluateRadiusOfGyration2("unit", {"beta_unit": 1.0}), 4/3
            )
            for integrand in ("q * t^2", "qAlias * t^2"):
                invalid = source.replace('integrand: "t^2"', f'integrand: "{integrand}"')
                invalid += '\nvariables: {qAlias: "q"}\n'
                path.write_text(invalid)
                with self.assertRaisesRegex(RuntimeError, "size expressions cannot depend on q"):
                    pyseb.load_subunit_definition(str(path))

    def test_thin_circle_radius_of_gyration_matches_guinier_expansion(self):
        world = pyseb.World("ring-guinier")
        world.Add("ThinCircle", "ring")
        for radius in (0.5, 1.0, 3.0):
            with self.subTest(radius=radius):
                parameters = {"R_ring": radius, "beta_ring": 1.0}
                reported = world.EvaluateRadiusOfGyration2("ring", parameters)
                self.assertAlmostEqual(reported, radius**2, places=12)
                q = 0.02 / radius
                estimates = [
                    3 * (1 - world.EvaluateFormFactor("ring", parameters, step)) / step**2
                    for step in (q, q / 2)
                ]
                extrapolated = (4 * estimates[1] - estimates[0]) / 3
                self.assertAlmostEqual(reported, extrapolated, delta=1e-7 * radius**2)

    def test_sphere_kernel_remains_stable_through_variable_aliases(self):
        source = (MODELS / "SolidSphere.pyseb.yaml").read_text()
        source = source.replace("pyseb/SolidSphere", "test/AliasedSphere").replace(
            "api_name: SolidSphere", "api_name: AliasedSphere"
        ).replace(
            '  x: "q * R"', '  x: "q * R"\n  shape: "sphere_amplitude(x)"'
        ).replace('  A: "sphere_amplitude(x)"', '  A: "shape"')
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sphere.pyseb.yaml"
            path.write_text(source)
            world = pyseb.World("aliased")
            world.register_subunit_file(str(path))
            world.Add("AliasedSphere", "sphere")
            parameters = {"R_sphere": 2.0, "beta_sphere": -4.0}
            self.assertAlmostEqual(
                world.EvaluateFormFactorUnnormalized("sphere", parameters, 1e-10), 16.0
            )

    def test_bundled_files_load_and_validate(self):
        paths = sorted(MODELS.glob("*.pyseb.yaml"))
        self.assertTrue(paths)
        world = pyseb.World("catalogue")
        catalogue = {
            model.id: model
            for model in world.list_subunit_models()
            if model.source.startswith("<bundled:")
        }
        self.assertEqual(len(catalogue), len(paths))
        for index, path in enumerate(paths):
            definition = pyseb.load_subunit_definition(str(path))
            report = pyseb.validate_subunit_file(str(path))
            self.assertTrue(definition.api_name[0].isalpha())
            self.assertTrue(definition.api_name.isalnum())
            self.assertTrue(report.ok, [failure.message for failure in report.failures])
            self.assertIn(definition.id, catalogue)
            self.assertEqual(catalogue[definition.id].api_name, definition.api_name)
            world.Add(definition.id, f"canonical{index}")
            world.Add(definition.api_name, f"alias{index}")

    def test_direct_registered_and_legacy_models_are_equivalent(self):
        path = str(MODELS / "GaussianPolymer.pyseb.yaml")
        worlds = [pyseb.World("legacy"), pyseb.World("canonical"), pyseb.World("direct")]
        worlds[0].Add("GaussianPolymer", "poly")
        worlds[1].Add("pyseb/GaussianPolymer", "poly")
        worlds[2].add_subunit_file(path, "poly")
        parameters = {"beta_poly": 2.0, "Rg_poly": 2.0}

        for q in (0.0, 1e-6, 0.1, 0.5, 1.0):
            values = [world.EvaluateFormFactor("poly", parameters, q) for world in worlds]
            self.assertTrue(all(math.isfinite(value) for value in values))
            self.assertAlmostEqual(values[0], values[1], places=12)
            self.assertAlmostEqual(values[0], values[2], places=12)

    def test_directory_registration_is_atomic(self):
        source = (MODELS / "GaussianPolymer.pyseb.yaml").read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            (directory / "a.pyseb.yaml").write_text(
                source.replace("pyseb/GaussianPolymer", "lab/First").replace(
                    "api_name: GaussianPolymer", "api_name: First"
                ), encoding="utf-8"
            )
            (directory / "b.pyseb.yaml").write_text(
                source.replace("pyseb/GaussianPolymer", "lab/Second").replace(
                    "api_name: GaussianPolymer", "api_name: Second"
                ).replace(
                    "form_factor: F", "form_factor: system(1)"
                ),
                encoding="utf-8",
            )
            world = pyseb.World("atomic")
            with self.assertRaises(RuntimeError):
                world.register_subunit_directory(str(directory))
            with self.assertRaises(RuntimeError):
                world.Add("lab/First", "first")

    def test_directory_registration_order_is_deterministic(self):
        source = (MODELS / "Point.pyseb.yaml").read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            (directory / "z.pyseb.yaml").write_text(
                source.replace("pyseb/Point", "lab/Last").replace(
                    "api_name: Point", "api_name: Last"
                ), encoding="utf-8"
            )
            (directory / "a.pyseb.yaml").write_text(
                source.replace("pyseb/Point", "lab/First").replace(
                    "api_name: Point", "api_name: First"
                ), encoding="utf-8"
            )
            registered = pyseb.World("ordered").register_subunit_directory(str(directory))
            self.assertEqual([model.id for model in registered], ["lab/First", "lab/Last"])

    def test_registration_and_collision_behavior(self):
        source = (MODELS / "GaussianPolymer.pyseb.yaml").read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "custom.pyseb.yaml"
            path.write_text(
                source.replace("pyseb/GaussianPolymer", "lab/Polymer").replace(
                    "api_name: GaussianPolymer", "api_name: LabPolymer"
                ),
                encoding="utf-8",
            )
            world = pyseb.World("registry")
            info = world.register_subunit_file(str(path))
            self.assertEqual(info.id, "lab/Polymer")
            self.assertEqual(info.api_name, "LabPolymer")
            world.Add("lab/Polymer", "poly")
            world.Add("LabPolymer", "polyalias")
            with self.assertRaises(RuntimeError):
                world.register_subunit_file(str(path))

    def test_registration_rejects_api_name_collisions_atomically(self):
        source = (MODELS / "Point.pyseb.yaml").read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            (directory / "a.pyseb.yaml").write_text(
                source.replace("pyseb/Point", "lab/First").replace(
                    "api_name: Point", "api_name: SharedName"
                ),
                encoding="utf-8",
            )
            (directory / "b.pyseb.yaml").write_text(
                source.replace("pyseb/Point", "lab/Second").replace(
                    "api_name: Point", "api_name: SharedName"
                ),
                encoding="utf-8",
            )
            world = pyseb.World("collision")
            with self.assertRaises(RuntimeError):
                world.register_subunit_directory(str(directory))
            with self.assertRaises(RuntimeError):
                world.Add("lab/First", "first")

            bundled_collision = directory / "bundled.pyseb.yaml"
            bundled_collision.write_text(
                source.replace("pyseb/Point", "lab/BundledCollision"),
                encoding="utf-8",
            )
            with self.assertRaises(RuntimeError):
                world.register_subunit_file(str(bundled_collision))

            builtin_collision = directory / "builtin.pyseb.yaml"
            builtin_collision.write_text(
                source.replace("pyseb/Point", "lab/BuiltinCollision").replace(
                    "api_name: Point", "api_name: ThinRod"
                ),
                encoding="utf-8",
            )
            with self.assertRaises(RuntimeError):
                world.register_subunit_file(str(builtin_collision))

    def test_removed_model_classes_use_registry_construction(self):
        import pyseb.subunits as subunits

        self.assertFalse(hasattr(subunits, "Point"))
        self.assertFalse(hasattr(subunits, "GaussianPolymer"))
        self.assertFalse(hasattr(subunits, "ThinDisk"))
        self.assertFalse(hasattr(subunits, "SolidCylinder"))

    def test_integral_metadata_is_read_only_and_discoverable(self):
        definition = pyseb.load_subunit_definition(
            str(MODELS / "SolidCylinder.pyseb.yaml")
        )
        self.assertEqual(
            definition.integration.method,
            pyseb.IntegrationMethod.CQUAD,
        )
        self.assertEqual(definition.integration.workspace_size, 1000)
        self.assertEqual(definition.integration.qag_rule, 61)
        self.assertEqual(len(definition.integrals), 18)
        self.assertEqual(
            set(definition.specific_references),
            {"center", "top", "bottom"},
        )
        variants = definition.visualization.references["ends"].variants
        self.assertEqual(set(variants), {"top", "bottom"})
        self.assertEqual(variants["top"].geometry, "ends_top")
        self.assertEqual(variants["top"].sampling, "surface_area")
        self.assertEqual(variants["bottom"].geometry, "ends_bottom")
        form_factor = definition.integrals["formFactor"]
        self.assertEqual(form_factor.variable, "theta")
        self.assertEqual(form_factor.lower, "0")
        self.assertEqual(form_factor.upper, "pi / 2")
        with self.assertRaises(AttributeError):
            form_factor.variable = "other"

        spheroid = pyseb.load_subunit_definition(str(MODELS / "Spheroid.pyseb.yaml"))
        self.assertEqual(
            set(spheroid.specific_references),
            {"center", "pole", "north", "south"},
        )
        surface_amplitude = spheroid.integrals["surfaceAmplitude"]
        self.assertEqual(
            [(dimension.variable, dimension.lower, dimension.upper)
             for dimension in surface_amplitude.dimensions],
            [("orientation", "0", "pi / 2"), ("theta", "0", "pi")],
        )
        self.assertEqual(
            len(spheroid.integrals["surfaceSurfacePhase"].dimensions), 3
        )

    def test_solid_cylinder_has_fixed_cap_center_references(self):
        world = pyseb.World("cylinder-cap-centers")
        world.add_subunit_file(
            str(MODELS / "SolidCylinder.pyseb.yaml"), "cylinder"
        )
        parameters = {
            "beta_cylinder": 1.0,
            "R_cylinder": 1.0,
            "L_cylinder": 1.5,
        }
        q = 0.2

        top = world.EvaluateFormFactorAmplitude(
            "cylinder.top", parameters, q
        )
        bottom = world.EvaluateFormFactorAmplitude(
            "cylinder.bottom", parameters, q
        )
        self.assertAlmostEqual(top, bottom, places=12)
        self.assertAlmostEqual(
            world.EvaluatePhaseFactor(
                "cylinder.top", "cylinder.bottom", parameters, q
            ),
            math.sin(q * 1.5) / (q * 1.5),
            places=12,
        )
        self.assertAlmostEqual(
            world.EvaluatePhaseFactor(
                "cylinder.center", "cylinder.top", parameters, q
            ),
            math.sin(q * 1.5 / 2.0) / (q * 1.5 / 2.0),
            places=12,
        )
        self.assertAlmostEqual(
            world.EvaluateSMSDRef2Scat("cylinder.top", parameters),
            1.25,
            places=12,
        )
        self.assertAlmostEqual(
            world.EvaluateSMSDRef2Ref(
                "cylinder.top", "cylinder.bottom", parameters
            ),
            2.25,
            places=12,
        )

    def test_schema_and_expression_errors_include_context(self):
        source = (MODELS / "Point.pyseb.yaml").read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "broken.pyseb.yaml"
            path.write_text(source + "unknown_field: true\n", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, r"unknown_field"):
                pyseb.load_subunit_definition(str(path))

        cylinder = (MODELS / "SolidCylinder.pyseb.yaml").read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "broken-variant.pyseb.yaml"
            path.write_text(
                cylinder.replace("geometry: ends_top", "geometry: missing_patch"),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(RuntimeError, r"missing_patch"):
                pyseb.load_subunit_definition(str(path))

    def test_schema_rejects_reference_errors_and_multiple_documents(self):
        point = (MODELS / "Point.pyseb.yaml").read_text(encoding="utf-8")
        polymer = (MODELS / "GaussianPolymer.pyseb.yaml").read_text(encoding="utf-8")
        invalid_sources = (
            point.replace("specific: [point]", "specific: [point, point]"),
            point.replace("specific: [point]", "specific: [bad_name]"),
            polymer.replace(
                '    - references: [contour, contour]\n      expression: F\n', "", 1
            ),
            point + "---\nformat: pyseb-subunit\n",
        )
        with tempfile.TemporaryDirectory() as directory:
            for index, source in enumerate(invalid_sources):
                path = Path(directory) / f"invalid-{index}.pyseb.yaml"
                path.write_text(source, encoding="utf-8")
                with self.subTest(index=index), self.assertRaises(RuntimeError):
                    pyseb.load_subunit_definition(str(path))

    def test_validation_reports_failure_and_valid_without_cases(self):
        source = (MODELS / "Point.pyseb.yaml").read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            failing = directory / "failing.pyseb.yaml"
            failing.write_text(
                source.replace("expected: 1.0", "expected: 2.0", 1), encoding="utf-8"
            )
            report = pyseb.validate_subunit_file(str(failing))
            self.assertFalse(report.ok)
            self.assertEqual(subunit_main(["validate", str(failing)]), 1)

            no_cases = directory / "no-cases.pyseb.yaml"
            no_cases.write_text(source.split("\nvalidation:\n", 1)[0] + "\n", encoding="utf-8")
            report = pyseb.validate_subunit_file(str(no_cases))
            self.assertTrue(report.ok)
            self.assertEqual(report.case_count, 0)
            self.assertTrue(report.warnings)
            self.assertEqual(subunit_main(["validate", str(no_cases)]), 0)

    def test_cli_success_and_failure_codes(self):
        self.assertEqual(
            subunit_main(["validate", str(MODELS / "GaussianPolymer.pyseb.yaml")]), 0
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.pyseb.yaml"
            path.write_text("format: wrong\n", encoding="utf-8")
            self.assertEqual(subunit_main(["validate", str(path)]), 2)


if __name__ == "__main__":
    unittest.main()
