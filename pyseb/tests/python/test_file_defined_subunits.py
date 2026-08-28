import math
from pathlib import Path
import tempfile
import unittest

import pyseb

from pyseb.cli import subunit_main


ROOT = Path(__file__).resolve().parents[3]
MODELS = ROOT / "models"


class TestFileDefinedSubunits(unittest.TestCase):
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
        self.assertEqual(len(definition.integrals), 14)
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
