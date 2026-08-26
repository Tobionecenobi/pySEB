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
        for name, model_id in (
            ("Point.pyseb.yaml", "pyseb/Point"),
            ("GaussianPolymer.pyseb.yaml", "pyseb/GaussianPolymer"),
        ):
            path = MODELS / name
            definition = pyseb.load_subunit_definition(str(path))
            report = pyseb.validate_subunit_file(str(path))
            self.assertEqual(definition.id, model_id)
            self.assertTrue(report.ok, [failure.message for failure in report.failures])
            self.assertGreater(report.case_count, 0)

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
                source.replace("pyseb/GaussianPolymer", "lab/First"), encoding="utf-8"
            )
            (directory / "b.pyseb.yaml").write_text(
                source.replace("pyseb/GaussianPolymer", "lab/Second").replace(
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
                source.replace("pyseb/Point", "lab/Last"), encoding="utf-8"
            )
            (directory / "a.pyseb.yaml").write_text(
                source.replace("pyseb/Point", "lab/First"), encoding="utf-8"
            )
            registered = pyseb.World("ordered").register_subunit_directory(str(directory))
            self.assertEqual([model.id for model in registered], ["lab/First", "lab/Last"])

    def test_registration_and_collision_behavior(self):
        source = (MODELS / "GaussianPolymer.pyseb.yaml").read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "custom.pyseb.yaml"
            path.write_text(source.replace("pyseb/GaussianPolymer", "lab/Polymer"), encoding="utf-8")
            world = pyseb.World("registry")
            info = world.register_subunit_file(str(path))
            self.assertEqual(info.id, "lab/Polymer")
            world.Add("lab/Polymer", "poly")
            with self.assertRaises(RuntimeError):
                world.register_subunit_file(str(path))

    def test_schema_and_expression_errors_include_context(self):
        source = (MODELS / "Point.pyseb.yaml").read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "broken.pyseb.yaml"
            path.write_text(source + "unknown_field: true\n", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, r"unknown_field"):
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
