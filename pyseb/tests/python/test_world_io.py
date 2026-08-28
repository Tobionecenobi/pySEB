import tempfile
import unittest
from pathlib import Path

import pyseb


WORLD_YAML = """\
format: pyseb-world
schema_version: 1
world:
  id: diblock-pair
subunits:
  - name: A
    model: pyseb/GaussianPolymer
    tag: A
  - name: B
    model: pyseb/GaussianPolymer
    tag: B
structures:
  - name: diblock1
    graph: diblock-graph
  - name: diblock2
    graph: diblock-graph
  - name: pair
    graph: pair-graph
graphs:
  - id: diblock-graph
    members: [A, B]
    links:
      - [A.end2, B.end1]
  - id: pair-graph
    members: [diblock1, diblock2]
    links:
      - [diblock1:A.end1, diblock2:A.end1]
  - id: root
    members: [pair]
    links: []
roots: [root]
"""


def nested_world_yaml(depth: int) -> str:
    lines = [
        "format: pyseb-world",
        "schema_version: 1",
        "world: {id: deep}",
        "subunits: [{name: A, model: pyseb/Point, tag: A}]",
        "structures:",
    ]
    lines.extend(
        f"  - {{name: S{index}, graph: g{index + 1}}}"
        for index in range(depth)
    )
    lines.append("graphs:")
    lines.extend(
        f"  - {{id: g{index}, members: [S{index}], links: []}}"
        for index in range(depth)
    )
    lines.extend(
        (
            f"  - {{id: g{depth}, members: [A], links: []}}",
            "roots: [g0]",
        )
    )
    return "\n".join(lines) + "\n"


class TestWorldIO(unittest.TestCase):
    def test_nested_world_import_and_round_trip(self):
        world = pyseb.world_from_yaml(WORLD_YAML)
        exported = world.to_yaml()
        reloaded = pyseb.world_from_yaml(exported)

        self.assertEqual(exported, reloaded.to_yaml())
        self.assertEqual(
            str(world.FormFactor("pair")),
            str(reloaded.FormFactor("pair")),
        )
        parameters = {
            "beta_A": 1.0,
            "beta_B": 2.0,
            "Rg_A": 3.0,
            "Rg_B": 4.0,
        }
        for q in (0.0, 0.1, 0.4):
            self.assertAlmostEqual(
                world.EvaluateFormFactor("pair", parameters, q),
                reloaded.EvaluateFormFactor("pair", parameters, q),
            )

    def test_save_and_load_world_file(self):
        world = pyseb.world_from_yaml(WORLD_YAML)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "pair.pyseb-world.yaml"
            pyseb.save_world(world, str(path))
            loaded = pyseb.load_world(str(path))
            self.assertEqual(world.to_yaml(), loaded.to_yaml())

    def test_schema_validation_rejects_unknown_fields_and_disconnected_graphs(self):
        with self.assertRaises(RuntimeError):
            pyseb.world_from_yaml(WORLD_YAML.replace("world:\n", "extra: true\nworld:\n"))

        disconnected = WORLD_YAML.replace(
            "links:\n      - [A.end2, B.end1]",
            "links: []",
        )
        with self.assertRaises(RuntimeError):
            pyseb.world_from_yaml(disconnected)

    def test_schema_uses_world_name_rules(self):
        underscored = WORLD_YAML.replace("name: A\n", "name: A_one\n").replace(
            "members: [A, B]", "members: [A_one, B]",
        ).replace("A.end2", "A_one.end2")
        with self.assertRaises(RuntimeError):
            pyseb.world_from_yaml(underscored)

    def test_deep_structure_nesting_does_not_use_the_call_stack(self):
        pyseb.world_from_yaml(nested_world_yaml(15_000))

    def test_control_characters_round_trip_in_yaml_strings(self):
        control_character = WORLD_YAML.replace("id: diblock-pair", "id: \"diblock-\\x01pair\"")
        world = pyseb.world_from_yaml(control_character)
        exported = world.to_yaml()
        self.assertIn('id: "diblock-\\u0001pair"', exported)
        self.assertEqual(exported, pyseb.world_from_yaml(exported).to_yaml())

    def test_import_rejects_models_without_world_schema_representation(self):
        symbolic = WORLD_YAML.replace("pyseb/GaussianPolymer", "SymbolicSubunit")
        with self.assertRaises(RuntimeError):
            pyseb.world_from_yaml(symbolic)

    def test_export_rejects_callback_subunits(self):
        numerical = pyseb.NumericalSubunit(pyseb.NormalizationMode.Unnormalized)
        numerical.addReferencePoint("center")
        numerical.setTotalBeta(1.0)
        numerical.setFormFactorFunction(lambda q, parameters: 1.0)
        world = pyseb.World()
        world.Add(numerical, "cloud")
        with self.assertRaises(RuntimeError):
            pyseb.world_to_yaml(world)

    def test_failed_save_does_not_truncate_existing_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "existing.pyseb-world.yaml"
            path.write_text("preserve this file")
            with self.assertRaises(RuntimeError):
                pyseb.World().save(str(path))
            self.assertEqual(path.read_text(), "preserve this file")


if __name__ == "__main__":
    unittest.main()
