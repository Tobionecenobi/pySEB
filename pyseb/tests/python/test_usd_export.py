import math
from pathlib import Path
import re
import tempfile
import unittest

import pyseb


def _vector(text):
    return tuple(float(value) for value in text.split(", "))


def _pose(document, name):
    pattern = re.compile(
        rf'    def Xform "{re.escape(name)}" \{{\n'
        rf'        double3 xformOp:translate = \(([^)]+)\)\n'
        rf'        quatd xformOp:orient = \(([^,]+), \(([^)]+)\)\)'
    )
    match = pattern.search(document)
    if not match:
        raise AssertionError(f"missing pose for {name}")
    return _vector(match.group(1)), (float(match.group(2)),) + _vector(match.group(3))


def _reference(document, instance, reference):
    start = document.index(f'    def Xform "{instance}"')
    next_instance = document.find('\n    def Xform "', start + 1)
    block = document[start:next_instance if next_instance >= 0 else len(document)]
    match = re.search(
        rf'def Xform "ref_{re.escape(reference)}".*?'
        rf'double3 xformOp:translate = \(([^)]+)\)',
        block,
        re.S,
    )
    if not match:
        raise AssertionError(f"missing reference {instance}.{reference}")
    return _vector(match.group(1))


def _rotate(quaternion, point):
    w, x, y, z = quaternion
    ux, uy, uz = x, y, z
    px, py, pz = point
    cross1 = (uy*pz-uz*py, uz*px-ux*pz, ux*py-uy*px)
    cross2 = (
        uy*cross1[2]-uz*cross1[1],
        uz*cross1[0]-ux*cross1[2],
        ux*cross1[1]-uy*cross1[0],
    )
    return tuple(point[i] + 2*w*cross1[i] + 2*cross2[i] for i in range(3))


class TestUSDExport(unittest.TestCase):
    def _sphere_rods(self, seed, output):
        world = pyseb.World("usd")
        graph = world.Add("SolidSphere", "aSphere")
        graph = world.Link("ThinRod", "rod1.end1", "aSphere.surface#anchor1")
        graph = world.Link("ThinRod", "rod2.end1", "aSphere.surface#anchor2")
        world.Add(graph, "assembly")
        options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
        options.seed = seed
        options.surface_samples = 24
        options.curve_samples = 8
        world.export_usd(
            "assembly",
            str(output),
            {"R_aSphere": 3.0, "L_rod1": 2.0, "L_rod2": 2.0},
            options,
        )
        return world

    def test_position_only_links_coincide_after_free_rotation(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "assembly.usda"
            self._sphere_rods(42, output)
            document = output.read_text(encoding="utf-8")
            sphere_translation, sphere_rotation = _pose(document, "aSphere")
            self.assertEqual(sphere_rotation, (1.0, 0.0, 0.0, 0.0))

            rotations = []
            for index in (1, 2):
                rod_translation, rod_rotation = _pose(document, f"rod{index}")
                rotations.append(rod_rotation)
                sphere_reference = _reference(
                    document, "aSphere", f"surface_anchor{index}"
                )
                rod_reference = _reference(document, f"rod{index}", "end1")
                sphere_world = tuple(
                    sphere_translation[i] + _rotate(sphere_rotation, sphere_reference)[i]
                    for i in range(3)
                )
                rod_world = tuple(
                    rod_translation[i] + _rotate(rod_rotation, rod_reference)[i]
                    for i in range(3)
                )
                for actual, expected in zip(rod_world, sphere_world):
                    self.assertAlmostEqual(actual, expected, places=12)
                self.assertAlmostEqual(
                    sum(value*value for value in rod_rotation), 1.0, places=12
                )
            self.assertNotEqual(rotations[0], rotations[1])
            self.assertIn(
                'pyseb:orientationSemantics = "representative_free_rotation"',
                document,
            )
            self.assertIn('pyseb:overlapPolicy = "permitted"', document)
            self.assertIn("rel pyseb:link_0", document)

    def test_seed_is_reproducible_and_does_not_change_scattering(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            first, second, different = (
                directory / "first.usda",
                directory / "second.usda",
                directory / "different.usda",
            )
            world = self._sphere_rods(7, first)
            parameters = {
                "R_aSphere": 3.0,
                "L_rod1": 2.0,
                "L_rod2": 2.0,
                "beta_aSphere": 1.0,
                "beta_rod1": 1.0,
                "beta_rod2": 1.0,
            }
            before = world.EvaluateFormFactor("assembly", parameters, 0.15)
            options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
            options.seed = 7
            options.surface_samples = 24
            options.curve_samples = 8
            world.export_usd("assembly", str(second), parameters, options)
            after = world.EvaluateFormFactor("assembly", parameters, 0.15)
            options.seed = 8
            world.export_usd("assembly", str(different), parameters, options)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            self.assertNotEqual(first.read_bytes(), different.read_bytes())
            self.assertAlmostEqual(before, after, places=14)


if __name__ == "__main__":
    unittest.main()
