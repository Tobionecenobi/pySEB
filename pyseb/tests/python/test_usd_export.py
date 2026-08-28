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


def _world_point(pose, point):
    translation, rotation = pose
    rotated = _rotate(rotation, point)
    return tuple(translation[i] + rotated[i] for i in range(3))


def _segment_distance(first, second):
    p1, q1 = first
    p2, q2 = second
    subtract = lambda a, b: tuple(a[i] - b[i] for i in range(3))
    dot = lambda a, b: sum(a[i] * b[i] for i in range(3))
    add_scaled = lambda a, b, scale: tuple(a[i] + scale * b[i] for i in range(3))
    u, v, w = subtract(q1, p1), subtract(q2, p2), subtract(p1, p2)
    a, b, c = dot(u, u), dot(u, v), dot(v, v)
    d, e = dot(u, w), dot(v, w)
    if a < 1e-15 and c < 1e-15:
        return math.dist(p1, p2)
    if a < 1e-15:
        return math.dist(p1, add_scaled(p2, v, max(0.0, min(1.0, e / c))))
    if c < 1e-15:
        return math.dist(add_scaled(p1, u, max(0.0, min(1.0, -d / a))), p2)
    denominator = a*c - b*b
    s = 0.0 if denominator < 1e-15 else max(0.0, min(1.0, (b*e-c*d)/denominator))
    t = max(0.0, min(1.0, (b*s+e)/c))
    s = max(0.0, min(1.0, (b*t-d)/a))
    return math.dist(add_scaled(p1, u, s), add_scaled(p2, v, t))


def _minimum_nonlinked_rod_distance(document, rod_count):
    segments = {}
    for index in range(1, rod_count + 1):
        pose = _pose(document, f"rod{index}")
        segments[index] = (
            _world_point(pose, (0.0, 0.0, -0.5)),
            _world_point(pose, (0.0, 0.0, 0.5)),
        )
    linked = {
        frozenset((int(first), int(second)))
        for first, second in re.findall(
            r'/rod(\d+)/ref_[^>]+>, </[^>]+/rod(\d+)/ref_', document
        )
    }
    return min(
        _segment_distance(segments[first], segments[second])
        for first in segments
        for second in segments
        if first < second and frozenset((first, second)) not in linked
    )


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

    def test_readable_layout_improves_clearance_and_preserves_links(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            world = pyseb.World("readable")
            graph = world.Add("ThinRod", "rod1", "rod")
            parents = ["rod1"]
            next_index = 2
            for generation in range(1, 4):
                children = []
                for parent in parents:
                    for branch in (1, 2):
                        child = f"rod{next_index}"
                        graph = world.Link(
                            "ThinRod",
                            f"{child}.contour#attachment",
                            f"{parent}.contour#g{generation}_branch{branch}",
                            "rod",
                        )
                        children.append(child)
                        next_index += 1
                parents = children
            world.Add(graph, "dendrimer")

            random_output = directory / "random.usda"
            greedy_output = directory / "greedy.usda"
            readable_output = directory / "readable.usda"
            repeated_output = directory / "repeated.usda"
            options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
            options.seed = 42
            options.curve_samples = 2
            world.export_usd("dendrimer", str(random_output), {"L_rod": 1.0}, options)
            options.layout_mode = pyseb.LayoutMode.Readable
            options.orientation_trials = 64
            options.minimum_clearance = 0.08
            options.relaxation_sweeps = 0
            world.export_usd("dendrimer", str(greedy_output), {"L_rod": 1.0}, options)
            options.relaxation_sweeps = 4
            world.export_usd("dendrimer", str(readable_output), {"L_rod": 1.0}, options)
            world.export_usd("dendrimer", str(repeated_output), {"L_rod": 1.0}, options)

            random_document = random_output.read_text(encoding="utf-8")
            greedy_document = greedy_output.read_text(encoding="utf-8")
            readable_document = readable_output.read_text(encoding="utf-8")
            self.assertEqual(readable_output.read_bytes(), repeated_output.read_bytes())
            self.assertIn('pyseb:layoutMode = "readable"', readable_document)
            self.assertIn('pyseb:overlapPolicy = "best_effort_minimized"', readable_document)
            self.assertGreater(
                _minimum_nonlinked_rod_distance(readable_document, 15),
                _minimum_nonlinked_rod_distance(greedy_document, 15),
            )
            self.assertGreater(
                _minimum_nonlinked_rod_distance(greedy_document, 15),
                _minimum_nonlinked_rod_distance(random_document, 15),
            )

            relationships = re.findall(
                r'rel pyseb:link_\d+ = \[</dendrimer/(rod\d+)/ref_([^>]+)>, '
                r'</dendrimer/(rod\d+)/ref_([^>]+)>\]',
                readable_document,
            )
            self.assertEqual(len(relationships), 14)
            for first, first_reference, second, second_reference in relationships:
                first_world = _world_point(
                    _pose(readable_document, first),
                    _reference(readable_document, first, first_reference),
                )
                second_world = _world_point(
                    _pose(readable_document, second),
                    _reference(readable_document, second, second_reference),
                )
                for actual, expected in zip(first_world, second_world):
                    self.assertAlmostEqual(actual, expected, places=12)

    def test_readable_layout_options_are_validated(self):
        with tempfile.TemporaryDirectory() as directory:
            world = pyseb.World("invalid_readable")
            world.Add("ThinRod", "rod")
            options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
            options.orientation_trials = 0
            with self.assertRaisesRegex(RuntimeError, "orientationTrials"):
                world.export_usd(
                    "rod", str(Path(directory) / "invalid.usda"), {"L_rod": 1.0}, options
                )

    def test_reference_variants_use_both_ends_of_cylinder_chain(self):
        with tempfile.TemporaryDirectory() as directory:
            world = pyseb.World("cylinder_ends")
            graph = world.Add("SolidCylinder", "cylinder1", "cylinder")
            for index in range(2, 21):
                graph = world.Link(
                    "SolidCylinder",
                    f"cylinder{index}.ends#bottom#joint{index - 1}",
                    f"cylinder{index - 1}.ends#top#joint{index - 1}",
                    "cylinder",
                )
            world.Add(graph, "chain")
            output = Path(directory) / "cylinders.usda"
            options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
            options.seed = 42
            options.surface_samples = 8
            world.export_usd(
                "chain",
                str(output),
                {"R_cylinder": 0.25, "L_cylinder": 3.0},
                options,
            )
            document = output.read_text(encoding="utf-8")
            for index in range(1, 21):
                if index > 1:
                    bottom = _reference(
                        document, f"cylinder{index}", f"ends_bottom_joint{index - 1}"
                    )
                    self.assertAlmostEqual(bottom[2], -1.5)
                if index < 20:
                    top = _reference(
                        document, f"cylinder{index}", f"ends_top_joint{index}"
                    )
                    self.assertAlmostEqual(top[2], 1.5)

            relationships = re.findall(
                r'rel pyseb:link_\d+ = \[</chain/(cylinder\d+)/ref_([^>]+)>, '
                r'</chain/(cylinder\d+)/ref_([^>]+)>\]',
                document,
            )
            self.assertEqual(len(relationships), 19)
            for first, first_reference, second, second_reference in relationships:
                first_world = _world_point(
                    _pose(document, first), _reference(document, first, first_reference)
                )
                second_world = _world_point(
                    _pose(document, second), _reference(document, second, second_reference)
                )
                for actual, expected in zip(first_world, second_world):
                    self.assertAlmostEqual(actual, expected, places=12)

    def test_unknown_reference_variant_fails_atomically(self):
        with tempfile.TemporaryDirectory() as directory:
            world = pyseb.World("invalid_variant")
            graph = world.Add("SolidCylinder", "cylinder1", "cylinder")
            graph = world.Link(
                "SolidCylinder",
                "cylinder2.ends#side#joint1",
                "cylinder1.ends#top#joint1",
                "cylinder",
            )
            world.Add(graph, "chain")
            output = Path(directory) / "invalid.usda"
            options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
            with self.assertRaisesRegex(RuntimeError, "unknown visualization reference variant 'side'"):
                world.export_usd(
                    "chain",
                    str(output),
                    {"R_cylinder": 0.25, "L_cylinder": 3.0},
                    options,
                )
            self.assertFalse(output.exists())

    def test_spheroid_chain_uses_north_and_south_poles(self):
        with tempfile.TemporaryDirectory() as directory:
            world = pyseb.World("spheroid_poles")
            graph = world.Add("Spheroid", "spheroid1", "spheroid")
            for index in range(2, 6):
                graph = world.Link(
                    "Spheroid",
                    f"spheroid{index}.south",
                    f"spheroid{index - 1}.north",
                    "spheroid",
                )
            world.Add(graph, "chain")
            output = Path(directory) / "spheroids.usda"
            options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
            options.seed = 42
            options.surface_samples = 8
            world.export_usd(
                "chain",
                str(output),
                {"a_spheroid": 0.4, "c_spheroid": 1.2},
                options,
            )
            document = output.read_text(encoding="utf-8")
            for index in range(1, 6):
                self.assertEqual(_reference(document, f"spheroid{index}", "north"), (0.0, 0.0, 1.2))
                self.assertEqual(_reference(document, f"spheroid{index}", "south"), (0.0, 0.0, -1.2))
                self.assertEqual(_reference(document, f"spheroid{index}", "pole"), (0.0, 0.0, 1.2))

            relationships = re.findall(
                r'rel pyseb:link_\d+ = \[</chain/(spheroid\d+)/ref_([^>]+)>, '
                r'</chain/(spheroid\d+)/ref_([^>]+)>\]',
                document,
            )
            self.assertEqual(len(relationships), 4)
            for first, first_reference, second, second_reference in relationships:
                first_world = _world_point(
                    _pose(document, first), _reference(document, first, first_reference)
                )
                second_world = _world_point(
                    _pose(document, second), _reference(document, second, second_reference)
                )
                for actual, expected in zip(first_world, second_world):
                    self.assertAlmostEqual(actual, expected, places=12)

    def test_debye_cloud_reference_attaches_to_analytic_subunit(self):
        with tempfile.TemporaryDirectory() as directory:
            cloud = pyseb.DebyeSphereCloud(
                [
                    pyseb.SphereScatterer(0.0, 0.0, 0.0, 0.4, 1.0),
                    pyseb.SphereScatterer(0.8, 0.2, 0.0, 0.3, 0.5),
                    pyseb.SphereScatterer(1.6, -0.1, 0.2, 0.0, -0.25),
                ]
            )
            cloud.addReferencePoint("join", 2.0, 0.0, 0.0)
            world = pyseb.World("cloud_attachment")
            graph = world.Add(cloud, "cloud")
            world.Link("ThinRod", "rod.end1", "cloud.join", "rod")
            world.Add(graph, "assembly")

            output = Path(directory) / "cloud.usda"
            options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
            options.seed = 42
            options.curve_samples = 8
            options.zero_radius_marker_size = 0.125
            options.color_overrides = {"cloud": (0.2, 0.6, 0.8)}
            parameters = {"L_rod": 2.0, "beta_rod": 1.0}
            before = world.EvaluateFormFactor("assembly", parameters, 0.2)
            world.export_usd("assembly", str(output), parameters, options)
            after = world.EvaluateFormFactor("assembly", parameters, 0.2)
            self.assertAlmostEqual(before, after, places=14)
            document = output.read_text(encoding="utf-8")

            self.assertEqual(_reference(document, "cloud", "join"), (2.0, 0.0, 0.0))
            cloud_world = _world_point(
                _pose(document, "cloud"), _reference(document, "cloud", "join")
            )
            rod_world = _world_point(
                _pose(document, "rod"), _reference(document, "rod", "end1")
            )
            for actual, expected in zip(cloud_world, rod_world):
                self.assertAlmostEqual(actual, expected, places=12)

            self.assertIn('def PointInstancer "cloud"', document)
            self.assertIn("custom float[] pyseb:radius = [0.40000000000000002, 0.29999999999999999, 0, ]", document)
            self.assertIn("custom float[] pyseb:beta = [1, 0.5, -0.25, ]", document)
            self.assertIn("custom int[] pyseb:index = [0, 1, 2, ]", document)
            self.assertIn("(0.125, 0.125, 0.125)", document)
            self.assertIn('custom string pyseb:model_id = "pyseb/DebyeSphereCloud"', document)
            self.assertIn("rel pyseb:link_0", document)


if __name__ == "__main__":
    unittest.main()
