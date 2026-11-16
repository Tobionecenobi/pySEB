import hashlib
import math
import os
from pathlib import Path
import tempfile
import unittest
import urllib.request

import numpy as np
import pyseb


class TestStructures(unittest.TestCase):
    def test_debye_sphere_cloud_matches_single_sphere(self):
        cloud = pyseb.DebyeSphereCloud(
            [pyseb.SphereScatterer(0.0, 0.0, 0.0, 2.0, 3.0)]
        )
        world = pyseb.World()
        world.Add(cloud, "cloud")
        world.Add("SolidSphere", "sphere")
        params = {"beta_sphere": 3.0, "R_sphere": 2.0}

        for q in (0.0, 0.05, 0.2, 0.7):
            self.assertAlmostEqual(
                world.EvaluateFormFactor("cloud", params, q),
                world.EvaluateFormFactor("sphere", params, q),
                places=10,
            )
            self.assertAlmostEqual(
                world.EvaluateFormFactorAmplitude("cloud.center", params, q),
                world.EvaluateFormFactorAmplitude("sphere.center", params, q),
                places=10,
            )

        self.assertAlmostEqual(
            world.EvaluateRadiusOfGyration2("cloud", params),
            12.0 / 5.0,
        )

    def test_debye_reference_distances_and_vector_evaluation(self):
        cloud = pyseb.DebyeSphereCloud(
            [
                pyseb.SphereScatterer(0.0, 0.0, 0.0, 0.0, 1.0),
                pyseb.SphereScatterer(2.0, 0.0, 0.0, 0.0, 1.0),
            ]
        )
        cloud.addReferencePoint("left", 0.0, 0.0, 0.0)
        cloud.addReferencePoint("right", 2.0, 0.0, 0.0)
        world = pyseb.World()
        world.Add(cloud, "cloud")

        q_values = [0.1, 0.2, 0.4]
        vector = world.EvaluateFormFactor("cloud", {}, q_values)
        scalar = [
            world.EvaluateFormFactor("cloud", {}, q)
            for q in q_values
        ]
        self.assertEqual(len(vector), len(q_values))
        for actual, expected in zip(vector, scalar):
            self.assertAlmostEqual(actual, expected)

        self.assertAlmostEqual(
            world.EvaluateSMSDRef2Ref("cloud.left", "cloud.right", {}),
            4.0,
        )

    def test_debye_zero_total_beta_only_allows_unnormalized_scattering(self):
        cloud = pyseb.DebyeSphereCloud(
            [
                pyseb.SphereScatterer(0.0, 0.0, 0.0, 0.0, 1.0),
                pyseb.SphereScatterer(1.0, 0.0, 0.0, 0.0, -1.0),
            ]
        )
        world = pyseb.World()
        world.Add(cloud, "cloud")

        self.assertTrue(
            math.isfinite(
                world.EvaluateFormFactorUnnormalized("cloud", {}, 0.2)
            )
        )
        with self.assertRaises(RuntimeError):
            world.EvaluateFormFactor("cloud", {}, 0.2)

    def test_pdb_parser_profile_and_cloud_builder_pipeline(self):
        pdb = (
            "ATOM      1  N   ALA A   1      11.104  13.207   9.234"
            "  0.50 20.00           N  \n"
            "ATOM      2  CA AALA A   1      12.560  13.100   9.500"
            "  1.00 21.00           C  \n"
            "ATOM      3  CA BALA A   1      12.700  13.200   9.600"
            "  1.00 22.00           C  \n"
        )

        parser = pyseb.PDBParser()
        atoms = parser.parse_string(pdb)
        self.assertEqual(len(atoms), 2)
        self.assertEqual(atoms[0].element, "N")
        self.assertEqual(atoms[1].alternate_location, "A")

        profile = pyseb.AtomParameterProfile()
        profile.set_element("N", radius=1.0, beta=2.0)
        profile.set_element("C", radius=1.5, beta=-1.0)
        profile.set_atom("ALA", "CA", radius=1.7, beta=3.0)

        build_options = pyseb.AtomCloudBuildOptions()
        build_options.reference_atom_serials = {"n_term": 1}
        cloud = pyseb.AtomCloudBuilder.build(
            atoms,
            profile,
            build_options,
        )

        world = pyseb.World()
        world.Add(cloud, "protein")
        self.assertEqual(cloud.scattererCount(), 2)
        self.assertAlmostEqual(
            world.EvaluateFormFactorUnnormalized("protein", {}, 0.0),
            16.0,
        )
        self.assertIn("n_term", cloud.getReferenceCoordinates())

        with tempfile.TemporaryDirectory() as tmpdir:
            filename = Path(tmpdir) / "protein.pdb"
            filename.write_text(pdb)
            loaded_cloud = pyseb.StructureCloudLoader.load_pdb(
                str(filename),
                profile,
                build_options=build_options,
            )
            loaded_world = pyseb.World()
            loaded_world.Add(loaded_cloud, "loaded")
            self.assertAlmostEqual(
                loaded_world.EvaluateFormFactorUnnormalized(
                    "loaded", {}, 0.0
                ),
                16.0,
            )

    @unittest.skipUnless(
        os.environ.get("PYSEB_RUN_NETWORK_TESTS") == "1",
        "Set PYSEB_RUN_NETWORK_TESTS=1 to run network integration tests",
    )

    def test_downloads_parses_and_removes_official_rcsb_pdb(self):
        url = "https://files.rcsb.org/download/1CRN.pdb"
        expected_sha256 = (
            "42199a30a0701864a2a5cc76cd7f35cc"
            "544cd0e65fbcf63e03c166543249b811"
        )

        downloaded_path = None
        with tempfile.TemporaryDirectory() as tmpdir:
            downloaded_path = Path(tmpdir) / "1CRN.pdb"
            with urllib.request.urlopen(url, timeout=30) as response:
                contents = response.read()
            downloaded_path.write_bytes(contents)

            self.assertEqual(
                hashlib.sha256(contents).hexdigest(),
                expected_sha256,
            )

            parser = pyseb.PDBParser()
            atoms = parser.parse_file(str(downloaded_path))
            self.assertEqual(len(atoms), 327)

            profile = pyseb.AtomParameterProfile()
            profile.set_element("C", radius=1.70, beta=1.0)
            profile.set_element("N", radius=1.55, beta=1.0)
            profile.set_element("O", radius=1.52, beta=1.0)
            profile.set_element("S", radius=1.80, beta=1.0)

            cloud = pyseb.AtomCloudBuilder.build(atoms, profile)
            world = pyseb.World()
            world.Add(cloud, "crambin")
            self.assertTrue(
                math.isfinite(
                    world.EvaluateFormFactor("crambin", {}, 0.1)
                )
            )

        self.assertIsNotNone(downloaded_path)
        self.assertFalse(downloaded_path.exists())

    def test_python_mixed_cloud_integrated_and_analytic_structure(self):
        cloud = pyseb.DebyeSphereCloud(
            [
                pyseb.SphereScatterer(0.0, 0.0, 0.0, 0.0, 1.0),
                pyseb.SphereScatterer(1.0, 0.0, 0.0, 0.0, 1.0),
            ]
        )
        cloud.addReferencePoint("join", 0.0, 0.0, 0.0)

        world = pyseb.World()
        graph = world.Add(cloud, "cloud")
        world.Link(
            "SolidCylinder",
            "cylinder.center",
            "cloud.join",
        )
        world.Link(
            "SolidSphere",
            "sphere.center",
            "cylinder.center",
        )
        world.Add(graph, "mixed")

        params = {
            "beta_cylinder": 3.0,
            "R_cylinder": 1.0,
            "L_cylinder": 2.0,
            "beta_sphere": 4.0,
            "R_sphere": 1.5,
        }
        values = pyseb.evaluate_form_factor(
            world,
            "mixed",
            params,
            np.array([0.0, 0.1, 0.3]),
        )

        self.assertEqual(values[0], 1.0)
        self.assertTrue(np.all(np.isfinite(values)))
        self.assertTrue(
            math.isfinite(
                world.EvaluateRadiusOfGyration2("mixed", params)
            )
        )


if __name__ == "__main__":
    unittest.main()
