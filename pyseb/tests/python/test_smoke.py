import unittest

import pyseb


class TestSmoke(unittest.TestCase):
    def test_import_and_world_creation(self):
        world = pyseb.World()
        self.assertIsNotNone(world)

    def test_build_diblock_structure_and_form_factor(self):
        world = pyseb.World()
        graph_id = world.Add("GaussianPolymer", "poly1")
        world.Link("GaussianPolymer", "poly2.end1", "poly1.end2")
        world.Add(graph_id, "diblockcopolymer")

        form_factor = world.FormFactor("diblockcopolymer")
        self.assertIsNotNone(form_factor)


if __name__ == "__main__":
    unittest.main()
