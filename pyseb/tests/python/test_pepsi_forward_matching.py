"""Check the signed forward-matching diagnostic without a Pepsi executable."""
import unittest

import pyseb
from Examples.compare_atomic_saxs_pepsi import (
    forward_amplitude, hydration_contrast_multiplier,
)


class TestPepsiForwardMatching(unittest.TestCase):
    def test_signed_amplitude_branches(self):
        # Include a hydration-induced sign crossing and the review counterexample.
        for dry, shell in [(-2., 12.), (2., 12.), (-2., -12.), (2., -12.)]:
            multiplier = hydration_contrast_multiplier(dry, dry + shell, 4.)
            adjusted = dry + multiplier * shell
            self.assertAlmostEqual(adjusted / dry, 2.)
            self.assertAlmostEqual(adjusted**2 / dry**2, 4.)
        self.assertAlmostEqual(hydration_contrast_multiplier(-2., 10., 4.), -1/6)

    def test_native_negative_dry_amplitude(self):
        # Fraser CH3 is a supported group whose solvent-subtracted amplitude is negative.
        atom = pyseb.AtomRecord()
        atom.element, atom.residue_name, atom.atom_name = 'C', 'ALA', 'CB'
        parameters = pyseb.AtomicSAXSParameterProfile.protein_defaults().resolve(atom)
        point = pyseb.AtomicSAXSScatterer()
        point.vacuum = parameters.vacuum
        point.implicit_hydrogens = parameters.implicit_hydrogens
        point.excluded_volume = parameters.excluded_volume
        options = pyseb.AtomicSAXSOptions()
        options.include_hydration = False
        cloud = pyseb.AtomicSAXSCloud([point], options)
        amplitude = forward_amplitude(cloud)
        self.assertLess(amplitude, 0.)
        world = pyseb.World('negative_forward')
        world.Add(cloud, 'cloud')
        intensity = world.EvaluateFormFactorUnnormalized('cloud', {}, [0.])[0]
        self.assertAlmostEqual(amplitude**2, intensity)

    def test_undefined_matching_rejected(self):
        for dry, wet, ratio in [(0., 1., 4.), (2., 2., 4.),
                                (2., 3., 0.), (2., 3., -1.),
                                (2., float('nan'), 4.)]:
            with self.assertRaises(ValueError):
                hydration_contrast_multiplier(dry, wet, ratio)


if __name__ == '__main__':
    unittest.main()
