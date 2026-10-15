"""Independent native-program regression fixtures; no external executables needed."""
import json
from pathlib import Path
import unittest

import numpy as np
import pyseb


FIXTURES = Path(__file__).resolve().parents[3]/'Tests/data/hydration_reference'


class TestHydrationReferences(unittest.TestCase):
    def test_pepsi_hydration_from_atoms_to_proteins(self):
        fixtures = sorted(FIXTURES.glob('*.json'))
        self.assertTrue(fixtures,'Native hydration reference fixtures are missing')
        for path in fixtures:
            with self.subTest(fixture=path.stem):
                reference = json.loads(path.read_text())
                atoms = pyseb.PDBParser().parse_string(reference['pdb'])
                options = pyseb.AtomicSAXSOptions()
                options.hydration_model = pyseb.AtomicSAXSHydrationModel.SpatialShell
                options.hydration_water_model = pyseb.AtomicSAXSHydrationWaterModel.BondedWater
                options.include_hydration = True
                options.hydration_shell_contrast = .0167
                build = pyseb.AtomicSAXSCloudBuildOptions()
                build.calculate_solvent_accessibility = False
                build.hydration_shell_convention = pyseb.AtomicSAXSHydrationShellConvention.Pepsi3
                full = pyseb.AtomicSAXSCloudBuilder.build(atoms,
                    pyseb.AtomicSAXSParameterProfile.protein_defaults(),options,build)
                shell = [p for p in full.get_scatterers() if p.hydration_volume>0]
                self.assertEqual(len(shell),reference['pepsi_shell_points'])
                cloud = pyseb.AtomicSAXSCloud(shell,options)
                world = pyseb.World(path.stem)
                world.Add(cloud,'water')
                ours = np.array(world.EvaluateFormFactorUnnormalized('water',{},reference['q']))
                expected = np.array(reference['pepsi_hydration'])
                # Four- vs five-Gaussian atomic water factors remain distinct.
                self.assertLess(np.sqrt(np.mean((ours/expected-1)**2)),.002)
                self.assertAlmostEqual(ours[0]/expected[0],1.,delta=.0002)

    def test_foxs_hydration_with_native_exposure_controls(self):
        fixtures = sorted(FIXTURES.glob('*.json'))
        self.assertTrue(fixtures,'Native hydration reference fixtures are missing')
        for path in fixtures:
            reference = json.loads(path.read_text())
            atoms = pyseb.PDBParser().parse_string(reference['pdb'])
            self.assertEqual(len(atoms),len(reference['foxs_exposure']))
            for mode in ('default_fast','reciprocal'):
                with self.subTest(fixture=path.stem,mode=mode):
                    scatterers = []
                    for atom,exposure in zip(atoms,reference['foxs_exposure']):
                        point = pyseb.AtomicSAXSScatterer()
                        point.center = atom.coordinate
                        point.solvent_accessibility = exposure
                        scatterers.append(point)
                    options = pyseb.AtomicSAXSOptions()
                    options.include_hydration = True
                    options.hydration_shell_contrast = reference['foxs_water_f0'][mode]/30.
                    options.hydration_water_model = getattr(pyseb.AtomicSAXSHydrationWaterModel,
                        'FoXSReciprocal' if mode=='reciprocal' else 'FoXSFast')
                    world = pyseb.World(path.stem+mode)
                    world.Add(pyseb.AtomicSAXSCloud(scatterers,options),'water')
                    ours = np.array(world.EvaluateFormFactorUnnormalized('water',{},reference['q']))
                    expected = np.array(reference['foxs_hydration'][mode])
                    tolerance = 1e-6 if mode=='reciprocal' else .001
                    self.assertLess(np.sqrt(np.mean((ours/expected-1)**2)),tolerance)


if __name__ == '__main__':
    unittest.main()
