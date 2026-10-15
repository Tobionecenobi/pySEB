import math
import unittest

import numpy as np

import pyseb


def carbon_parameters():
    atom = pyseb.AtomRecord()
    atom.serial = 1
    atom.element = "C"
    atom.residue_name = "ALA"
    atom.atom_name = "C"
    return pyseb.AtomicSAXSParameterProfile.element_defaults().resolve(atom)


def scatterer(x, parameters, occupancy=1.0):
    value = pyseb.AtomicSAXSScatterer()
    value.center = pyseb.CartesianPoint3D(x, 0.0, 0.0)
    value.vacuum = parameters.vacuum
    value.excluded_volume = parameters.excluded_volume
    value.occupancy = occupancy
    value.element = "C"
    return value


def atom_record(
    serial,
    x=0.0,
    y=0.0,
    z=0.0,
    occupancy=1.0,
    residue="ALA",
    atom_name="CA",
    element="C",
):
    atom = pyseb.AtomRecord()
    atom.serial = serial
    atom.coordinate = pyseb.CartesianPoint3D(x, y, z)
    atom.occupancy = occupancy
    atom.residue_name = residue
    atom.atom_name = atom_name
    atom.element = element
    return atom


class TestAtomicSAXSCloud(unittest.TestCase):
    def test_profile_surface_radii_and_probe_match_two_sphere_exposure(self):
        atoms = [atom_record(1),atom_record(2,x=2.)]
        profile = pyseb.AtomicSAXSParameterProfile.protein_defaults()
        p = profile.resolve(atoms[0])
        p.sasa_radius = 4.
        profile.set_atom('ALA','CA',p)
        build = pyseb.AtomicSAXSCloudBuildOptions()
        build.sasa_use_profile_radii = True
        build.sasa_probe_radius = 1.8
        build.sasa_point_count = 4000
        cloud = pyseb.AtomicSAXSCloudBuilder.build(atoms,profile,pyseb.AtomicSAXSOptions(),build)
        # Two equal expanded spheres: exposed fraction is 1/2+d/(4R).
        expected = .5+2/(4*(4+1.8))
        for scatterer in cloud.get_scatterers():
            self.assertAlmostEqual(scatterer.solvent_accessibility,expected,delta=.003)

    def test_pepsi_shell_convention_has_fixed_grid_and_weight_mapping(self):
        atom = atom_record(1, atom_name='C')
        options = pyseb.AtomicSAXSOptions()
        options.hydration_model = pyseb.AtomicSAXSHydrationModel.SpatialShell
        options.hydration_water_model = pyseb.AtomicSAXSHydrationWaterModel.BondedWater
        options.include_hydration = True
        options.hydration_shell_contrast = .0167
        build = pyseb.AtomicSAXSCloudBuildOptions()
        build.hydration_shell_convention = pyseb.AtomicSAXSHydrationShellConvention.Pepsi3
        profile = pyseb.AtomicSAXSParameterProfile.protein_defaults()
        for atom_list, expected_count in [([atom],24),([atom,atom_record(2,x=2.,atom_name='C')],32)]:
            cloud = pyseb.AtomicSAXSCloudBuilder.build(atom_list,profile,options,build)
            shell = [p for p in cloud.get_scatterers() if p.hydration_volume > 0]
            self.assertEqual(len(shell),expected_count)
            for point in shell:
                self.assertAlmostEqual(point.hydration_volume,27*.334,places=12)
        # All isolated shell nodes are the 24 permutations of (+/-4.5,
        # +/-1.5,+/-1.5), a geometry independent of any fitted intensity.
        cloud = pyseb.AtomicSAXSCloudBuilder.build([atom],profile,options,build)
        shell = [p for p in cloud.get_scatterers() if p.hydration_volume > 0]
        for p in shell:
            self.assertEqual(sorted([abs(p.center.x),abs(p.center.y),abs(p.center.z)]),
                             [1.5,1.5,4.5])
        zero = profile.resolve(atom)
        zero.sasa_radius = 0.
        profile.set_atom('ALA','C',zero)
        with self.assertRaises(RuntimeError):
            pyseb.AtomicSAXSCloudBuilder.build([atom],profile,options,build)

    def test_hydration_water_conventions_and_guinier_moments(self):
        point = pyseb.AtomicSAXSScatterer()
        point.hydration_volume = 1.0
        options = pyseb.AtomicSAXSOptions()
        options.hydration_model = pyseb.AtomicSAXSHydrationModel.SpatialShell
        options.include_hydration = True
        options.hydration_shell_contrast = 1.0
        oxygen = pyseb.AtomicSAXSParameterProfile.element_defaults().resolve(
            atom_record(1, element='O', atom_name='O')).vacuum
        h = pyseb.AtomicSAXSCloud.hydrogen_coefficients()
        for model in ('NeutralWater','BondedWater','FoXSReciprocal','FoXSFast'):
            options.hydration_water_model = getattr(pyseb.AtomicSAXSHydrationWaterModel,model)
            cloud = pyseb.AtomicSAXSCloud([point], options)
            world = pyseb.World(model)
            world.Add(cloud,'water')
            for q in (.0,.2,.5):
                beta = oxygen.value_at_zero()+2*h.value_at_zero()
                expected = oxygen.evaluate(q)+2*h.evaluate(q)*(
                    float(np.sinc(q*.967/np.pi)) if model=='BondedWater' else 1.)
                if model=='FoXSReciprocal':
                    expected -= .334*(9.13*math.exp(-q*q*9.13**(2/3)/(16*math.pi))+
                                     2*5.15*math.exp(-q*q*5.15**(2/3)/(16*math.pi)))
                    beta -= .334*(9.13+2*5.15)
                expected = math.exp(-.115*q*q) if model=='FoXSFast' else expected/beta
                self.assertAlmostEqual(cloud.atomic_amplitude(0,q),expected,places=13)
            step = .002
            slope = lambda q: 3*(1-world.EvaluateFormFactor('water',{},q))/q**2
            self.assertAlmostEqual(world.EvaluateRadiusOfGyration2('water',{}),
                                   (4*slope(step)-slope(2*step))/3, delta=1e-7)

    def test_spatial_points_include_atom_shell_interference_and_references(self):
        atom = pyseb.AtomicSAXSScatterer()
        atom.vacuum = pyseb.CromerMannCoefficients([0.]*4, [0.]*4, 2.)
        left, right = pyseb.AtomicSAXSScatterer(), pyseb.AtomicSAXSScatterer()
        left.center = pyseb.CartesianPoint3D(-2., 0., 0.)
        right.center = pyseb.CartesianPoint3D(2., 0., 0.)
        left.hydration_volume = right.hydration_volume = 2.
        options = pyseb.AtomicSAXSOptions()
        options.hydration_model = pyseb.AtomicSAXSHydrationModel.SpatialShell
        options.include_hydration = True
        options.hydration_shell_contrast = .5
        cloud = pyseb.AtomicSAXSCloud([atom, left, right], options)
        cloud.add_reference_point('end', 2., 0., 0.)
        world = pyseb.World('mixed-shell')
        world.Add(cloud, 'cloud')
        for contrast in (.5, -.25):
            cloud.set_hydration(True, contrast, 999.)
            q = .4
            # Obtain the water factor independently from the elemental profile.
            oxygen = pyseb.AtomicSAXSParameterProfile.element_defaults().resolve(
                atom_record(2, element='O', atom_name='O')).vacuum
            h = pyseb.AtomicSAXSCloud.hydrogen_coefficients()
            water = (oxygen.evaluate(q)+2*h.evaluate(q))/(oxygen.value_at_zero()+2*h.value_at_zero())
            b = 2*contrast*water
            sinc = lambda x: math.sin(x)/x
            expected = 4+2*b*b*(1+sinc(4*q))+8*b*sinc(2*q)
            self.assertAlmostEqual(world.EvaluateFormFactorUnnormalized('cloud', {}, q), expected, places=12)
            self.assertAlmostEqual(world.EvaluateFormFactorAmplitude('cloud.center', {}, q),
                                   (2+2*b*sinc(2*q))/(2+4*contrast), places=12)
            self.assertAlmostEqual(world.EvaluatePhaseFactor('cloud.center', 'cloud.end', {}, q),
                                   sinc(2*q), places=12)

    def test_spatial_shell_sphere_and_attachment_amplitude(self):
        # A zero-scattering central atom leaves only a uniform 3-to-6 A shell.
        profile = pyseb.AtomicSAXSParameterProfile.protein_defaults()
        parameters = profile.resolve(atom_record(1))
        parameters.vacuum = pyseb.CromerMannCoefficients()
        parameters.excluded_volume = 0
        parameters.implicit_hydrogens = 0
        profile.set_atom('ALA', 'CA', parameters)
        options = pyseb.AtomicSAXSOptions()
        options.hydration_model = pyseb.AtomicSAXSHydrationModel.SpatialShell
        options.include_hydration = True
        options.hydration_shell_contrast = 0.02
        build = pyseb.AtomicSAXSCloudBuildOptions()
        build.hydration_shell_spacing = 1.0
        cloud = pyseb.AtomicSAXSCloudBuilder.build([atom_record(1)], profile, options, build)
        points = [p for p in cloud.get_scatterers() if p.hydration_volume > 0]
        volume = sum(p.hydration_volume for p in points)
        exact_volume = 4*math.pi/3*(6**3-3**3)
        self.assertLess(abs(volume/exact_volume-1), 0.02)
        world = pyseb.World('shell')
        world.Add(cloud, 'shell')
        oxygen = pyseb.AtomicSAXSParameterProfile.element_defaults().resolve(
            atom_record(2, element='O', atom_name='O')).vacuum
        h = pyseb.AtomicSAXSCloud.hydrogen_coefficients()
        def water(q):
            return (oxygen.evaluate(q)+2*h.evaluate(q))/(oxygen.value_at_zero()+2*h.value_at_zero())
        def sphere(x):
            return 3*(math.sin(x)-x*math.cos(x))/x**3
        for q in (0.1, 0.3, 0.5):
            amplitude = world.EvaluateFormFactorAmplitude('shell.center', {}, q)
            expected = water(q)*(6**3*sphere(6*q)-3**3*sphere(3*q))/(6**3-3**3)
            self.assertAlmostEqual(amplitude, expected, delta=0.014)
            self.assertAlmostEqual(world.EvaluateFormFactor('shell', {}, q),
                                   expected**2, delta=0.025)
        q = 1e-3
        rg2 = world.EvaluateRadiusOfGyration2('shell', {})
        slope = lambda x: 3*(1-world.EvaluateFormFactor('shell', {}, x))/x**2
        self.assertAlmostEqual(rg2, (4*slope(q)-slope(2*q))/3, delta=1e-6)
        self.assertAlmostEqual(world.EvaluateFormFactorUnnormalized('shell', {}, 0),
                               (volume*0.02)**2, delta=1e-8)
        build.hydration_shell_spacing = 2.0
        coarse = pyseb.AtomicSAXSCloudBuilder.build([atom_record(1)], profile, options, build)
        coarse_volume = sum(p.hydration_volume for p in coarse.get_scatterers())
        self.assertLess(abs(volume-exact_volume), abs(coarse_volume-exact_volume))
        coarse_world = pyseb.World('coarse')
        coarse_world.Add(coarse, 'shell')
        expected = water(0.3)*(6**3*sphere(1.8)-3**3*sphere(0.9))/(6**3-3**3)
        self.assertLess(abs(world.EvaluateFormFactorAmplitude('shell.center', {}, .3)-expected),
                        abs(coarse_world.EvaluateFormFactorAmplitude('shell.center', {}, .3)-expected))

    def test_spatial_shell_geometry_switches_and_translation(self):
        profile = pyseb.AtomicSAXSParameterProfile.protein_defaults()
        options = pyseb.AtomicSAXSOptions()
        options.hydration_model = pyseb.AtomicSAXSHydrationModel.SpatialShell
        options.include_hydration = True
        options.hydration_shell_contrast = 0.0167
        build = pyseb.AtomicSAXSCloudBuildOptions()
        atoms = [atom_record(1), atom_record(2, x=2.0)]
        cloud = pyseb.AtomicSAXSCloudBuilder.build(atoms, profile, options, build)
        points = [p for p in cloud.get_scatterers() if p.hydration_volume > 0]
        self.assertTrue(points)
        for p in points:
            distances = [math.dist((p.center.x,p.center.y,p.center.z),
                                  (a.coordinate.x,a.coordinate.y,a.coordinate.z)) for a in atoms]
            self.assertGreaterEqual(min(distances), 3.0-1e-12)
            self.assertLess(min(distances), 6.0)
        world = pyseb.World('switches')
        world.Add(cloud, 'protein')
        dry_options = cloud.get_options()
        dry_options.include_hydration = False
        cloud.set_options(dry_options)
        dry = world.EvaluateFormFactor('protein', {}, 0.3)
        baseline = pyseb.AtomicSAXSCloudBuilder.build(atoms, profile)
        baseline_world = pyseb.World('baseline')
        baseline_world.Add(baseline, 'protein')
        self.assertAlmostEqual(dry, baseline_world.EvaluateFormFactor('protein', {}, 0.3), places=13)
        cloud.set_options(options)
        wet = world.EvaluateFormFactor('protein', {}, 0.3)
        self.assertNotAlmostEqual(dry, wet)
        for atom in atoms:
            atom.coordinate = pyseb.CartesianPoint3D(atom.coordinate.x+17,
                                                     atom.coordinate.y-8, atom.coordinate.z+3)
        translated = pyseb.AtomicSAXSCloudBuilder.build(atoms, profile, options, build)
        other_world = pyseb.World('translated')
        other_world.Add(translated, 'protein')
        self.assertAlmostEqual(wet, other_world.EvaluateFormFactor('protein', {}, 0.3), places=12)
        self.assertAlmostEqual(translated.get_reference_coordinates()['center'].x, 18.)
        for bad in (0., -1., float('nan')):
            build.hydration_shell_spacing = bad
            with self.assertRaises(RuntimeError):
                pyseb.AtomicSAXSCloudBuilder.build(atoms, profile, options, build)
        with self.assertRaises(RuntimeError):
            baseline.set_options(options)

    def test_bonded_hydrogen_shell_matches_explicit_orientation_average(self):
        # Average one C-H pair over bond directions, with the heavy atom at
        # the reference origin. The amplitude is C(q) + H(q)*sinc(q*d).
        profile = pyseb.AtomicSAXSParameterProfile.protein_defaults()
        record = atom_record(1)
        parameters = profile.resolve(record)
        self.assertAlmostEqual(parameters.implicit_hydrogen_distance, 1.099)
        options = pyseb.AtomicSAXSOptions()
        options.include_excluded_solvent = False
        cloud = pyseb.AtomicSAXSCloudBuilder.build([record], profile, options)
        hydrogen = pyseb.AtomicSAXSCloud.hydrogen_coefficients()
        nodes, weights = np.polynomial.legendre.leggauss(64)
        for q in (0.0, 0.1, 0.5, 1.0):
            averaged = np.dot(weights, np.cos(q * 1.099 * nodes)) / 2
            expected = parameters.vacuum.evaluate(q) + hydrogen.evaluate(q) * averaged
            self.assertAlmostEqual(cloud.atomic_amplitude(0, q), expected, places=12)
        world = pyseb.World('bonded-hydrogen')
        world.Add(cloud, 'group')
        q = 1e-4
        numeric_rg2 = 3 * (1 - world.EvaluateFormFactor('group', {}, q)) / q**2
        self.assertAlmostEqual(world.EvaluateRadiusOfGyration2('group', {}),
                               numeric_rg2, delta=1e-6)
        for distance in (-1.0, float('nan')):
            parameters.implicit_hydrogen_distance = distance
            with self.assertRaises(RuntimeError):
                profile.set_atom('ALA', 'CA', parameters)

    def test_literature_volume_profiles_and_custom_sasa_radius(self):
        alanine_ca = atom_record(1)
        cysteine_sg = atom_record(
            2, residue="CYS", atom_name="SG", element="S"
        )
        serine_og = atom_record(
            3, residue="SER", atom_name="OG", element="O"
        )

        bondi = pyseb.AtomicSAXSParameterProfile.protein_defaults(
            pyseb.AtomicSAXSVolumeModel.BondiElement
        )
        fraser = pyseb.AtomicSAXSParameterProfile.protein_defaults()
        pontius = pyseb.AtomicSAXSParameterProfile.protein_defaults(
            pyseb.AtomicSAXSVolumeModel.Pontius1996
        )

        self.assertAlmostEqual(
            bondi.resolve(alanine_ca).excluded_volume,
            4.0 * math.pi * 1.7**3 / 3.0,
        )
        self.assertAlmostEqual(
            fraser.resolve(alanine_ca).excluded_volume, 21.59
        )
        self.assertAlmostEqual(
            fraser.resolve(cysteine_sg).excluded_volume, 25.01
        )
        fraser_elements = {
            ("PRO", "N", "N"): 2.49,
            ("ARG", "CZ", "C"): 16.44,
            ("ASP", "OD1", "O"): 9.13,
            ("MET", "SD", "S"): 19.86,
            ("LIG", "P", "P"): 5.73,
            ("LIG", "H1", "H"): 5.15,
        }
        for serial, ((residue, name, element), expected) in enumerate(
            fraser_elements.items(), start=10
        ):
            self.assertAlmostEqual(
                fraser.resolve(
                    atom_record(
                        serial,
                        residue=residue,
                        atom_name=name,
                        element=element,
                    )
                ).excluded_volume,
                expected,
            )
        self.assertAlmostEqual(
            pontius.resolve(alanine_ca).excluded_volume, 11.8
        )
        self.assertAlmostEqual(
            pontius.resolve(serine_og).excluded_volume, 23.9
        )

        custom = fraser.resolve(alanine_ca)
        custom.sasa_radius = 2.25
        pontius.set_atom("LIG", "QX", custom)
        resolved = pontius.resolve(
            atom_record(4, residue="LIG", atom_name="QX", element="C")
        )
        self.assertAlmostEqual(resolved.sasa_radius, 2.25)

        unsupported = atom_record(
            5, residue="LIG", atom_name="QX", element="X"
        )
        with self.assertRaisesRegex(RuntimeError, r"LIG:QX.*element X"):
            pyseb.AtomicSAXSParameterProfile().resolve(unsupported)

    def test_standard_residue_group_coverage_and_terminal_atoms(self):
        # Expected implicit hydrogens for every standard amino-acid heavy atom.
        # Zero-hydrogen atoms intentionally resolve through explicit element
        # rules when no more-specific group rule is necessary.
        standard_atoms = {
            "ALA": {"N": 1, "CA": 1, "C": 0, "O": 0, "CB": 3},
            "ARG": {"N": 1, "CA": 1, "C": 0, "O": 0, "CB": 2,
                    "CG": 2, "CD": 2, "NE": 1, "CZ": 0, "NH1": 2,
                    "NH2": 2},
            "ASN": {"N": 1, "CA": 1, "C": 0, "O": 0, "CB": 2,
                    "CG": 0, "OD1": 0, "ND2": 2},
            "ASP": {"N": 1, "CA": 1, "C": 0, "O": 0, "CB": 2,
                    "CG": 0, "OD1": 0, "OD2": 0},
            "CYS": {"N": 1, "CA": 1, "C": 0, "O": 0, "CB": 2,
                    "SG": 1},
            "GLN": {"N": 1, "CA": 1, "C": 0, "O": 0, "CB": 2,
                    "CG": 2, "CD": 0, "OE1": 0, "NE2": 2},
            "GLU": {"N": 1, "CA": 1, "C": 0, "O": 0, "CB": 2,
                    "CG": 2, "CD": 0, "OE1": 0, "OE2": 0},
            "GLY": {"N": 1, "CA": 2, "C": 0, "O": 0},
            "HIS": {"N": 1, "CA": 1, "C": 0, "O": 0, "CB": 2,
                    "CG": 0, "ND1": 1, "CD2": 1, "CE1": 1, "NE2": 0},
            "ILE": {"N": 1, "CA": 1, "C": 0, "O": 0, "CB": 1,
                    "CG1": 2, "CG2": 3, "CD1": 3},
            "LEU": {"N": 1, "CA": 1, "C": 0, "O": 0, "CB": 2,
                    "CG": 1, "CD1": 3, "CD2": 3},
            "LYS": {"N": 1, "CA": 1, "C": 0, "O": 0, "CB": 2,
                    "CG": 2, "CD": 2, "CE": 2, "NZ": 3},
            "MET": {"N": 1, "CA": 1, "C": 0, "O": 0, "CB": 2,
                    "CG": 2, "SD": 0, "CE": 3},
            "PHE": {"N": 1, "CA": 1, "C": 0, "O": 0, "CB": 2,
                    "CG": 0, "CD1": 1, "CD2": 1, "CE1": 1,
                    "CE2": 1, "CZ": 1},
            "PRO": {"N": 0, "CA": 1, "C": 0, "O": 0, "CB": 2,
                    "CG": 2, "CD": 2},
            "SER": {"N": 1, "CA": 1, "C": 0, "O": 0, "CB": 2,
                    "OG": 1},
            "THR": {"N": 1, "CA": 1, "C": 0, "O": 0, "CB": 1,
                    "OG1": 1, "CG2": 3},
            "TRP": {"N": 1, "CA": 1, "C": 0, "O": 0, "CB": 2,
                    "CG": 0, "CD1": 1, "CD2": 0, "NE1": 1,
                    "CE2": 0, "CE3": 1, "CZ2": 1, "CZ3": 1,
                    "CH2": 1},
            "TYR": {"N": 1, "CA": 1, "C": 0, "O": 0, "CB": 2,
                    "CG": 0, "CD1": 1, "CD2": 1, "CE1": 1,
                    "CE2": 1, "CZ": 0, "OH": 1},
            "VAL": {"N": 1, "CA": 1, "C": 0, "O": 0, "CB": 1,
                    "CG1": 3, "CG2": 3},
        }
        elements = {"N": "N", "O": "O", "OXT": "O", "SG": "S",
                    "SD": "S"}
        profile = pyseb.AtomicSAXSParameterProfile.protein_defaults()
        serial = 1
        for residue, atoms in standard_atoms.items():
            for atom_name, hydrogen_count in atoms.items():
                element = elements.get(atom_name, atom_name[0])
                atom = atom_record(
                    serial, residue=residue, atom_name=atom_name,
                    element=element,
                )
                self.assertTrue(profile.has_parameters(atom))
                self.assertEqual(
                    profile.resolve(atom).implicit_hydrogens,
                    hydrogen_count,
                    f"{residue}:{atom_name}",
                )
                serial += 1

        oxt = atom_record(
            serial, residue="ALA", atom_name="OXT", element="O"
        )
        self.assertEqual(profile.resolve(oxt).implicit_hydrogens, 0)
        self.assertAlmostEqual(profile.resolve(oxt).excluded_volume, 9.13)

    def test_pontius_published_group_values(self):
        profile = pyseb.AtomicSAXSParameterProfile.protein_defaults(
            pyseb.AtomicSAXSVolumeModel.Pontius1996
        )
        groups = [
            ("PRO", "N", "N", 8.8),
            ("ASP", "OD1", "O", 22.3),
            ("SER", "OG", "O", 23.9),
            ("ARG", "NE", "N", 14.1),
            ("VAL", "CB", "C", 11.8),
            ("ARG", "CB", "C", 20.9),
            ("ALA", "CB", "C", 33.9),
        ]
        for serial, (residue, name, element, expected) in enumerate(
            groups, start=1
        ):
            value = profile.resolve(
                atom_record(
                    serial, residue=residue, atom_name=name, element=element
                )
            ).excluded_volume
            self.assertAlmostEqual(value, expected)

    def test_freesasa_exposure_golden_and_occupancy(self):
        profile = pyseb.AtomicSAXSParameterProfile.protein_defaults()
        cloud = pyseb.AtomicSAXSCloudBuilder.build(
            [atom_record(1, x=-1.0), atom_record(2, x=1.0)], profile
        )
        exposure = [
            atom.solvent_accessibility for atom in cloud.get_scatterers()
        ]
        # FreeSASA 2.1.3, Shrake--Rupley, 1.4 A probe, 960 points.
        self.assertAlmostEqual(exposure[0], 0.6541666666666667, places=12)
        self.assertAlmostEqual(exposure[1], 0.65, places=12)

        translated = pyseb.AtomicSAXSCloudBuilder.build(
            [
                atom_record(1, x=9.0, y=5.0, z=-2.0),
                atom_record(2, x=11.0, y=5.0, z=-2.0),
            ],
            profile,
        )
        translated_exposure = [
            atom.solvent_accessibility
            for atom in translated.get_scatterers()
        ]
        self.assertEqual(exposure, translated_exposure)

        rotated = pyseb.AtomicSAXSCloudBuilder.build(
            [atom_record(1, y=-1.0), atom_record(2, y=1.0)], profile
        )
        rotated_exposure = [
            atom.solvent_accessibility for atom in rotated.get_scatterers()
        ]
        for before, after in zip(exposure, rotated_exposure):
            self.assertAlmostEqual(before, after, delta=2.0 / 960.0)

        occupancy = pyseb.AtomicSAXSCloudBuilder.build(
            [atom_record(1), atom_record(2, occupancy=0.0)], profile
        )
        self.assertAlmostEqual(
            occupancy.get_scatterers()[0].solvent_accessibility, 1.0
        )
        self.assertAlmostEqual(
            occupancy.get_scatterers()[1].solvent_accessibility, 0.0
        )

        fractional = pyseb.AtomicSAXSCloudBuilder.build(
            [
                atom_record(1, x=-1.0, occupancy=0.1),
                atom_record(2, x=1.0, occupancy=0.2),
            ],
            profile,
        )
        self.assertEqual(
            exposure,
            [atom.solvent_accessibility
             for atom in fractional.get_scatterers()],
        )

    def test_freesasa_custom_radius_and_option_validation(self):
        profile = pyseb.AtomicSAXSParameterProfile()
        carbon = carbon_parameters()
        carbon.sasa_radius = 2.25
        profile.set_atom("LIG", "QX", carbon)
        ligand = atom_record(
            1, residue="LIG", atom_name="QX", element="C"
        )
        cloud = pyseb.AtomicSAXSCloudBuilder.build([ligand], profile)
        self.assertAlmostEqual(
            cloud.get_scatterers()[0].solvent_accessibility, 1.0
        )

        inner = carbon_parameters()
        inner.sasa_radius = 1.0
        outer = carbon_parameters()
        outer.sasa_radius = 5.0
        profile.set_atom("LIG", "INNER", inner)
        profile.set_atom("LIG", "OUTER", outer)
        buried = pyseb.AtomicSAXSCloudBuilder.build(
            [
                atom_record(
                    2, residue="LIG", atom_name="INNER", element="C"
                ),
                atom_record(
                    3, residue="LIG", atom_name="OUTER", element="C"
                ),
            ],
            profile,
        )
        self.assertEqual(
            [atom.solvent_accessibility for atom in buried.get_scatterers()],
            [0.0, 1.0],
        )

        build_options = pyseb.AtomicSAXSCloudBuildOptions()
        build_options.sasa_probe_radius = 0.0
        with self.assertRaises(RuntimeError):
            pyseb.AtomicSAXSCloudBuilder.build(
                [ligand], profile, pyseb.AtomicSAXSOptions(), build_options
            )

        build_options.sasa_probe_radius = 1.4
        build_options.sasa_point_count = 0
        with self.assertRaises(RuntimeError):
            pyseb.AtomicSAXSCloudBuilder.build(
                [ligand], profile, pyseb.AtomicSAXSOptions(), build_options
            )

        build_options.sasa_point_count = 960
        build_options.calculate_solvent_accessibility = False
        build_options.sasa_probe_radius = float("nan")
        with self.assertRaises(RuntimeError):
            pyseb.AtomicSAXSCloudBuilder.build(
                [ligand], profile, pyseb.AtomicSAXSOptions(), build_options
            )

        build_options.sasa_probe_radius = 1.4
        malformed = atom_record(
            2, x=float("nan"), residue="LIG", atom_name="QX", element="C"
        )
        with self.assertRaises(RuntimeError):
            pyseb.AtomicSAXSCloudBuilder.build(
                [malformed], profile, pyseb.AtomicSAXSOptions(), build_options
            )

    def test_cromer_mann_carbon_coefficients(self):
        parameters = carbon_parameters()
        coefficients = parameters.vacuum

        self.assertAlmostEqual(coefficients.value_at_zero(), 5.9992, places=12)
        q = 0.37
        s2 = (q / (4.0 * math.pi)) ** 2
        expected = coefficients.c + sum(
            a * math.exp(-b * s2)
            for a, b in zip(coefficients.a, coefficients.b)
        )
        self.assertAlmostEqual(coefficients.evaluate(q), expected, places=12)

    def test_single_atom_vacuum_profile_is_q_dependent(self):
        parameters = carbon_parameters()
        options = pyseb.AtomicSAXSOptions()
        options.include_excluded_solvent = False
        cloud = pyseb.AtomicSAXSCloud(
            [scatterer(0.0, parameters)], options
        )
        world = pyseb.World("atomic")
        world.Add(cloud, "carbon")

        q = 0.4
        f0 = parameters.vacuum.value_at_zero()
        fq = parameters.vacuum.evaluate(q)
        self.assertAlmostEqual(
            world.EvaluateFormFactorUnnormalized("carbon", {}, q),
            fq * fq,
            places=12,
        )
        self.assertAlmostEqual(
            world.EvaluateFormFactor("carbon", {}, q),
            (fq / f0) ** 2,
            places=12,
        )
        self.assertLess(world.EvaluateFormFactor("carbon", {}, q), 1.0)

    def test_excluded_solvent_uses_fraser_gaussian(self):
        parameters = carbon_parameters()
        options = pyseb.AtomicSAXSOptions()
        options.solvent_electron_density = 0.334
        options.excluded_volume_scale = 0.9
        cloud = pyseb.AtomicSAXSCloud(
            [scatterer(0.0, parameters)], options
        )

        q = 0.25
        volume = 0.9 * parameters.excluded_volume
        expected = parameters.vacuum.evaluate(q) - (
            0.334
            * volume
            * math.exp(-q * q * volume ** (2.0 / 3.0) / (4.0 * math.pi))
        )
        self.assertAlmostEqual(cloud.atomic_amplitude(0, q), expected, places=12)

    def test_two_atoms_and_references_use_q_dependent_debye_terms(self):
        parameters = carbon_parameters()
        options = pyseb.AtomicSAXSOptions()
        options.include_excluded_solvent = False
        cloud = pyseb.AtomicSAXSCloud(
            [scatterer(0.0, parameters), scatterer(2.0, parameters)],
            options,
        )
        cloud.add_reference_point("left", 0.0, 0.0, 0.0)
        cloud.add_reference_point("right", 2.0, 0.0, 0.0)
        world = pyseb.World("atomic")
        world.Add(cloud, "pair")

        q = 0.3
        fq = parameters.vacuum.evaluate(q)
        f0 = parameters.vacuum.value_at_zero()
        phase = math.sin(2.0 * q) / (2.0 * q)
        expected_raw = fq * fq * (2.0 + 2.0 * phase)
        self.assertAlmostEqual(
            world.EvaluateFormFactorUnnormalized("pair", {}, q),
            expected_raw,
            places=11,
        )
        self.assertAlmostEqual(
            world.EvaluateFormFactorAmplitude("pair.left", {}, q),
            fq * (1.0 + phase) / (2.0 * f0),
            places=11,
        )
        self.assertAlmostEqual(
            world.EvaluatePhaseFactor("pair.left", "pair.right", {}, q),
            phase,
            places=12,
        )

    def test_analytic_rg_matches_low_q_profile(self):
        parameters = carbon_parameters()
        options = pyseb.AtomicSAXSOptions()
        options.include_excluded_solvent = False
        cloud = pyseb.AtomicSAXSCloud(
            [scatterer(0.0, parameters), scatterer(3.0, parameters)],
            options,
        )
        world = pyseb.World("atomic")
        world.Add(cloud, "pair")

        rg2 = world.EvaluateRadiusOfGyration2("pair", {})
        q = 1.0e-4
        form_factor = world.EvaluateFormFactor("pair", {}, q)
        finite_difference = 3.0 * (1.0 - form_factor) / (q * q)
        self.assertAlmostEqual(rg2, finite_difference, places=5)

    def test_pdb_builder_applies_protein_hydrogens_and_references(self):
        pdb = (
            "ATOM      1  N   ALA A   1      11.104  13.207   9.234  0.50 20.00           N  \n"
            "ATOM      2  CA  ALA A   1      12.560  13.100   9.500  1.00 21.00           C  \n"
            "ATOM      3  CB  ALA A   1      13.100  14.400   9.900  1.00 22.00           C  \n"
            "ATOM      4  C   ALA A   1      13.100  11.800   9.900  1.00 23.00           C  \n"
            "ATOM      5  O   ALA A   1      14.200  11.600  10.100  1.00 24.00           O  \n"
        )
        atoms = pyseb.PDBParser().parse_string(pdb)
        profile = pyseb.AtomicSAXSParameterProfile.protein_defaults()
        self.assertEqual(profile.resolve(atoms[2]).implicit_hydrogens, 3)

        build_options = pyseb.AtomicSAXSCloudBuildOptions()
        build_options.reference_atom_serials = {
            "nterminus": 1,
            "cterminus": 4,
        }
        cloud = pyseb.AtomicSAXSCloudBuilder.build(
            atoms,
            profile,
            pyseb.AtomicSAXSOptions(),
            build_options,
        )
        self.assertEqual(cloud.scatterer_count(), 5)
        references = cloud.get_reference_coordinates()
        self.assertAlmostEqual(references["nterminus"].x, 11.104)
        self.assertAlmostEqual(references["cterminus"].x, 13.100)

        expected_electrons = (
            0.5 * (7.0 + 1.0)
            + 6.0 + 1.0
            + 6.0 + 3.0
            + 6.0
            + 8.0
        )
        options = cloud.get_options()
        options.include_excluded_solvent = False
        cloud.set_options(options)
        world = pyseb.World("atomic")
        world.Add(cloud, "ala")
        # Cromer--Mann coefficient sums differ from integer electron counts by
        # less than 0.02 electrons for this small residue.
        self.assertAlmostEqual(
            math.sqrt(world.EvaluateFormFactorUnnormalized("ala", {}, 0.0)),
            expected_electrons,
            delta=0.02,
        )

        build_options.use_implicit_hydrogens = False
        explicit_hydrogen_cloud = pyseb.AtomicSAXSCloudBuilder.build(
            atoms,
            profile,
            pyseb.AtomicSAXSOptions(),
            build_options,
        )
        self.assertTrue(
            all(
                atom.implicit_hydrogens == 0
                for atom in explicit_hydrogen_cloud.get_scatterers()
            )
        )

    def test_vector_evaluation_is_finite_after_solvent_change(self):
        parameters = carbon_parameters()
        cloud = pyseb.AtomicSAXSCloud(
            [scatterer(0.0, parameters), scatterer(2.0, parameters)]
        )
        world = pyseb.World("atomic")
        world.Add(cloud, "pair")
        before = np.asarray(
            world.EvaluateFormFactor("pair", {}, [0.0, 0.1, 0.2])
        )
        cloud.set_solvent_electron_density(0.30)
        after = np.asarray(
            world.EvaluateFormFactor("pair", {}, [0.0, 0.1, 0.2])
        )
        self.assertTrue(np.all(np.isfinite(after)))
        self.assertAlmostEqual(after[0], 1.0, places=12)
        self.assertFalse(np.allclose(before[1:], after[1:]))

    def test_hydration_is_opt_in_and_requires_supplied_exposure(self):
        parameters = carbon_parameters()
        atom = scatterer(0.0, parameters)
        atom.solvent_accessibility = 0.5
        options = pyseb.AtomicSAXSOptions()
        options.include_excluded_solvent = False
        options.include_hydration = True
        options.hydration_shell_contrast = 0.02
        options.hydration_volume_per_exposed_atom = 30.0
        cloud = pyseb.AtomicSAXSCloud([atom], options)

        expected = parameters.vacuum.value_at_zero() + 0.02 * 30.0 * 0.5
        self.assertAlmostEqual(cloud.atomic_amplitude(0, 0.0), expected, places=12)

        q = 0.3
        hydrogen = pyseb.AtomicSAXSCloud.hydrogen_coefficients()
        oxygen_atom = atom_record(1, element="O", atom_name="O")
        oxygen = (
            pyseb.AtomicSAXSParameterProfile.element_defaults()
            .resolve(oxygen_atom)
            .vacuum
        )
        water_zero = oxygen.value_at_zero() + 2.0 * hydrogen.value_at_zero()
        normalized_water = (
            oxygen.evaluate(q) + 2.0 * hydrogen.evaluate(q)
        ) / water_zero
        expected_q = (
            parameters.vacuum.evaluate(q)
            + 0.02 * 30.0 * 0.5 * normalized_water
        )
        self.assertAlmostEqual(cloud.atomic_amplitude(0, q), expected_q, places=12)

        options.include_hydration = False
        disabled = pyseb.AtomicSAXSCloud([atom], options)
        self.assertAlmostEqual(
            disabled.atomic_amplitude(0, q),
            parameters.vacuum.evaluate(q),
            places=12,
        )

        atom.solvent_accessibility = 1.1
        invalid = pyseb.AtomicSAXSCloud([atom], options)
        with self.assertRaises(RuntimeError):
            invalid.atomic_amplitude(0, 0.0)

    def test_zero_solution_contrast_keeps_only_raw_scattering(self):
        parameters = carbon_parameters()
        options = pyseb.AtomicSAXSOptions()
        options.solvent_electron_density = (
            parameters.vacuum.value_at_zero() / parameters.excluded_volume
        )
        cloud = pyseb.AtomicSAXSCloud(
            [scatterer(0.0, parameters)], options
        )
        world = pyseb.World("match-point")
        world.Add(cloud, "carbon")

        self.assertTrue(
            math.isfinite(
                world.EvaluateFormFactorUnnormalized("carbon", {}, 0.2)
            )
        )
        with self.assertRaises(RuntimeError):
            world.EvaluateFormFactor("carbon", {}, 0.2)
        with self.assertRaises(RuntimeError):
            world.EvaluateRadiusOfGyration2("carbon", {})

    def test_atomic_cloud_composes_with_gaussian_linker(self):
        parameters = carbon_parameters()
        domain1 = pyseb.AtomicSAXSCloud(
            [scatterer(0.0, parameters), scatterer(2.0, parameters)]
        )
        domain1.add_reference_point("end", 2.0, 0.0, 0.0)
        domain2 = pyseb.AtomicSAXSCloud(
            [scatterer(0.0, parameters), scatterer(2.0, parameters)]
        )
        domain2.add_reference_point("end", 0.0, 0.0, 0.0)

        world = pyseb.World("linked-atomic")
        graph = world.Add(domain1, "domain1")
        world.Link(
            "GaussianPolymer", "linker.end1", "domain1.end"
        )
        world.Link(domain2, "domain2.end", "linker.end2")
        world.Add(graph, "assembly")

        values = {
            "Rg_linker": 4.0,
            "beta_linker": 3.0,
        }
        profile = np.asarray(
            world.EvaluateFormFactor(
                "assembly", values, [0.0, 0.05, 0.2]
            )
        )
        self.assertTrue(np.all(np.isfinite(profile)))
        self.assertAlmostEqual(profile[0], 1.0, places=12)


if __name__ == "__main__":
    unittest.main()
