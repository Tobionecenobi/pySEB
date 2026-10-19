"""Correctness guards for future numerical-evaluation optimizations.

References deliberately use symbolic substitution or an independent double
Debye sum, so scalar/vector implementations cannot merely agree on one bug.
No timing assertions or new evaluation algorithm are introduced here.
"""
import math
import unittest

import numpy as np
import pyseb
import sympy as sp


class TestNumericalOptimization(unittest.TestCase):
    def assert_close(self, actual, expected):
        np.testing.assert_allclose(actual, expected, rtol=2e-11, atol=2e-11)

    def test_numpy_combined_evaluation_preserves_independent_intensity_and_amplitude(self):
        world = pyseb.World()
        graph = world.Add('GaussianPolymer', 'a')
        world.Link('GaussianPolymer', 'b.end1', 'a.end2')
        world.Add(graph, 'particle')
        reference = 'particle:a.end1'
        # A strided, signed grid exercises buffer conversion and sample order.
        q = np.array([.2, 9., 0., 9., -.2, 9., .03, 9.])[::2]
        for beta in (-2., 3., -2.):
            p = dict(Rg_a=3., Rg_b=5., beta_a=beta, beta_b=2.)
            intensity, amplitude = world._evaluate_form_factor_numpy('particle', p, q, reference)
            self.assert_close(intensity, world.EvaluateFormFactorUnnormalized('particle', p, q.tolist()))
            self.assert_close(amplitude, world.EvaluateFormFactorAmplitudeUnnormalized(reference, p, q.tolist()))
            self.assert_close(world._evaluate_form_factor_numpy('particle', p, q), intensity)
            # Contrast matching does not make the finite-q particle invisible.
            self.assertGreater(np.max(np.abs(intensity-amplitude**2)), 1e-4)
            self.assertEqual(intensity.dtype, np.dtype('float64'))
            saved = intensity.copy()
            # Later calls must not overwrite returned arrays or retain parameters.
            world._evaluate_form_factor_numpy('particle', dict(p, beta_a=8.), q, reference)
            self.assert_close(intensity, saved)
        del world
        self.assert_close(intensity, saved)

    def test_numpy_evaluation_empty_grids_and_error_recovery(self):
        world = pyseb.World()
        world.Add('SolidSphere', 'sphere')
        p = dict(R_sphere=3., beta_sphere=-2.)
        for reference in ('', 'sphere.center'):
            empty = world._evaluate_form_factor_numpy('sphere', p, np.array([]), reference)
            for array in (empty if reference else (empty,)):
                self.assertEqual(array.shape, (0,))
            for q in (np.array([np.nan]), np.array([np.inf])):
                with self.assertRaises(RuntimeError):
                    world._evaluate_form_factor_numpy('sphere', p, q, reference)
            with self.assertRaises(ValueError):
                world._evaluate_form_factor_numpy('sphere', p, np.zeros((2, 2)), reference)
            with self.assertRaises(RuntimeError):
                world._evaluate_form_factor_numpy('sphere', dict(beta_sphere=-2.), [.1], reference)
            actual = world._evaluate_form_factor_numpy('sphere', p, np.array([0.], dtype='float32'), reference)
            self.assert_close(actual[0] if reference else actual, [4.])
            if reference:
                self.assert_close(actual[1], [-2.])
            for dtype in ('float32', '>f8'):
                q = np.array([.25, 0., .125], dtype=dtype)
                actual = world._evaluate_form_factor_numpy('sphere', p, q, reference)
                self.assert_close(actual[0] if reference else actual,
                                  world.EvaluateFormFactorUnnormalized('sphere', p, q.tolist()))
        with self.assertRaises(RuntimeError):
            world._evaluate_form_factor_numpy('sphere', p, [.1], 'other.center')

    def test_analytical_values_match_symbolic_substitution(self):
        # These exercise arithmetic, powers, trigonometric functions, exp,
        # and the rod's sine integral through the actual scattering API.
        for kind, dimensions, ref in [
            ('SolidSphere', {'R_unit': 3.}, 'center'),
            ('GaussianPolymer', {'Rg_unit': 2.}, 'end1'),
            ('ThinRod', {'L_unit': 4.}, 'end1'),
        ]:
            world = pyseb.World()
            world.Add(kind, 'unit')
            for beta in (2., -3.):
                parameters = dict(dimensions, beta_unit=beta)
                for expression, method, args in [
                    (world.FormFactor_Unnormalized('unit'),
                     world.EvaluateFormFactorUnnormalized, ('unit',)),
                    (world.FormFactorAmplitude('unit.'+ref),
                     world.EvaluateFormFactorAmplitude, ('unit.'+ref,)),
                ]:
                    # Freeze the substitution-based reference before adding a
                    # direct tree evaluator; evaluate at high precision here.
                    expression = pyseb.to_sympy(expression)
                    for q in (.03, .17, .8):
                        expected = float(expression.evalf(50, subs=dict(parameters, q=q)))
                        with self.subTest(kind=kind, beta=beta, q=q, args=args):
                            self.assert_close(method(*args, parameters, q), expected)

    def test_missing_parameters_still_fail_after_successful_evaluation(self):
        world = pyseb.World()
        world.Add('SolidSphere', 'sphere')
        complete = {'R_sphere': 3., 'beta_sphere': -2.}
        # A warmed evaluator must not retain a previously supplied parameter.
        for missing in complete:
            world.EvaluateFormFactorUnnormalized('sphere', complete, [.1, .3])
            incomplete = {k: v for k, v in complete.items() if k != missing}
            for q in (.1, [.3, .1, .3]):
                with self.subTest(missing=missing, q=q):
                    with self.assertRaises(RuntimeError):
                        world.EvaluateFormFactorUnnormalized('sphere', incomplete, q)
            self.assert_close(world.EvaluateFormFactorUnnormalized('sphere', complete, 0.), 4.)

    def test_chain_parameter_updates_match_a_fresh_world_and_scalar_calls(self):
        def build():
            world = pyseb.World()
            graph = world.Add('GaussianPolymer', 'a')
            world.Link('GaussianPolymer', 'b.end1', 'a.end2')
            world.Link('GaussianPolymer', 'c.end1', 'b.end2')
            world.Add(graph, 'particle')
            return world

        world = build()
        q = [.2, 0., 1e-6, -.2, .03, .2]
        # Return to the original settings to expose stale state in both
        # directions. Different names prevent accidental sharing by type.
        for size, beta in [(2., 3.), (5., -2.), (2., 3.)]:
            p = dict(Rg_a=size, Rg_b=3., Rg_c=4., beta_a=beta, beta_b=-1., beta_c=2.)
            fresh = build()
            for name, args in [
                ('EvaluateFormFactorUnnormalized', ('particle',)),
                ('EvaluateFormFactor', ('particle',)),
                ('EvaluateFormFactorAmplitudeUnnormalized', ('particle:a.end1',)),
                ('EvaluatePhaseFactor', ('particle:a.end1', 'particle:c.end2')),
            ]:
                method = getattr(world, name)
                with self.subTest(size=size, beta=beta, method=name):
                    values = method(*args, p, q)
                    self.assert_close(values, getattr(fresh, name)(*args, p, q))
                    self.assert_close(values, [method(*args, p, x) for x in q])

    def test_numerical_callbacks_follow_parameter_and_external_state_changes(self):
        # Numerical subunits can depend on data outside ParameterList. Reuse
        # must stop at the evaluation boundary even when q/parameters repeat.
        state = {'width': 2.}
        unit = pyseb.NumericalSubunit(pyseb.NormalizationMode.Normalized)
        unit.addReferencePoint('join')
        unit.setTotalBeta(3.)
        unit.setFormFactorFunction(lambda q, p: math.exp(-q*q*state['width']**2))
        unit.setFormFactorAmplitudeFunction(
            'join', lambda q, p: math.exp(-.5*q*q*state['width']**2))
        world = pyseb.World()
        graph = world.Add(unit, 'numeric')
        for i in range(3):
            world.Link('GaussianPolymer', f'p{i}.end1', 'numeric.join')
        world.Add(graph, 'particle')
        q = [.2, .03, .2, 0., -.2]
        for width, contrast in [(2., 1.), (5., -2.), (2., 1.)]:
            state['width'] = width
            p = {f'Rg_p{i}': 2.+i for i in range(3)}
            p.update({f'beta_p{i}': contrast*(i+1) for i in range(3)})
            for name, args in [
                ('EvaluateFormFactorUnnormalized', ('particle',)),
                ('EvaluateFormFactorAmplitudeUnnormalized', ('particle:p0.end2',)),
            ]:
                method = getattr(world, name)
                with self.subTest(width=width, contrast=contrast, method=name):
                    self.assert_close(method(*args, p, q), [method(*args, p, x) for x in q])
            # Errors must not poison a subsequent successful evaluation.
            with self.assertRaises(RuntimeError):
                world.EvaluateFormFactorUnnormalized('particle', p, [.2, math.nan])
            self.assert_close(world.EvaluateFormFactorUnnormalized('particle', p, q),
                              [world.EvaluateFormFactorUnnormalized('particle', p, x) for x in q])

    def test_common_atomic_factors_match_individually_evaluated_atoms(self):
        # Each public atomic_amplitude call starts independently. Comparing
        # with their Debye sum checks lazy reuse inside the whole-cloud path.
        points = []
        coordinates = [(0., 0., 0.), (2., -1., 0.), (-1., 3., 4.)]
        for i, xyz in enumerate(coordinates):
            atom = pyseb.AtomicSAXSScatterer()
            atom.center = pyseb.CartesianPoint3D(*xyz)
            atom.vacuum = pyseb.CromerMannCoefficients(
                [2.+i, 1., .5, .2], [2., 4., 8., 16.], .1+i)
            atom.excluded_volume = 10.+20*i
            atom.occupancy = 1.-.2*i
            atom.implicit_hydrogens = i
            atom.implicit_hydrogen_distance = .9 if i else 0.
            atom.temperature_factor = 5.*i
            atom.solvent_accessibility = .4*i
            points.append(atom)
        cloud = pyseb.AtomicSAXSCloud(points)
        cloud.add_reference_point('join', 0., 0., 0.)
        world = pyseb.World()
        world.Add(cloud, 'cloud')
        q = [.3, 0., 1e-10, -.3, .03, .3]
        models = pyseb.AtomicSAXSHydrationWaterModel
        for water in (models.NeutralWater, models.BondedWater,
                      models.FoXSReciprocal, models.FoXSFast, models.NeutralWater):
            for density, hydrated in ((.2, False), (.334, True), (.4, True)):
                options = pyseb.AtomicSAXSOptions()
                options.solvent_electron_density = density
                options.hydration_water_model = water
                options.include_hydration = hydrated
                options.hydration_shell_contrast = -.02
                options.use_temperature_factors = True
                cloud.set_options(options)
                intensities, amplitudes = [], []
                for x in q:
                    f = [cloud.atomic_amplitude(i, x) for i in range(len(points))]
                    intensities.append(math.fsum(
                        f[i]*f[j]*float(np.sinc(x*math.dist(a, b)/math.pi))
                        for i, a in enumerate(coordinates) for j, b in enumerate(coordinates)))
                    amplitudes.append(math.fsum(
                        value*float(np.sinc(x*math.dist(xyz, (0., 0., 0.))/math.pi))
                        for value, xyz in zip(f, coordinates)))
                with self.subTest(water=water, density=density, hydrated=hydrated):
                    self.assert_close(world.EvaluateFormFactorUnnormalized('cloud', {}, q), intensities)
                    self.assert_close(world.EvaluateFormFactorAmplitudeUnnormalized('cloud.join', {}, q), amplitudes)

    def test_large_cloud_vector_lengths_and_hybrid_interference(self):
        rng = np.random.default_rng(2026)
        atoms = []
        for i, xyz in enumerate(rng.normal(0., 10., (257, 3))):
            atom = pyseb.AtomicSAXSScatterer()
            atom.center = pyseb.CartesianPoint3D(*xyz)
            atom.vacuum = pyseb.CromerMannCoefficients([2., 1., 1., 1.], [2., 4., 8., 16.], 1.)
            atom.excluded_volume = 12. if i % 2 else 30.
            atom.solvent_accessibility = .5
            atom.occupancy = .8 if i % 3 else 1.
            atoms.append(atom)
        cloud = pyseb.AtomicSAXSCloud(atoms)
        cloud.add_reference_point('join', 0., 0., 0.)
        world = pyseb.World()
        graph = world.Add(cloud, 'cloud')
        world.Link('SolidSphere', 'sphere.center', 'cloud.join')
        world.Add(graph, 'hybrid')
        p = dict(R_sphere=3., beta_sphere=-20.)
        # Repeated/signed q and lengths around possible tile boundaries guard
        # future vector kernels while preserving the scalar interference.
        grid = [0., -.2, .2, .03, .2]+np.geomspace(.001, .7, 28).tolist()
        for density, hydration in ((.334, .02), (.2, 0.), (.4, -.01)):
            cloud.set_solvent_electron_density(density)
            cloud.set_hydration(hydration != 0., hydration, 30.)
            for count in (1, 3, 4, 15, 16, 17, 31, 32, 33):
                q = grid[:count]
                for target in ('cloud', 'hybrid'):
                    with self.subTest(density=density, hydration=hydration, count=count, target=target):
                        actual = world.EvaluateFormFactorUnnormalized(target, p, q)
                        scalar = [world.EvaluateFormFactorUnnormalized(target, p, x) for x in q]
                        np.testing.assert_array_equal(actual, scalar)
                        beta2 = world.EvaluateFormFactorUnnormalized(target, p, 0.)
                        self.assert_close(world.EvaluateFormFactor(target, p, q), np.array(scalar)/beta2)
            # The combined fitting path must retain hydration and hybrid cross
            # terms when both intensity and a signed amplitude are requested.
            intensity, amplitude = world._evaluate_form_factor_numpy(
                'hybrid', p, np.array(grid), 'hybrid:cloud.join')
            self.assert_close(intensity, world.EvaluateFormFactorUnnormalized('hybrid', p, grid))
            self.assert_close(amplitude, world.EvaluateFormFactorAmplitudeUnnormalized(
                'hybrid:cloud.join', p, grid))
            with self.assertRaises(RuntimeError):
                world.EvaluateFormFactorUnnormalized('hybrid', p, [.1, .2, math.nan, .3])
        # Inherited callback replacement is public API and must also be
        # respected in vector curves and when embedded in a hybrid.
        cloud.setFormFactorFunction(lambda q, parameters: 7.+q*q)
        q = [.1, .2, .3, .4, .5]
        self.assert_close(world.EvaluateFormFactorUnnormalized('cloud', {}, q), 7.+np.array(q)**2)
        self.assert_close(world.EvaluateFormFactorUnnormalized('hybrid', p, q),
                          [world.EvaluateFormFactorUnnormalized('hybrid', p, x) for x in q])

    @staticmethod
    def make_cloud(density, hydration):
        # Unequal, non-collinear sites with occupancies and signed solution
        # contrasts defeat shortcuts based on equal amplitudes or |amplitude|.
        sites = [(0., 0., 0., 6., 12., 1., .5),
                 (3., 1., 0., 8., 40., .6, 1.),
                 (-1., 2., 4., 9., 8., .8, .2)]
        atoms = []
        for x, y, z, electrons, volume, occupancy, exposure in sites:
            atom = pyseb.AtomicSAXSScatterer()
            atom.center = pyseb.CartesianPoint3D(x, y, z)
            atom.vacuum = pyseb.CromerMannCoefficients([0.]*4, [0.]*4, electrons)
            atom.excluded_volume = volume
            atom.occupancy = occupancy
            atom.solvent_accessibility = exposure
            atoms.append(atom)
        options = pyseb.AtomicSAXSOptions()
        options.solvent_electron_density = density
        options.include_hydration = hydration != 0.
        options.hydration_shell_contrast = hydration
        options.hydration_water_model = pyseb.AtomicSAXSHydrationWaterModel.FoXSFast
        cloud = pyseb.AtomicSAXSCloud(atoms, options)
        cloud.add_reference_point('join', 0., 0., 0.)
        return cloud, sites

    @staticmethod
    def reference(sites, q, density, hydration):
        # Independent full i,j sum rather than pySEB's cached i<j pair loop.
        amplitudes = [occ*(e-density*v*math.exp(-q*q*v**(2/3)/(4*math.pi))
                          + hydration*30*exposure*math.exp(-.115*q*q))
                      for x, y, z, e, v, occ, exposure in sites]
        sinc = lambda x: float(np.sinc(x/math.pi))
        intensity = math.fsum(amplitudes[i]*amplitudes[j]*sinc(q*math.dist(a[:3], b[:3]))
                              for i, a in enumerate(sites) for j, b in enumerate(sites))
        amplitude = math.fsum(f*sinc(q*math.dist(a[:3], (0., 0., 0.)))
                              for f, a in zip(amplitudes, sites))
        return intensity, amplitude

    def test_contrast_matched_cloud_has_no_normalized_intensity(self):
        # At matching, signed atomic amplitudes cancel at q=0 even though
        # finite-q intensity remains. Taking absolute amplitudes breaks this.
        density = 18./42.4
        cloud, sites = self.make_cloud(density, 0.)
        world = pyseb.World()
        world.Add(cloud, 'cloud')
        for q in (0., 1e-10, .2):
            intensity, amplitude = self.reference(sites, q, density, 0.)
            self.assert_close(world.EvaluateFormFactorUnnormalized('cloud', {}, q), intensity)
            self.assert_close(world.EvaluateFormFactorAmplitudeUnnormalized('cloud.join', {}, q), amplitude)
        with self.assertRaises(RuntimeError):
            world.EvaluateFormFactor('cloud', {}, [.1, .2])

    def test_sphere_tiny_nonzero_q_matches_forward_limit(self):
        world = pyseb.World()
        world.Add('SolidSphere', 'sphere')
        p = dict(R_sphere=2., beta_sphere=-4.)
        self.assert_close(world.EvaluateFormFactorUnnormalized('sphere', p, 1e-10), 16.)

    def test_sphere_amplitudes_match_high_precision_across_small_q_cutoff(self):
        # The reference uses the original mathematical expression at 80-digit
        # precision, independently of the numerical small-argument series.
        x = sp.Symbol('x')
        exact = 3 * (sp.sin(x) - x * sp.cos(x)) / x**3
        arguments = np.array([
            0., 1e-20, 1e-10, 1e-4, 1e-3,
            np.nextafter(.1, 0.), .1, np.nextafter(.1, np.inf), .2, 1., 5., 10.,
        ])
        arguments = np.concatenate((arguments, -arguments[1:]))
        world = pyseb.World()
        world.Add('SolidSphere', 'sphere')
        for radius, beta in [(2., -4.), (.25, 3.), (20., -2.)]:
            q = arguments / radius
            qr = q * radius
            amplitude = np.array([
                1. if value == 0 else float(exact.evalf(80, subs={x: value}))
                for value in qr
            ])
            surface = amplitude * np.sinc(qr / np.pi)
            parameters = dict(R_sphere=radius, beta_sphere=beta)
            with self.subTest(radius=radius, beta=beta):
                self.assert_close(world.EvaluateFormFactor('sphere', parameters, q.tolist()),
                                  amplitude**2)
                self.assert_close(world.EvaluateFormFactorUnnormalized('sphere', parameters, q.tolist()),
                                  beta**2 * amplitude**2)
                self.assert_close(world.EvaluateFormFactorAmplitude('sphere.center', parameters, q.tolist()),
                                  amplitude)
                self.assert_close(world.EvaluateFormFactorAmplitudeUnnormalized('sphere.surface', parameters, q.tolist()),
                                  beta * surface)
                self.assert_close([world.EvaluateFormFactorUnnormalized('sphere', parameters, value)
                                   for value in q], beta**2 * amplitude**2)
                intensity, center = world._evaluate_form_factor_numpy(
                    'sphere', parameters, q, 'sphere.center')
                self.assert_close(intensity, beta**2 * amplitude**2)
                self.assert_close(center, beta * amplitude)
                self.assert_close([pyseb.DebyeSphereCloud.sphereAmplitude(value) for value in qr],
                                  amplitude)

    def test_cloud_tiny_nonzero_q_matches_independent_debye_sum(self):
        cloud, sites = self.make_cloud(.334, .02)
        world = pyseb.World()
        world.Add(cloud, 'cloud')
        for q in (0., 1e-10, -1e-10):
            intensity, amplitude = self.reference(sites, q, .334, .02)
            self.assert_close(world.EvaluateFormFactorUnnormalized('cloud', {}, q), intensity)
            self.assert_close(world.EvaluateFormFactorAmplitudeUnnormalized('cloud.join', {}, q), amplitude)

    def test_cloud_and_hybrid_preserve_debye_sum_and_signed_cross_terms(self):
        cloud, sites = self.make_cloud(.334, .02)
        world = pyseb.World()
        graph = world.Add(cloud, 'cloud')
        # A center attachment gives an independent closed-form cross term:
        # I_hybrid = I_cloud + (beta*f_sphere)^2 + 2*A_cloud*beta*f_sphere.
        world.Link('SolidSphere', 'sphere.center', 'cloud.join')
        world.Add(graph, 'hybrid')
        q = [.2, 0., .003, -.2, .7, .2]
        for density, hydration, beta, radius in [
            (.334, .02, -4., 2.), (.4, 0., 3., 4.),
            (.2, -.01, -2., 3.), (.334, .02, -4., 2.),
        ]:
            cloud.set_solvent_electron_density(density)
            cloud.set_hydration(hydration != 0., hydration, 30.)
            p = dict(beta_sphere=beta, R_sphere=radius)
            expected_cloud, expected_amplitude, expected_hybrid = [], [], []
            for x in q:
                intensity, amplitude = self.reference(sites, x, density, hydration)
                z = x*radius
                # Avoid importing the old double-precision cancellation into
                # the independent hybrid reference near contrast matching.
                if z == 0:
                    sphere = 1.0
                else:
                    high_precision_z = sp.Float(z, 80)
                    sphere = float(3 * (sp.sin(high_precision_z)
                        - high_precision_z * sp.cos(high_precision_z)) / high_precision_z**3)
                expected_cloud.append(intensity)
                expected_amplitude.append(amplitude)
                expected_hybrid.append(intensity+(beta*sphere)**2+2*amplitude*beta*sphere)
            with self.subTest(density=density, hydration=hydration, beta=beta):
                self.assert_close(world.EvaluateFormFactorUnnormalized('cloud', {}, q), expected_cloud)
                self.assert_close(world.EvaluateFormFactorAmplitudeUnnormalized('cloud.join', {}, q), expected_amplitude)
                self.assert_close(world.EvaluateFormFactorUnnormalized('hybrid', p, q), expected_hybrid)
                self.assert_close([world.EvaluateFormFactorUnnormalized('hybrid', p, x) for x in q], expected_hybrid)
                # A fresh cloud provides a separate guard against stale mutable
                # solvent/hydration state, including returning to the first state.
                fresh, _ = self.make_cloud(density, hydration)
                other = pyseb.World()
                other.Add(fresh, 'cloud')
                self.assert_close(other.EvaluateFormFactorUnnormalized('cloud', {}, q), expected_cloud)


if __name__ == '__main__':
    unittest.main()
