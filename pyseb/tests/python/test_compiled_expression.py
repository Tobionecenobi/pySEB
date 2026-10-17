"""Reusable analytical-equation evaluation must preserve mathematical values."""
import math
import unittest
import warnings
from unittest.mock import patch

import numpy as np
import sympy as sp
import pyseb


class TestCompiledExpression(unittest.TestCase):
    def test_prepared_equations_match_native_and_existing_helper(self):
        q = np.array([[.03, .17], [.4, .03]])
        for kind, dimensions in [('SolidSphere', {'R_unit': 3.}),
                                 ('GaussianPolymer', {'Rg_unit': 2.}),
                                 ('ThinRod', {'L_unit': 4.})]:
            if kind == 'ThinRod':
                try:
                    import scipy.special
                except ImportError:
                    continue  # Special-function support is optional.
            world = pyseb.World()
            world.Add(kind, 'unit')
            equation = world.FormFactor_Unnormalized('unit')
            evaluate = pyseb.compile_expression(equation)
            for beta in (2., -3., 2.):
                parameters = dict(dimensions, beta_unit=beta)
                with self.subTest(kind=kind, beta=beta):
                    actual = evaluate(parameters, q)
                    native = np.array(world.EvaluateFormFactorUnnormalized(
                        'unit', parameters, q.ravel().tolist())).reshape(q.shape)
                    np.testing.assert_allclose(actual, native, rtol=2e-9, atol=1e-10)
                    # The old NumPy-only helper cannot evaluate the rod's Si.
                    # Use high-precision substitution as its independent oracle.
                    if kind == 'ThinRod':
                        expr = pyseb.to_sympy(equation)
                        legacy = np.array([float(expr.evalf(50, subs=dict(parameters, q=x)))
                                           for x in q.ravel()]).reshape(q.shape)
                    else:
                        legacy = pyseb.evaluate_expression(world, equation, parameters, q)
                    np.testing.assert_allclose(actual, legacy, rtol=2e-9, atol=1e-10)

    def test_connected_chain_with_signed_contrasts_and_parameter_changes(self):
        world = pyseb.World()
        graph = world.Add('GaussianPolymer', 'a')
        world.Link('GaussianPolymer', 'b.end1', 'a.end2')
        world.Add(graph, 'particle')
        equation = world.FormFactor_Unnormalized('particle')
        evaluate = pyseb.compile_expression(equation)
        q = [.03, .2, -.2, .03]
        for beta, size in ((2., 3.), (-2., 5.), (2., 3.)):
            p = dict(beta_a=beta, beta_b=-1., Rg_a=size, Rg_b=4.)
            np.testing.assert_allclose(evaluate(p, q),
                world.EvaluateFormFactorUnnormalized('particle', p, q), rtol=2e-9, atol=1e-10)

    def test_shapes_constants_and_empty_arrays(self):
        q, a = sp.symbols('q a')
        for expression in (a, sp.Integer(7), a+q):
            evaluate = pyseb.compile_expression(expression)
            self.assertIsInstance(evaluate({'a': 3.}, 0.), float)
            for grid in (np.zeros((2, 3)), np.empty((2, 0, 3))):
                original = grid.copy()
                actual = evaluate({'a': 3.}, grid)
                self.assertEqual(actual.shape, grid.shape)
                np.testing.assert_array_equal(grid, original)
        self.assertEqual(pyseb.compile_expression(a)({'a': 3.}, 0.), 3.)

    def test_parameters_are_rebound_without_regenerating_code(self):
        a, q = sp.symbols('a q')
        evaluate = pyseb.compile_expression((a+q)**2+sp.sin(a+q))
        self.assertEqual(evaluate.parameter_names, ('a',))
        # Code generation must finish at preparation, not happen inside a fit.
        with patch('sympy.lambdify', side_effect=AssertionError('recompiled')):
            p = {'a': 2., 'unused': 99., 'q': 100.}
            for value in (2., -3., 2.):
                p['a'] = value
                self.assertAlmostEqual(evaluate(p, .1), (value+.1)**2+math.sin(value+.1))
            with self.assertRaisesRegex(ValueError, 'Missing parameter: a'):
                evaluate({}, .1)
            self.assertTrue(math.isfinite(evaluate({'a': 2.}, .1)))

    def test_complex_results_are_checked_before_float_conversion(self):
        q = sp.Symbol('q')
        evaluate = pyseb.compile_expression(sp.I*q)
        # Reject even very small imaginary values, without emitting a cast
        # warning or accepting the real part as a different fitting curve.
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            for grid in (1., [1., 2.], np.array([[0., 1e-100]])):
                with self.subTest(grid=grid):
                    with self.assertRaisesRegex(ValueError, 'complex value'):
                        evaluate({}, grid)
            # Check each evaluation rather than rejecting all complex-typed
            # equations: these inputs yield exactly real values.
            self.assertEqual(evaluate({}, 0.), 0.)
            real = pyseb.compile_expression(sp.exp(sp.I*q))
            self.assertEqual(real({}, 0.), 1.)
            np.testing.assert_array_equal(real({}, np.zeros((2, 3))), np.ones((2, 3)))
            self.assertEqual(real({}, np.empty((2, 0))).shape, (2, 0))

    def test_parameter_names_and_symbol_assumptions_are_unambiguous(self):
        name = 'beta:domain.end'
        unusual = sp.Symbol(name)
        evaluate = pyseb.compile_expression(unusual+sp.Symbol('q'))
        self.assertEqual(evaluate.parameter_names, (name,))
        self.assertAlmostEqual(evaluate({name: -2.}, .3), -1.7)
        # Name-based dictionaries cannot distinguish differently assumed
        # symbols with the same spelling. Reject rather than bind incorrectly.
        a = sp.Symbol('a')
        positive_a = sp.Symbol('a', positive=True)
        with self.assertRaises(ValueError):
            pyseb.compile_expression(a+positive_a)

    def test_special_functions_and_near_zero_q(self):
        try:
            from scipy import special
        except ImportError:
            self.skipTest('Optional SciPy special functions unavailable')
        q, a = sp.symbols('q a')
        expression = sp.erf(a*q)+sp.besselj(0, a*q)+sp.Si(a*q)
        evaluate = pyseb.compile_expression(expression)
        grid = np.array([0., 1e-10, -.2, .2])
        expected = special.erf(2*grid)+special.j0(2*grid)+special.sici(2*grid)[0]
        np.testing.assert_allclose(evaluate({'a': 2.}, grid), expected, rtol=1e-12, atol=1e-12)

    def test_numpy_equations_do_not_require_scipy(self):
        import sys
        q = sp.Symbol('q')
        with patch.dict(sys.modules, {'scipy': None, 'scipy.special': None}):
            evaluate = pyseb.compile_expression(sp.sin(q)+2)
            self.assertAlmostEqual(evaluate({}, .3), math.sin(.3)+2)

    def test_rejects_unresolved_operations_and_invalid_inputs(self):
        q, a = sp.symbols('q a')
        for expression in (sp.Integral(sp.sin(q), q), sp.Function('unknown')(q),
                           sp.Derivative(sp.Function('unknown')(q), q)):
            with self.subTest(expression=expression):
                with self.assertRaises(ValueError):
                    pyseb.compile_expression(expression)
        evaluate = pyseb.compile_expression(a+q)
        for value in (math.nan, math.inf):
            with self.assertRaises(ValueError):
                evaluate({'a': 1.}, [0., value])
            with self.assertRaises(ValueError):
                evaluate({'a': value}, .1)
        # Generic equations have no structural metadata to resolve a removable
        # q=0 singularity. Report it instead of returning non-finite fit values.
        with np.errstate(invalid='ignore', divide='ignore'):
            with self.assertRaises(ValueError):
                pyseb.compile_expression(sp.sin(q)/q)({}, 0.)


if __name__ == '__main__':
    unittest.main()
