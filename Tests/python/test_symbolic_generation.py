"""Mathematical regression tests for analytical equation generation.

Compare generated equations with the independent numerical traversal rather
than fixing one printed layout. The large-chain case also has a closed-form
reference and guards against silently returning an unparsed expression symbol.
"""
import contextlib
import io
import unittest

import numpy as np
import sympy as sp
import pyseb


def chain(count, shared=False):
    world = pyseb.World('symbolic_chain')
    graph = world.Add('GaussianPolymer', 'p0', 'poly' if shared else '')
    for i in range(1, count):
        world.Link('GaussianPolymer', f'p{i}.end1', f'p{i-1}.end2', 'poly' if shared else '')
    world.Add(graph, 'particle')
    return world


class TestSymbolicGeneration(unittest.TestCase):
    def setUp(self):
        self.backend = pyseb.get_backend()
        pyseb.set_backend('portable')
        self.q = np.array([.02, .08, .2, .5])

    def tearDown(self):
        pyseb.set_backend(self.backend)

    def assert_curve(self, expression, parameters, expected):
        expr = pyseb.to_sympy(expression)
        # A parse failure must not masquerade as a valid symbolic equation.
        self.assertTrue(expr.free_symbols <= {sp.Symbol(k) for k in [*parameters, 'q']},
                        'Generated equation contains unexpected or unparsed symbols')
        bound = expr.subs(parameters)
        actual = sp.lambdify(sp.Symbol('q'), bound, modules=['numpy', 'scipy'])(self.q)
        np.testing.assert_allclose(actual, expected, rtol=2e-9, atol=1e-9)

    def test_chain_intensity_amplitudes_and_phase_with_signed_contrasts(self):
        world = chain(4)
        intensity = world.FormFactor_Unnormalized('particle')
        normalized = world.FormFactor('particle')
        ref = 'particle:p0.end1'
        amplitude = world.FormFactorAmplitude(ref)
        last = 'particle:p3.end2'
        phase = world.PhaseFactor(ref, last)
        for sign in (1., -1.):
            p = {f'Rg_p{i}': 3.+i for i in range(4)}
            p.update({f'beta_p{i}': sign*(i+1) if i % 2 else 2. for i in range(4)})
            for equation, method, refs in [
                (intensity, world.EvaluateFormFactorUnnormalized, ('particle',)),
                (normalized, world.EvaluateFormFactor, ('particle',)),
                (amplitude, world.EvaluateFormFactorAmplitude, (ref,)),
                (phase, world.EvaluatePhaseFactor, (ref, last)),
            ]:
                self.assert_curve(equation, p, method(*refs, p, self.q.tolist()))

    def test_distributed_references_and_repeated_parameter_tags(self):
        world = pyseb.World('symbolic_star')
        graph = world.Add('SolidSphere', 'core')
        for i in range(4):
            # Two arms share a sampled anchor; the other anchors are independent.
            # Memoization must preserve that distinction in the cross terms.
            anchor = 'shared' if i < 2 else f'arm{i}'
            world.Link('GaussianPolymer', f'p{i}.end1', f'core.surface#{anchor}', 'poly')
        world.Add(graph, 'particle')
        equation = world.FormFactor_Unnormalized('particle')
        p = dict(R_core=8., beta_core=-10., Rg_poly=3., beta_poly=2.)
        self.assert_curve(equation, p, world.EvaluateFormFactorUnnormalized('particle', p, self.q.tolist()))
        self.assertEqual(set(world.getParams()), set(p))

    def test_nested_reused_graphs_and_depth(self):
        world = pyseb.World('symbolic_nested')
        block = world.Add('GaussianPolymer', 'a')
        world.Link('GaussianPolymer', 'b.end1', 'a.end2')
        graph = world.Add(block, 'left')
        world.Link(block, 'right:a.end1', 'left:b.end2')
        world.Add(graph, 'particle')
        p = dict(Rg_a=3., Rg_b=5., beta_a=2., beta_b=-1.)
        equation = world.FormFactor_Unnormalized('particle')
        self.assert_curve(equation, p, world.EvaluateFormFactorUnnormalized('particle', p, self.q.tolist()))
        ref, last = 'particle:left:a.end1', 'particle:right:b.end2'
        self.assert_curve(world.FormFactorAmplitude(ref), p,
                          world.EvaluateFormFactorAmplitude(ref, p, self.q.tolist()))
        self.assert_curve(world.PhaseFactor(ref, last), p,
                          world.EvaluatePhaseFactor(ref, last, p, self.q.tolist()))
        # Shallow output must remain generic and must not inherit a deeper call's cache.
        shallow = pyseb.to_sympy(world.FormFactor_Unnormalized('particle', depth=0))
        self.assertTrue(sp.Symbol('F_particle') in shallow.free_symbols)
        self.assertNotIn(sp.Symbol('q'), shallow.free_symbols)

    def test_contrast_matching_and_zero_q_limit(self):
        world = chain(2)
        p = dict(Rg_p0=3., Rg_p1=5., beta_p0=2., beta_p1=-2.)
        equation = world.FormFactor_Unnormalized('particle')
        self.assert_curve(equation, p, world.EvaluateFormFactorUnnormalized('particle', p, self.q.tolist()))
        self.assertEqual(sp.limit(pyseb.to_sympy(equation).subs(p), sp.Symbol('q'), 0), 0)
        self.assertEqual(world.EvaluateFormFactorUnnormalized('particle', p, 0.), 0.)

    def test_generation_after_world_edit(self):
        world = chain(8, shared=True)
        old = world.FormFactor_Unnormalized('particle')
        world.Link('GaussianPolymer', 'p8.end1', 'p7.end2', 'poly')
        new = world.FormFactor_Unnormalized('particle')
        p = dict(Rg_poly=3., beta_poly=2.)
        self.assert_curve(new, p, world.EvaluateFormFactorUnnormalized('particle', p, self.q.tolist()))
        before = sp.lambdify(sp.Symbol('q'), pyseb.to_sympy(old).subs(p), 'numpy')(self.q)
        after = sp.lambdify(sp.Symbol('q'), pyseb.to_sympy(new).subs(p), 'numpy')(self.q)
        self.assertGreater(np.max(np.abs(before-after)), 1.)

    def test_expression_modes_and_parameter_tracking(self):
        world = chain(2)
        p = dict(Rg_p0=3., Rg_p1=5., beta_p0=2., beta_p1=-1.)
        beta_mode = pyseb.to_sympy(world.FormFactor_Unnormalized('particle', varForm=4))
        self.assertEqual(float(beta_mode.subs(p)), 1.)
        self.assertEqual(set(world.getParams()), {'beta_p0', 'beta_p1'})
        count_mode = pyseb.to_sympy(world.FormFactor_Unnormalized('particle', varForm=7))
        self.assertEqual(count_mode, 4)
        self.assertEqual(world.getParams(), {})
        generic = pyseb.to_sympy(world.FormFactor_Unnormalized('particle', varForm=1))
        self.assertNotIn(sp.Symbol('q'), generic.free_symbols)
        guinier = pyseb.to_sympy(world.FormFactor_Unnormalized('particle', varForm=5)).subs(p)
        q = sp.Symbol('q')
        rg2 = -3.*float(sp.diff(guinier, q, 2).subs(q, 0))/2.
        self.assertAlmostEqual(rg2, world.EvaluateRadiusOfGyration2('particle', p), places=10)
        self.assert_curve(world.FormFactor_Unnormalized('particle'), p,
                          world.EvaluateFormFactorUnnormalized('particle', p, self.q.tolist()))
        self.assertEqual(set(world.getParams()), set(p))

    def test_deep_native_expression_conversion(self):
        # This case cannot be repaired just by balancing World's generated sums:
        # it exercises direct conversion of a user's existing native tree.
        x = pyseb.symbol('deep_x')
        expression = pyseb.constant(0.)
        for _ in range(300):
            expression = expression+x
        with contextlib.redirect_stdout(io.StringIO()):
            converted = pyseb.to_sympy(expression)
        self.assertEqual(converted, 300*sp.Symbol('deep_x'))

    def test_shallow_native_copies_preserve_operation_and_operand_order(self):
        x, y = pyseb.symbol('copy_x'), pyseb.symbol('copy_y')
        # Arithmetic copies a root while retaining its children. Distinct
        # operations on those same children must never share a converted value.
        a, b = x.sin()+y, x.cos()+2.
        expressions = [a+b, a-b, b-a, a*b, a/b, b/a, a**b,
                       a.sin()+a.cos()+a.exp(), -(a*b), (a*b)+(a*b)]
        for expression in expressions:
            converted = pyseb.to_sympy(expression)
            for xv, yv in ((.4, .7), (.8, 1.2)):
                expected = expression.subs('copy_x', xv).subs('copy_y', yv).eval()
                actual = float(converted.evalf(subs={'copy_x': xv, 'copy_y': yv}))
                self.assertAlmostEqual(actual, expected, places=11)

    def test_independently_parameterized_nested_blocks(self):
        world = pyseb.World('independent_blocks')
        p = {}
        for i in range(6):
            inner = world.Add('GaussianPolymer', f'a{i}')
            world.Link('GaussianPolymer', f'b{i}.end1', f'a{i}.end2')
            if i == 0:
                graph = world.Add(inner, 'block0')
            else:
                world.Link(inner, f'block{i}:a{i}.end1', f'block{i-1}:b{i-1}.end2')
            p.update({f'Rg_a{i}': 3.+i, f'Rg_b{i}': 4.+i,
                      f'beta_a{i}': 2., f'beta_b{i}': -1.})
        world.Add(graph, 'particle')
        first, last = 'particle:block0:a0.end1', 'particle:block5:b5.end2'
        intensity = world.FormFactor_Unnormalized('particle')
        amplitude = world.FormFactorAmplitude(first)
        phase = world.PhaseFactor(first, last)
        self.assertEqual(pyseb.to_sympy(intensity).free_symbols,
                         {sp.Symbol(name) for name in [*p, 'q']})
        for parameters in (p, dict(p, beta_a0=-2., Rg_b1=9.)):
            self.assert_curve(intensity, parameters,
                world.EvaluateFormFactorUnnormalized('particle', parameters, self.q.tolist()))
            self.assert_curve(amplitude, parameters,
                world.EvaluateFormFactorAmplitude(first, parameters, self.q.tolist()))
            self.assert_curve(phase, parameters,
                world.EvaluatePhaseFactor(first, last, parameters, self.q.tolist()))

    def test_large_chain_sympy_conversion_and_closed_form(self):
        count = 32
        world = chain(count, shared=True)
        # Suppress the current parser's enormous warning, while still failing
        # on the actual bug: a returned Symbol instead of the complete equation.
        with contextlib.redirect_stdout(io.StringIO()):
            equation = world.FormFactor_Unnormalized('particle')
        p = dict(Rg_poly=3., beta_poly=2.)
        native = world.EvaluateFormFactorUnnormalized('particle', p, self.q.tolist())
        self.assert_curve(equation, p, native)
        self.assert_curve(world.FormFactor('particle'), p,
                          world.EvaluateFormFactor('particle', p, self.q.tolist()))
        ref, last = 'particle:p0.end1', f'particle:p{count-1}.end2'
        self.assert_curve(world.FormFactorAmplitude(ref), p,
                          world.EvaluateFormFactorAmplitude(ref, p, self.q.tolist()))
        self.assert_curve(world.PhaseFactor(ref, last), p,
                          world.EvaluatePhaseFactor(ref, last, p, self.q.tolist()))
        # Identical Gaussian segments concatenate into one Gaussian chain.
        x = self.q**2 * count * p['Rg_poly']**2
        closed = (count*p['beta_poly'])**2 * 2.*(np.expm1(-x)+x)/x**2
        np.testing.assert_allclose(native, closed, rtol=2e-9)


if __name__ == '__main__':
    unittest.main()
