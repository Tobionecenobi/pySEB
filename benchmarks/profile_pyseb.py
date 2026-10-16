"""Measure native pySEB workloads; save curves and timings for before/after checks.

Run from the repository root with its Python environment, for example:
  .venv-sympy/bin/python benchmarks/profile_pyseb.py --output /tmp/pyseb-before
Timings exclude construction and one warm-up call. cProfile identifies Python
and binding costs; use a native profiler to resolve work inside C++ calls.
"""
import argparse
import cProfile
import json
from pathlib import Path
import platform
import statistics
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pyseb


def workloads(q):
    # Analytical shapes exercise both expression evaluation and graph traversal.
    sphere = pyseb.World('profile_sphere')
    sphere.Add('SolidSphere', 'sphere')
    yield 'sphere', lambda: sphere.EvaluateFormFactorUnnormalized(
        'sphere', dict(R_sphere=20., beta_sphere=100.), q)
    chain = pyseb.World('profile_chain')
    graph = chain.Add('GaussianPolymer', 'p0', 'chain')
    for i in range(1, 8):
        chain.Link('GaussianPolymer', f'p{i}.end1', f'p{i-1}.end2', 'chain')
    chain.Add(graph, 'particle')
    yield 'chain_8', lambda: chain.EvaluateFormFactorUnnormalized(
        'particle', dict(Rg_chain=8., beta_chain=10.), q)

    # Fixed seeded geometry makes results reproducible without external assets.
    for count in (128, 512):
        points = []
        for xyz in np.random.default_rng(2026).normal(0., 15., (count, 3)):
            atom = pyseb.AtomicSAXSScatterer()
            atom.center = pyseb.CartesianPoint3D(*xyz)
            atom.vacuum = pyseb.CromerMannCoefficients([2., 1., 1., 1.], [2., 4., 8., 16.], 1.)
            atom.excluded_volume = 12.
            atom.solvent_accessibility = .5
            points.append(atom)
        options = pyseb.AtomicSAXSOptions()
        options.include_hydration = True
        options.hydration_shell_contrast = .02
        cloud = pyseb.AtomicSAXSCloud(points, options)
        world = pyseb.World(f'profile_cloud_{count}')
        world.Add(cloud, 'cloud')
        # Bind loop variables in defaults so every closure keeps its own world.
        yield f'cloud_{count}_intensity', lambda w=world: w.EvaluateFormFactorUnnormalized('cloud', {}, q)
        yield f'cloud_{count}_rg2', lambda w=world: w.EvaluateRadiusOfGyration2('cloud', {})

    # Reuse the maintained hybrid builder, bypassing sasmodels for native timing.
    from Examples.sasmodels_hybrid import make_model
    spec = make_model()
    ctx = spec.world_factory()
    p = spec.bind(ctx, {v.name: v.default for v in spec.parameters})
    yield 'hybrid', lambda: ctx.world.EvaluateFormFactorUnnormalized(spec.structure, p, q)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=7)
    parser.add_argument('--q-count', type=int, default=101)
    parser.add_argument('--compare', type=Path, help='Earlier results.json for parity and speed comparison')
    args = parser.parse_args()
    if args.repeats < 1 or args.q_count < 1:
        parser.error('repeats and q-count must be positive')
    args.output.mkdir(parents=True, exist_ok=True)
    q = np.geomspace(.001, .5, args.q_count).tolist()
    report = dict(python=platform.python_version(), platform=platform.platform(),
                  extension=pyseb._pyseb.__file__, q_count=len(q), repeats=args.repeats, cases={})
    # Save the same call sequence as a Python profile, separately from timings.
    profile = cProfile.Profile()
    for name, fn in workloads(q):
        expected = fn()
        elapsed = []
        for _ in range(args.repeats):
            start = time.perf_counter()
            actual = fn()
            elapsed.append(time.perf_counter() - start)
            np.testing.assert_allclose(actual, expected, rtol=1e-12)
        profile.runcall(fn)
        report['cases'][name] = dict(median_seconds=statistics.median(elapsed),
                                    min_seconds=min(elapsed), values=np.asarray(expected).tolist())
        print(f'{name}: {statistics.median(elapsed)*1000:.3f} ms', flush=True)
    profile.dump_stats(str(args.output / 'python.prof'))
    if args.compare:
        # Fail on changed outputs before interpreting a shorter runtime as a win.
        before = json.loads(args.compare.read_text())
        if before['q_count'] != len(q) or before['cases'].keys() != report['cases'].keys():
            raise ValueError('Comparison requires the same workloads and q grid')
        for name, case in report['cases'].items():
            baseline = before['cases'][name]
            np.testing.assert_allclose(case['values'], baseline['values'], rtol=1e-12, atol=1e-12)
            case['speedup'] = baseline['median_seconds'] / case['median_seconds']
            print(f'{name}: {case["speedup"]:.2f}x baseline speed', flush=True)
    (args.output / 'results.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
