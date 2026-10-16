"""Measure analytical graph traversal as connected structures grow.

Construction is measured separately from repeated curve evaluation. All
calculations use pySEB's C++ evaluator, without symbolic conversion, fitting,
or importing sasmodels. Timings are measured after one warm-up call.
"""
import argparse
import json
from pathlib import Path
import statistics
import sys
import time

import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pyseb


def make_world(topology, size):
    world = pyseb.World(f'profile_{topology}_{size}')
    if topology == 'chain':
        graph = world.Add('GaussianPolymer', 'p0', 'poly')
        for i in range(1, size):
            world.Link('GaussianPolymer', f'p{i}.end1', f'p{i-1}.end2', 'poly')
        p = dict(Rg_poly=8., beta_poly=10.)
    else:
        graph = world.Add('SolidSphere', 'core')
        for i in range(size):
            world.Link('GaussianPolymer', f'p{i}.end1', f'core.surface#arm{i}', 'poly')
        p = dict(R_core=20., beta_core=100., Rg_poly=8., beta_poly=10.)
    world.Add(graph, 'particle')
    return world, p


def timed(fn, repeats):
    fn()  # Warm symbolic machinery and any numerical preparation.
    samples = []
    for _ in range(repeats):
        start = time.perf_counter()
        result = fn()
        samples.append(time.perf_counter()-start)
    return result, statistics.median(samples)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sizes', type=int, nargs='+', default=[8, 16, 32])
    parser.add_argument('--q-count', type=int, default=31)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--compare', type=Path, help='Earlier results.json for numerical parity')
    parser.add_argument('--scalar-reference', action='store_true',
                        help='Also time the uncached scalar loop on the same build')
    args = parser.parse_args()
    if min(args.sizes+[args.q_count, args.repeats]) < 1:
        parser.error('sizes, q-count and repeats must be positive')
    args.output.mkdir(parents=True, exist_ok=True)
    q = np.geomspace(.001, .5, args.q_count).tolist()
    report = dict(q_count=len(q), repeats=args.repeats, cases=[])
    baseline = json.loads(args.compare.read_text()) if args.compare else None
    if baseline and baseline['q_count'] != len(q):
        raise ValueError('Comparison requires the same q grid')
    for topology in ('chain', 'star'):
        for size in args.sizes:
            start = time.perf_counter()
            world, p = make_world(topology, size)
            build_seconds = time.perf_counter()-start
            direct, native_seconds = timed(
                lambda: world.EvaluateFormFactorUnnormalized('particle', p, q), args.repeats)
            # Preserve output values for comparison after future optimizations.
            case = dict(topology=topology, size=size, construction_seconds=build_seconds,
                        native_seconds=native_seconds, values=direct)
            if args.scalar_reference:
                # Scalar traversal is unchanged and bypasses the vector cache.
                # This approximates the old vector loop, with added Python call
                # overhead, and provides an independent same-build reference.
                scalar, seconds = timed(lambda: [world.EvaluateFormFactorUnnormalized(
                    'particle', p, value) for value in q], args.repeats)
                np.testing.assert_allclose(direct, scalar, rtol=1e-12, atol=1e-12)
                case['scalar_loop_seconds'] = seconds
                case['scalar_loop_over_vector'] = seconds/native_seconds
            if baseline:
                old = next(c for c in baseline['cases']
                           if c['topology'] == topology and c['size'] == size)
                np.testing.assert_allclose(direct, old['values'], rtol=1e-12, atol=1e-12)
                case['speedup'] = old['native_seconds']/native_seconds
            report['cases'].append(case)
            print(f'{topology} {size}: native {native_seconds*1000:.2f} ms', flush=True)
            (args.output/'results.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
