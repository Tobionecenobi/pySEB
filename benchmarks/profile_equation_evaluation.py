"""Separate analytical generation, evaluator preparation and repeated curves.

The prepared NumPy/SymPy path evaluates exported equations; the native vector
path evaluates the connected World directly. No timing thresholds are tests.
"""
import argparse
import json
from pathlib import Path
import statistics
import sys
import time

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pyseb
from profile_analytical import make_world


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--sizes', type=int, nargs='+', default=[8, 16])
    parser.add_argument('--repeats', type=int, default=7)
    args = parser.parse_args()
    if min(args.sizes+[args.repeats]) < 1:
        parser.error('sizes and repeats must be positive')
    q = np.geomspace(.03, .5, 101)
    report = dict(backend=pyseb.get_backend(), q=q.tolist(), cases=[])
    # Exclude cold optional-library imports from per-equation preparation.
    pyseb.compile_expression(sp.Integer(1))
    for topology in ('chain', 'star'):
        for count in args.sizes:
            world, parameters = make_world(topology, count)
            start = time.perf_counter()
            equation = world.FormFactor_Unnormalized('particle')
            generation = time.perf_counter()-start
            start = time.perf_counter()
            evaluate = pyseb.compile_expression(equation)
            preparation = time.perf_counter()-start
            samples = {name: [] for name in ('existing_helper', 'prepared', 'native_vector')}
            values = []
            # Change bindings like an optimizer; compilation happens only once.
            for repeat in range(args.repeats+1):
                p = dict(parameters, Rg_poly=6.+repeat/4,
                         beta_poly=(-1 if repeat % 2 else 1)*10.)
                functions = dict(
                    existing_helper=lambda: pyseb.evaluate_expression(world, equation, p, q),
                    prepared=lambda: evaluate(p, q),
                    native_vector=lambda: world.EvaluateFormFactorUnnormalized('particle', p, q.tolist()),
                )
                curves = {}
                for name, function in functions.items():
                    start = time.perf_counter()
                    curves[name] = np.asarray(function())
                    elapsed = time.perf_counter()-start
                    if repeat:
                        samples[name].append(elapsed)
                # Validate each changed parameter set, not just one checksum.
                for curve in curves.values():
                    np.testing.assert_allclose(curve, curves['native_vector'], rtol=2e-8, atol=1e-8)
                values.append(curves['prepared'].tolist())
            case = dict(topology=topology, count=count, generation_seconds=generation,
                        preparation_seconds=preparation,
                        median_seconds={name: statistics.median(times) for name, times in samples.items()},
                        samples=samples, prepared_values=values)
            report['cases'].append(case)
            print(json.dumps({key: value for key, value in case.items()
                              if key not in ('samples', 'prepared_values')}), flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
