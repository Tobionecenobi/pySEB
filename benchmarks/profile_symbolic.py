"""Measure native generation plus Python conversion, including conversion failures.

Use profile_symbolic.cpp to isolate native generation and rendering. This script
measures the public Python API, then checks each successfully converted equation
against numerical traversal. No timings are enforced as test thresholds.
"""
import argparse
import contextlib
import io
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
    parser.add_argument('--sizes', type=int, nargs='+', default=[8, 16, 32])
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    if min(args.sizes+[args.repeats]) < 1:
        parser.error('sizes and repeats must be positive')
    args.output.mkdir(parents=True, exist_ok=True)
    report = dict(backend=pyseb.get_backend(), repeats=args.repeats, cases=[])
    q = np.array([.02, .08, .2, .5])
    for topology in ('chain', 'star'):
        for size in args.sizes:
            world, p = make_world(topology, size)
            samples = []
            for repeat in range(args.repeats+1):
                # The old string parser prints the entire failed expression.
                # Capture it without hiding the failure in the recorded report.
                with contextlib.redirect_stdout(io.StringIO()):
                    start = time.perf_counter()
                    expression = world.FormFactor_Unnormalized('particle')
                    elapsed = time.perf_counter()-start
                if repeat:
                    samples.append(elapsed)
            expr = pyseb.to_sympy(expression)
            converted = expr.free_symbols <= {sp.Symbol(k) for k in [*p, 'q']}
            case = dict(topology=topology, size=size, seconds=statistics.median(samples),
                        converted=bool(converted))
            if converted:
                actual = sp.lambdify(sp.Symbol('q'), expr.subs(p), ['numpy', 'scipy'])(q)
                expected = world.EvaluateFormFactorUnnormalized('particle', p, q.tolist())
                np.testing.assert_allclose(actual, expected, rtol=2e-8, atol=1e-8)
                case['values'] = np.asarray(actual).tolist()
            report['cases'].append(case)
            print(f'{topology} {size}: {case["seconds"]*1000:.2f} ms; converted={converted}', flush=True)
            (args.output/'results.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
