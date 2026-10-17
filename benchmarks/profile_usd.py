"""Compare native USD exporters with interleaved timings and exact scene checks.

Build each executable with build_usd_profile.py, selecting an original source
snapshot for --before. All generated scenes and JSON are written to --output.
Only Python's standard library is required; model construction and warm-up are
handled by the native driver and excluded from its reported export duration.
"""
import argparse
import hashlib
import json
from pathlib import Path
import statistics
import subprocess


CASES = [
    ('analytical', 64, 'random', 'points'),
    ('analytical', 16, 'readable', 'points'),
    ('debye', 512, 'random', 'envelope'),
    ('debye', 512, 'random', 'points'),
    ('atomic', 4096, 'random', 'points'),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', type=Path, required=True)
    parser.add_argument('--after', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=5)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('repeats must be positive')
    args.output.mkdir(parents=True, exist_ok=True)
    report = {'repeats': args.repeats, 'cases': []}
    for kind, count, layout, envelope in CASES:
        name = f'{kind}-{count}-{layout}-{envelope}'
        samples = {'before': [], 'after': []}
        digest = None
        for repeat in range(args.repeats):
            # Alternate order to reduce systematic effects from CPU/filesystem
            # warm-up. Fresh processes also make peak RSS comparable on Linux.
            order = ('before', 'after') if repeat % 2 == 0 else ('after', 'before')
            for version in order:
                scene = args.output/f'{name}-{version}.usda'
                result = subprocess.run([
                    str(getattr(args, version).resolve()), kind, str(count),
                    layout, envelope, '1', str(scene.resolve()),
                ], capture_output=True, text=True, check=True)
                row = json.loads(result.stdout)
                samples[version].append(row)
                hasher = hashlib.sha256()
                with scene.open('rb') as stream:
                    for chunk in iter(lambda: stream.read(1024*1024), b''):
                        hasher.update(chunk)
                actual = hasher.hexdigest()
                if digest is not None and actual != digest:
                    raise AssertionError(f'Exported scene changed: {name} {version}')
                digest = actual
        # Preserve individual stage measurements instead of reporting only the
        # total. Independent stage medians need not add to the median total.
        medians = {}
        for version, rows in samples.items():
            medians[version] = {
                key: statistics.median(row[key] for row in rows)
                for key in rows[0] if key.endswith('_seconds') or key == 'peak_rss_kib'
            }
        speedup = medians['before']['total_seconds']/medians['after']['total_seconds']
        report['cases'].append(dict(name=name, sha256=digest, samples=samples,
                                    medians=medians, speedup=speedup))
        (args.output/'results.json').write_text(json.dumps(report, indent=2)+'\n')
        print(f'{name}: {speedup:.2f}x; identical bytes', flush=True)


if __name__ == '__main__':
    main()
