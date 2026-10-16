"""Isolate atom-group and pair differences using the Pepsi comparison harness."""
import argparse
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

from compare_atomic_saxs_pepsi import ROOT
import pyseb


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pepsi', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    groups = [('C', 'ALA', 'C'), ('CH', 'ALA', 'CA'),
              ('CH2', 'GLY', 'CA'), ('CH3', 'ALA', 'CB'),
              ('N', 'PRO', 'N'), ('NH', 'ALA', 'N'),
              ('O', 'ALA', 'O'), ('OH', 'SER', 'OG'),
              ('S', 'MET', 'SD'), ('SH', 'CYS', 'SG')]
    report = {}
    for label, residue, name in groups:
        directory = args.output / label
        directory.mkdir(exist_ok=True)
        element = name[0]
        pdb = directory / 'single.pdb'
        pdb.write_text(f'ATOM      1 {name:>4} {residue} A   1    '
                       f'{0.:8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00          {element:>2}\nEND\n')
        command = [sys.executable, str(ROOT / 'Examples/compare_atomic_saxs_pepsi.py'),
                   '--pepsi', str(args.pepsi.resolve()), '--pdb', str(pdb.resolve()),
                   '--output', str(directory.resolve())]
        subprocess.run(command, check=True, capture_output=True, text=True)
        data = np.loadtxt(directory / 'hydration_0pct.dat')
        log = (directory / 'hydration_0pct.stdout.txt').read_text()
        scaling = float(re.search(r'Intensity scaling\.+ : (\S+)', log)[1])
        atom = pyseb.AtomRecord()
        atom.residue_name, atom.atom_name, atom.element = residue, name, element
        p = pyseb.AtomicSAXSParameterProfile.protein_defaults().resolve(atom)
        ours_f0 = p.vacuum.value_at_zero() + p.implicit_hydrogens * pyseb.AtomicSAXSCloud.hydrogen_coefficients().value_at_zero()
        pepsi_f0 = np.sqrt(data[0, 2] / scaling)
        pepsi_volume = np.sqrt(data[0, 3] / scaling) / 0.334
        report[label] = {'pyseb_f0': ours_f0, 'pepsi_f0': float(pepsi_f0),
                         'pyseb_volume': p.excluded_volume,
                         'pepsi_effective_volume': float(pepsi_volume),
                         'comparison': json.loads((directory / 'comparison.json').read_text())}
    (args.output / 'groups.json').write_text(json.dumps(report, indent=2) + '\n')
    for name, values in report.items():
        print(name, 'f0', round(values['pyseb_f0'], 4), round(values['pepsi_f0'], 4),
              'volume', round(values['pyseb_volume'], 4), round(values['pepsi_effective_volume'], 4),
              'vacuum shape max', values['comparison']['cases']['hydration_0pct']['vacuum_normalized_relative_max'])

    pairs = {}
    for second, distance in [('C', 2.0), ('C', 10.0), ('O', 2.0), ('O', 10.0)]:
        label = f'C{second}_{distance:g}A'
        directory = args.output / label
        directory.mkdir(exist_ok=True)
        pdb = directory / 'pair.pdb'
        lines = []
        for serial, (element, x) in enumerate([('C', -distance / 2),
                                                (second, distance / 2)], 1):
            lines.append(f'ATOM  {serial:5d} {element:>4} ALA A{serial:4d}    '
                         f'{x:8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00          {element:>2}')
        pdb.write_text('\n'.join(lines) + '\nEND\n')
        subprocess.run([sys.executable, str(ROOT / 'Examples/compare_atomic_saxs_pepsi.py'),
                        '--pepsi', str(args.pepsi.resolve()), '--pdb', str(pdb.resolve()),
                        '--output', str(directory.resolve())],
                       check=True, capture_output=True, text=True)
        data = np.loadtxt(directory / 'hydration_0pct.dat')
        log = (directory / 'hydration_0pct.stdout.txt').read_text()
        scaling = float(re.search(r'Intensity scaling\.+ : (\S+)', log)[1])
        def single_amplitude(group, column):
            single = np.loadtxt(args.output / group / 'hydration_0pct.dat')
            single_log = (args.output / group / 'hydration_0pct.stdout.txt').read_text()
            scale = float(re.search(r'Intensity scaling\.+ : (\S+)', single_log)[1])
            return np.sqrt(single[:, column] / scale)
        checks = {}
        for name, column in [('vacuum', 2), ('excluded', 3), ('solution', 1)]:
            a, b = single_amplitude('C', column), single_amplitude(second, column)
            predicted = a*a + b*b + 2*a*b*np.sinc(data[:, 0] * distance / np.pi)
            actual = data[:, column] / scaling
            # Normalize to I0: Pepsi log scale is printed with limited precision.
            error = (actual / actual[0]) / (predicted / predicted[0]) - 1
            checks[name + '_pair_identity_max_relative'] = float(np.max(np.abs(error)))
        pairs[label] = checks
        print(label, checks)
    (args.output / 'pairs.json').write_text(json.dumps(pairs, indent=2) + '\n')


if __name__ == '__main__':
    main()
