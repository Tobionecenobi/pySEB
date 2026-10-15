"""Compare fixed-coordinate SAXS with an externally installed Pepsi-SAXS 3.0.

Run with .venv-sympy/bin/python Examples/compare_atomic_saxs_pepsi.py
  --pepsi /path/to/Pepsi-SAXS --output /tmp/atomic-pepsi
No experimental fitting, background adjustment, or downloaded binary is bundled.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import pyseb


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pepsi', type=Path, required=True)
    parser.add_argument('--pdb', type=Path,
                        default=ROOT / 'docs/tutorials/data/1TIU-heavy-ensemble.pdb')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--spatial-shell', action='store_true')
    parser.add_argument('--shell-spacing', type=float, default=3.0)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    # Both calculators receive exactly the same first model and ATOM records.
    records = []
    for line in args.pdb.read_text().splitlines():
        if line.startswith('ENDMDL'):
            break
        if line.startswith('ATOM  '):
            records.append(line)
    pdb = args.output / 'input.pdb'
    pdb.write_text('\n'.join(records) + '\nEND\n')
    result = {'input_sha256': hashlib.sha256(pdb.read_bytes()).hexdigest(),
              'pepsi_sha256': hashlib.sha256(args.pepsi.read_bytes()).hexdigest(),
              'atom_records': len(records), 'cases': {}}
    curves = {}
    for contrast in (0, 5):
        label = f'hydration_{contrast}pct'
        target = args.output / f'{label}.dat'
        command = [str(args.pepsi.resolve()), str(pdb.resolve()), '-o',
                   str(target.resolve()), '--dro', str(contrast), '--noSmearing',
                   '--esmodel', '0', '-ns', '101', '-ms', '0.5', '-n', '30']
        completed = subprocess.run(command, check=True, capture_output=True, text=True)
        (args.output / f'{label}.stdout.txt').write_text(completed.stdout)
        data = np.loadtxt(target)
        q, pepsi = data[:, 0], data[:, 1]
        options = pyseb.AtomicSAXSOptions()
        options.include_hydration = contrast != 0
        options.hydration_shell_contrast = 0.334 * contrast / 100
        build = pyseb.AtomicSAXSCloudBuildOptions()
        if args.spatial_shell:
            options.hydration_model = pyseb.AtomicSAXSHydrationModel.SpatialShell
            build.hydration_shell_spacing = args.shell_spacing
        cloud = pyseb.AtomicSAXSStructureLoader.load_pdb(
            str(pdb), pyseb.AtomicSAXSParameterProfile.protein_defaults(),
            pyseb.StructureParseOptions(), options,
            build)
        world = pyseb.World(label)
        world.Add(cloud, 'protein')
        ours = np.asarray(world.EvaluateFormFactorUnnormalized('protein', {}, q))
        vacuum_options = pyseb.AtomicSAXSOptions()
        vacuum_options.include_excluded_solvent = False
        cloud.set_options(vacuum_options)
        vacuum = np.asarray(world.EvaluateFormFactorUnnormalized('protein', {}, q))
        # Iat is Pepsi's vacuum component. Compare shapes to remove unit scaling.
        normalized = ours / ours[0]
        reference = pepsi / pepsi[0]
        rel = normalized / reference - 1
        vacuum_rel = (vacuum / vacuum[0]) / (data[:, 2] / data[0, 2]) - 1
        result['cases'][label] = {
            'shell_contrast_e_per_A3': options.hydration_shell_contrast,
            'normalized_curve_relative_rms': float(np.sqrt(np.mean(rel**2))),
            'normalized_curve_relative_max': float(np.max(np.abs(rel))),
            'vacuum_normalized_relative_rms': float(np.sqrt(np.mean(vacuum_rel**2))),
            'vacuum_normalized_relative_max': float(np.max(np.abs(vacuum_rel))),
            'pyseb_I0_electrons_squared': float(ours[0]),
            'pepsi_I0_cm_inverse': float(pepsi[0]),
        }
        if args.spatial_shell:
            shell = [p for p in cloud.get_scatterers() if p.hydration_volume > 0]
            result['cases'][label]['shell_points'] = len(shell)
            result['cases'][label]['shell_volume_A3'] = sum(p.hydration_volume for p in shell)
            result['cases'][label]['shell_spacing_A'] = args.shell_spacing
        curves[label] = (ours, pepsi)
        np.savetxt(args.output / f'{label}.csv',
                   np.column_stack((q, ours, pepsi, normalized, reference, rel,
                                    vacuum / vacuum[0], data[:, 2] / data[0, 2])),
                   delimiter=',', header='q,pyseb,pepsi,pyseb_normalized,pepsi_normalized,relative_difference,pyseb_vacuum_normalized,pepsi_vacuum_normalized')
    dry, wet = curves['hydration_0pct'], curves['hydration_5pct']
    result['hydration_I0_ratio'] = {'pyseb': float(wet[0][0] / dry[0][0]),
                                    'pepsi': float(wet[1][0] / dry[1][0])}
    if args.spatial_shell:
        # Sensitivity diagnostic: match the hydration-induced I0 ratio to
        # isolate geometry from each program's shell-volume convention.
        # This is sample-specific and is never a production default.
        wanted = np.sqrt(dry[0][0]) * (np.sqrt(wet[1][0]/dry[1][0])-1)
        present = np.sqrt(wet[0][0])-np.sqrt(dry[0][0])
        options.hydration_shell_contrast *= wanted/present
        cloud.set_options(options)
        matched = np.asarray(world.EvaluateFormFactorUnnormalized('protein', {}, q))
        difference = (matched/matched[0])/(wet[1]/wet[1][0])-1
        result['forward_hydration_ratio_matched'] = {
            'contrast_e_per_A3': options.hydration_shell_contrast,
            'normalized_relative_rms': float(np.sqrt(np.mean(difference**2))),
            'normalized_relative_max': float(np.max(np.abs(difference)))}
        np.savetxt(args.output/'forward-matched.csv', np.column_stack((q,matched,wet[1])),
                   delimiter=',',header='q,pyseb_forward_matched,pepsi')
    (args.output / 'comparison.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
