"""Validate hydration components from isolated groups through proteins.

Run with the pySEB Python and --pepsi PATH --imp-python PATH --output PATH.
External programs are used only to produce independent reference outputs.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

from hydration_reference_clouds import pepsi_shell, foxs_hydration, foxs_geometry_exposure, pyseb

ROOT = Path(__file__).resolve().parents[1]
_evaluations = {}


def intensity(cloud, q, label):
    # World owns its added numerical subunit. Retain one World per cloud
    # while varying options, rather than returning a dangling Python handle.
    key = id(cloud)
    if key not in _evaluations:
        world = pyseb.World(label)
        world.Add(cloud, 'sample')
        _evaluations[key] = (world,cloud)
    world = _evaluations[key][0]
    return np.asarray(world.EvaluateFormFactorUnnormalized('sample', {}, q))


def metrics(ours, reference):
    return {'forward_ratio': float(ours[0]/reference[0]),
            'shape_rms_absolute': float(np.sqrt(np.mean((ours/ours[0]-reference/reference[0])**2))),
            'normalized_relative_rms': float(np.sqrt(np.mean(((ours/ours[0])/(reference/reference[0])-1)**2))),
            'relative_rms': float(np.sqrt(np.mean((ours/reference-1)**2))),
            'relative_max': float(np.max(np.abs(ours/reference-1)))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pepsi', type=Path, required=True)
    parser.add_argument('--imp-python', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--proteins', action='store_true', help='Run proteins after the small controls')
    parser.add_argument('--extra-protein',type=Path,action='append',default=[])
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/'environment.json').write_text(json.dumps({
        'pepsi_sha256': hashlib.sha256(args.pepsi.read_bytes()).hexdigest(),
        'pyseb_extension_sha256': hashlib.sha256(Path(pyseb._pyseb.__file__).read_bytes()).hexdigest(),
        'contrast_e_per_A3': .0167,
        'pepsi_point_weight': '3^3 * 0.334, before multiplying by contrast',
        'foxs_mapping': 'unit c2, native water forward amplitude per unit exposure',
        'molecule_specific_fitting': False},indent=2)+'\n')
    fixtures = {}
    for label, names, coordinates in [
        ('C', ['C'], [(0,0,0)]), ('O', ['O'], [(0,0,0)]),
        ('CH3', ['CB'], [(0,0,0)]),
        ('CC2', ['C','C'], [(-1,0,0),(1,0,0)]),
        ('CO2', ['C','O'], [(-1,0,0),(1,0,0)]),
        ('CC10', ['C','C'], [(-5,0,0),(5,0,0)]),
        ('triad', ['C','N','O'], [(0,0,0),(2,0,0),(0,2,0)])]:
        lines = [f'ATOM  {i:5d} {name:>4} ALA A{i:4d}    '
                 f'{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {name[0]:>2}'
                 for i,(name,(x,y,z)) in enumerate(zip(names,coordinates),1)]
        fixtures[label] = '\n'.join(lines)+'\nEND\n'
    fixtures['glycine'] = (ROOT/'Tests/data/minimal_glycine.pdb').read_text()
    if args.proteins:
        fixtures['1TIU'] = (ROOT/'docs/tutorials/data/1TIU-heavy-ensemble.pdb').read_text().split('ENDMDL')[0]+'\nEND\n'
        for path in args.extra_protein:
            fixtures[path.stem] = path.read_text().split('ENDMDL')[0]+'\nEND\n'
    results = {}
    for label, pdb_text in fixtures.items():
        if label not in ('C','O','CH3','CC2','CO2','CC10','triad','glycine'):
            for control in results.values():
                if (control['pepsi']['relative_rms'] > .01 or
                    control['foxs']['reciprocal']['relative_rms'] > 1e-6 or
                    control['foxs']['default_fast']['relative_rms'] > .001):
                    raise RuntimeError('Small-control discrepancy must be resolved before proteins')
        directory = args.output/label
        directory.mkdir(exist_ok=True)
        pdb = directory/'input.pdb'
        pdb.write_text(pdb_text)
        reference_path = directory/'pepsi.dat'
        command = [str(args.pepsi.resolve()), str(pdb.resolve()), '-o', str(reference_path.resolve()),
                   '--dro','5','--noSmearing','--esmodel','0','-ns','101','-ms','0.5','-n','30']
        run = subprocess.run(command, check=True, capture_output=True, text=True)
        (directory/'pepsi.stdout').write_text(run.stdout+run.stderr)
        data = np.loadtxt(reference_path)
        q = data[:,0]
        scale = float(re.search(r'Intensity scaling\.+ : (\S+)',run.stdout)[1])
        count = int(re.search(r'Number of points in the hydration shell\.+ : (\d+)',run.stdout)[1])
        rg = float(re.search(r'Radius of gyration of the hydration shell\.+ : (\S+)',run.stdout)[1])
        atoms = pyseb.PDBParser().parse_string(pdb_text)
        cloud, points, full = pepsi_shell(atoms, pyseb.AtomicSAXSParameterProfile.protein_defaults())
        ours = intensity(cloud, q, label+'-pepsi')
        reference = data[:,4]/scale
        full_intensity = intensity(full,q,label+'-total')
        comparison = {'input_sha256': hashlib.sha256(pdb.read_bytes()).hexdigest(),
                      'atoms': len(atoms), 'pepsi': metrics(ours, reference),
                      'pepsi_total': metrics(full_intensity,data[:,1]/scale),
                      'pepsi_shell_points': count, 'pyseb_shell_points': len(points),
                      'pepsi_shell_rg': rg,
                      'pyseb_shell_rg': float(np.sqrt(np.mean(np.sum((points-points.mean(0))**2,axis=1))))}
        comparison['pepsi_contrast_transfer'] = {}
        if label in ('C','CC2','glycine','1TIU'):
            for percent in (2,10):
                changed = command.copy()
                changed[changed.index('--dro')+1] = str(percent)
                target = directory/f'pepsi-{percent}pct.dat'
                changed[changed.index('-o')+1] = str(target.resolve())
                native_run = subprocess.run(changed,check=True,capture_output=True,text=True)
                (directory/f'pepsi-{percent}pct.stdout').write_text(native_run.stdout)
                changed_data = np.loadtxt(target)
                changed_options = cloud.get_options()
                changed_options.hydration_shell_contrast = .334*percent/100
                cloud.set_options(changed_options)
                comparison['pepsi_contrast_transfer'][str(percent)] = metrics(
                    intensity(cloud,q,label+str(percent)),changed_data[:,4]/scale)
        np.savetxt(directory/'pepsi-hydration.csv', np.column_stack((q,ours,reference)),
                   delimiter=',',header='q,pyseb_hydration_only,pepsi_hydration_only')
        np.savetxt(directory/'pepsi-total.csv',np.column_stack((q,full_intensity,data[:,1]/scale)),
                   delimiter=',',header='q,pyseb_total,pepsi_total')
        worker_output = directory/'foxs.json'
        run = subprocess.run([str(args.imp_python.resolve()),
            str(ROOT/'Examples/compare_atomic_saxs_foxs.py'),'--imp-worker',
            '--pdb',str(pdb.resolve()),'--output',str(worker_output.resolve())],
            check=True,capture_output=True,text=True)
        (directory/'foxs.stdout').write_text(run.stdout+run.stderr)
        native = json.loads(worker_output.read_text())
        (directory/'reference.json').write_text(json.dumps({
            'pdb': pdb_text, 'q': q.tolist(), 'pepsi_hydration': reference.tolist(),
            'pepsi_version': '3.0 May 2020', 'imp_version': native['imp_version'],
            'pepsi_shell_points': count, 'foxs_exposure': native['surface'],
            'foxs_water_f0': native['water_f0'],
            'foxs_hydration': native['hydration_only']},indent=2)+'\n')
        if len(atoms) != native['atoms']:
            raise ValueError('Different atom selections in FoXS')
        comparison['foxs'] = {}
        independent_exposure = foxs_geometry_exposure(atoms)
        comparison['foxs_independent_exposure'] = {}
        for mode in ('default_fast','reciprocal'):
            cloud = foxs_hydration(atoms,native['surface'],mode,native['water_f0'][mode])
            ours = intensity(cloud,q,label+'-'+mode)
            reference = np.array(native['hydration_only'][mode])
            comparison['foxs'][mode] = metrics(ours,reference)
            independent = intensity(foxs_hydration(atoms,independent_exposure,mode,
                native['water_f0'][mode]),q,label+'-independent-'+mode)
            comparison['foxs_independent_exposure'][mode] = metrics(independent,reference)
            np.savetxt(directory/f'foxs-{mode}-hydration.csv',np.column_stack((q,ours,reference)),
                       delimiter=',',header='q,pyseb_hydration_only,foxs_hydration_only')
        results[label] = comparison
        (args.output/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
        print(label, json.dumps(comparison), flush=True)


if __name__ == '__main__':
    main()
