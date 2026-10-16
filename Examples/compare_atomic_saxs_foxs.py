"""Compare pySEB with fixed-parameter IMP/FoXS in a separate Python process.

Run with the pySEB Python, --imp-python /path/to/IMP/python,
--pdb input.pdb and --output /tmp/comparison. No experimental fitting is used.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np


def imp_worker(pdb, output):
    import IMP
    import IMP.atom
    import IMP.core
    import IMP.saxs
    model = IMP.Model()
    hierarchy = IMP.atom.read_pdb(str(pdb), model,
                                 IMP.atom.NonWaterNonHydrogenPDBSelector())
    particles = IMP.atom.get_by_type(hierarchy, IMP.atom.ATOM_TYPE)
    default_table = IMP.saxs.get_default_form_factor_table()
    radii = [default_table.get_radius(p) for p in particles]
    for particle, radius in zip(particles, radii):
        if IMP.core.XYZR.get_is_setup(particle):
            IMP.core.XYZR(particle).set_radius(radius)
        else:
            IMP.core.XYZR.setup_particle(particle, radius)
    surface = IMP.saxs.SolventAccessibleSurface().get_solvent_accessibility(
        IMP.core.XYZRs(particles))
    # Equal hydration amplitude per unit exposure at q=0, not a fitted mapping
    # of physical shell properties. The native SASA algorithms still differ.
    c2 = 0.0167 * 30 / default_table.get_water_form_factor()
    result = {'imp_version': IMP.get_module_version(), 'atoms': len(particles),
              'foxs_water_f0': default_table.get_water_form_factor(),
              'c2_forward_amplitude_matched': c2,
              'surface': list(surface), 'radii': radii, 'modes': {},
              'hydration_only': {}, 'water_f0': {}}
    for mode in ('default_fast', 'reciprocal'):
        profile = IMP.saxs.Profile(0.0, 0.5, 0.005)
        profile.set_average_radius(sum(radii) / len(radii))
        if mode == 'reciprocal':
            table = IMP.saxs.FormFactorTable(IMP.saxs.get_data_path(
                'formfactors-int_tab_solvation.lib'), 0.0, 0.51, 0.005)
            profile.set_ff_table(table)
            profile.calculate_profile_reciprocal_partial(particles, surface)
        else:
            profile.calculate_profile_partial(particles, surface)
        water_f0 = (table.get_water_form_factors()[0] if mode == 'reciprocal'
                    else default_table.get_water_form_factor())
        result['water_f0'][mode] = water_f0
        values = {}
        for name, c1_value, c2_value in [('vacuum', 0, 0), ('solution', 1, 0),
                                        ('hydrated', 1, c2)]:
            profile.sum_partial_profiles(c1_value, c2_value, False)
            values[name] = [profile.get_intensity(i) for i in range(profile.size())]
        result['modes'][mode] = values
        # Isolate hydration's quadratic term without a baseline/model fit.
        hydration = []
        for c2_value in (1.0, -1.0, 0.0):
            profile.sum_partial_profiles(0.0, c2_value, False)
            hydration.append(np.array([profile.get_intensity(i) for i in range(profile.size())]))
        result['hydration_only'][mode] = ((hydration[0]+hydration[1])/2-hydration[2]).tolist()
        result['q'] = [profile.get_q(i) for i in range(profile.size())]
    output.write_text(json.dumps(result, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--imp-python', type=Path)
    parser.add_argument('--pdb', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--imp-worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.imp_worker:
        imp_worker(args.pdb, args.output)
        return
    if args.imp_python is None:
        parser.error('--imp-python is required')
    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root))
    import pyseb
    args.output.mkdir(parents=True, exist_ok=True)
    lines = []
    for line in args.pdb.read_text().splitlines():
        if line.startswith('ENDMDL'):
            break
        if line.startswith('ATOM  '):
            lines.append(line)
    pdb = args.output / 'input.pdb'
    pdb.write_text('\n'.join(lines) + '\nEND\n')
    worker_output = args.output / 'foxs.json'
    command = [str(args.imp_python.resolve()), str(Path(__file__).resolve()),
               '--imp-worker', '--pdb', str(pdb.resolve()),
               '--output', str(worker_output.resolve())]
    process = subprocess.run(command, capture_output=True, text=True)
    (args.output / 'foxs.stdout.txt').write_text(process.stdout + process.stderr)
    process.check_returncode()
    reference = json.loads(worker_output.read_text())
    q = np.array(reference['q'])
    curves = {}
    for name in ('vacuum', 'solution', 'hydrated'):
        options = pyseb.AtomicSAXSOptions()
        options.include_excluded_solvent = name != 'vacuum'
        options.include_hydration = name == 'hydrated'
        options.hydration_shell_contrast = 0.0167
        cloud = pyseb.AtomicSAXSStructureLoader.load_pdb(
            str(pdb), pyseb.AtomicSAXSParameterProfile.protein_defaults(),
            pyseb.StructureParseOptions(), options,
            pyseb.AtomicSAXSCloudBuildOptions())
        if cloud.scatterer_count() != reference['atoms']:
            raise ValueError('Different atom selections')
        world = pyseb.World(name)
        world.Add(cloud, 'protein')
        curves[name] = np.asarray(world.EvaluateFormFactorUnnormalized('protein', {}, q))
    report = {'imp_version': reference['imp_version'], 'atoms': reference['atoms'],
              'input_sha256': hashlib.sha256(pdb.read_bytes()).hexdigest(),
              'foxs_c2': reference['c2_forward_amplitude_matched'],
              'foxs_water_f0': reference['foxs_water_f0'],
              'hydration_mapping': 'equal forward amplitude per unit exposure; native SASA retained',
              'modes': {}}
    for mode, values in reference['modes'].items():
        report['modes'][mode] = {}
        for name, theirs in values.items():
            ours = curves[name]
            theirs = np.array(theirs)
            if not np.all(np.isfinite(theirs)) or np.any(theirs <= 0):
                raise ValueError(f'Invalid reference intensity: {mode}/{name}')
            error = (ours / ours[0]) / (theirs / theirs[0]) - 1
            report['modes'][mode][name] = {
                'normalized_relative_rms': float(np.sqrt(np.mean(error**2))),
                'normalized_relative_max': float(np.max(np.abs(error))),
                'pyseb_I0': float(ours[0]), 'foxs_I0': float(theirs[0])}
            np.savetxt(args.output / f'{mode}-{name}.csv',
                       np.column_stack((q, ours, theirs, error)), delimiter=',',
                       header='q,pyseb,foxs,normalized_relative_difference')
    atoms = pyseb.PDBParser().parse_string(pdb.read_text())
    if len(atoms) == 1:
        parameters = pyseb.AtomicSAXSParameterProfile.protein_defaults().resolve(atoms[0])
        if parameters.implicit_hydrogens == 0:
            vacuum = np.array([parameters.vacuum.evaluate(float(x)) for x in q])
            volume = parameters.excluded_volume
            predicted = (vacuum - 0.334 * volume * np.exp(
                -q*q * volume**(2/3) / (16*np.pi)))**2
            actual = np.array(reference['modes']['reciprocal']['solution'])
            report['single_group_foxs_16pi_formula_max_relative_error'] = float(
                np.max(np.abs(predicted / actual - 1)))
    (args.output / 'comparison.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
