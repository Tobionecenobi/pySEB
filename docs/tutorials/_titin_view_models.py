"""Rebuild tutorial display worlds from saved fit values, without rerunning fits.

Use the topology and cloud-builder definitions in the tagged tutorial source.
The saved scientific bundle remains unchanged; both display layouts are newly
exported by structure_figure from those worlds.
"""
import ast
from pathlib import Path

import numpy as np
import pyseb
from _titin_sources import generation_cells


class FittedTitinViews:
    def __init__(self, tutorial_directory, results):
        self.directory = Path(tutorial_directory)
        self.results = results
        self.namespace = {'np': np, 'pyseb': pyseb,
                          'data_directory': self.directory / 'data'}
        self.worlds = {}
        self.modules = [ast.parse(cell) for cell in generation_cells(self.directory)]
        constants = {'n_tail_sequence', 'linker_sequence', 'c_tail_sequence',
                     'domain_beta', 'linker_beta', 'n_tail_beta', 'c_tail_beta'}
        functions = {'construct_tail_parameters', 'construct_polymer_parameters',
                     'spheroid_parameters', 'build_pdb_world', 'atomic_visual_world'}
        nodes = [node for module in self.modules for node in module.body
                 if (isinstance(node, ast.FunctionDef) and node.name in functions)
                 or (isinstance(node, ast.Assign) and any(
                     isinstance(target, ast.Name) and target.id in constants
                     for target in node.targets))]
        self._execute(nodes)

    def _execute(self, nodes):
        exec(compile(ast.Module(body=nodes, type_ignores=[]),
                     '<titin display builders>', 'exec'), self.namespace)

    def _setup_cell(self, variable):
        matches = [module for module in self.modules if any(
            isinstance(node, ast.Assign) and any(
                isinstance(target, ast.Name) and target.id == variable
                for target in node.targets) for node in module.body)]
        if len(matches) != 1:
            raise ValueError(f'Expected one display setup for {variable}')
        # These setup cells create a topology or parse one PDB model, with no fits.
        self._execute([node for node in matches[0].body if not (
            isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
            and isinstance(node.value.func, ast.Name) and node.value.func.id == 'print')])

    def model(self, name):
        kind, fit_key = {
            'analytical-collapsed': ('analytic', 'analytical_collapsed'),
            'analytical-constrained': ('analytic', 'analytical_constrained'),
            'debye-domains': ('debye', 'debye'),
            'atomic-displaced-solvent': ('atomic-dry', 'fixed_fraser_excluded'),
            'atomic-hydrated': ('atomic-hydrated', 'fixed_fraser_complete'),
            'atomic-weighted-hybrid': ('atomic-hydrated', 'weighted_hybrid'),
        }[name]
        ns = self.namespace
        if kind not in self.worlds:
            if kind == 'analytic':
                self._setup_cell('spheroid_world')
                world = ns['spheroid_world']
            else:
                if 'unscaled_atoms' not in ns:
                    self._setup_cell('pdb_path')
                if kind == 'debye':
                    world = ns['build_pdb_world']()
                else:
                    ns['parsed_ensemble'] = {1: ns['unscaled_atoms']}
                    ns['atomic_profiles'] = {
                        'fraser': pyseb.AtomicSAXSParameterProfile.protein_defaults(
                            pyseb.AtomicSAXSVolumeModel.Fraser1978)}
                    world = ns['atomic_visual_world']('fraser', kind == 'atomic-hydrated')
            self.worlds[kind] = world
        values = self.results['fits'][fit_key]['x']
        parameters = (ns['spheroid_parameters'](values) if kind == 'analytic'
                      else ns['construct_polymer_parameters'](*values[:2]))
        if name == 'atomic-weighted-hybrid':
            for part in ('ntail', 'linker', 'ctail'):
                parameters[f'beta_{part}'] *= values[4]
        options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
        options.seed = 42
        options.layout_mode = pyseb.LayoutMode.Readable
        options.orientation_trials = 32
        options.relaxation_sweeps = 2
        options.zero_radius_marker_size = .8
        options.color_overrides = {
            'domain1': (.10, .55, .78), 'domain2': (.92, .43, .12),
            'linker': (.2, .7, .3), 'ntail': (.6, .6, .6), 'ctail': (.6, .6, .6)}
        return self.worlds[kind], parameters, options
