"""Tutorial layout coverage and the saved-fit/compact display adapters."""
import ast
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import sys
import unittest
from unittest.mock import patch

import pyseb

TUTORIALS = Path(__file__).resolve().parents[3] / 'docs/tutorials'
HAS_VIEWER = bool(importlib.util.find_spec('plotly') and importlib.util.find_spec('pxr'))


def load_helper(name):
    with patch.object(sys, 'path', [str(TUTORIALS), *sys.path]):
        spec = importlib.util.spec_from_file_location(name, TUTORIALS / f'{name}.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    return module


class TestTutorialLayoutRequests(unittest.TestCase):
    def test_every_tutorial_structure_call_requests_both_layouts(self):
        checked = 0
        for path in TUTORIALS.glob('*.qmd'):
            for block in re.findall(r'```(?:\{python\}|python)\n(.*?)\n```', path.read_text(), re.S):
                for node in ast.walk(ast.parse(block)):
                    if not isinstance(node, ast.Call):
                        continue
                    name = (node.func.attr if isinstance(node.func, ast.Attribute)
                            else node.func.id if isinstance(node.func, ast.Name) else '')
                    if name not in ('visualize', 'structure_figure'):
                        continue
                    keyword = next((k for k in node.keywords if k.arg == 'layouts'), None)
                    self.assertIsNotNone(keyword, f'{path.name}:{node.lineno}')
                    self.assertEqual(set(ast.literal_eval(keyword.value)), {'random', 'readable'})
                    checked += 1
        self.assertGreater(checked, 30)

    def test_titin_provenance_accepts_only_the_recorded_display_migration(self):
        source = load_helper('_titin_sources')
        manifest = json.loads((TUTORIALS / 'data/titin-results/manifest.json').read_text())
        digest = hashlib.sha256('\n'.join(source.generation_cells(TUTORIALS)).encode()).hexdigest()
        self.assertEqual(digest, manifest['viewer_migration']['generation_sha256'])
        for name in ('curves.npz', 'results.json'):
            self.assertEqual(hashlib.sha256((TUTORIALS / 'data/titin-results' / name).read_bytes()).hexdigest(),
                             manifest['files'][name])


@unittest.skipUnless(HAS_VIEWER, 'optional visualization dependencies not installed')
class TestTutorialViewAdapters(unittest.TestCase):
    def test_compact_gallery_preserves_both_geometries_and_connectivity(self):
        compact = load_helper('_compact_view')
        world = pyseb.World('CompactTest')
        graph = world.Add('SolidSphere', 'core')
        previous = 'core.surface#anchor'
        for index in range(3):
            world.Link('ThinRod', f'arm{index}.end1', previous, 'segment')
            previous = f'arm{index}.end2'
        world.Add(graph, 'micelle')
        parameters = {'R_core': 2., 'L_segment': 3., 'beta_core': 1., 'beta_segment': .5}
        options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
        options.seed = 7
        options.surface_samples = 8
        options.orientation_trials = 2
        options.relaxation_sweeps = 0
        source = pyseb.structure_figure(world, 'micelle', parameters, options,
                                        layouts=['random', 'readable'])
        figure = compact.compact_structure_view(world, parameters, options, 'Test')
        self.assertLess(len(figure.data), len(source.data))
        self.assertEqual({v['layout'] for v in figure._pyseb_variants}, {'random', 'readable'})

        def vertices(data):
            return Counter((float(x), float(y), float(z))
                           for trace in data if trace['meta']['pyseb']['role'] == 'geometry'
                           for x, y, z in zip(trace['x'], trace['y'], trace['z']) if x is not None)

        for original, reduced in zip(source._pyseb_variants, figure._pyseb_variants):
            self.assertEqual(vertices(original['data']), vertices(reduced['data']))
            graph = {t['meta']['pyseb']['instance']: t['meta']['pyseb']['neighbors']
                     for t in reduced['data']}
            self.assertTrue(any(graph.values()))
            for node, neighbors in graph.items():
                for neighbor in neighbors:
                    self.assertIn(node, graph[neighbor])
        self.assertIn('pyseb-fold', figure.to_html(full_html=False, include_plotlyjs=False))

    def test_every_titin_fit_rebuilds_a_paired_view_without_optimization(self):
        module = load_helper('_titin_view_models')
        results = json.loads((TUTORIALS / 'data/titin-results/results.json').read_text())
        views = module.FittedTitinViews(TUTORIALS, results)
        for name in ('analytical-collapsed', 'analytical-constrained', 'debye-domains',
                     'atomic-displaced-solvent', 'atomic-hydrated', 'atomic-weighted-hybrid'):
            with self.subTest(name=name):
                world, parameters, options = views.model(name)
                options.surface_samples = 8
                options.curve_samples = 16
                options.orientation_trials = 2
                options.relaxation_sweeps = 0
                figure = pyseb.structure_figure(world, 'titin', parameters, options,
                                                layouts=['random', 'readable'])
                self.assertEqual({v['layout'] for v in figure._pyseb_variants}, {'random', 'readable'})
                self.assertEqual({v['seed'] for v in figure._pyseb_variants}, {'42'})
                self.assertTrue(figure.data)
        self.assertNotIn('least_squares', views.namespace)
