"""Public viewer API, optional imports and real OpenUSD conversion."""
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from urllib.parse import unquote, urlsplit

import pyseb


class TestVisualization(unittest.TestCase):
    @unittest.skipUnless(importlib.util.find_spec('plotly'), 'Plotly not installed')
    def test_browser_page_survives_script_exit(self):
        code = '''
import webbrowser
from pyseb._viewer import figure_class
webbrowser.open = lambda url, *args, **kwargs: print(url)
figure_class()().show(renderer='browser')
'''
        with tempfile.TemporaryDirectory() as directory:
            result = subprocess.run(
                [sys.executable, '-c', code], check=True, capture_output=True,
                text=True, env={**os.environ, 'TMPDIR': directory, 'TEMP': directory, 'TMP': directory},
            )
            url = urlsplit(result.stdout.strip())
            self.assertEqual(url.scheme, 'file')
            path = Path(unquote(url.path))
            self.assertTrue(path.is_file(), 'browser HTML disappeared when Python exited')
            self.assertIn('pyseb-settings', path.read_text())

    def test_import_without_optional_dependencies(self):
        code = '''
import builtins
original = builtins.__import__
def blocked(name, *args, **kwargs):
    if name.split('.')[0] in {'plotly', 'pxr'}:
        raise ImportError('blocked optional dependency')
    return original(name, *args, **kwargs)
builtins.__import__ = blocked
import pyseb
try:
    pyseb.structure_figure(None, 'assembly')
except ImportError as exc:
    assert "pip install 'pyseb[visualization]'" in str(exc)
else:
    raise AssertionError('expected actionable missing-dependency error')
'''
        subprocess.run([sys.executable, '-c', code], check=True)

    @unittest.skipUnless(importlib.util.find_spec('plotly') and importlib.util.find_spec('pxr'),
                         'optional visualization dependencies not installed')
    def test_real_structure_and_display(self):
        world = pyseb.World('viewer')
        graph = world.Add('SolidSphere', 'sphere')
        graph = world.Link('ThinRod', 'rod.end1', 'sphere.surface#anchor')
        world.Add(graph, 'assembly')
        options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
        options.surface_samples = 24
        options.curve_samples = 8
        paths = []

        class RecordingWorld:
            def export_usd(self, structure, path, parameters, export_options):
                paths.append(Path(path))
                world.export_usd(structure, path, parameters, export_options)

        parameters = {'R_sphere': 3.0, 'L_rod': 2.0}
        figure = pyseb.structure_figure(RecordingWorld(), 'assembly', parameters, options, height=420)
        self.assertTrue(figure.data)
        self.assertIn('mesh3d', {trace.type for trace in figure.data})
        self.assertIn('scatter3d', {trace.type for trace in figure.data})
        self.assertEqual(figure.layout.height, 420)
        self.assertFalse(paths[0].exists())
        with patch('plotly.graph_objects.Figure.show') as show:
            result = pyseb.visualize(world, 'assembly', parameters, options, renderer='json')
        self.assertTrue(result.data)
        show.assert_called_once_with(renderer='json', config={'displaylogo': False, 'scrollZoom': True})


    @unittest.skipUnless(importlib.util.find_spec('plotly') and importlib.util.find_spec('pxr'),
                         'optional visualization dependencies not installed')
    def test_hidden_geometry_stays_hidden(self):
        from pxr import Usd, UsdGeom
        from pyseb.visualization import usd_to_plotly

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'visibility.usda'
            stage = Usd.Stage.CreateNew(str(path))
            parent = UsdGeom.Xform.Define(stage, '/hidden')
            parent.CreateVisibilityAttr(UsdGeom.Tokens.invisible)
            for name in ('/visible', '/hidden/mesh', '/transparent'):
                mesh = UsdGeom.Mesh.Define(stage, name)
                mesh.CreatePointsAttr([(0, 0, 0), (1, 0, 0), (0, 1, 0)])
                mesh.CreateFaceVertexCountsAttr([3])
                mesh.CreateFaceVertexIndicesAttr([0, 1, 2])
                if name == "/transparent":
                    mesh.CreateDisplayOpacityAttr([0.0])
            stage.GetRootLayer().Save()
            figure = usd_to_plotly(path)
        self.assertEqual([trace.name for trace in figure.data], ['/visible'])


    @unittest.skipUnless(importlib.util.find_spec('plotly') and importlib.util.find_spec('pxr'),
                         'optional visualization dependencies not installed')
    def test_realizations_and_reference_metadata(self):
        world = pyseb.World('realizations')
        graph = world.Add('SolidSphere', 'core')
        world.Link('ThinRod', 'rod.end1', 'core.surface#anchor')
        world.Add(graph, 'assembly')
        options = pyseb.USDExportOptions(pyseb.LengthUnit.Nanometer)
        options.seed = 7
        options.orientation_trials = 2
        options.relaxation_sweeps = 0
        options.surface_samples = 8
        snapshots = []

        class RecordingWorld:
            def export_usd(self, name, path, parameters, selected):
                snapshots.append((selected.seed, int(selected.layout_mode), selected.meters_per_unit))
                world.export_usd(name, path, parameters, selected)

        figure = pyseb.structure_figure(
            RecordingWorld(), 'assembly', {'R_core': 3.0, 'L_rod': 4.0}, options,
            seeds=[7, 8], layouts=['random', 'readable', 'readable'],
        )
        self.assertEqual([(seed, mode) for seed, mode, _ in snapshots], [
            (7, int(pyseb.LayoutMode.Random)), (8, int(pyseb.LayoutMode.Random)),
            (7, int(pyseb.LayoutMode.Readable)), (8, int(pyseb.LayoutMode.Readable)),
        ])
        self.assertTrue(all(unit == options.meters_per_unit for _, _, unit in snapshots))
        self.assertEqual(options.seed, 7)
        self.assertEqual(options.layout_mode, pyseb.LayoutMode.Random)
        self.assertEqual(len(figure._pyseb_variants), 4)
        reference_traces = [t for t in figure.data if t.meta['pyseb']['role'] == 'reference']
        self.assertTrue(reference_traces)
        self.assertTrue(all(t.visible is False for t in reference_traces))
        self.assertEqual({t.meta['pyseb']['instance'] for t in reference_traces}, {'core', 'rod'})
        self.assertNotEqual(figure._pyseb_variants[0]['data'][1]['x'],
                            figure._pyseb_variants[1]['data'][1]['x'])
        plain = pyseb.structure_figure(world, 'assembly', {'R_core': 3., 'L_rod': 4.},
                                       options, controls=False)
        import plotly.graph_objects as go
        self.assertIs(type(plain), go.Figure)
        self.assertFalse(any(t.meta['pyseb']['role'] == 'reference' for t in plain.data))
        with self.assertRaisesRegex(ValueError, 'layouts must contain'):
            pyseb.structure_figure(world, 'assembly', layouts=['invalid'])
        with self.assertRaisesRegex(ValueError, 'require controls'):
            pyseb.structure_figure(world, 'assembly', controls=False, seeds=[1])

    @unittest.skipUnless(importlib.util.find_spec('plotly') and importlib.util.find_spec('pxr'),
                         'optional visualization dependencies not installed')
    @unittest.skipUnless(importlib.util.find_spec('IPython'), 'IPython not installed')
    def test_notebook_html_and_user_post_script(self):
        from pyseb.visualization import usd_to_plotly
        from pxr import Usd, UsdGeom

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'mesh.usda'
            stage = Usd.Stage.CreateNew(str(path))
            mesh = UsdGeom.Mesh.Define(stage, '/mesh')
            mesh.CreatePointsAttr([(0, 0, 0), (1, 0, 0), (0, 1, 0)])
            mesh.CreateFaceVertexCountsAttr([3])
            mesh.CreateFaceVertexIndicesAttr([0, 1, 2])
            stage.GetRootLayer().Save()
            figure = usd_to_plotly(path)
        notebook = type('Notebook', (), {'kernel': object()})()
        with patch('IPython.get_ipython', return_value=notebook), patch('IPython.display.display') as display:
            figure.show(renderer='notebook', config={'scrollZoom': True})
        html = display.call_args.args[0]['text/html']
        self.assertIn('pyseb-settings', html)
        self.assertIn('plotly.js', html)
        self.assertEqual(display.call_args.kwargs, {'raw': True})
        exported = figure.to_html(post_script='window.pysebUserScript = true;', full_html=False)
        self.assertIn('window.pysebUserScript = true;', exported)
        self.assertIn('pyseb-settings', exported)


    @unittest.skipUnless(importlib.util.find_spec('plotly') and importlib.util.find_spec('pxr'),
                         'optional visualization dependencies not installed')
    def test_scattering_metadata_uses_tags_and_material_volumes(self):
        import math

        world = pyseb.World('density')
        graph = world.Add('SolidSphere', 'first', 'shared')
        world.Link('GaussianPolymer', 'coil.end1', 'first.surface#anchor')
        world.Add(graph, 'sample')
        params = {'R_shared': 2.0, 'beta_shared': -3.0, 'Rg_coil': 4.0, 'beta_coil': 0.0}
        figure = pyseb.structure_figure(world, 'sample', params)
        infos = {t.meta['pyseb']['instance']: t.meta['pyseb'] for t in figure.data
                 if t.meta['pyseb']['role'] == 'geometry'}
        self.assertEqual(infos['first']['beta'], -3.0)
        self.assertAlmostEqual(infos['first']['volume'], 4 * math.pi * 8 / 3)
        self.assertEqual(infos['coil']['beta'], 0.0)
        self.assertIsNone(infos['coil']['volume'])
        explicit = pyseb.structure_figure(world, 'sample', params, volumes={'shared': 12., 'coil': 20.})
        infos = {t.meta['pyseb']['instance']: t.meta['pyseb'] for t in explicit.data
                 if t.meta['pyseb']['role'] == 'geometry'}
        self.assertEqual(infos['first']['volume'], 12.)
        self.assertEqual(infos['coil']['volume'], 20.)
        self.assertEqual(infos['coil']['volume_source'], 'supplied material volume')
        for invalid in (0, -1, float('nan'), float('inf')):
            with self.assertRaisesRegex(ValueError, 'volumes must be finite and positive'):
                pyseb.structure_figure(world, 'sample', params, volumes={'coil': invalid})

    @unittest.skipUnless(importlib.util.find_spec('plotly') and importlib.util.find_spec('pxr'),
                         'optional visualization dependencies not installed')
    def test_analytic_material_volumes_and_cloud_beta(self):
        import math

        cases = [
            ('SolidSphere', {'R_shape': 2.}, 32 * math.pi / 3),
            ('SolidSphericalShell', {'Ri_shape': 1., 'Ro_shape': 2.}, 28 * math.pi / 3),
            ('SolidCylinder', {'R_shape': 2., 'L_shape': 3.}, 12 * math.pi),
            ('Spheroid', {'a_shape': 2., 'c_shape': 3.}, 16 * math.pi),
        ]
        for model, dims, expected in cases:
            with self.subTest(model=model):
                world = pyseb.World('volume')
                world.Add(model, 'shape')
                figure = pyseb.structure_figure(world, 'shape', {**dims, 'beta_shape': 5.})
                self.assertAlmostEqual(figure.data[0].meta['pyseb']['volume'], expected)
        cloud = pyseb.DebyeSphereCloud([
            pyseb.SphereScatterer(0, 0, 0, 1., 3.),
            pyseb.SphereScatterer(1, 0, 0, 1., -1.),
        ])
        world = pyseb.World('cloud_beta')
        world.Add(cloud, 'cloud')
        figure = pyseb.structure_figure(world, 'cloud')
        info = next(t.meta['pyseb'] for t in figure.data if t.meta['pyseb']['role'] == 'geometry')
        self.assertEqual(info['beta'], 2.)
        self.assertIsNone(info['volume'], 'overlapping sphere volumes must not be summed as material volume')


if __name__ == '__main__':
    unittest.main()
