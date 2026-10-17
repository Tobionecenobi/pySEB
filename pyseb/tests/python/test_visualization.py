"""Public viewer API, optional imports and real OpenUSD conversion."""
import importlib.util
from pathlib import Path
import subprocess
import sys
import unittest
from unittest.mock import patch

import pyseb


class TestVisualization(unittest.TestCase):
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


if __name__ == '__main__':
    unittest.main()
