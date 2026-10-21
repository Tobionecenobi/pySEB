"""Scene caches survive equation edits but reject stale or incomplete exports."""
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location('scene_cache', ROOT / 'docs/tutorials/_scene_cache.py')
cache = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cache)


class TestSceneCache(unittest.TestCase):
    def test_key_tracks_fit_and_display_settings(self):
        options = SimpleNamespace(seed=127, opacity_overrides={'anchor': .25})
        with patch.object(cache, 'implementation_key', return_value='implementation'):
            def key(parameters=None, arms=96):
                return cache.scene_key(ROOT, ('cylinder', 1), arms,
                                       parameters or {'R_core': 42.}, options, 'Micelle')
            original = key()
            self.assertEqual(original, key())
            self.assertNotEqual(original, key({'R_core': 43.}))
            self.assertNotEqual(original, key(arms=97))
            options.seed = 128
            self.assertNotEqual(original, key())

    def test_equation_edits_preserve_scene_key_but_model_edits_invalidate_it(self):
        inputs = [('docs/tutorials/compare-f127-coronas.qmd', json.dumps(['old equation', 'scene_key(ROOT)']).encode()),
                  ('docs/tutorials/_f127-model.qmd', b'original graph')]
        cache.implementation_key.cache_clear()
        self.addCleanup(cache.implementation_key.cache_clear)
        with patch.object(cache, 'source_inputs', side_effect=lambda root: iter(inputs)):
            original = cache.implementation_key(ROOT)
            inputs[0] = (inputs[0][0], json.dumps(['new equation', 'scene_key(ROOT)']).encode())
            cache.implementation_key.cache_clear()
            self.assertEqual(original, cache.implementation_key(ROOT))
            inputs[1] = (inputs[1][0], b'changed graph')
            cache.implementation_key.cache_clear()
            self.assertNotEqual(original, cache.implementation_key(ROOT))

    def test_missing_changed_and_partial_outputs_are_regenerated(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'micelle.usda'
            outputs = [path, path.with_suffix('.html'), path.parent / 'plotly.min.js']
            self.assertFalse(cache.reuse_scene(path, 'key'))
            for output in outputs:
                output.write_text('complete')
            cache.record_scene(path, 'key')
            self.assertTrue(cache.reuse_scene(path, 'key'))
            self.assertFalse(cache.reuse_scene(path, 'changed settings'))
            path.write_text('interrupted export')
            self.assertFalse(cache.reuse_scene(path, 'key'))
            path.write_text('complete')
            outputs[2].unlink()
            self.assertFalse(cache.reuse_scene(path, 'key'))
