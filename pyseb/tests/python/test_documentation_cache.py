"""External-input invalidation for cached executable documentation."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location('docs_cache', ROOT / 'docs/_prepare_cache.py')
cache = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cache)


class TestDocumentationCache(unittest.TestCase):
    def setUp(self):
        self.workspace = tempfile.TemporaryDirectory()
        self.addCleanup(self.workspace.cleanup)
        self.root = Path(self.workspace.name)
        (self.root / 'docs/tutorials').mkdir(parents=True)

    def write(self, name, text):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    def key(self):
        return cache.fingerprint(self.root, 'fixed test runtime')

    def test_prose_reuses_results_but_code_changes_refresh_the_cache(self):
        page = self.write('docs/tutorials/example.qmd',
                          'First explanation\n```{python}\nanswer = 42\n```\n')
        before = self.key()
        page.write_text(page.read_text().replace('First explanation', 'Revised explanation'))
        self.assertEqual(before, self.key())
        page.write_text(page.read_text().replace('answer = 42', 'answer = 43'))
        self.assertNotEqual(before, self.key())

    def test_cross_page_titin_generation_inputs_invalidate_the_cache(self):
        page = self.write('docs/tutorials/titin.qmd',
                          'Explanation\n```{python}\n#| tags: [titin-generation]\nvalue = 1\n```\n')
        before = self.key()
        page.write_text(page.read_text().replace('Explanation', 'Updated explanation'))
        self.assertEqual(before, self.key())
        page.write_text(page.read_text().replace('value = 1', 'value = 2'))
        self.assertNotEqual(before, self.key())

    def test_external_inputs_invalidate_even_when_notebook_cells_do_not_change(self):
        self.write('docs/tutorials/example.qmd', '```{python}\nrun_model()\n```\n')
        for name in ('seb/Model.hpp', 'models/Shape.pyseb.yaml',
                     'docs/tutorials/models/Custom.pyseb.yaml',
                     'docs/tutorials/_helper.py', 'docs/tutorials/data/measurement.csv',
                     'docs/tutorials/cpp/example.cpp', 'pyseb/evaluation.py',
                     'cmake/Backend.cmake.in', '.github/workflows/documentation.yml'):
            with self.subTest(input=name):
                path = self.write(name, 'original calculation input')
                before = self.key()
                path.write_text('updated calculation input')
                self.assertNotEqual(before, self.key())

    def test_missing_generated_downloads_refresh_cached_notebooks(self):
        self.write('docs/tutorials/example.qmd', '```{python}\nexport_scene()\n```\n')
        key = self.key()
        self.assertTrue(cache.prepare(self.root, key))
        root_store = self.write('docs/.jupyter_cache/store.txt', 'cached outputs')
        tutorial_store = self.write('docs/tutorials/.jupyter_cache/store.txt', 'cached outputs')
        scene = self.write('docs/tutorials/generated/example.html', 'interactive scene')
        cache.record_outputs(self.root)
        self.assertFalse(cache.prepare(self.root, key))
        self.assertTrue(root_store.exists())
        self.assertTrue(tutorial_store.exists())
        scene.unlink()
        self.assertTrue(cache.prepare(self.root, key))
        self.assertFalse(root_store.exists())
        self.assertFalse(tutorial_store.exists())

    def test_input_change_refreshes_both_notebook_cache_locations(self):
        self.write('docs/tutorials/example.qmd', '```{python}\nrun_model()\n```\n')
        cache.prepare(self.root, self.key())
        root_store = self.write('docs/.jupyter_cache/store.txt', 'cached outputs')
        tutorial_store = self.write('docs/tutorials/.jupyter_cache/store.txt', 'cached outputs')
        self.write('docs/tutorials/data/measurement.csv', 'new data')
        self.assertTrue(cache.prepare(self.root, self.key()))
        self.assertFalse(root_store.exists())
        self.assertFalse(tutorial_store.exists())
