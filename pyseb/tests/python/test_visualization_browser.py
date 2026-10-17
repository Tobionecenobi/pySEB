"""Optional Chromium integration tests for the in-frame structure controls.

Run with Playwright installed and a Chromium browser available. No browser or
Playwright dependency is needed for ordinary pyseb imports or unit tests.
"""
import importlib.util
from pathlib import Path
import tempfile
import unittest

import pyseb

try:
    _HAS_PLAYWRIGHT = importlib.util.find_spec("playwright.sync_api") is not None
except ModuleNotFoundError:
    _HAS_PLAYWRIGHT = False


@unittest.skipUnless(_HAS_PLAYWRIGHT, 'Playwright not installed')
class TestViewerBrowser(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from playwright.sync_api import sync_playwright

        cls.workspace = tempfile.TemporaryDirectory()
        world = pyseb.World('ControlsPreview')
        graph = world.Add('SolidSphere', 'core')
        world.Link('ThinRod', 'rod.end1', 'core.surface#anchor')
        world.Add(graph, 'assembly')
        options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
        options.surface_samples = 12
        options.curve_samples = 8
        options.orientation_trials = 2
        options.relaxation_sweeps = 0
        options.color_overrides = {'core': (0.1, 0.55, 0.8), 'rod': (0.95, 0.45, 0.1)}
        figure = pyseb.structure_figure(
            world, 'assembly', {'R_core': 3.0, 'L_rod': 5.0}, options,
            seeds=[1], layouts=['readable'], height=480,
        )
        path = Path(cls.workspace.name) / 'viewer.html'
        figure.write_html(path)
        cls.url = path.as_uri()

        nested = pyseb.World('NestedControls')
        diblock = nested.Add('GaussianPolymer', 'A')
        nested.Link('GaussianPolymer', 'B.end1', 'A.end2')
        star = nested.Add(diblock, 'diblock1')
        for index in range(2, 5):
            nested.Link(diblock, f'diblock{index}:A.end1', 'diblock1:A.end1')
        chain = nested.Add(star, 'star4')
        nested.Link(star, 'star40:diblock1:B.end2', 'star4:diblock4:B.end2')
        nested.Add(chain, 'chain')
        nested_figure = pyseb.structure_figure(
            nested, 'chain', {'Rg_A': 3.0, 'Rg_B': 4.0, 'beta_A': 1.0, 'beta_B': 0.8},
            options, seeds=[1], height=540,
        )
        nested_path = Path(cls.workspace.name) / 'nested.html'
        nested_figure.write_html(nested_path)
        cls.nested_url = nested_path.as_uri()
        contrast_world = pyseb.World('ContrastControls')
        contrast_graph = contrast_world.Add('SolidSphere', 'positive')
        contrast_world.Link('SolidSphere', 'negative.surface#in', 'positive.surface#out')
        contrast_world.Link('SolidSphere', 'zero.surface#in', 'negative.surface#out')
        contrast_world.Link('GaussianPolymer', 'coil.end1', 'zero.surface#coil')
        contrast_world.Add(contrast_graph, 'sample')
        contrast = pyseb.structure_figure(contrast_world, 'sample', {
            'R_positive': 1., 'beta_positive': 2., 'R_negative': 2., 'beta_negative': -2.,
            'R_zero': 1., 'beta_zero': 0., 'Rg_coil': 2., 'beta_coil': 1.,
        }, options, layouts=['readable'], height=540)
        contrast_path = Path(cls.workspace.name) / 'contrast.html'
        contrast.write_html(contrast_path)
        cls.contrast_url = contrast_path.as_uri()
        cls.playwright = sync_playwright().start()
        if not Path(cls.playwright.chromium.executable_path).is_file():
            cls.playwright.stop()
            cls.workspace.cleanup()
            raise unittest.SkipTest("Playwright Chromium not installed")
        cls.browser = cls.playwright.chromium.launch(args=[
            '--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader',
        ])

    @classmethod
    def tearDownClass(cls):
        cls.browser.close()
        cls.playwright.stop()
        cls.workspace.cleanup()

    def setUp(self):
        self.page = self.browser.new_page(viewport={'width': 1000, 'height': 760})
        self.errors = []
        self.page.on('pageerror', lambda error: self.errors.append(str(error)))
        self.page.goto(self.url)
        self.page.wait_for_selector('.pyseb-settings')

    def tearDown(self):
        self.assertEqual(self.errors, [])
        self.page.close()

    def test_panel_changes_geometry_and_preserves_camera(self):
        page = self.page
        page.mouse.move(950, 730)
        page.wait_for_timeout(200)
        toggle = page.get_by_role('button', name='Structure settings', exact=True)
        self.assertEqual(toggle.evaluate('(el)=>getComputedStyle(el).opacity'), '0')
        page.locator('.plotly-graph-div').hover()
        toggle.click()
        self.assertEqual(toggle.get_attribute('aria-expanded'), 'true')
        page.mouse.move(950, 730)
        self.assertTrue(page.locator('.pyseb-panel').is_visible())

        page.get_by_role('checkbox', name='Show core', exact=True).uncheck()
        page.wait_for_function("""() => document.querySelector('.plotly-graph-div').data
            .filter(t=>t.meta.pyseb.instance==='core'&&t.meta.pyseb.role==='geometry')
            .every(t=>t.visible===false)""")
        page.get_by_role('checkbox', name='Reference points', exact=True).check()
        page.wait_for_function("""() => document.querySelector('.plotly-graph-div').data
            .some(t=>t.meta.pyseb.instance==='rod'&&t.meta.pyseb.role==='reference'&&t.visible===true)""")
        page.get_by_label('Color rod', exact=True).evaluate("""el=>{
            el.value='#ff00ff';el.dispatchEvent(new Event('input',{bubbles:true}));} """)
        page.get_by_role('slider', name='Opacity rod', exact=True).fill('0.35')
        page.wait_for_function("""() => document.querySelector('.plotly-graph-div').data
            .some(t=>t.meta.pyseb.instance==='rod'&&t.meta.pyseb.role==='geometry'
                &&t.line.color==='#ff00ff'&&t.opacity===0.35)""")
        page.evaluate("""()=>Plotly.relayout(document.querySelector('.plotly-graph-div'),
            {'scene.camera.eye':{x:4,y:2,z:1}})""")
        page.get_by_role('combobox', name='Realization', exact=True).select_option('3')
        page.wait_for_function("""() => !document.querySelector('.pyseb-panel select').disabled""")
        self.assertEqual(page.evaluate("document.querySelector('.plotly-graph-div').layout.scene.camera.eye.x"), 4)
        page.wait_for_function("""() => document.querySelector('.plotly-graph-div').data
            .some(t=>t.meta.pyseb.instance==='rod'&&t.line&&t.line.color==='#ff00ff'&&t.opacity===0.35)""")
        page.get_by_role('button', name='Reset camera', exact=True).click()
        page.wait_for_function("""() => document.querySelector('.plotly-graph-div').layout.scene.camera.eye.x===1.5""")
        page.get_by_role('button', name='Reset style', exact=True).click()
        page.wait_for_function("""() => document.querySelector('.plotly-graph-div').data
            .filter(t=>t.meta.pyseb.role==='geometry').every(t=>t.visible!==false&&t.opacity===1)""")
        page.get_by_role('button', name='Close structure settings', exact=True).click()
        self.assertFalse(page.locator('.pyseb-panel').is_visible())
        toggle.press('Enter')
        self.assertTrue(page.locator('.pyseb-panel').is_visible())
        page.get_by_role('slider', name='Opacity rod', exact=True).press('Escape')
        self.assertFalse(page.locator('.pyseb-panel').is_visible())
        # Plotly's original camera toolbar is still available in its own corner.
        self.assertTrue(page.locator('.modebar').count())

    def test_model_unfolded_preserves_seed_camera_and_display_settings(self):
        page = self.page
        page.locator('.plotly-graph-div').hover()
        page.get_by_role('button', name='Structure settings', exact=True).click()
        selector = page.get_by_role('combobox', name='Realization', exact=True)
        selector.select_option('1')  # random layout, seed 1
        page.wait_for_function("() => !document.querySelector('.pyseb-panel select').disabled")
        coordinates = """() => document.querySelector('.plotly-graph-div').data
            .filter(t=>t.meta.pyseb.role==='geometry').map(t=>[t.x,t.y,t.z])"""
        original = page.evaluate(coordinates)
        page.get_by_role('checkbox', name='Show core', exact=True).uncheck()
        page.get_by_role('slider', name='Opacity rod', exact=True).fill('0.45')
        page.get_by_role('checkbox', name='Reference points', exact=True).check()
        page.evaluate("""()=>Plotly.relayout(document.querySelector('.plotly-graph-div'),
            {'scene.camera.eye':{x:3,y:2,z:1}})""")
        page.get_by_role('button', name='See unfolded', exact=True).click()
        page.get_by_role('button', name='See model', exact=True).wait_for()
        self.assertEqual(selector.input_value(), '3')  # readable layout, same seed
        self.assertEqual(page.evaluate("document.querySelector('.plotly-graph-div').layout.scene.camera.eye.x"), 3)
        page.get_by_role('button', name='See model', exact=True).click()
        page.get_by_role('button', name='See unfolded', exact=True).wait_for()
        self.assertEqual(selector.input_value(), '1')
        self.assertEqual(page.evaluate(coordinates), original)
        self.assertTrue(page.get_by_role('checkbox', name='Reference points', exact=True).is_checked())
        self.assertFalse(page.get_by_role('checkbox', name='Show core', exact=True).is_checked())
        self.assertTrue(page.evaluate("""() => document.querySelector('.plotly-graph-div').data
            .some(t=>t.meta.pyseb.instance==='rod'&&t.meta.pyseb.role==='geometry'&&t.opacity===0.45)"""))
        self.assertEqual(page.evaluate("document.querySelector('.plotly-graph-div').layout.scene.camera.eye.x"), 3)

    def test_hierarchy_cascades_with_leaf_overrides_and_mixed_parents(self):
        page = self.page
        page.goto(self.nested_url)
        page.wait_for_selector('.pyseb-settings')
        page.locator('.plotly-graph-div').hover()
        page.get_by_role('button', name='Structure settings', exact=True).click()
        page.get_by_role('button', name='Expand star4', exact=True).click()
        page.get_by_role('button', name='Expand star4:diblock4', exact=True).click()
        leaf = page.get_by_role('checkbox', name='Show star4:diblock4:A', exact=True)
        self.assertTrue(leaf.is_visible())
        self.assertFalse(page.get_by_role('checkbox', name='Show star40:diblock1:A', exact=True).is_visible())

        star = page.get_by_role('checkbox', name='Show star4', exact=True)
        star.uncheck()
        page.wait_for_function("""() => document.querySelector('.plotly-graph-div').data
            .filter(t=>t.meta.pyseb.role==='geometry'&&t.meta.pyseb.instance.replaceAll('/',':').startsWith('star4:'))
            .every(t=>t.visible===false)""")
        self.assertTrue(page.evaluate("""() => document.querySelector('.plotly-graph-div').data
            .filter(t=>t.meta.pyseb.role==='geometry'&&t.meta.pyseb.instance.replaceAll('/',':').startsWith('star40:'))
            .every(t=>t.visible!==false)"""))
        leaf.check()
        self.assertTrue(star.evaluate('(el)=>el.indeterminate'))
        self.assertTrue(page.get_by_role('checkbox', name='Show star4:diblock4', exact=True)
                        .evaluate('(el)=>el.indeterminate'))
        page.wait_for_function("""() => document.querySelector('.plotly-graph-div').data
            .filter(t=>t.meta.pyseb.role==='geometry'&&t.meta.pyseb.instance.replaceAll('/',':').startsWith('star4:'))
            .every(t=>t.visible===(t.meta.pyseb.instance.replaceAll('/',':')==='star4:diblock4:A'))""")
        page.get_by_role('checkbox', name='Reference points', exact=True).check()
        page.wait_for_function("""() => document.querySelector('.plotly-graph-div').data
            .some(t=>t.meta.pyseb.role==='reference'&&t.meta.pyseb.instance.replaceAll('/',':')==='star4:diblock4:A'&&t.visible===true)""")
        self.assertTrue(page.evaluate("""() => document.querySelector('.plotly-graph-div').data
            .filter(t=>t.meta.pyseb.role==='reference'&&t.meta.pyseb.instance.replaceAll('/',':')==='star4:diblock4:B')
            .every(t=>t.visible===false)"""))

        page.get_by_label('Color star4', exact=True).evaluate("""el=>{
            el.value='#dd2244';el.dispatchEvent(new Event('input',{bubbles:true}));}""")
        page.get_by_role('slider', name='Opacity star4:diblock4', exact=True).fill('0.4')
        page.get_by_label('Color star4:diblock4:A', exact=True).evaluate("""el=>{
            el.value='#22aadd';el.dispatchEvent(new Event('input',{bubbles:true}));}""")
        page.get_by_role('slider', name='Opacity star4:diblock4:A', exact=True).fill('0.7')
        page.wait_for_function("""() => document.querySelector('.plotly-graph-div').data
            .filter(t=>t.meta.pyseb.role==='geometry'&&t.meta.pyseb.instance.replaceAll('/',':')==='star4:diblock4:B')
            .every(t=>t.line.color==='#dd2244'&&t.opacity===0.4)""")
        page.wait_for_function("""() => document.querySelector('.plotly-graph-div').data
            .filter(t=>t.meta.pyseb.role==='geometry'&&t.meta.pyseb.instance.replaceAll('/',':')==='star4:diblock4:A')
            .every(t=>t.line.color==='#22aadd'&&t.opacity===0.7)""")
        self.assertEqual(page.get_by_role('slider', name='Opacity star4:diblock4', exact=True)
                         .get_attribute('aria-valuetext'), 'Mixed opacity')
        camera = page.evaluate("JSON.stringify(document.querySelector('.plotly-graph-div').layout.scene.camera)")
        page.get_by_role('combobox', name='Realization', exact=True).select_option('1')
        page.wait_for_function("() => !document.querySelector('.pyseb-panel select').disabled")
        self.assertEqual(page.evaluate("JSON.stringify(document.querySelector('.plotly-graph-div').layout.scene.camera)"), camera)
        page.wait_for_function("""() => document.querySelector('.plotly-graph-div').data
            .filter(t=>t.meta.pyseb.role==='geometry'&&t.meta.pyseb.instance.replaceAll('/',':')==='star4:diblock4:A')
            .every(t=>t.line.color==='#22aadd'&&t.opacity===0.7&&t.visible===true)""")
        page.get_by_role('checkbox', name='Show All subunits', exact=True).check()
        page.wait_for_function("""() => document.querySelector('.plotly-graph-div').data
            .filter(t=>t.meta.pyseb.role==='geometry').every(t=>t.visible===true)""")
        page.get_by_role('button', name='Reset style', exact=True).click()
        page.wait_for_function("""() => document.querySelector('.plotly-graph-div').data
            .filter(t=>t.meta.pyseb.role==='geometry').every(t=>t.visible!==false&&t.opacity===1)""")
        self.assertTrue(page.get_by_role('checkbox', name='Show star4:diblock4:A', exact=True).is_visible())

    def test_scattering_opacity_density_zero_sign_and_manual_restore(self):
        page = self.page
        page.goto(self.contrast_url)
        page.wait_for_selector('.pyseb-settings')
        page.locator('.plotly-graph-div').hover()
        page.get_by_role('button', name='Structure settings', exact=True).click()
        page.get_by_role('slider', name='Opacity positive', exact=True).fill('0.4')
        checkbox = page.get_by_role('checkbox', name='Opacity by scattering', exact=True)
        checkbox.check()
        def wait_opacity(name, opacity):
            page.wait_for_function("""([name,expected])=>document.querySelector('.plotly-graph-div').data
                .filter(t=>t.meta.pyseb.role==='geometry'&&t.meta.pyseb.instance===name)
                .every(t=>Math.abs(t.opacity-expected)<1e-9)""", arg=[name, opacity])
        wait_opacity('positive', 0.4)
        wait_opacity('negative', 1.0)
        wait_opacity('zero', 0.0)
        wait_opacity('coil', 0.5)
        page.get_by_role('combobox', name='Opacity weighting', exact=True).select_option('density')
        wait_opacity('positive', 0.4)
        wait_opacity('negative', 0.125)
        wait_opacity('zero', 0.0)
        wait_opacity('coil', 1.0)  # no invented volume for the Gaussian polymer
        self.assertIn('1 subunit has no beta/material volume', page.locator('.pyseb-panel p[role="status"]').inner_text())
        page.get_by_role('button', name='See unfolded', exact=True).click()
        page.get_by_role('button', name='See model', exact=True).wait_for()
        wait_opacity('negative', 0.125)
        self.assertTrue(checkbox.is_checked())
        checkbox.uncheck()
        wait_opacity('positive', 0.4)
        wait_opacity('negative', 1.)
        wait_opacity('zero', 1.)
        wait_opacity('coil', 1.)
        page.evaluate("""()=>document.querySelector('.plotly-graph-div').data.forEach(t=>{
            if(t.meta.pyseb.role==='geometry')t.meta.pyseb.beta=0;
        })""")
        checkbox.check()
        wait_opacity('coil', 0.)  # zero contrast stays transparent even without a material volume
        page.get_by_role('combobox', name='Opacity weighting', exact=True).select_option('beta')
        for name in ('positive', 'negative', 'zero', 'coil'):
            wait_opacity(name, 0.)

    def test_touch_and_overlay_wheel_do_not_move_camera(self):
        page = self.browser.new_page(viewport={'width': 390, 'height': 650}, has_touch=True, is_mobile=True)
        try:
            page.goto(self.url)
            page.wait_for_selector('.pyseb-settings')
            toggle = page.get_by_role('button', name='Structure settings', exact=True)
            toggle.tap()
            self.assertTrue(page.locator('.pyseb-panel').is_visible())
            frame = page.locator('.plotly-graph-div').bounding_box()
            panel = page.locator('.pyseb-panel').bounding_box()
            self.assertLessEqual(panel['x'] + panel['width'], frame['x'] + frame['width'])
            self.assertLessEqual(panel['y'] + panel['height'], frame['y'] + frame['height'])
            camera = page.evaluate("JSON.stringify(document.querySelector('.plotly-graph-div').layout.scene.camera)")
            page.locator('.pyseb-panel').hover()
            page.mouse.wheel(0, 120)
            page.wait_for_timeout(150)
            self.assertEqual(page.evaluate("JSON.stringify(document.querySelector('.plotly-graph-div').layout.scene.camera)"), camera)
        finally:
            page.close()


if __name__ == '__main__':
    unittest.main()
