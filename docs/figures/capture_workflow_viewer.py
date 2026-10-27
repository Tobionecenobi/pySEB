"""Capture the real pySEB Plotly viewer with a mouse-hover label.

Requires pyseb[visualization], NumPy, Playwright and its Chromium browser.
Run from the repository root with the Python environment containing those tools.
Set PLAYWRIGHT_BROWSERS_PATH if Chromium is installed outside its default cache.
"""
from pathlib import Path
import sys
import tempfile
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np
import pyseb
from playwright.sync_api import sync_playwright

out=Path(__file__).with_name('workflow-viewer-hover.png')
rng=np.random.default_rng(42)
positions=rng.uniform(-.3,.3,(20,3))
cloud=pyseb.DebyeSphereCloud([
    pyseb.SphereScatterer(float(x),float(y),float(z),0.,.1)
    for x,y,z in positions])
cloud.addReferencePoint('join',*positions[positions[:,2].argmax()].tolist())
world=pyseb.World('WorkflowViewer')
graph=world.Add(cloud,'cloud')
world.Link('SolidCylinder','cylinder.bottom','cloud.join')
world.Link('SolidSphere','sphere.center','cylinder.top')
world.Add(graph,'assembly')
params=dict(beta_cylinder=3.,R_cylinder=1.,L_cylinder=4.,beta_sphere=4.,R_sphere=1.5)
options=pyseb.USDExportOptions(pyseb.LengthUnit.Nanometer)
options.zero_radius_marker_size=.08
options.layout_mode=pyseb.LayoutMode.Readable
options.surface_samples=40
options.orientation_trials=30
options.relaxation_sweeps=10
fig=pyseb.structure_figure(world,'assembly',params,options,height=440,layouts=['readable'])
fig.update_layout(
    width=568, paper_bgcolor='white', hoverlabel=dict(font_size=20),
    scene_camera=dict(eye=dict(x=1.25,y=1.25,z=.9)),
    title_font_size=18,
)
with tempfile.TemporaryDirectory(prefix='pyseb-viewer-capture-') as directory:
    html=Path(directory)/'viewer.html'
    fig.write_html(html,include_plotlyjs=True,config={'displaylogo':False,'scrollZoom':True})
    with sync_playwright() as p:
        browser=p.chromium.launch(args=['--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader'])
        page=browser.new_page(viewport={'width':568,'height':440},device_scale_factor=2)
        errors=[]
        page.on('pageerror',lambda error: errors.append(str(error)))
        page.goto(html.as_uri())
        page.add_style_tag(content='html,body{margin:0;padding:0;background:white;overflow:hidden}')
        page.wait_for_selector('.pyseb-settings',timeout=30000)
        page.wait_for_function("document.querySelector('.plotly-graph-div')._fullLayout.scene._scene.glplot !== undefined")
        page.wait_for_timeout(1000)
        found=False
        for x,y in [(x,y) for y in (240,200,280,160,320) for x in range(140,441,30)]:
            page.mouse.move(x,y)
            page.wait_for_timeout(90)
            tooltip=page.locator('.hovertext')
            if tooltip.count() and tooltip.first.is_visible():
                label=tooltip.first.text_content()
                if 'β' in label or 'sphere' in label or 'cylinder' in label:
                    found=True
                    print('Captured hover:',repr(label),flush=True)
                    break
        assert found,'No real geometry hover label found'
        assert not errors,errors
        page.screenshot(path=str(out))
        browser.close()
print(out,flush=True)
