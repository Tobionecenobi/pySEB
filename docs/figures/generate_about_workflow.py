"""Generate the initial workflow draft with real pySEB parameter sweeps.

Run from the repository root: venv/bin/python docs/figures/generate_about_workflow.py
The mixed assembly follows tutorials/mixed-model-evaluation.qmd (lengths in nm).
Writes about-workflow-draft.svg; the approved about-workflow.svg is maintained
separately and includes subsequent user-supplied vector artwork and edits.
"""
from html import escape
from pathlib import Path
from io import BytesIO
import math
import os
import random
import sys
import xml.etree.ElementTree as ET

os.environ.setdefault("MPLCONFIGDIR", "/tmp/pyseb-figure-matplotlib")
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from matplotlib.mathtext import math_to_image
import numpy as np
import pyseb

rng = np.random.default_rng(42)
positions = rng.uniform(-0.3, 0.3, size=(20, 3))
cloud = pyseb.DebyeSphereCloud([
    pyseb.SphereScatterer(float(x), float(y), float(z), 0.0, 0.1)
    for x, y, z in positions
])
cloud.addReferencePoint("join", *positions[positions[:, 2].argmax()].tolist())
world = pyseb.World("AboutWorkflow")
graph = world.Add(cloud, "cloud")
world.Link("SolidCylinder", "cylinder.bottom", "cloud.join")
world.Link("SolidSphere", "sphere.center", "cylinder.top")
world.Add(graph, "assembly")
expression = world.FormFactor("assembly")
parameters = {
    "beta_cylinder": 3.0, "R_cylinder": 1.0, "L_cylinder": 4.0,
    "beta_sphere": 4.0, "R_sphere": 1.5,
}
q = np.geomspace(0.02, 8.0, 220)
curves = []
for radius in (0.8, 1.5, 2.5):
    parameters["R_sphere"] = radius
    values = np.asarray(world.EvaluateFormFactor("assembly", parameters, q.tolist()))
    assert np.isfinite(values).all() and (values > 0).all()
    curves.append(values)

svg = ['''<svg xmlns="http://www.w3.org/2000/svg" width="1440" height="1010" viewBox="0 0 1440 1010" role="img" aria-labelledby="title desc">
<title id="title">Connect subunits in Python, then choose expressions or direct evaluation</title>
<desc id="desc">Point clouds and analytical shapes snap together at named reference points using World.Add and World.Link. The assembled model branches to a symbolic expression or directly to numerical evaluation. Evaluation produces a scattering plot. Changing the sphere radius gives three curves for the same connected cloud, cylinder, and sphere.</desc>
<defs>
 <marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="8" markerHeight="8" orient="auto-start-reverse"><path d="M0 0L10 5L0 10Z" fill="#64748b"/></marker>
 <linearGradient id="sphere"><stop stop-color="#d1f1ec"/><stop offset="1" stop-color="#35a99b"/></linearGradient>
 <linearGradient id="rod" x2="0" y2="1"><stop stop-color="#dbe7fb"/><stop offset="1" stop-color="#719be1"/></linearGradient>
 <style>
 text {font-family:DejaVu Sans,Arial,sans-serif;fill:#192b40}
 .heading {font-size:25px;font-weight:700}
 .body {font-size:20px}
 .small {font-size:18px;fill:#526477}
 .code {font-family:DejaVu Sans Mono,monospace;font-size:18px}
 .route {fill:none;stroke:#64748b;stroke-width:2.5;marker-end:url(#arrow)}
 </style>
</defs>
<rect width="1440" height="1010" rx="22" fill="#f4f7fb"/>
''']
def text(x, y, value, cls="body", **attrs):
    extra = " ".join(f'{k.replace("_", "-")}="{v}"' for k, v in attrs.items())
    svg.append(f'<text x="{x}" y="{y}" class="{cls}" {extra}>{escape(value)}</text>')
def rect(x, y, w, h, fill="white", stroke="#dce4ee", rx=18):
    svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}"/>')
def path(d, cls="route"):
    svg.append(f'<path d="{d}" class="{cls}"/>')
def circle(x, y, r, fill, stroke="none"):
    svg.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="3"/>')
def link_symbol(x, y):
    svg.append(f'<g transform="translate({x} {y})" fill="none" '
               'stroke="#526477" stroke-width="2.5" stroke-linecap="round" '
               'stroke-linejoin="round" aria-label="linked to">'
               '<path d="M12 3 H7 a5 5 0 0 0 0 10 H10 '
               'M10 13 H15 a5 5 0 0 0 0 -10 H14 M8 8 H14"/></g>')
def cloud_shape(cx, cy, scale=1):
    randomizer = random.Random(12)
    for _ in range(65):
        x, y = randomizer.gauss(0, 25), randomizer.gauss(0, 21)
        if x*x/1800 + y*y/1200 < 1:
            circle(round(cx+x*scale, 1), round(cy+y*scale, 1), 3.1*scale, "#be8051")
def cylinder(cx, cy, width=85, height=40):
    rect(cx-width/2, cy-height/2, width, height, "url(#rod)", "#6084bc", 3)
    svg.append(f'<ellipse cx="{cx-width/2}" cy="{cy}" rx="8" ry="{height/2}" fill="#b9cef0" stroke="#6084bc"/>')
    svg.append(f'<ellipse cx="{cx+width/2}" cy="{cy}" rx="8" ry="{height/2}" fill="#dbe7fb" stroke="#6084bc"/>')
def code(x, y, lines):
    for i, line in enumerate(lines):
        text(x, y+29*i, line, "code")

def latex_equation(x, y, width, height):
    # A compact schematic of the normalized composite form factor. Render
    # mathematical glyphs as vector paths so no browser math engine is needed.
    formula = (r'$F(q)=\frac{\sum_i\beta_i^2F_i(q)+'
               r'2\sum_{i<j}\beta_i\beta_j C_{ij}(q)}'
               r'{(\sum_i\beta_i)^2}$')
    buffer = BytesIO()
    math_to_image(formula, buffer, format="svg", color="#192b40")
    equation = ET.fromstring(buffer.getvalue())
    equation.set("x", str(x))
    equation.set("y", str(y))
    equation.set("width", str(width))
    equation.set("height", str(height))
    svg.append(ET.tostring(equation, encoding="unicode"))

text(30, 42, "Build once. Explore many scattering curves.", "heading")
rect(24, 68, 308, 375)
text(46, 110, "1  Choose subunits", "heading")
text(46, 148, "Point clouds", "body")
cloud_shape(112, 204)
text(177, 197, "Atoms or", "small")
text(177, 223, "scatterers", "small")
text(46, 284, "Analytical models", "body")
cylinder(104, 333)
circle(234, 333, 31, "url(#sphere)", "#31988f")
text(46, 410, "Mix both in one structure.", "small")
path("M340 260 H377")
rect(390, 68, 1026, 375)
text(414, 110, "2  Click together at reference points", "heading")
text(414, 147, "Connect named attachment points in Python.", "body")
# A schematic of exactly the topology used for the plotted curves.
cloud_shape(458, 258, 0.88)
cylinder(548, 258, 94, 43)
circle(634, 258, 38, "url(#sphere)", "#31988f")
for x in (501, 595):
    circle(x, 258, 6, "white", "#da9445")
text(423, 331, "cloud", "small")
link_symbol(484, 317)
text(516, 331, "cylinder", "small")
link_symbol(597, 317)
text(629, 331, "sphere", "small")
text(423, 370, "○  connection / reference point", "small")
rect(716, 174, 676, 238, "#f3f6fa", "#e4eaf2", 10)
code(734, 208, [
    "world = pyseb.World()",
    'graph = world.Add(cloud, "cloud")',
    'world.Link("SolidCylinder",',
    '    "cylinder.bottom", "cloud.join")',
    'world.Link("SolidSphere",',
    '    "sphere.center", "cylinder.top")',
    'world.Add(graph, "assembly")',
])
# Explicit alternatives: the direct route bypasses expression construction.
path("M903 445 V475 H224 V535")
path("M903 475 H668 V535")
rect(24, 548, 400, 350)
rect(446, 548, 360, 350)
rect(828, 548, 588, 350)
text(46, 590, "3a  Get an expression", "heading")
rect(40, 612, 368, 78, "#f3f6fa", "#e4eaf2", 10)
code(46, 638, ['expr = world.FormFactor(', '    "assembly")'])
text(46, 726, "LaTeX output (schematic)", "small")
latex_equation(46, 745, 352, 86)
text(46, 867, "Cloud terms remain numerical.", "small")
text(468, 590, "3b  Evaluate directly", "heading")
rect(462, 612, 328, 107, "#f3f6fa", "#e4eaf2", 10)
code(468, 638, ['F = world.EvaluateFormFactor(', '    "assembly", params,', '    q.tolist())'])
text(468, 763, "Supply q, sizes, and contrasts.", "small")
text(468, 848, "No symbolic step required.", "small")
path("M426 721 H444")
path("M808 721 H826")
text(850, 590, "4  Plot & vary parameters", "heading")
colors = ("#2878b5", "#d27a32", "#279487")
for i, (radius, color) in enumerate(zip((0.8, 1.5, 2.5), colors)):
    x = 866+i*174
    svg.append(f'<path d="M{x} 620 h20" stroke="{color}" stroke-width="3"/>')
    text(x+26, 626, f"R = {radius}", "small")
# Real pySEB curves on logarithmic axes.
x0, y0, pw, ph = 906, 650, 472, 199
for exponent in (0, -2, -4):
    y = y0-exponent/4*ph
    svg.append(f'<path d="M{x0} {y} h{pw}" stroke="#e1e8ef"/>')
    text(x0-12, y+6, {0:"1", -2:"10⁻²", -4:"10⁻⁴"}[exponent], "small", text_anchor="end")
for value, label in ((0.02,"0.02"), (0.1,"0.1"), (1,"1"), (8,"8")):
    x = x0+math.log(value/0.02)/math.log(8/0.02)*pw
    text(x, y0+ph+25, label, "small", text_anchor="middle")
svg.append(f'<path d="M{x0} {y0} V{y0+ph} H{x0+pw}" fill="none" stroke="#8092a5" stroke-width="1.5"/>')
text(862, 684, "F(q)", "small", transform="rotate(-90 862 684)")
text(1142, 887, "q (nm⁻¹)", "small", text_anchor="middle")
for values, color in zip(curves, colors):
    coords = []
    for qi, fi in zip(q, values):
        x = x0+math.log(qi/0.02)/math.log(8/0.02)*pw
        y = y0-math.log10(max(fi, 1e-4))/4*ph
        coords.append(f"{x:.2f},{y:.2f}")
    svg.append(f'<polyline points="{" ".join(coords)}" fill="none" stroke="{color}" stroke-width="2.7"/>')
svg.append('</svg>')
Path(__file__).with_name("about-workflow-draft.svg").write_text("\n".join(svg)+"\n")
