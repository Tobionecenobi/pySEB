"""Build the expanded workflow from the approved overview's vector artwork.

Run: venv/bin/python docs/figures/generate_about_workflow_full.py
Only about-workflow-full.svg is written; the About-page overview is preserved.
"""
import base64
import math
from copy import deepcopy
from pathlib import Path
import keyword
import random
import re
import xml.etree.ElementTree as E

NS = 'http://www.w3.org/2000/svg'
E.register_namespace('', NS)
E.register_namespace('xlink', 'http://www.w3.org/1999/xlink')
BASE = Path(__file__).with_name('about-workflow.svg')
OUT = BASE.with_name('about-workflow-full.svg')
source = E.parse(BASE).getroot()
ns = {'s': NS}
def tag(name): return '{'+NS+'}'+name
root = E.Element(tag('svg'), dict(width='1440', height='1852', viewBox='0 0 1440 1852', role='img', **{'aria-labelledby':'title desc'}))
def add(name, text=None, parent=root, **attrs):
    node = E.SubElement(parent, tag(name), {k.replace('_','-'):str(v) for k,v in attrs.items()})
    node.text = text
    return node
add('title', 'pySEB: from analytical and PDB subunits to expressions, visualization, and fitted scattering', id='title')
add('desc', 'Install pySEB; choose analytical or custom models, parse PDB files into numerical clouds, and optionally include excluded solvent and hydration. Connect named reference points in Python, nest structures, and save topology. The assembly branches to symbolic expressions and SymPy, direct numerical evaluation, 3D visualization through OpenUSD and Blender, or the SasView ecosystem through the Python sasmodels adapter. Direct evaluation feeds both plotting and parameter fitting. The plotted curves use the same cloud-cylinder-sphere assembly at three sphere radii.', id='desc')
root.append(deepcopy(source.find('s:defs', ns)))
add('rect', width=1440,height=1852,rx=22,fill='#f1f5fa')

def text(x,y,s,cls='small',**attrs): return add('text',s,x=x,y=y,**{'class':cls},**attrs)
def rect(x,y,w,h,fill='#fff',stroke='#dde5ef',rx=18):
    return add('rect',x=x,y=y,width=w,height=h,rx=rx,fill=fill,stroke=stroke)
def card(x,y,w,h,number,title):
    add('rect',x=x,y=y+3,width=w,height=h,rx=18,fill='#1b2f4a',opacity='.06')
    rect(x,y,w,h)
    if number is None:
        text(x+22,y+45,title,'heading')
        return
    badgewidth=44 if len(str(number))>1 else 32
    rect(x+22,y+20,badgewidth,32,'#1b2f4a','#1b2f4a',10)
    text(x+22+badgewidth/2,y+42,str(number),'badge')
    text(x+34+badgewidth,y+45,title,'heading')
def route(d): return add('path',d=d,**{'class':'route'})
def highlight(line, node):
    """Use the approved terminal palette without changing code or spacing."""
    pattern = r'''("[^"\\]*(?:\\.[^"\\]*)*"|'[^'\\]*(?:\\.[^'\\]*)*'|\#[^\n]*|[A-Za-z_]\w*|\d+(?:\.\d+)?|\s+|.)'''
    for match in re.finditer(pattern, line):
        token = match.group()
        remainder = line[match.end():].lstrip()
        cls = None
        if token.startswith(('"', "'")) or token[0].isdigit(): cls = 's'
        elif token.startswith('#'): cls = 'p'
        elif keyword.iskeyword(token) or token == 'pip': cls = 'k'
        elif token in ('pyseb', 'sympy', 'plt', 'np'): cls = 'm'
        elif token.isidentifier() and (remainder.startswith('(') or token == 'install'): cls = 'f'
        elif token in '.=(),[]:+-*/': cls = 'p'
        attrs = {'class': cls} if cls else {}
        add('tspan', token, parent=node, **attrs)

def code(x,y,w,lines,size=16):
    rect(x,y,w,24+29*len(lines),'#16233a','#16233a',12)
    for i,line in enumerate(lines):
        node = text(x+16,y+31+i*29,None,'code',style=f'font-size:{size}px',**{'xml:space':'preserve'})
        highlight(line, node)
def section(y,label,detail=''):
    text(24,y,label,'label')
    if detail: text(450,y,detail,'small')

def relocate_since(start, dx):
    """Move a complete panel, including nested vector artwork, as one group."""
    children=list(root)[start:]
    for child in children: root.remove(child)
    group=add('g',transform=f'translate({dx} 0)')
    for child in children: group.append(child)

def dots(cx,cy,hydration=False):
    rng=random.Random(8)
    if hydration:
        for i in range(70):
            import math
            a=rng.uniform(0,2*math.pi);r=rng.uniform(50,66)
            add('circle',cx=round(cx+r*math.cos(a),2),cy=round(cy+r*math.sin(a),2),r=2.8,fill='#6cb8da',opacity='.65')
    for i in range(65):
        dx,dy=rng.gauss(0,19),rng.gauss(0,17)
        if dx*dx+dy*dy<1600:
            add('circle',cx=round(cx+dx,2),cy=round(cy+dy,2),r=3,fill='#c48a5a')

# Installation is a prerequisite outside the four main sections.
card(24,24,1392,105,None,'Install pySEB')
code(280,43,354,['pip install pyseb'],17.5)
text(670,70,'Optional extras: 3D visualization and SasView tools','small')
text(670,99,'pip install "pyseb[visualization,sasmodels]"','mono')

# 1: alternative sources; hydration is optional, and only on the cloud path.
section(165,'1  Prepare building blocks','Analytical models and numerical clouds can be mixed.')
card(24,190,390,260,'1a','Analytical & custom')
text(46,277,'Spheres, cylinders, polymers, …','body')
add('rect',x=70,y=307,width=92,height=42,fill='#94b5e8',stroke='#5f84c0')
add('rect',x=74,y=312,width=84,height=4,rx=2,fill='#e4edfc',opacity='.65')
add('ellipse',cx=70,cy=328,rx=8,ry=21,fill='#86a9df',stroke='#5f84c0')
add('ellipse',cx=162,cy=328,rx=8,ry=21,fill='#bfd3f3',stroke='#5f84c0')
add('circle',cx=249,cy=328,r=34,fill='url(#sphere)',stroke='#31988f')
text(46,394,'YAML equations + reference points','small')
text(46,424,'Python callbacks for numerical subunits','small')
card(450,190,380,260,'1b','PDB & point clouds')
code(470,260,340,['parser = pyseb.PDBParser()', 'atoms = parser.parse_file(', '    "protein.pdb")'],15)
text(472,398,'Atom profiles → numerical cloud','small')
text(472,426,'Or supply coordinates and weights.','small')
card(866,190,550,260,'1c','Optional solvent & hydration')
dots(954,336,True)
text(1043,278,'Atomic SAXS cloud','label')
text(1043,314,'Displaced-solvent contribution','small')
text(1043,344,'Solvent-accessible hydration shell','small')
text(1043,382,'Set solvent density and','small')
text(1043,408,'hydration contrast.','small')
route('M838 320 H858')
# Each source can reach assembly. The cloud bypass is explicitly optional.
route('M219 452 V470 H515 V504')
route('M640 452 V504')
route('M1141 452 V479 H790 V504')
section(495,'2  Assemble in Python')
text(866,495,'Use clouds with or without hydration.','tiny')

# 2: restore the approved full-size assembly artwork, with current colours.
assembly = next(g for g in source.findall('s:g',ns) if g.get('transform')=='translate(238 0)')
parts=list(assembly)
panel=add('g',transform='translate(-348 418)')
for element in parts:
    item=deepcopy(element)
    if item.tag==tag('ellipse') and item.get('cx')=='686' and item.get('cy')=='329':
        continue  # No drop shadow beneath the assembly sphere.
    if item.tag==tag('rect') and item.get('x')=='372':item.set('width','1044')
    if item.tag==tag('text') and item.get('class')=='badge' and item.text=='3':item.text='2'
    panel.append(item)
codegroup = next(g for g in source.findall('s:g',ns) if g.get('transform')=='translate(218 0)')
panelcode=add('g',transform='translate(-348 418)')
for element in codegroup:
    item=deepcopy(element)
    if item.tag==tag('rect') and item.get('x')=='788':item.set('width','604')
    if item.tag==tag('text') and item.get('class')=='code':item.attrib.pop('style',None)
    if item.tag==tag('text') and ''.join(item.itertext())=='python':item.set('x','1374')
    panelcode.append(item)
card(1100,510,316,390,'2+','Reuse the model')
text(1122,600,'Nest complete assemblies','label')
text(1122,632,'Chains, stars, and hierarchies','small')
text(1122,684,'Save / load topology','label')
code(1120,704,276,['world.save("model.yaml")'],14)
text(1122,799,'Reload with pyseb.load_world.','tiny')
text(1122,835,'Store parameter values','small')
text(1122,861,'separately from topology.','small')
route('M1074 705 H1094')

# A shared bus makes all four outputs peers of the same assembled model.
add('path',d='M542 900 V950 H1254',fill='none',stroke='#7b8da3',stroke_width=2.5)
add('path',d='M542 950 H186',fill='none',stroke='#7b8da3',stroke_width=2.5)
for cx in (186,542,898,1254):route(f'M{cx} 950 V984')
text(24,932,'3  Choose outputs','label')
text(1092,932,'4  Visualize the structure','label')

card(24,990,324,330,'3a','SymPy & LaTeX')
code(44,1058,284,['expr = world.FormFactor(', '    "assembly")', 's = pyseb.to_sympy(expr)', 'sympy.latex(s)'],14)
equation=deepcopy(source.find('s:svg',ns))
equation.set('x','55');equation.set('y','1208');equation.set('width','262');equation.set('height','64')
root.append(equation)
text(46,1298,'Schematic; cloud terms numerical.','tiny')
card(380,990,324,330,'3b','Evaluate directly')
code(400,1058,284,['F = world.EvaluateFormFactor(', '    "assembly", params,', '    q.tolist())'],13)
text(402,1210,'q grid + dimensions + contrasts','small')
text(402,1245,'Normalized F(q) or intensity','small')
text(402,1298,'No symbolic step required.','tiny')
visual_start=len(root)
card(736,990,324,330,'4','3D visualization')
text(758,1067,'OpenUSD export (.usda)','small')
code(756,1084,284,['world.export_usd(', '    "assembly", "model.usda",', '    params, options)'],13.5)
# Official secondary Blender branding, embedded as vectors for a portable SVG.
# Source: https://download.blender.org/branding/blender_logo_socket.svg
# Usage guidance: https://www.blender.org/about/logo/
blender_logo=E.parse(BASE.with_name('blender-logo.svg')).getroot()
blender_logo.set('x','788');blender_logo.set('y','1207')
blender_logo.set('width','220');blender_logo.set('height','67.2')
add('title','Blender: view and render the exported 3D model',parent=blender_logo)
root.append(blender_logo)
text(758,1300,'View and render the model in 3D.','tiny')
relocate_since(visual_start,356)
adapter_start=len(root)
card(1092,990,324,330,'3c','SasView ecosystem')
root[-1].set('style','font-size:21px')
code(1112,1058,284,['spec = ModelSpec(...)', 'info = make_model_info(spec)'],13)
text(1114,1170,'Python adapter via sasmodels','small',style='font-size:16px')
text(1114,1196,"SasView’s scattering engine",'tiny')
text(1114,1230,'Declared size distributions','small')
text(1114,1260,'Interactions and resolution','small')
text(1114,1298,'SciPy / Bumps fitting','small')
relocate_since(adapter_start,-356)

# Calculation feeds compact plot/fit panels; visualization also has a Python path.
route('M350 1157 H374')
route('M542 1322 V1350 H239 V1374')
add('path',d='M542 1350 H898',fill='none',stroke='#7b8da3',stroke_width=2.5,stroke_linejoin='round')
route('M898 1322 V1374')
route('M1254 1322 V1374')

card(24,1380,430,440,'3d','Plot & vary parameters')
root[-1].set('style','font-size:22px')
text(46,1471,'Sphere radius R (nm)','small')
colors=('#2878b5','#d27a32','#279487')
for i,(radius,color) in enumerate(zip((.8,1.5,2.5),colors)):
    x=55+i*128
    add('path',d=f'M{x} 1492 h22',stroke=color,stroke_width=3.5,stroke_linecap='round')
    text(x+30,1498,str(radius),'small')
x0,y0,pw,ph=86,1516,340,224
for exponent,label in ((0,'1'),(-2,'10⁻²'),(-4,'10⁻⁴')):
    y=y0-exponent*ph/4
    add('path',d=f'M{x0} {y} h{pw}',stroke='#e1e8f0',stroke_width=1)
    text(x0-10,y+5,label,'tiny',text_anchor='end')
for value,label in ((.02,'0.02'),(.1,'0.1'),(1.,'1'),(8.,'8')):
    x=x0+math.log(value/.02)/math.log(8/.02)*pw
    text(x,y0+ph+23,label,'tiny',text_anchor='middle')
add('path',d=f'M{x0} {y0} V{y0+ph} H{x0+pw}',fill='none',stroke='#8a9bb0',stroke_width=1.5)
text(34,1628,'F(q)','tiny',transform='rotate(-90 34 1628)',text_anchor='middle')
text(256,1796,'q (nm⁻¹)','tiny',text_anchor='middle')
# Re-map the approved curve coordinates; no scattering data are approximated.
for curve in source.findall('.//s:polyline',ns):
    points=[]
    for pair in curve.get('points').split():
        x,y=map(float,pair.split(','))
        points.append(f'{x0+(x-954)/436*pw:.2f},{y0+(y-686)/168*ph:.2f}')
    add('polyline',points=' '.join(points),fill='none',stroke=curve.get('stroke'),stroke_width=3,stroke_linecap='round',stroke_linejoin='round')

card(490,1380,566,440,'3+','Explore & fit data')
# Three short actions, with recognizable tools rather than implementation detail.
text(512,1471,'Vary parameters','label')
text(512,1502,'Sizes, contrast, hydration.','small',style='font-size:17px')
text(512,1528,'Reuse the same model.','small',style='font-size:17px')
text(512,1570,'Fit measured scattering','label',style='font-size:18px')
# Official SciPy logo: https://github.com/scipy/scipy/blob/main/doc/source/_static/logo.svg
scipy_logo=E.parse(BASE.with_name('scipy-logo.svg')).getroot()
scipy_logo.set('x','512');scipy_logo.set('y','1591')
scipy_logo.set('width','40');scipy_logo.set('height','40')
root.append(scipy_logo)
text(564,1617,'SciPy / Bumps','small',style='font-size:17px')
text(512,1660,'Extend the model','label',style='font-size:18px')
text(512,1690,'Distributions, interactions,','small',style='font-size:16px')
text(512,1714,'and instrument resolution.','small',style='font-size:16px')
# Official website logo: https://www.sasview.org/img/sasview_logo.png
sasview_logo=BASE.with_name('sasview-logo.png').read_bytes()
add('image',x=512,y=1730,width=180,height=52,
    href='data:image/png;base64,'+base64.b64encode(sasview_logo).decode('ascii'))
text(512,1801,'via sasmodels (Python)','tiny')

# Actual SciPy fit of simulated data, kept separate from experimental results.
fit_plot=E.parse(BASE.with_name('workflow-fit-example.svg')).getroot()
fit_plot.set('x','760');fit_plot.set('y','1463')
fit_plot.set('width','274');fit_plot.set('height','328.8')
root.append(fit_plot)

card(1092,1380,324,440,'4+','Python 3D viewer')
code(1112,1448,284,['pyseb.visualize(world,', '    "assembly", params,', '    options)'],16)
# Embed the real browser capture so the SVG remains portable and self-contained.
screenshot=BASE.with_name('workflow-viewer-hover.png').read_bytes()
add('image',x=1112,y=1570,width=284,height=220,
    href='data:image/png;base64,'+base64.b64encode(screenshot).decode('ascii'))
rect(1112,1570,284,220,'none','#e5ebf3',12)
text(1114,1810,'Rotate • zoom • hover','small',style='font-size:16px')

OUT.write_text(E.tostring(root,encoding='unicode')+'\n')
print(f'Wrote {OUT}')
