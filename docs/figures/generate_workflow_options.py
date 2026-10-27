"""Ten editable figure concepts. Run with venv/bin/python from the repo root.

Produces SVGs, PNG previews, a contact sheet and a standalone comparison gallery.
All plotted curves are pySEB evaluations of the depicted mixed assembly.
"""
from pathlib import Path
import os
import sys
from html import escape

os.environ.setdefault('MPLCONFIGDIR', '/tmp/pyseb-figure-matplotlib')
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Ellipse, FancyArrowPatch
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import pyseb

OUT = Path(__file__).with_name('workflow-options')
OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'svg.fonttype': 'none'})
rng = np.random.default_rng(42)
positions = rng.uniform(-0.3, 0.3, size=(20, 3))
cloud = pyseb.DebyeSphereCloud([
    pyseb.SphereScatterer(float(x), float(y), float(z), 0., .1)
    for x, y, z in positions])
cloud.addReferencePoint('join', *positions[positions[:, 2].argmax()].tolist())
world = pyseb.World('WorkflowOptions')
graph = world.Add(cloud, 'cloud')
world.Link('SolidCylinder', 'cylinder.bottom', 'cloud.join')
world.Link('SolidSphere', 'sphere.center', 'cylinder.top')
world.Add(graph, 'assembly')
expr = world.FormFactor('assembly')
params = dict(beta_cylinder=3., R_cylinder=1., L_cylinder=4., beta_sphere=4., R_sphere=1.5)
q = np.geomspace(.02, 8., 240)
radii = (.8, 1.5, 2.5)
curves = [np.asarray(world.EvaluateFormFactor('assembly', dict(params, R_sphere=r), q.tolist())) for r in radii]
assert all(np.isfinite(f).all() and (f > 0).all() for f in curves)
FORMULA = r'$F(q)=\frac{\sum_i\beta_i^2 F_i(q)+2\sum_{i<j}\beta_i\beta_j C_{ij}(q)}{(\sum_i\beta_i)^2}$'
CODE = 'world = pyseb.World()\ngraph = world.Add(cloud, "cloud")\nworld.Link("SolidCylinder",\n    "cylinder.bottom", "cloud.join")\nworld.Link("SolidSphere",\n    "sphere.center", "cylinder.top")\nworld.Add(graph, "assembly")'
PALETTES = {
 'light': dict(bg='#f5f7fa', card='#ffffff', ink='#172d42', muted='#5b7083', line='#dbe4ed', accent='#287ba0', cloud='#bd794e', rod='#6492d4', sphere='#269b8c'),
 'paper': dict(bg='#fffdf8', card='#fffdf8', ink='#242b32', muted='#68717a', line='#d7d5ce', accent='#306b85', cloud='#b07b52', rod='#668bb3', sphere='#3f8d83'),
 'dark': dict(bg='#0c1524', card='#142237', ink='#eef5fe', muted='#a1b5cd', line='#2c415a', accent='#82c8ff', cloud='#e4a371', rod='#8eaeff', sphere='#55d6b8'),
 'blue': dict(bg='#102b46', card='#153753', ink='#f1f8ff', muted='#afcce1', line='#3a5e7b', accent='#83d9ff', cloud='#ffca93', rod='#8cbfff', sphere='#78e1ca'),
}

class Design:
 def __init__(self, number, title, subtitle, palette='light'):
    self.n, self.title, self.subtitle = number, title, subtitle
    self.p = PALETTES[palette]
    self.fig = plt.figure(figsize=(16, 11.2), facecolor=self.p['bg'])
    self.ax = self.fig.add_axes([0, 0, 1, 1], xlim=(0,100), ylim=(70,0))
    self.ax.set_axis_off()
    self.text(4, 4, f'{number:02d}  /  pySEB', 11, color=self.p['accent'], weight='bold')
    self.text(4, 8.2, title, 27, weight='bold')
    self.text(4, 12, subtitle, 12, color=self.p['muted'])
 def text(self,x,y,s,size=13,color=None,**kw):
    kw.setdefault('va', 'center')
    return self.ax.text(x,y,s,fontsize=size,color=color or self.p['ink'],**kw)
 def box(self,x,y,w,h,fill=None,edge=None,r=1.4):
    self.ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle=f'round,pad=0,rounding_size={r}',facecolor=fill or self.p['card'],edgecolor=edge or self.p['line'],linewidth=1))
 def arrow(self,a,b,label=None,rad=0):
    self.ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=15,connectionstyle=f'arc3,rad={rad}',color=self.p['muted'],linewidth=1.5,shrinkA=3,shrinkB=3))
    if label:
      self.text((a[0]+b[0])/2,(a[1]+b[1])/2-1.3,label,10,color=self.p['muted'],ha='center',bbox=dict(facecolor=self.p['bg'],edgecolor='none',pad=2))
 def label(self,x,y,s):
    self.text(x,y,s,11,color=self.p['accent'],weight='bold')
 def assembly(self,x,y,w=28,separate=False,labels=True):
    cy=y
    pts=np.random.default_rng(13).normal(size=(70,2))
    pts=pts[(pts**2).sum(axis=1)<4]
    centers=(x+w*.14,x+w*.51,x+w*.83)
    cx,rx,sx=centers
    self.ax.scatter(cx+pts[:,0]*w*.045,cy+pts[:,1]*w*.04,s=w*1.2,color=self.p['cloud'],alpha=.82,zorder=4,linewidth=0)
    rodw=w*.23 if separate else w*.36
    self.box(rx-rodw/2,cy-w*.055,rodw,w*.11,self.p['rod'],self.p['rod'],.35)
    self.ax.add_patch(Ellipse((rx-rodw/2,cy),w*.025,w*.11,facecolor=self.p['rod'],edgecolor=self.p['ink'],lw=.4))
    self.ax.add_patch(Circle((sx,cy),w*.11,facecolor=self.p['sphere'],edgecolor=self.p['sphere'],lw=1))
    self.ax.add_patch(Circle((sx-w*.035,cy-w*.035),w*.033,facecolor='white',alpha=.18,lw=0))
    if separate:
      for left,right in ((cx+w*.12,rx-rodw/2-.3),(rx+rodw/2+.3,sx-w*.12)):
        self.arrow((left,cy),(right,cy))
    else:
      for a,b in ((cx+w*.08,rx-rodw/2),(rx+rodw/2,sx)):
        self.ax.plot([a,b],[cy,cy],color=self.p['muted'],lw=1.5,zorder=2)
      for px in (rx-rodw/2,rx+rodw/2):
        self.ax.add_patch(Circle((px,cy),w*.017,facecolor=self.p['card'],edgecolor='#e6a344',lw=1.8,zorder=6))
    if labels:
      for px,name in zip(centers,('Cloud','Cylinder','Sphere')):
        self.text(px,cy+w*.19,name,10,color=self.p['muted'],ha='center')
 def library(self,x,y,w,h,title='Choose building blocks'):
    self.box(x,y,w,h)
    self.label(x+2,y+2.7,title.upper())
    self.assembly(x+w*.06,y+h*.46,w*.88,True)
 def python(self,x,y,w,h,compact=False):
    self.box(x,y,w,h)
    self.label(x+2,y+2.7,'CONNECT IN PYTHON')
    snippet=CODE if not compact else 'graph = world.Add(cloud, "cloud")\nworld.Link("SolidCylinder",\n    "cylinder.bottom", "cloud.join")\n# Attach the sphere, name the assembly.'
    self.text(x+2,y+5.4,snippet,11 if w<43 else 12,fontfamily='DejaVu Sans Mono',va='top',linespacing=1.55)
 def equation(self,x,y,w,h,full=True):
    self.box(x,y,w,h)
    self.label(x+2,y+2.7,'EXPRESSION → LATEX')
    if h >= 14:
      self.text(x+2,y+6,'expr = world.FormFactor("assembly")',9 if w<36 else 10,fontfamily='DejaVu Sans Mono')
    self.text(x+w/2,y+h*.59,FORMULA if full else r'$F(q;\boldsymbol{\theta})$',19 if w>55 else 15 if full else 25,ha='center')
    if h >= 14:
      self.text(x+2,y+h-2.2,'Schematic; cloud terms numerical.',8 if w<36 else 9,color=self.p['muted'])
 def evaluate(self,x,y,w,h):
    self.box(x,y,w,h)
    self.label(x+2,y+2.7,'DIRECT NUMERICAL EVALUATION')
    self.text(x+2,y+6,'F = world.EvaluateFormFactor(\n    "assembly", params, q.tolist())',11,fontfamily='DejaVu Sans Mono',va='top',linespacing=1.6)
    if h>17:self.text(x+2,y+h-2.4,'Supply sizes, contrasts, and a q grid.',10,color=self.p['muted'])
 def plot(self,x,y,w,h,controls=True):
    self.box(x,y,w,h)
    self.label(x+2,y+2.7,'CHANGE PARAMETERS · PLOT AGAIN')
    left=x+5; top=y+6.3; pw=w-7; ph=h-12.5-(3 if controls else 0)
    ax=self.fig.add_axes([left/100,1-(top+ph)/70,pw/100,ph/70],facecolor=self.p['card'])
    for f,c,r in zip(curves,(self.p['rod'],self.p['cloud'],self.p['sphere']),radii):
      ax.loglog(q,f,color=c,lw=2,label=f'{r:g}')
    ax.set_xlim(.02,8); ax.set_ylim(.0007,1.2)
    ax.set_xlabel(r'$q$ (nm$^{-1}$)',fontsize=10,color=self.p['muted'],labelpad=2)
    ax.set_ylabel('$F(q)$',fontsize=10,color=self.p['muted'],labelpad=2)
    ax.tick_params(which='both',colors=self.p['muted'],labelsize=8)
    for s in ax.spines.values():s.set_color(self.p['line'])
    ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    ax.grid(which='major',color=self.p['line'],lw=.5)
    ax.legend(title='Sphere radius R (nm)',fontsize=8,title_fontsize=8,frameon=False,labelcolor=self.p['muted'],loc='lower left')
    ax.get_legend().get_title().set_color(self.p['muted'])
    if controls:
      yy=y+h-2.2
      self.text(x+2,yy,'R',10,color=self.p['muted'])
      self.ax.plot([x+6,x+w*.53],[yy,yy],color=self.p['line'],lw=3,solid_capstyle='round')
      self.ax.add_patch(Circle((x+w*.3,yy),.65,fc=self.p['accent'],lw=0))
      self.text(x+w*.58,yy,'vary → re-evaluate',9,color=self.p['muted'])
 def footer(self,s='One topology · expression or direct evaluation · many parameter choices'):
    self.text(4,67.5,s,10,color=self.p['muted'])
 def save(self,slug):
    base=OUT/f'{self.n:02d}-{slug}'
    self.fig.savefig(base.with_suffix('.svg'),facecolor=self.p['bg'])
    self.fig.savefig(base.with_suffix('.png'),dpi=100,facecolor=self.p['bg'])
    plt.close(self.fig)
    return dict(n=self.n,title=self.title,subtitle=self.subtitle,stem=base.name)

items=[]
# 01: familiar horizontal process with genuinely separate computational routes.
d=Design(1,'From building blocks to scattering','A clean overview with a clear choice between symbolic and direct calculation.')
d.library(4,17,27,19);d.python(35,17,61,23)
d.arrow((31,26),(35,26));d.assembly(6,46,23);d.label(6,40,'ONE CONNECTED ASSEMBLY')
d.equation(35,44,29,18,False);d.plot(69,44,27,20,False)
d.arrow((65,40),(65,42));d.text(68,41.8,'evaluate directly',10,color=d.p['muted'])
d.arrow((66,42),(80,44));d.arrow((65,42),(50,44));d.arrow((64,53),(69,53),'evaluate')
d.footer();items.append(d.save('clean-overview'))
# 02: spacious, publication-minded two-row figure.
d=Design(2,'A structure becomes an equation—and a curve','A restrained scientific figure with a prominent LaTeX result.','paper')
d.library(4,17,27,18);d.python(35,17,61,23);d.arrow((31,26),(35,26))
d.equation(4,44,56,19);d.plot(65,44,31,19,False)
d.arrow((52,40),(32,44),'symbolic');d.arrow((85,40),(85,44),'direct')
d.arrow((60,53),(65,53));d.footer('Analytical models and point clouds share the same structure graph.');items.append(d.save('publication'))
# 03: central model is the hero, outputs branch to right.
d=Design(3,'The assembly is the model','A large central structure, with Python below and two outputs to the right.')
d.library(4,17,28,19);d.box(36,17,29,22);d.label(38,20,'CLICK AT REFERENCE POINTS');d.assembly(38,29,25)
d.arrow((32,26),(36,26));d.python(4,41,61,22)
d.equation(70,17,26,19,False);d.plot(70,40,26,23,False)
d.arrow((65,25),(70,25));d.arrow((65,34),(70,47),'direct');d.arrow((83,36),(83,40),'evaluate')
d.footer();items.append(d.save('assembly-centre'))
# 04: dark computing aesthetic, large plot.
d=Design(4,'Compose. Calculate. Explore.','A dark scientific dashboard with a visible parameter sweep.','dark')
d.library(4,17,27,19);d.python(35,17,61,23);d.arrow((31,26),(35,26))
d.equation(4,43,40,20);d.plot(49,43,47,20,False)
d.arrow((52,40),(24,43),'expression');d.arrow((82,40),(82,43),'direct');d.arrow((44,53),(49,53))
d.footer('Shown: R = 0.8, 1.5, 2.5 nm. Same cloud–cylinder–sphere connections.');items.append(d.save('dark-laboratory'))
# 05: vertical story with aligned sequence numbers.
d=Design(5,'One model, four steps','A vertical reading order for an explanatory About section.')
d.library(10,17,34,18);d.python(48,17,48,23);d.arrow((44,26),(48,26))
d.equation(10,44,26,19,False);d.evaluate(39,44,27,19);d.plot(69,44,27,19,False)
d.arrow((72,40),(53,44),'direct');d.arrow((51,40),(24,44),'expression');d.arrow((36,53),(39,53));d.arrow((66,53),(69,53))
# This variant puts a small preview between construction and outputs.
d.assembly(10,39,24,labels=False);d.text(4,25,'01',14,color=d.p['accent'],weight='bold');d.text(4,53,'02',14,color=d.p['accent'],weight='bold')
d.footer('Choose and connect above; calculate and plot below. Repeat with new parameters.');
items.append(d.save('vertical-explainer'))
# 06: notebook style with real code and equation output.
d=Design(6,'A Python notebook, made visual','Inputs on the left. Equations and curves on the right.','paper')
d.library(4,17,39,18);d.python(4,39,47,24)
d.equation(56,17,40,19);d.plot(56,40,40,23)
d.arrow((51,43),(56,26),'expression');d.arrow((51,51),(56,51),'direct')
d.footer('Choose blocks → connect in Python → inspect an expression or plot directly.');items.append(d.save('notebook'))
# 07: generous plot dominates, compact inputs down the left.
d=Design(7,'See what your parameters change','The scattering curves take centre stage; the construction stays concise.')
d.library(4,17,37,17);d.python(4,37,37,20,True)
d.plot(46,17,50,33);d.equation(46,53,50,11,False)
d.arrow((41,43),(46,35),'direct');d.arrow((41,53),(46,58),'expression');d.arrow((74,53),(74,50),'evaluate')
d.footer('Radius controls illustrate re-evaluation; the exported figure is static.');items.append(d.save('plot-first'))
# 08: precise blueprint style with linked nodes and technical panels.
d=Design(8,'A connected graph. Two calculation routes.','A technical blueprint that makes the named attachment points explicit.','blue')
for x in range(4,98,4):d.ax.plot([x,x],[16,65],color=d.p['line'],lw=.3,alpha=.35,zorder=0)
for y in range(17,66,4):d.ax.plot([4,96],[y,y],color=d.p['line'],lw=.3,alpha=.35,zorder=0)
d.box(4,17,40,20);d.label(6,20,'cloud.join ↔ cylinder.bottom');d.assembly(7,27,32)
d.python(49,17,47,23);d.arrow((44,27),(49,27))
d.equation(4,44,40,19);d.plot(49,44,47,19,False)
d.arrow((60,40),(28,44),'expression');d.arrow((84,40),(84,44),'direct');d.arrow((44,53),(49,53))
d.footer('Named references connect geometry; the calculation includes interference.');items.append(d.save('blueprint'))
# 09: split palette makes construction and analysis visually distinct.
d=Design(9,'Build on the left. Explore on the right.','Two clear halves connected by one reusable assembly.')
d.box(3,16,44,49,fill='#eaf2f8',edge='#eaf2f8');d.box(50,16,47,49)
d.library(5,18,40,17);d.python(5,39,40,23)
d.equation(52,18,43,16);d.plot(52,38,43,25)
d.arrow((45,43),(52,26),'expression');d.arrow((45,51),(52,51),'direct');d.arrow((76,34),(76,38))
d.footer();items.append(d.save('split-story'))
# 10: image-led structure, code as compact subtitle, large equation and plot.
d=Design(10,'Click the pieces together','An illustration-led opening with a simple mathematical and numerical finish.','paper')
d.label(6,18,'CLOUDS + ANALYTICAL SUBUNITS');d.assembly(5,25,34,True)
d.arrow((41,25),(48,25));d.assembly(50,25,39,labels=True)
d.text(51,18,'ONE CONNECTED ASSEMBLY',11,color=d.p['accent'],weight='bold')
d.box(4,35,92,7);d.text(6,38.5,'Python:  world.Add(cloud, "cloud")    →    world.Link("SolidCylinder", ...)',12,fontfamily='DejaVu Sans Mono')
d.equation(4,45,43,19);d.plot(54,45,42,19,False)
d.arrow((49,42),(26,45),'expression');d.arrow((80,42),(80,45),'direct');d.arrow((47,55),(54,55),'evaluate')
d.footer('Connection markers indicate named reference points; curves are real pySEB evaluations.');items.append(d.save('illustrated-story'))

# Standalone gallery; PNGs work without extra libraries or a server.
styles='''body{margin:0;background:#eef2f6;color:#172d42;font:16px system-ui,sans-serif}header,main{max-width:1500px;margin:auto;padding:32px}header p{max-width:850px;line-height:1.6}nav{display:flex;flex-wrap:wrap;gap:8px}a{color:#226e91}nav a,.links a{display:inline-block;padding:8px 12px;border:1px solid #cbd7e2;border-radius:8px;text-decoration:none;background:white}main{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:24px;padding-top:0}article{background:white;padding:18px;border-radius:16px;box-shadow:0 5px 25px #19364a0c}h2{font-size:20px;margin:0 0 8px}article p{line-height:1.5;color:#526477;min-height:48px}img{width:100%;height:auto;border-radius:8px}.links{display:flex;gap:8px;margin-top:12px}@media(max-width:850px){main{grid-template-columns:1fr}header,main{padding:18px}}'''
html=['<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>pySEB · 10 figure directions</title><style>'+styles+'</style><header><h1>10 directions for the About figure</h1><p>Compare the numbered layouts, then choose a number or combine elements. The equation is schematic; plots use the same real cloud–cylinder–sphere calculation with sphere radii 0.8, 1.5 and 2.5 nm. Slider drawings illustrate the workflow; these are static figures.</p><nav>']
html.extend(f'<a href="#option-{v["n"]}">{v["n"]:02d}</a>' for v in items)
html.append('</nav><p><a href="contact-sheet.jpg">Open the contact sheet</a> · <a href="workflow-options.zip">Download all 10 figures</a></p></header><main>')
for v in items:
 html.append(f'<article id="option-{v["n"]}"><h2>{v["n"]:02d} · {escape(v["title"])}</h2><p>{escape(v["subtitle"])}</p><a href="{v["stem"]}.png"><img loading="lazy" src="{v["stem"]}.png" alt="Option {v["n"]}: {escape(v["title"])}"></a><div class="links"><a href="{v["stem"]}.png">Full size PNG</a><a href="{v["stem"]}.svg">Vector SVG</a></div></article>')
html.append('</main></html>')
(OUT/'index.html').write_text('\n'.join(html))
# Contact sheet with legible headings; full figures are available in the gallery.
thumbw, thumbh=640,448
sheet=Image.new('RGB',(2*thumbw+72,5*(thumbh+48)+48),'#e9eef3')
draw=ImageDraw.Draw(sheet)
font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf',19)
for i,v in enumerate(items):
 x=24+(i%2)*(thumbw+24);y=24+(i//2)*(thumbh+48)
 draw.text((x,y),f'{v["n"]:02d}  {v["title"]}',fill='#172d42',font=font)
 with Image.open(OUT/(v['stem']+'.png')) as im:
  im.thumbnail((thumbw,thumbh),Image.Resampling.LANCZOS)
  sheet.paste(im,(x,y+30))
sheet.save(OUT/'contact-sheet.jpg',quality=92)
from zipfile import ZipFile, ZIP_DEFLATED
with ZipFile(OUT/'workflow-options.zip','w',ZIP_DEFLATED) as z:
 for f in sorted(OUT.iterdir()):
  if f.suffix in ('.png','.svg','.html','.jpg'):z.write(f,f.name)
 z.write(Path(__file__), 'generate_workflow_options.py')
print(f'Created {len(items)} figure variants in {OUT}')
