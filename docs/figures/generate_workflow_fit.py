"""Generate a small, explicitly simulated-data fitting example for the full figure.

Run with a Python environment containing pySEB, NumPy, SciPy and Matplotlib.
Fits the sphere radius in the same cloud-cylinder-sphere assembly as the overview.
"""
from pathlib import Path
import os
import sys
os.environ.setdefault('MPLCONFIGDIR','/tmp/pyseb-figure-matplotlib')
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares
import pyseb

rng=np.random.default_rng(42)
positions=rng.uniform(-.3,.3,(20,3))
cloud=pyseb.DebyeSphereCloud([pyseb.SphereScatterer(float(x),float(y),float(z),0.,.1) for x,y,z in positions])
cloud.addReferencePoint('join',*positions[positions[:,2].argmax()].tolist())
w=pyseb.World('FitIllustration');g=w.Add(cloud,'cloud')
w.Link('SolidCylinder','cylinder.bottom','cloud.join')
w.Link('SolidSphere','sphere.center','cylinder.top');w.Add(g,'assembly')
params=dict(beta_cylinder=3.,R_cylinder=1.,L_cylinder=4.,beta_sphere=4.,R_sphere=1.5)
q=np.geomspace(.05,3.,28)
def evaluate(radius,grid):
 return np.asarray(w.EvaluateFormFactor('assembly',dict(params,R_sphere=float(radius)),grid.tolist()))
truth=evaluate(1.5,q)
sigma=.045*truth
observed=truth+np.random.default_rng(7).normal(size=q.size)*sigma
residual=lambda x:(evaluate(x[0],q)-observed)/sigma
fit=least_squares(residual,[.85],bounds=([.3],[3.]))
assert fit.success and abs(fit.x[0]-1.5)<.12
assert np.linalg.norm(fit.fun)<np.linalg.norm(residual([.85]))
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(2.5,3.0),layout='constrained',facecolor='#f6f8fc')
ax.set_facecolor('#f6f8fc')
ax.errorbar(q,observed,yerr=sigma,fmt='o',ms=3,color='#596e87',ecolor='#aab8c8',elinewidth=.7,capsize=1.3,label='Data')
fine=np.geomspace(q[0],q[-1],250)
ax.plot(fine,evaluate(fit.x[0],fine),color='#d27a32',lw=2,label='Fit')
ax.set_xscale('log');ax.set_yscale('log')
ax.set_xlabel(r'$q$ (nm$^{-1}$)',fontsize=11,color='#526477')
ax.set_ylabel('$F(q)$',fontsize=11,color='#526477')
ax.set_title('Simulated data + fit',fontsize=11,color='#192b40',pad=10)
ax.text(.04,.06,f'R = {fit.x[0]:.2f} nm',transform=ax.transAxes,fontsize=10,color='#526477')
ax.tick_params(labelsize=9,colors='#526477')
ax.spines[['top','right']].set_visible(False)
for s in ax.spines.values():s.set_color('#dce4ee')
ax.grid(color='#e1e8ef',lw=.7)
ax.legend(loc='upper right',fontsize=9,frameon=False,labelcolor='#526477')
path=Path(__file__).with_name('workflow-fit-example.svg')
fig.savefig(path)
print(f'{path}: fitted R={fit.x[0]:.4f} nm (simulated R=1.5 nm)')
