"""Time prepared sasmodels kernels while changing fit parameters on each call."""
import argparse
import json
import pathlib
import sys
import time
import cProfile
import importlib.metadata
import platform

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('repo',type=pathlib.Path)
parser.add_argument('output',type=pathlib.Path)
parser.add_argument('--profile',action='store_true')
parser.add_argument('--fit',action='store_true')
args=parser.parse_args()
args.output.parent.mkdir(parents=True,exist_ok=True)
repo=args.repo.resolve()
sys.path[:0]=[str(repo/'build-python/pyseb'),str(repo)]
import numpy as np
import pyseb
assert pathlib.Path(pyseb._pyseb.__file__).resolve().is_relative_to(repo), pyseb._pyseb.__file__
from sasmodels.core import build_model
from sasmodels.direct_model import call_kernel
from sasmodels.resolution import Pinhole1D
from pyseb.sasmodels import make_model_info, with_hardsphere
from Examples.sasmodels_aot import make_model as aot
from Examples.sasmodels_hybrid import make_model as hybrid

# Kernel and resolution construction are deliberately outside timed fitting calls.
cases=[]
for name,spec,interaction,distributed in [
    ('hybrid',hybrid(),False,False),('hybrid_pd21',hybrid(),False,True),
    ('aot',aot(),True,False),('aot_pd11',aot(),True,True)]:
    q=np.geomspace(.003,.4,101)
    resolution=Pinhole1D(q,.05*q)
    model=build_model(with_hardsphere(spec) if interaction else make_model_info(spec))
    kernel=model.make_kernel([resolution.q_calc])
    parameters={v.name:v.default for v in spec.parameters}
    parameters.update(background=0.)
    if interaction:
        parameters.update(volfraction=.05,structure_factor_mode=1,radius_effective_mode=1)
    if distributed:
        parameters.update(radius_pd=.12,radius_pd_n=11 if interaction else 21)
    def evaluate(index,kernel=kernel,parameters=parameters,resolution=resolution):
        # Change radius on every call, as an optimizer does; no repeated-value cache.
        p=dict(parameters,radius=parameters['radius']*(1.+.01*(index%7-3)))
        return resolution.apply(call_kernel(kernel,p))
    cases.append((name,model,kernel,evaluate,int(resolution.q_calc.size)))

report=dict(extension=pyseb._pyseb.__file__,python=platform.python_version(),
            versions={name:importlib.metadata.version(name) for name in ('numpy','scipy','sasmodels')},
            cases={})
try:
    for name,model,kernel,evaluate,q_count in cases:
        references=[evaluate(i).tolist() for i in range(7)]
        samples=[]
        for repeat in range(3):
            started=time.perf_counter();count=0
            while time.perf_counter()-started < .5:
                evaluate(count);count+=1
            samples.append((time.perf_counter()-started)/count)
        report['cases'][name]=dict(q_count=q_count,
                                  seconds=samples,curves=references)
        print(name,samples,flush=True)
    if args.fit:
        from scipy.optimize import least_squares
        # Reuse a prepared interacting/disperse kernel throughout the fit.
        name,model,kernel,evaluate,q_count=cases[-1]
        q=np.geomspace(.003,.4,101)
        resolution=Pinhole1D(q,.05*q)
        p={v.name:v.default for v in aot().parameters}
        p.update(background=0.,volfraction=.05,structure_factor_mode=1,
                 radius_effective_mode=1,radius_pd=.12,radius_pd_n=11)
        def predict(x):
            return resolution.apply(call_kernel(kernel,dict(p,radius=x[0],thickness=x[1])))
        truth=np.array([29.,6.5])
        measured=predict(truth)
        sigma=.02*measured+.002*measured.max()
        fits=[]
        for repeat in range(3):
            started=time.perf_counter()
            fit=least_squares(lambda x:(predict(x)-measured)/sigma,[25.,8.],
                              bounds=([10.,1.],[50.,20.]),x_scale='jac')
            fits.append(dict(seconds=time.perf_counter()-started,x=fit.x.tolist(),
                             nfev=fit.nfev,cost=fit.cost,success=fit.success))
            if not fit.success:
                raise RuntimeError(fit.message)
            np.testing.assert_allclose(fit.x,truth,rtol=1e-6)
        report['fits']=fits
        print('fits',fits,flush=True)
    if args.profile:
        profiler=cProfile.Profile()
        profiler.enable()
        for name,model,kernel,evaluate,q_count in cases:
            for index in range(50):evaluate(index)
        profiler.disable()
        profiler.dump_stats(str(args.output.with_suffix('.prof')))
finally:
    for name,model,kernel,evaluate,q_count in cases:kernel.release();model.release()
args.output.write_text(json.dumps(report,indent=2)+'\n')
