"""Integration contracts for the optional sasmodels bridge."""
import dataclasses
import importlib.util
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
import unittest

import numpy as np
import pyseb

HAS_SASMODELS = importlib.util.find_spec("sasmodels") is not None


@unittest.skipUnless(HAS_SASMODELS, "install pyseb[sasmodels]")
class TestSasmodels(unittest.TestCase):
    def setUp(self):
        from pyseb.sasmodels import ModelSpec, Parameter
        self.ModelSpec, self.Parameter = ModelSpec, Parameter
        self.q = np.array([0., 1e-5, .01, .1, .25])

    def calculate(self, spec, values=None, backend="runtime", q=None):
        from pyseb.sasmodels import make_model_info
        from sasmodels.core import build_model
        from sasmodels.direct_model import call_kernel
        model = build_model(make_model_info(spec, backend), platform="dll", dtype="double")
        kernel = model.make_kernel([self.q if q is None else q])
        try:
            return call_kernel(kernel, dict(values or {}, background=0.))
        finally:
            kernel.release()
            model.release()

    def simple(self, shape="SolidSphere", connected=False):
        def factory():
            w=pyseb.World("fixture")
            graph=w.Add(shape,"a")
            if connected:
                w.Link(shape,"b.end1","a.end2")
                w.Add(graph,"particle")
            return w
        defaults={"SolidSphere":dict(R_a=3.,beta_a=2.),
                  "SolidCylinder":dict(R_a=2.,L_a=5.,beta_a=-2.),
                  "ThinRod":dict(L_a=5.,beta_a=2.),
                  "GaussianPolymer":dict(Rg_a=3.,beta_a=2.)}[shape]
        if connected:
            defaults.update({k.replace("_a","_b"):v for k,v in list(defaults.items())})
        return self.ModelSpec("fixture",[self.Parameter(k,v) for k,v in defaults.items()],
                              factory,"particle" if connected else "a",lambda c,p:p,
                              intensity_factor=2.)

    def test_native_parity_across_shapes_and_links(self):
        for shape,connected in [("SolidSphere",False),("SolidCylinder",False),
                                 ("ThinRod",True),("GaussianPolymer",True)]:
            spec=self.simple(shape,connected)
            for beta in (2.,-2.):
                p={v.name:v.default for v in spec.parameters}
                p["beta_a"]=beta
                native=np.asarray(spec.world_factory().EvaluateFormFactorUnnormalized(
                    spec.structure,p,self.q.tolist()))*spec.intensity_factor
                np.testing.assert_allclose(self.calculate(spec,p),native,rtol=1e-11,atol=1e-12)

    def test_cloud_and_numerical_callback(self):
        from pyseb.sasmodels import SymbolicSpec, make_model_info
        def cloud_factory():
            world = pyseb.World("cloud")
            cloud = pyseb.DebyeSphereCloud([
                pyseb.SphereScatterer(0., 0., 0., 1., 2.),
                pyseb.SphereScatterer(3., 0., 0., 1., -2.),
            ])
            world.Add(cloud, "cloud")
            return world
        spec = self.ModelSpec("cloud", [], cloud_factory, "cloud", lambda c,p: {}, 1.)
        expected = cloud_factory().EvaluateFormFactorUnnormalized("cloud", {}, self.q.tolist())
        np.testing.assert_allclose(self.calculate(spec), expected, rtol=1e-12)
        def callback_factory():
            world = pyseb.World("callback")
            subunit = pyseb.NumericalSubunit(pyseb.NormalizationMode.Unnormalized)
            subunit.setTotalBeta(2.)
            subunit.setFormFactorFunction(lambda q,p: 4.*np.exp(-q*q))
            world.Add(subunit, "callback")
            return world
        spec = self.ModelSpec("callback", [], callback_factory, "callback", lambda c,p: {}, 1.)
        np.testing.assert_allclose(self.calculate(spec), 4.*np.exp(-self.q**2), rtol=1e-12)
        with self.assertRaises(ValueError):
            make_model_info(dataclasses.replace(spec, symbolic=SymbolicSpec({})), "compiled")

    def test_compiled_intensity_without_interactions(self):
        import sympy as sp
        from pyseb.sasmodels import SymbolicSpec, make_model_info
        spec=self.simple()
        spec=dataclasses.replace(spec,symbolic=SymbolicSpec(
            dict(R_a=sp.Symbol("R_a"),beta_a=sp.Symbol("beta_a"))))
        np.testing.assert_allclose(self.calculate(spec),self.calculate(spec,backend="compiled"),
                                   rtol=1e-6,atol=1e-10)
        bad=dataclasses.replace(spec,symbolic=SymbolicSpec(
            dict(R_a=sp.Symbol("undeclared"),beta_a=sp.Symbol("beta_a"))))
        with self.assertRaisesRegex(ValueError,"Unsupported"):
            make_model_info(bad,"compiled")

    def test_core_shell_backends_contrasts_and_distribution(self):
        from Examples.sasmodels_aot import make_model
        spec=make_model()
        for contrast in [(6.39,.62,-.28),(-.55,.62,6.68),(6.39,.62,6.68)]:
            p=dict(zip(("sld_core","sld_shell","sld_solvent"),contrast))
            for width in (0.,.12):
                pars=dict(p,radius_pd=width,radius_pd_n=11,thickness_pd=.15,thickness_pd_n=11)
                np.testing.assert_allclose(self.calculate(spec,pars),
                                           self.calculate(spec,pars,"compiled"),rtol=2e-7,atol=1e-8)
        # Total forward contrast match must not invoke normalized scattering.
        p={v.name:v.default for v in spec.parameters}
        p["sld_core"]=2.;p["sld_shell"]=1.
        vc=4*np.pi*p["radius"]**3/3;vo=spec.volume(p)
        p["sld_solvent"]=(2*vc+(vo-vc))/vo
        np.testing.assert_allclose(self.calculate(spec,p,q=self.q[2:]),
                                   self.calculate(spec,p,"compiled",q=self.q[2:]),rtol=1e-7,atol=1e-8)

    def test_distributions_and_absolute_volume(self):
        from Examples.sasmodels_aot import make_model
        from pyseb.sasmodels import positive_gaussian
        spec=make_model(); q=np.geomspace(.001,.2,31)
        coarse=dict(radius_pd=.12,radius_pd_n=41,radius_pd_type="lognormal",radius_pd_nsigma=8.,
                    thickness_pd=.36,thickness_pd_n=41,thickness_pd_type=positive_gaussian(),
                    thickness_pd_nsigma=8.)
        fine=dict(coarse,radius_pd_n=81,thickness_pd_n=81)
        np.testing.assert_allclose(self.calculate(spec,coarse,"compiled",q),
                                   self.calculate(spec,fine,"compiled",q),rtol=1e-5)
        p={v.name:v.default for v in spec.parameters}
        w=spec.world_factory(); native=np.asarray(w.EvaluateFormFactorUnnormalized(
            spec.structure,spec.bind(w,p),self.q.tolist()))*spec.intensity_factor/spec.volume(p)
        np.testing.assert_allclose(self.calculate(spec),native,rtol=1e-11)

    def test_hardsphere_and_limits(self):
        from Examples.sasmodels_aot import make_model
        from pyseb.sasmodels import with_hardsphere
        from sasmodels.core import build_model,load_model
        from sasmodels.direct_model import call_kernel,call_Fq
        spec=make_model(); q=self.q[2:]
        from pyseb.sasmodels import make_model_info
        self.assertEqual(make_model_info(spec, "compiled").parameters.nmagnetic, 0)
        self.assertEqual(with_hardsphere(spec, "compiled").parameters.nmagnetic, 0)
        for backend in ("runtime","compiled"):
            model=build_model(with_hardsphere(spec,backend),platform="dll",dtype="double")
            kernel=model.make_kernel([q])
            ref=load_model("core_shell_sphere@hardsphere",platform="dll",dtype="double")
            rk=ref.make_kernel([q])
            pars={v.name:v.default for v in spec.parameters}
            pars.update(volfraction=.06,structure_factor_mode=1,radius_effective_mode=1,
                        radius_pd=.12,radius_pd_n=15,background=0.)
            try:
                np.testing.assert_allclose(call_kernel(kernel,pars),call_kernel(rk,pars),rtol=1e-7)
                tiny=dict(pars,volfraction=1e-10,scale=.06/1e-10)
                plain={k:v for k,v in pars.items() if k not in (
                    "volfraction","structure_factor_mode","radius_effective_mode")}
                plain["scale"]=.06
                np.testing.assert_allclose(call_kernel(kernel,tiny),
                    self.calculate(spec,plain,backend,q),rtol=1e-8)
            finally:
                kernel.release();rk.release();model.release();ref.release()

    def test_hybrid_hydration_cross_terms_and_reset(self):
        from Examples.sasmodels_hybrid import make_model
        from pyseb.sasmodels import make_model_info,with_hardsphere
        from sasmodels.core import build_model
        from sasmodels.direct_model import call_kernel
        spec=make_model();ctx=spec.world_factory()
        model=build_model(make_model_info(spec));k=model.make_kernel([self.q])
        defaults={p.name:p.default for p in spec.parameters}
        try:
            outputs=[]
            for changes in ({},dict(hydration=0.),dict(solvent_density=.2),{}):
                p=dict(defaults,**changes)
                native=np.asarray(ctx.world.EvaluateFormFactorUnnormalized(
                    spec.structure,spec.bind(ctx,p),self.q.tolist()))*spec.intensity_factor
                actual=call_kernel(k,dict(p,background=0.))
                np.testing.assert_allclose(actual,native,rtol=1e-12)
                outputs.append(actual)
            np.testing.assert_allclose(outputs[0],outputs[-1],rtol=1e-12)
            self.assertGreater(np.max(np.abs(outputs[0]/outputs[1]-1)),.001)
            p=spec.bind(ctx,defaults)
            particle=np.asarray(ctx.world.EvaluateFormFactorUnnormalized("hybrid",p,self.q.tolist()))
            separate=np.asarray(ctx.world.EvaluateFormFactorUnnormalized("cloud",p,self.q.tolist()))
            separate+=ctx.world.EvaluateFormFactorUnnormalized("sphere",p,self.q.tolist())
            self.assertGreater(np.max(np.abs(particle-separate)),1.)
            with self.assertRaisesRegex(ValueError,"requires declared"):
                with_hardsphere(spec)
            with self.assertRaisesRegex(ValueError,"SymbolicSpec"):
                make_model_info(spec,"compiled")
        finally:k.release();model.release()

    def test_hybrid_scipy_workflow(self):
        from Examples.sasmodels_hybrid import make_model, fit_synthetic_data
        from pyseb.sasmodels import make_model_info
        from sasmodels.core import build_model
        spec = make_model()
        model = build_model(make_model_info(spec))
        parameters = {p.name: p.default for p in spec.parameters}
        try:
            report, measured, sigma, truth, initial, prediction = fit_synthetic_data(
                model, parameters, np.geomspace(.001, .5, 101)
            )
            self.assertTrue(report["success"])
            self.assertLess(abs(report["fitted"]["radius"]-report["truth"]["radius"]), .15)
            self.assertLess(abs(report["fitted"]["hydration"]-report["truth"]["hydration"]), .015)
            self.assertLess(abs(report["fitted"]["background"]-report["truth"]["background"]), .0003)
            self.assertGreater(report["reduced_chi_squared"], .5)
            self.assertLess(report["reduced_chi_squared"], 1.5)
            self.assertTrue(np.all(sigma > 0))
            self.assertTrue(np.all(np.isfinite(prediction)))
            self.assertLess(np.sum(((prediction-measured)/sigma)**2),
                            .1*np.sum(((initial-measured)/sigma)**2))
        finally:
            model.release()

    def test_thread_isolation(self):
        from pyseb.sasmodels import make_model_info
        from sasmodels.core import build_model
        from sasmodels.direct_model import call_kernel
        spec=self.simple();count=[];factory=spec.world_factory;barrier=threading.Barrier(2)
        def create():
            count.append(threading.get_ident())
            return factory()
        spec=dataclasses.replace(spec,world_factory=create)
        model=build_model(make_model_info(spec))
        def run(beta):
            k=model.make_kernel([self.q]);barrier.wait()
            try:return call_kernel(k,dict(beta_a=beta,background=0.))
            finally:k.release()
        with ThreadPoolExecutor(max_workers=2) as pool:
            a,b=list(pool.map(run,(2.,4.)))
        np.testing.assert_allclose(b,4*a)
        self.assertEqual(len(set(count)),2)
        model.release()

    def test_resolution_and_synthetic_fit(self):
        from scipy.optimize import least_squares
        from sasmodels.resolution import Perfect1D,Pinhole1D
        q=np.geomspace(.001,.5,51);spec=self.simple()
        perfect=Perfect1D(q);np.testing.assert_equal(perfect.apply(q*q),q*q)
        smear=Pinhole1D(q,.05*q)
        np.testing.assert_allclose(smear.apply(np.ones_like(smear.q_calc)),1.)
        truth=self.calculate(spec,dict(R_a=4.2),q=q)
        fit=least_squares(lambda p:(self.calculate(spec,dict(R_a=p[0]),q=q)-truth)/.01,
                          [3.],bounds=([1.],[8.]))
        self.assertTrue(fit.success);self.assertAlmostEqual(fit.x[0],4.2,places=6)

    def test_optional_dependency_and_capability_errors(self):
        script='''import sys, importlib.abc
class Block(importlib.abc.MetaPathFinder):
 def find_spec(self,name,*args):
  if name == "sasmodels" or name.startswith("sasmodels."): raise ModuleNotFoundError(name)
sys.meta_path.insert(0,Block())
import pyseb
import pyseb.sasmodels as bridge
try: bridge._sasmodels()
except ImportError as e: assert 'pyseb[sasmodels]' in str(e)
else: raise AssertionError('missing dependency not detected')
'''
        subprocess.run([sys.executable,"-c",script],check=True)
        from pyseb.sasmodels import make_model_info
        spec=self.simple()
        with self.assertRaises(ValueError):make_model_info(dataclasses.replace(spec,intensity_factor=0))
        with self.assertRaises(ValueError):make_model_info(dataclasses.replace(spec,name="bad name"))
        with self.assertRaises(ValueError):make_model_info(spec,"automatic")


if __name__=="__main__":unittest.main()
