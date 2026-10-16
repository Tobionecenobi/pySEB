"""Proof of concept: let sasmodels average a linked pySEB core-shell particle.

Run with the same environment as pySEB:
    python -m pip install 'sasmodels==1.1.0'
    python Examples/sasmodels_pyseb_adapter.py

This is an investigation, not a public adapter API. No saved fits are loaded.
"""
from types import ModuleType

import numpy as np
import pyseb
import sympy as sp
from sasmodels.core import build_model, load_model
from sasmodels.direct_model import call_Fq, call_kernel
from sasmodels.modelinfo import make_model_info
from sasmodels.weights import get_weights


def make_adapter():
    world = pyseb.World("SasmodelsAdapter")
    graph = world.Add("SolidSphere", "core")
    world.Link("SolidSphericalShell", "shell.center", "core.center")
    world.Add(graph, "droplet")

    q, radius, inner, outer, bc, bs = sp.symbols(
        "q R_core Ri_shell Ro_shell beta_core beta_shell"
    )
    args = (q, radius, inner, outer, bc, bs)
    amplitude = sp.sympify(str(world.FormFactorAmplitude("droplet:core.center")))
    intensity = sp.sympify(str(world.FormFactor("droplet")))
    # Cancel the normalization algebraically, including contrast-match cases.
    amplitude = sp.cancel(amplitude * (bc + bs))
    intensity = sp.cancel(intensity * (bc + bs)**2)
    np.testing.assert_equal(sp.cancel(amplitude**2 - intensity), 0)
    evaluate_amplitude = sp.lambdify(args, amplitude, "numpy")

    def volumes(radius, thickness):
        core = 4 * np.pi * radius**3 / 3
        total = 4 * np.pi * (radius + thickness)**3 / 3
        return core, total

    def betas(radius, thickness, sld_core, sld_shell, sld_solvent):
        core, total = volumes(radius, thickness)
        return ((sld_core - sld_solvent) * 1e-6 * core,
                (sld_shell - sld_solvent) * 1e-6 * (total - core))

    def Fq(q, radius, thickness, sld_core, sld_shell, sld_solvent):
        bc_value, bs_value = betas(radius, thickness, sld_core, sld_shell, sld_solvent)
        arguments = (q, radius, radius, radius + thickness, bc_value, bs_value)
        # sasmodels uses F in sqrt(cm^-1 * Angstrom^3), F^2 in cm^-1 * Angstrom^3.
        f = 1e4 * evaluate_amplitude(*arguments)
        # For this spherical particle only, F^2 equals the squared center
        # amplitude (proved symbolically above). This is not general for
        # anisotropic or flexible pySEB structures.
        return f, f**2

    Fq.vectorized = True
    module = ModuleType("pyseb_core_shell")
    module.__file__ = __file__
    module.name = "Linked pySEB core-shell adapter"
    module.parameters = [
        ["radius", "Ang", 25., [0.01, np.inf], "volume", "Core radius"],
        ["thickness", "Ang", 8., [0.01, np.inf], "volume", "Shell thickness"],
        ["sld_core", "1e-6/Ang^2", 6.39, [-np.inf, np.inf], "sld", "Core SLD"],
        ["sld_shell", "1e-6/Ang^2", .62, [-np.inf, np.inf], "sld", "Shell SLD"],
        ["sld_solvent", "1e-6/Ang^2", -.28, [-np.inf, np.inf], "sld", "Solvent SLD"],
    ]
    module.Fq = Fq
    module.form_volume = lambda radius, thickness: volumes(radius, thickness)[1]
    return build_model(make_model_info(module)), Fq, volumes


def check_adapter():
    model, Fq, volumes = make_adapter()
    q = np.geomspace(.001, .3, 81)
    kernel = model.make_kernel([q])
    reference = load_model("core_shell_sphere", platform="dll", dtype="double")
    reference_kernel = reference.make_kernel([q])
    phi, background = .06, .02
    # sasmodels lognormal center is the median; width is sigma in log space.
    # Convert the tutorial's arithmetic mean and coefficient of variation.
    mean_radius, radius_cv = 25., .155
    sigma_log = np.sqrt(np.log1p(radius_cv**2))
    median_radius = mean_radius * np.exp(-sigma_log**2 / 2)
    dispersity = dict(radius=median_radius, radius_pd=sigma_log,
                      radius_pd_type="lognormal", radius_pd_n=81,
                      radius_pd_nsigma=8., thickness=8., thickness_pd=.36,
                      thickness_pd_type="gaussian", thickness_pd_n=81,
                      thickness_pd_nsigma=8.)
    radii, rw = get_weights("lognormal", 81, sigma_log, 8., median_radius,
                           (.01, np.inf), True)
    thicknesses, tw = get_weights("gaussian", 81, .36, 8., 8., (.01, np.inf), True)
    contrasts = [(6.39, .62, -.28), (-.55, .62, 6.68), (6.39, .62, 6.68)]
    try:
        for rc, rs, solvent in contrasts:
            pars = dict(dispersity, sld_core=rc, sld_shell=rs, sld_solvent=solvent)
            f, f2, _, volume, ratio = call_Fq(kernel, pars.copy())
            # Independent explicit averaging on the library's identical nodes
            # isolates adapter normalization from differences in quadrature.
            r = radii[:, None, None]
            t = thicknesses[None, :, None]
            w = rw[:, None, None] * tw[None, :, None]
            expected_f, expected_f2 = Fq(q[None, None, :], r, t, rc, rs, solvent)
            expected_volume = np.sum(w * volumes(r, t)[1])
            np.testing.assert_allclose(f, np.sum(w * expected_f, axis=(0, 1)), rtol=1e-10, atol=1e-10)
            np.testing.assert_allclose(f2, np.sum(w * expected_f2, axis=(0, 1)), rtol=1e-10)
            np.testing.assert_allclose(volume, expected_volume, rtol=1e-12)
            np.testing.assert_allclose(ratio, 1.)
            result = call_kernel(kernel, dict(pars, scale=phi, background=background))
            np.testing.assert_allclose(result, phi * f2 / volume + background, rtol=1e-12)
            # A separate established particle implementation also checks
            # units, volume normalization, and contrast handling.
            reference_result = call_kernel(
                reference_kernel, dict(pars, scale=phi, background=background)
            )
            np.testing.assert_allclose(result, reference_result, rtol=1e-8, atol=1e-9)
            # Zero dispersity recovers the single pySEB particle.
            mono = dict(pars, radius_pd=0., thickness_pd=0.)
            mf, mf2, _, mv, _ = call_Fq(kernel, mono)
            ef, ef2 = Fq(q, median_radius, 8., rc, rs, solvent)
            np.testing.assert_allclose(mf, ef, rtol=1e-12, atol=1e-10)
            np.testing.assert_allclose(mf2, ef2, rtol=1e-12)
            np.testing.assert_allclose(mv, volumes(median_radius, 8.)[1])
            assert np.any(f < 0), "Signed amplitudes must survive averaging"
        print("PASS: three contrasts; signed moments; mean volume; absolute intensity; zero dispersity")
        print("PASS: agreement with independent sasmodels core_shell_sphere (rtol=1e-8)")
        print("Python Fq adapter works. Distribution grids and conventions differ from the tutorial.")
    finally:
        kernel.release()
        model.release()
        reference_kernel.release()
        reference.release()


if __name__ == "__main__":
    check_adapter()
