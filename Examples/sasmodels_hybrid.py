"""Atomic-cloud/analytical hybrid for the Python sasmodels fitting workflow.

No experimental or saved model data are needed. Atomic amplitudes and the
analytical sphere beta are in electrons. The classical electron-radius factor
converts the *complete* hybrid scattering, preserving interference. The script
also creates explicitly synthetic measurements and fits them with SciPy.
"""
from types import SimpleNamespace
from pathlib import Path
import argparse
import json
import time

import numpy as np
import pyseb
from pyseb.sasmodels import ModelSpec, Parameter


def make_model():
    # Retain both the world and the cloud handle: fitting solvent or hydration
    # parameters must update the actual cloud used by the connected structure.
    def factory():
        points = []
        # Use two simple six-electron scatterers to make this example completely
        # reproducible without a PDB file. Excluded volume and solvent exposure
        # let the native atomic-cloud model include solvent and hydration effects.
        for x in (-2., 2.):
            point = pyseb.AtomicSAXSScatterer()
            point.center = pyseb.CartesianPoint3D(x,0.,0.)
            point.vacuum = pyseb.CromerMannCoefficients([0.]*4,[0.]*4,6.)
            point.excluded_volume = 12.
            point.solvent_accessibility = 1.
            points.append(point)
        options = pyseb.AtomicSAXSOptions()
        cloud = pyseb.AtomicSAXSCloud(points, options)
        # Attach a sampled point on the analytical sphere surface to a fixed
        # cloud reference. #join names that distributed surface reference.
        # pySEB calculates the interference of the resulting connected particle.
        cloud.add_reference_point("join", 2.,0.,0.)
        world = pyseb.World("HybridAdapter")
        graph = world.Add(cloud,"cloud")
        world.Link("SolidSphere","sphere.surface#join","cloud.join")
        world.Add(graph,"hybrid")
        return SimpleNamespace(world=world,cloud=cloud)

    # Apply every adjustable cloud setting on every evaluation. This prevents
    # a dry or changed-solvent calculation from affecting the next parameter set.
    def bind(ctx,p):
        ctx.cloud.set_solvent_electron_density(p["solvent_density"])
        ctx.cloud.set_hydration(p["hydration"] != 0., p["hydration"], 30.)
        return dict(R_sphere=p["radius"], beta_sphere=p["beta"])

    # Keep this a per-particle model: no particle-volume normalization or
    # structure-factor coupling is declared. Squaring the electron radius in
    # angstroms and multiplying by 1e-16 converts electron-squared scattering
    # into cm^2 for the complete hybrid, including its cross terms.
    return ModelSpec(
        name="pyseb_hybrid",title="Atomic cloud linked to an analytical sphere",
        description="Relative per-particle hybrid intensity in cm^2; no volume-fraction or interaction claim.",
        parameters=[Parameter("radius",4.,(1.,10.),"Ang","Sphere radius",True),
                    Parameter("beta",20.,(-100.,100.),"electron","Sphere excess amplitude"),
                    Parameter("hydration",.02,(-.1,.1),"electron/Ang^3","Hydration contrast"),
                    Parameter("solvent_density",.334,(0.,.5),"electron/Ang^3","Solvent density")],
        world_factory=factory,structure="hybrid",bind=bind,
        intensity_factor=(2.8179403262e-5)**2 * 1e-16,
    )


def fit_synthetic_data(model, parameters, q):
    """Demonstrate native scattering -> resolution -> data -> bounded fit.

    Return the numerical fit and the data needed to inspect it. These are toy
    measurements, not experimental evidence for a hydration model.
    """
    from sasmodels.direct_model import call_kernel
    from sasmodels.resolution import Pinhole1D
    from scipy.optimize import least_squares

    # A known number density converts the per-particle cm^2 scattering into
    # intensity in cm^-1. Hold it, beta, and solvent density fixed so this small
    # example concentrates on radius, hydration contrast, and flat background.
    number_density = 1e20  # particles per cm^3, chosen for this synthetic sample
    resolution = Pinhole1D(q, .05*q)
    fit_kernel = model.make_kernel([resolution.q_calc])
    truth = np.array([4.8, .035, .0002])
    initial = np.array([3.5, .01, .0001])
    names = ("radius", "hydration", "background")

    def predict(values):
        # sasmodels calls pySEB for the complete hybrid at its resolution grid.
        # Hydration changes the cloud amplitude before hybrid interference is
        # calculated. Add the constant background after smearing the intensity.
        radius, hydration, background = values
        theory = call_kernel(fit_kernel, dict(
            parameters, radius=radius, hydration=hydration,
            scale=number_density, background=0.,
        ))
        return resolution.apply(theory) + background

    try:
        expected = predict(truth)
        # Use a fixed random seed and known uncertainties for reproducibility.
        # The small absolute noise term avoids enormous weights near minima.
        uncertainty = .02*expected + .002*expected.max()
        measured = expected + np.random.default_rng(2026).normal(0., uncertainty)

        def residual(values):
            # SciPy sees dimensionless residuals weighted by measurement error.
            return (predict(values) - measured) / uncertainty

        fit = least_squares(
            residual, initial,
            bounds=([1., 0., 0.], [10., .1, .01]), x_scale="jac",
        )
        if not fit.success:
            raise RuntimeError(f"Hybrid fit did not converge: {fit.message}")
        prediction = predict(fit.x)
        start_prediction = predict(initial)
        chi_squared = float(np.dot(fit.fun, fit.fun))
        report = dict(
            data_source="synthetic, fixed seed 2026; not experimental data",
            truth=dict(zip(names, truth.tolist())),
            initial=dict(zip(names, initial.tolist())),
            fitted=dict(zip(names, fit.x.tolist())),
            fixed_number_density_cm3=number_density,
            fixed_beta_electrons=parameters["beta"],
            fixed_solvent_density=parameters["solvent_density"],
            assumed_relative_q_resolution=.05,
            chi_squared=chi_squared,
            degrees_of_freedom=int(fit.fun.size-fit.x.size),
            reduced_chi_squared=chi_squared/(fit.fun.size-fit.x.size),
            success=bool(fit.success), message=fit.message,
        )
        return report, measured, uncertainty, expected, start_prediction, prediction
    finally:
        fit_kernel.release()


def main():
    from sasmodels.core import build_model
    from sasmodels.direct_model import call_kernel
    from sasmodels.resolution import Pinhole1D
    from pyseb.sasmodels import make_model_info
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=Path("/tmp/pyseb-hybrid"))
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    # Build the sasmodels runtime model. Its particle evaluations call pySEB
    # in C++; sasmodels supplies parameter handling and population averaging.
    spec=make_model()
    model=build_model(make_model_info(spec))
    q=np.geomspace(.001,.5,101)
    # This illustrative resolution width is assumed, not measured. Smearing
    # requests its own q_calc grid, so prepare a separate kernel for that grid.
    resolution=Pinhole1D(q,.05*q)
    kernel=model.make_kernel([q])
    smear_kernel=model.make_kernel([resolution.q_calc])
    p={v.name:v.default for v in spec.parameters}
    curves={}
    # Change one feature at a time. The distribution case asks sasmodels to
    # average 21 radius samples; pySEB still evaluates each complete hybrid.
    for name,changes in [("hydrated",{}),("dry",dict(hydration=0.)),
                         ("larger sphere",dict(radius=6.)),
                         ("radius distribution",dict(radius_pd=.15,radius_pd_n=21))]:
        curves[name]=call_kernel(kernel,dict(p,**changes,background=0.))
    # Smear the complete hydrated intensity, rather than separately smearing
    # the cloud and sphere and losing their interference contribution.
    curves["resolution"]=resolution.apply(call_kernel(smear_kernel,dict(p,background=0.)))
    # Return to the starting parameters after the variants and verify that
    # the adapter agrees with a separate direct pySEB evaluation.
    np.testing.assert_allclose(call_kernel(kernel,dict(p,background=0.)),curves["hydrated"])
    ctx=spec.world_factory()
    parameters=spec.bind(ctx,p)
    reference=np.asarray(ctx.world.EvaluateFormFactorUnnormalized("hybrid",parameters,q.tolist()))*spec.intensity_factor
    np.testing.assert_allclose(curves["hydrated"],reference,rtol=1e-12)
    # Fit the same hybrid through SciPy, now including the measurement model.
    # Keep the truth, noisy data, initial guess, and fitted curve available for
    # comparison rather than judging success from the optimizer flag alone.
    fit_report, measured, uncertainty, expected, start_prediction, prediction = (
        fit_synthetic_data(model, p, q)
    )
    (args.output/"fit.json").write_text(json.dumps(fit_report, indent=2)+"\n")
    np.savetxt(
        args.output/"fit.csv",
        np.column_stack([q, measured, uncertainty, expected, start_prediction, prediction]),
        delimiter=",", header="q,synthetic_intensity,uncertainty,truth,initial,fit",
    )
    # Warm each calculation once, then report the mean time over repeated
    # calls. Compare the ordinary adapter overhead separately from the extra
    # work of averaging a population; these timings are machine dependent.
    repeats=20
    def timed(fn):
        fn()
        start=time.perf_counter()
        for _ in range(repeats):fn()
        return (time.perf_counter()-start)/repeats
    native=timed(lambda:ctx.world.EvaluateFormFactorUnnormalized("hybrid",parameters,q.tolist()))
    runtime=timed(lambda:call_kernel(kernel,dict(p,background=0.)))
    averaged=timed(lambda:call_kernel(kernel,dict(p,radius_pd=.15,radius_pd_n=21,background=0.)))
    report=dict(native_seconds=native,runtime_seconds=runtime,
                runtime_over_native=runtime/native,averaged_21_seconds=averaged,
                q_count=len(q),repeats=repeats)
    # Save reproducible numeric curves, a small comparison plot, and timing
    # evidence in the requested output directory, outside the source examples.
    (args.output/"benchmark.json").write_text(json.dumps(report,indent=2)+"\n")
    np.savetxt(args.output/"curves.csv",np.column_stack([q,*curves.values()]),
               delimiter=",",header="q,"+",".join(curves))
    import matplotlib.pyplot as plt
    for name,curve in curves.items():plt.loglog(q,curve,label=name)
    plt.xlabel("q (1/Angstrom)");plt.ylabel("Per-particle intensity (cm^2)")
    plt.legend();plt.tight_layout();plt.savefig(args.output/"hybrid.png");plt.close()
    # Show the fit and its uncertainty-weighted residuals. A visually close
    # curve alone can hide systematic disagreement, especially on a log scale.
    figure, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True,
                                gridspec_kw={"height_ratios": [3, 1]})
    axes[0].errorbar(q, measured, yerr=uncertainty, fmt=".", markersize=3,
                     label="Synthetic measurements")
    axes[0].plot(q, expected, ":", label="Generating model")
    axes[0].plot(q, start_prediction, "--", label="Initial guess")
    axes[0].plot(q, prediction, label="SciPy fit")
    axes[0].set(xscale="log", yscale="log", ylabel="Intensity (cm$^{-1}$)",
                title=f"Synthetic hybrid fit: reduced chi-squared = {fit_report['reduced_chi_squared']:.2f}")
    axes[0].legend()
    axes[1].semilogx(q, (prediction-measured)/uncertainty, ".", markersize=3)
    axes[1].axhline(0., color="black", linewidth=.8)
    axes[1].set(xlabel="q (1/Angstrom)", ylabel="Residual / sigma")
    figure.tight_layout();figure.savefig(args.output/"fit.png");plt.close(figure)
    # Release sasmodels kernel resources when all evaluations are finished.
    kernel.release();smear_kernel.release();model.release()
    print("Synthetic-data fit (not experimental data):")
    print(json.dumps(fit_report, indent=2))
    print("Benchmark:")
    print(json.dumps(report,indent=2))


if __name__=="__main__":main()
