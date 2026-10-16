"""AOT core-shell model builder for fitting in Python with sasmodels."""
import numpy as np
import sympy as sp
import pyseb
from pyseb.sasmodels import Parameter, ModelSpec, InteractionSpec, SymbolicSpec


def make_model():
    # Build a fresh particle for each evaluation thread. The shared center
    # connects the core and shell without making their material volumes overlap.
    def world_factory():
        world = pyseb.World("AOTAdapter")
        graph = world.Add("SolidSphere", "core")
        world.Link("SolidSphericalShell", "shell.center", "core.center")
        world.Add(graph, "droplet")
        return world

    # Translate the fit parameters into native pySEB dimensions and signed
    # excess scattering lengths: beta = (component SLD - solvent SLD) * volume.
    # The 1e-6 factor converts the supplied SLD units into inverse angstroms squared.
    def bind(world, p):
        vc = 4*np.pi*p["radius"]**3/3
        vo = 4*np.pi*(p["radius"]+p["thickness"])**3/3
        return dict(R_core=p["radius"], Ri_shell=p["radius"],
                    Ro_shell=p["radius"]+p["thickness"],
                    beta_core=(p["sld_core"]-p["sld_solvent"])*1e-6*vc,
                    beta_shell=(p["sld_shell"]-p["sld_solvent"])*1e-6*(vo-vc))

    # Describe the same mapping symbolically for the optional compiled backend.
    # These symbols are model parameters, not a second implementation of scattering.
    r,t,rc,rs,solvent = sp.symbols("radius thickness sld_core sld_shell sld_solvent")
    vc = 4*sp.pi*r**3/3
    vo = 4*sp.pi*(r+t)**3/3
    world = world_factory()
    bc,bs = sp.symbols("beta_core beta_shell")
    # Remove the normalization from the center amplitude so its sign survives
    # contrast variation. sasmodels needs this amplitude for decoupling.
    amplitude = sp.cancel(sp.sympify(str(world.FormFactorAmplitude("droplet:core.center")))*(bc+bs))
    # Declare the public fitting interface. Only radius and thickness have
    # dispersity enabled (the final True argument to Parameter). The 1e8
    # conversion and total particle volume make scale a volume fraction.
    # InteractionSpec opts into an approximation; it does not establish that
    # these droplets physically obey a hard-sphere interaction model.
    return ModelSpec(
        name="pyseb_aot", title="Linked pySEB core-shell droplet",
        description="Connected sphere and shell; betas are excess scattering lengths in angstroms.",
        parameters=[
            Parameter("radius",25.,(0.,np.inf),"Ang","Core radius",True),
            Parameter("thickness",8.,(0.,np.inf),"Ang","Shell thickness",True),
            Parameter("sld_core",6.39,units="1e-6/Ang^2"),
            Parameter("sld_shell",.62,units="1e-6/Ang^2"),
            Parameter("sld_solvent",-.28,units="1e-6/Ang^2"),
        ], world_factory=world_factory, structure="droplet", bind=bind,
        intensity_factor=1e8,
        volume=lambda p: 4*np.pi*(p["radius"]+p["thickness"])**3/3,
        volume_parameters=("radius","thickness"),
        interaction=InteractionSpec("droplet:core.center",lambda p:p["radius"]+p["thickness"]),
        symbolic=SymbolicSpec(
            parameters=dict(R_core=r,Ri_shell=r,Ro_shell=r+t,
                            beta_core=(rc-solvent)*sp.Rational(1,10**6)*vc,
                            beta_shell=(rs-solvent)*sp.Rational(1,10**6)*(vo-vc)),
            volume=vo, amplitude=amplitude, effective_radius=r+t),
    )
