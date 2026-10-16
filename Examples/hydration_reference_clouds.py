"""Explicit reference conventions for the hydration component benchmarks.

These adapters keep the reference's point-weight normalization separate from
the physical volume-density default. They require pySEB, not IMP or Pepsi.
"""
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pyseb


def pepsi_shell(atoms, profile, contrast=0.0167):
    """Pepsi 3-style grid; atomic radii still come from the selected profile.

    The reference's bulk-density prefactor is retained explicitly in point
    weights. This maps a reference convention; it is not a density fit.
    """
    options = pyseb.AtomicSAXSOptions()
    options.hydration_model = pyseb.AtomicSAXSHydrationModel.SpatialShell
    options.hydration_water_model = pyseb.AtomicSAXSHydrationWaterModel.BondedWater
    options.include_hydration = True
    options.hydration_shell_contrast = contrast
    build = pyseb.AtomicSAXSCloudBuildOptions()
    build.calculate_solvent_accessibility = False
    build.hydration_shell_convention = pyseb.AtomicSAXSHydrationShellConvention.Pepsi3
    full = pyseb.AtomicSAXSCloudBuilder.build(atoms, profile, options, build)
    scatterers = [p for p in full.get_scatterers() if p.hydration_volume > 0]
    points = np.array([[p.center.x,p.center.y,p.center.z] for p in scatterers])
    return pyseb.AtomicSAXSCloud(scatterers, options), points, full


def foxs_hydration(atoms, exposure, mode, water_f0):
    """Hydration-only cloud with explicit native exposures as a control.

    Unit c2 maps to water_f0 electrons per unit exposure. Nothing is fitted
    to a molecular intensity or its forward value.
    """
    if len(atoms) != len(exposure):
        raise ValueError('One exposure value is required for each atom')
    scatterers = []
    for atom, fraction in zip(atoms, exposure):
        point = pyseb.AtomicSAXSScatterer()
        point.center = atom.coordinate
        point.solvent_accessibility = float(fraction)
        scatterers.append(point)
    options = pyseb.AtomicSAXSOptions()
    options.include_hydration = True
    options.hydration_shell_contrast = water_f0/30.0
    options.hydration_water_model = getattr(pyseb.AtomicSAXSHydrationWaterModel,
        'FoXSReciprocal' if mode == 'reciprocal' else 'FoXSFast')
    return pyseb.AtomicSAXSCloud(scatterers, options)


def foxs_geometry_exposure(atoms):
    """Independent FreeSASA estimate with FoXS radii and a 1.8 A probe.

    This does not import native exposures and retains a different surface
    quadrature, making its residual distinct from the common-exposure check.
    """
    profile = pyseb.AtomicSAXSParameterProfile.protein_defaults()
    for atom in atoms:
        p = profile.resolve(atom)
        p.sasa_radius = (3*p.excluded_volume/(4*np.pi))**(1/3)
        profile.set_atom(atom.residue_name,atom.atom_name,p)
    build = pyseb.AtomicSAXSCloudBuildOptions()
    build.sasa_use_profile_radii = True
    build.sasa_probe_radius = 1.8
    build.sasa_point_count = 2000
    cloud = pyseb.AtomicSAXSCloudBuilder.build(atoms,profile,pyseb.AtomicSAXSOptions(),build)
    return [p.solvent_accessibility for p in cloud.get_scatterers()]
