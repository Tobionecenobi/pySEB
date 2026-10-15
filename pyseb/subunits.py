"""
Subunits module for the SEB Python bindings.

This module provides Python wrappers for the SEB subunits.
"""

from . import _pyseb as _extension

_SymbolicSubunit = _extension.SymbolicSubunit
_NumericalSubunit = _extension.NumericalSubunit
_FileDefinedSubunit = _extension.FileDefinedSubunit
_DebyeSphereCloud = _extension.DebyeSphereCloud
_AtomicSAXSCloud = _extension.AtomicSAXSCloud
SphereScatterer = _extension.SphereScatterer
AtomicSAXSScatterer = _extension.AtomicSAXSScatterer
CromerMannCoefficients = _extension.CromerMannCoefficients
AtomicSAXSOptions = _extension.AtomicSAXSOptions
AtomicSAXSHydrationModel = _extension.AtomicSAXSHydrationModel
AtomicSAXSHydrationWaterModel = _extension.AtomicSAXSHydrationWaterModel
CartesianPoint3D = _extension.CartesianPoint3D
NormalizationMode = _extension.NormalizationMode


class SymbolicSubunit(_SymbolicSubunit):
    """A subunit represented by symbolic F, A and Psi placeholders."""
    pass


class NumericalSubunit(_NumericalSubunit):
    """A callback-driven numerical subunit."""
    pass


class FileDefinedSubunit(_FileDefinedSubunit):
    """An analytic subunit created from a validated .pyseb.yaml definition."""
    pass


class DebyeSphereCloud(_DebyeSphereCloud):
    """A Debye sphere cloud with finite-radius sphere scatterers."""
    pass


class AtomicSAXSCloud(_AtomicSAXSCloud):
    """An atomic solution-SAXS cloud with q-dependent form factors."""
    pass


__all__ = [
    "SymbolicSubunit",
    "NumericalSubunit",
    "FileDefinedSubunit",
    "DebyeSphereCloud",
    "AtomicSAXSCloud",
    "SphereScatterer",
    "AtomicSAXSScatterer",
    "CromerMannCoefficients",
    "AtomicSAXSOptions",
    "AtomicSAXSHydrationModel",
    "AtomicSAXSHydrationWaterModel",
    "CartesianPoint3D",
    "NormalizationMode",
]
