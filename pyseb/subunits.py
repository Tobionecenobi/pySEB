"""
Subunits module for the SEB Python bindings.

This module provides Python wrappers for the SEB subunits.
"""

from . import _pyseb as _extension

_SymbolicSubunit = _extension.SymbolicSubunit
_NumericalSubunit = _extension.NumericalSubunit
_FileDefinedSubunit = _extension.FileDefinedSubunit
_DebyeSphereCloud = _extension.DebyeSphereCloud
SphereScatterer = _extension.SphereScatterer
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


__all__ = [
    "SymbolicSubunit",
    "NumericalSubunit",
    "FileDefinedSubunit",
    "DebyeSphereCloud",
    "SphereScatterer",
    "CartesianPoint3D",
    "NormalizationMode",
]
