"""
Subunits module for the SEB Python bindings.

This module provides Python wrappers for the SEB subunits.
"""

from . import _pyseb as _extension

_GaussianPolymer = _extension.GaussianPolymer
_GaussianLoop = _extension.GaussianLoop
_ThinCircle = _extension.ThinCircle
_ThinRod = _extension.ThinRod
_ThinDisk = _extension.ThinDisk
_ThinSphericalShell = _extension.ThinSphericalShell
_SolidSphere = _extension.SolidSphere
_SolidSphericalShell = _extension.SolidSphericalShell
_SolidCylinder = _extension.SolidCylinder
_Point = _extension.Point
_SymbolicSubunit = _extension.SymbolicSubunit
_NumericalSubunit = _extension.NumericalSubunit
_DebyeSphereCloud = _extension.DebyeSphereCloud
SphereScatterer = _extension.SphereScatterer
CartesianPoint3D = _extension.CartesianPoint3D
NormalizationMode = _extension.NormalizationMode


class GaussianPolymer(_GaussianPolymer):
    """
    A Gaussian polymer subunit.
    
    This represents a polymer chain with Gaussian statistics.
    """
    pass


class GaussianLoop(_GaussianLoop):
    """
    A Gaussian loop subunit.
    
    This represents a closed loop polymer with Gaussian statistics.
    """
    pass


class ThinCircle(_ThinCircle):
    """
    A thin circle subunit.
    
    This represents a thin circular ring.
    """
    pass


class ThinRod(_ThinRod):
    """
    A thin rod subunit.
    
    This represents a thin rod with uniform scattering density.
    """
    pass


class ThinDisk(_ThinDisk):
    """
    A thin disk subunit.
    
    This represents a thin disk with uniform scattering density.
    """
    pass


class ThinSphericalShell(_ThinSphericalShell):
    """
    A thin spherical shell subunit.
    
    This represents a thin spherical shell with uniform scattering density.
    """
    pass


class SolidSphere(_SolidSphere):
    """
    A solid sphere subunit.
    
    This represents a solid sphere with uniform scattering density.
    """
    pass


class SolidSphericalShell(_SolidSphericalShell):
    """
    A solid spherical shell subunit.
    
    This represents a solid spherical shell with uniform scattering density.
    """
    pass


class SolidCylinder(_SolidCylinder):
    """
    A solid cylinder subunit.
    
    This represents a solid cylinder with uniform scattering density.
    """
    pass


class Point(_Point):
    """
    A point subunit.
    
    This represents a point scatterer.
    """
    pass


class SymbolicSubunit(_SymbolicSubunit):
    """A subunit represented by symbolic F, A and Psi placeholders."""
    pass


class NumericalSubunit(_NumericalSubunit):
    """A callback-driven numerical subunit."""
    pass


class DebyeSphereCloud(_DebyeSphereCloud):
    """A Debye sphere cloud with finite-radius sphere scatterers."""
    pass


__all__ = [
    "GaussianPolymer",
    "GaussianLoop",
    "ThinCircle",
    "ThinRod",
    "ThinDisk",
    "ThinSphericalShell",
    "SolidSphere",
    "SolidSphericalShell",
    "SolidCylinder",
    "Point",
    "SymbolicSubunit",
    "NumericalSubunit",
    "DebyeSphereCloud",
    "SphereScatterer",
    "CartesianPoint3D",
    "NormalizationMode",
]
