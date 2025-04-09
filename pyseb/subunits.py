"""
Subunits module for the SEB Python bindings.

This module provides Python wrappers for the SEB subunits.
"""

from ._pyseb import (
    GaussianPolymer as _GaussianPolymer,
    GaussianLoop as _GaussianLoop,
    ThinCircle as _ThinCircle,
    ThinRod as _ThinRod,
    ThinDisk as _ThinDisk,
    ThinSphericalShell as _ThinSphericalShell,
    SolidSphere as _SolidSphere,
    SolidSphericalShell as _SolidSphericalShell,
    SolidCylinder as _SolidCylinder,
    Point as _Point
)


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
