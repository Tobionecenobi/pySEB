"""
Python bindings for the Scattering Equation Builder (SEB) library.

SEB is a C++ library that allows you to build structures from sub-units such as spheres, 
polymers, rods, etc., and obtain their scattering properties symbolically.
"""

from ._pyseb import World, GraphID
from .subunits import (
    GaussianPolymer,
    GaussianLoop,
    ThinCircle,
    ThinRod,
    ThinDisk,
    ThinSphericalShell,
    SolidSphere,
    SolidSphericalShell,
    SolidCylinder,
    Point
)
from .utils import to_python, to_latex, to_cform, evaluate_expression

__version__ = "0.1.0"
