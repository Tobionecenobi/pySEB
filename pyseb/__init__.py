"""
Python bindings for the Scattering Equation Builder (SEB) library.
"""
# Import the extension module
import os
import sys
import importlib
import importlib.util

_pyseb = None

# First try normal package import (works for regular wheel installs)
try:
    _pyseb = importlib.import_module("_pyseb")
except Exception:
    pass

if _pyseb is None:
    # Fallback for editable/developer builds where extension is in a build folder
    dir_path = os.path.dirname(os.path.abspath(__file__))
    root_path = os.path.dirname(dir_path)

    candidate_dirs = [
        dir_path,
        os.path.join(dir_path, "pyseb"),
        os.path.join(root_path, "build-ninja-ginac", "pyseb"),
        os.path.join(root_path, "build-ninja", "pyseb"),
        os.path.join(root_path, "build", "pyseb"),
    ]

    candidate_paths = []
    for candidate_dir in candidate_dirs:
        if not os.path.isdir(candidate_dir):
            continue
        for file in os.listdir(candidate_dir):
            if file.startswith('_pyseb') and (file.endswith('.so') or file.endswith('.pyd')):
                candidate_paths.append(os.path.join(candidate_dir, file))

    last_error = None
    for module_path in candidate_paths:
        try:
            spec = importlib.util.spec_from_file_location('_pyseb', module_path)
            module = importlib.util.module_from_spec(spec)
            sys.modules['_pyseb'] = module
            spec.loader.exec_module(module)
            _pyseb = module
            break
        except Exception as exc:
            sys.modules.pop('_pyseb', None)
            last_error = exc

    if _pyseb is None and last_error is not None:
        raise ImportError(f"Could not load _pyseb extension module from candidates: {last_error}") from last_error

if _pyseb is None:
    raise ImportError("Could not find _pyseb extension module")

World = _pyseb.World
Expression = _pyseb.Expression
GraphID = getattr(_pyseb, "GraphID", int)
available_backends = _pyseb.available_backends
get_backend = _pyseb.get_backend
set_backend = _pyseb.set_backend
symbol = _pyseb.symbol
constant = _pyseb.constant
pi = _pyseb.pi
e = _pyseb.e
i = _pyseb.i
pow = _pyseb.pow
sin = _pyseb.sin
cos = _pyseb.cos
tan = _pyseb.tan
exp = _pyseb.exp
log = _pyseb.log
sqrt = _pyseb.sqrt
abs = _pyseb.abs

from .symbolic import (
    SymPyExpression,
    SymPyFactory,
    get_factory,
    set_factory,
    to_latex,
    to_python
)
from .utils import evaluate_expression

# Set SymPy as the default symbolic engine for Python
set_factory(SymPyFactory())

def to_sympy(expr):
    if isinstance(expr, SymPyExpression):
        return expr.expr
    if hasattr(_pyseb, "Expression") and isinstance(expr, _pyseb.Expression):
        converted = _pyseb.to_sympy(expr)
        if isinstance(converted, SymPyExpression):
            return converted.expr
        return converted
    if hasattr(expr, "to_python"):
        return SymPyExpression(expr.to_python()).expr
    return SymPyExpression(str(expr)).expr

__all__ = [
    'World',
    'Expression',
    'GraphID',
    'available_backends',
    'get_backend',
    'set_backend',
    'symbol',
    'constant',
    'pi',
    'e',
    'i',
    'pow',
    'sin',
    'cos',
    'tan',
    'exp',
    'log',
    'sqrt',
    'abs',
    'SymPyExpression',
    'SymPyFactory',
    'get_factory',
    'set_factory',
    'to_latex',
    'to_python',
    'to_sympy',
    'evaluate_expression'
]
