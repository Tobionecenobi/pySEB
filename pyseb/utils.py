"""
Utility functions for the SEB Python bindings.
"""

import numpy as np
import sympy
import _pyseb as _ext
from .symbolic import SymPyExpression

to_string_python = getattr(_ext, "to_string_python", None)
to_string_latex = getattr(_ext, "to_string_latex", None)
to_string_cform = getattr(_ext, "to_string_cform", None)


def to_python(expression):
    """
    Convert a GiNaC expression to a Python string.
    
    Parameters
    ----------
    expression : GiNaC expression
        The expression to convert.
        
    Returns
    -------
    str
        The expression as a Python string.
    """
    if to_string_python is not None:
        return to_string_python(expression)
    return SymPyExpression(str(expression)).to_python()


def to_latex(expression):
    """
    Convert a GiNaC expression to a LaTeX string.
    
    Parameters
    ----------
    expression : GiNaC expression
        The expression to convert.
        
    Returns
    -------
    str
        The expression as a LaTeX string.
    """
    if to_string_latex is not None:
        return to_string_latex(expression)
    return SymPyExpression(str(expression)).to_latex()


def to_cform(expression):
    """
    Convert a GiNaC expression to a C/C++ string.
    
    Parameters
    ----------
    expression : GiNaC expression
        The expression to convert.
        
    Returns
    -------
    str
        The expression as a C/C++ string.
    """
    if to_string_cform is not None:
        return to_string_cform(expression)
    return SymPyExpression(str(expression)).to_cform()


def evaluate_expression(world, expression, parameters, q_values=None):
    """
    Evaluate a GiNaC expression with the given parameters.
    
    Parameters
    ----------
    world : World
        The SEB World object.
    expression : GiNaC expression
        The expression to evaluate.
    parameters : dict
        Dictionary of parameter names and values.
    q_values : array_like, optional
        If provided, evaluate the expression at each q value.
        
    Returns
    -------
    float or ndarray
        The evaluated expression. If q_values is provided, returns an array
        of the same shape as q_values.
    """
    sympy_expr = SymPyExpression(str(expression)).expr

    for name, value in parameters.items():
        sympy_expr = sympy_expr.subs(sympy.Symbol(name), value)

    sympy_expr = sympy_expr.doit()

    if q_values is None:
        evaluated = sympy_expr.evalf()
        return float(evaluated)

    q_array = np.asarray(q_values)
    q_symbol = sympy.Symbol('q')
    eval_fn = sympy.lambdify(q_symbol, sympy_expr, modules=['numpy'])

    if q_array.ndim == 0:
        return float(eval_fn(float(q_array)))

    return np.asarray(eval_fn(q_array), dtype=float)
