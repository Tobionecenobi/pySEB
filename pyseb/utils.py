"""
Utility functions for the SEB Python bindings.
"""

import numpy as np
from ._pyseb import to_string_python, to_string_latex, to_string_cform


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
    return to_string_python(expression)


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
    return to_string_latex(expression)


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
    return to_string_cform(expression)


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
    param_list = {}
    for name, value in parameters.items():
        world.setParameter(param_list, name, value)
    
    if q_values is None:
        return world.Evaluate(expression, param_list)
    else:
        q_array = np.asarray(q_values)
        if q_array.ndim == 0:
            return world.Evaluate(expression, param_list, float(q_array))
        else:
            return world.Evaluate(expression, param_list, q_array.tolist())
