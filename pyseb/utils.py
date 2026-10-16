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


def compile_expression(expression):
    """Prepare a reusable NumPy/SciPy evaluator for an analytical equation.

    Call the result as ``evaluate(parameters, q_values)``. All free symbols
    except ``q`` require scalar parameters on every call; extra parameters are
    ignored and q comes from q_values. Scalar inputs return float; arrays keep
    their shape, including constants and empty arrays. No World is retained.

    SciPy is optional and enables additional special functions.
    Common subexpressions are calculated once per call. This prepares Python
    numerical code, not a machine-code kernel. Unresolved integrals and
    symbolic functions are rejected: use native World evaluation for numerical
    quadrature, clouds, and hybrids. Equation singularities are not repaired.
    Only compile trusted expressions, since SymPy lambdify generates code.
    """
    # Import locally to avoid the package's initialization cycle and retain its
    # native-expression conversion instead of reparsing a printed equation.
    from . import to_sympy
    from sympy.core.function import AppliedUndef

    expr = expression if isinstance(expression, sympy.Expr) else to_sympy(expression)
    if expr.has(sympy.Integral, sympy.Derivative, sympy.Sum, sympy.Product, AppliedUndef):
        raise ValueError("Unresolved symbolic operation; use native World evaluation")
    symbols = sorted(expr.free_symbols, key=str)
    names = [str(symbol) for symbol in symbols]
    if len(set(names)) != len(names):
        raise ValueError("Expression has distinct symbols with the same parameter name")
    modules = ['numpy']
    try:
        import scipy.special  # Optional: extend the available special functions.
    except ImportError:
        pass
    else:
        modules.append('scipy')
    try:
        function = sympy.lambdify(symbols, expr, modules=modules,
                                 cse=True, dummify=True, docstring_limit=0)
    except (NotImplementedError, TypeError) as error:
        raise ValueError("Unsupported analytical expression") from error

    def evaluate(parameters, q_values):
        q_array = np.asarray(q_values, dtype=float)
        if not np.all(np.isfinite(q_array)):
            raise ValueError("q must be finite")
        # Bind every parameter again, allowing fitting to change its values
        # without retaining mutable parameter or model state in the evaluator.
        arguments = []
        for name in names:
            if name == 'q':
                arguments.append(q_array)
            else:
                if name not in parameters:
                    raise ValueError(f"Missing parameter: {name}")
                value = float(parameters[name])
                if not np.isfinite(value):
                    raise ValueError(f"Parameter must be finite: {name}")
                arguments.append(value)
        try:
            result = np.asarray(function(*arguments), dtype=float)
        except NameError as error:
            raise ValueError("Unsupported numerical function; some functions require SciPy") from error
        result = np.broadcast_to(result, q_array.shape)
        if not np.all(np.isfinite(result)):
            raise ValueError("Equation produced a non-finite value")
        return float(result) if q_array.ndim == 0 else result.copy()

    evaluate.parameter_names = tuple(name for name in names if name != 'q')
    return evaluate


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


def _evaluate_q_input(q, scalar_evaluator, vector_evaluator):
    """Evaluate scalar or array-like q while preserving the input shape."""
    q_array = np.asarray(q, dtype=float)
    if q_array.ndim == 0:
        return float(scalar_evaluator(float(q_array)))

    flat_q = q_array.reshape(-1).tolist()
    values = np.asarray(vector_evaluator(flat_q), dtype=float)
    return values.reshape(q_array.shape)


def evaluate_form_factor(world, name, parameters, q, normalized=True):
    """Evaluate a named subunit or structure form factor over q."""
    if normalized:
        scalar = lambda value: world.EvaluateFormFactor(
            name, parameters, value
        )
        vector = lambda values: world.EvaluateFormFactor(
            name, parameters, values
        )
    else:
        scalar = lambda value: world.EvaluateFormFactorUnnormalized(
            name, parameters, value
        )
        vector = lambda values: world.EvaluateFormFactorUnnormalized(
            name, parameters, values
        )
    return _evaluate_q_input(q, scalar, vector)


def evaluate_form_factor_amplitude(
    world,
    reference,
    parameters,
    q,
    normalized=True,
):
    """Evaluate a form-factor amplitude over scalar or array-like q."""
    if normalized:
        scalar = lambda value: world.EvaluateFormFactorAmplitude(
            reference, parameters, value
        )
        vector = lambda values: world.EvaluateFormFactorAmplitude(
            reference, parameters, values
        )
    else:
        scalar = lambda value: world.EvaluateFormFactorAmplitudeUnnormalized(
            reference, parameters, value
        )
        vector = lambda values: world.EvaluateFormFactorAmplitudeUnnormalized(
            reference, parameters, values
        )
    return _evaluate_q_input(q, scalar, vector)


def evaluate_phase_factor(
    world,
    reference1,
    reference2,
    parameters,
    q,
):
    """Evaluate a phase factor over scalar or array-like q."""
    scalar = lambda value: world.EvaluatePhaseFactor(
        reference1, reference2, parameters, value
    )
    vector = lambda values: world.EvaluatePhaseFactor(
        reference1, reference2, parameters, values
    )
    return _evaluate_q_input(q, scalar, vector)
