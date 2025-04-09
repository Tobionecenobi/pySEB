"""
SymPy implementation of the symbolic expression interface.
"""

import sympy
import numpy as np
from typing import Dict, List, Union, Optional, Any


class SymPyExpression:
    """
    SymPy implementation of symbolic expressions.
    """

    def __init__(self, expr):
        """
        Initialize with a SymPy expression.

        Parameters
        ----------
        expr : sympy.Expr
            The SymPy expression.
        """
        self._expr = expr

    @property
    def expr(self):
        """Get the underlying SymPy expression."""
        return self._expr

    def add(self, other):
        """Add another expression."""
        return SymPyExpression(self._expr + other.expr)

    def subtract(self, other):
        """Subtract another expression."""
        return SymPyExpression(self._expr - other.expr)

    def multiply(self, other):
        """Multiply by another expression."""
        return SymPyExpression(self._expr * other.expr)

    def divide(self, other):
        """Divide by another expression."""
        return SymPyExpression(self._expr / other.expr)

    def negate(self):
        """Negate the expression."""
        return SymPyExpression(-self._expr)

    def pow(self, exponent):
        """Raise to a power."""
        return SymPyExpression(self._expr ** exponent.expr)

    def exp(self):
        """Exponential function."""
        return SymPyExpression(sympy.exp(self._expr))

    def log(self):
        """Natural logarithm."""
        return SymPyExpression(sympy.log(self._expr))

    def sin(self):
        """Sine function."""
        return SymPyExpression(sympy.sin(self._expr))

    def cos(self):
        """Cosine function."""
        return SymPyExpression(sympy.cos(self._expr))

    def tan(self):
        """Tangent function."""
        return SymPyExpression(sympy.tan(self._expr))

    def asin(self):
        """Inverse sine function."""
        return SymPyExpression(sympy.asin(self._expr))

    def acos(self):
        """Inverse cosine function."""
        return SymPyExpression(sympy.acos(self._expr))

    def atan(self):
        """Inverse tangent function."""
        return SymPyExpression(sympy.atan(self._expr))

    def sinh(self):
        """Hyperbolic sine function."""
        return SymPyExpression(sympy.sinh(self._expr))

    def cosh(self):
        """Hyperbolic cosine function."""
        return SymPyExpression(sympy.cosh(self._expr))

    def tanh(self):
        """Hyperbolic tangent function."""
        return SymPyExpression(sympy.tanh(self._expr))

    def sqrt(self):
        """Square root function."""
        return SymPyExpression(sympy.sqrt(self._expr))

    def abs(self):
        """Absolute value function."""
        return SymPyExpression(sympy.Abs(self._expr))

    def bessel_j0(self):
        """Bessel function of the first kind, order 0."""
        return SymPyExpression(sympy.besselj(0, self._expr))

    def bessel_j1(self):
        """Bessel function of the first kind, order 1."""
        return SymPyExpression(sympy.besselj(1, self._expr))

    def dawson(self):
        """Dawson function."""
        # SymPy doesn't have a built-in Dawson function, so we use the definition
        # dawson(x) = (sqrt(pi)/2) * exp(-x^2) * erfi(x)
        x = self._expr
        result = sympy.sqrt(sympy.pi) / 2 * sympy.exp(-x**2) * sympy.erfi(x)
        return SymPyExpression(result)

    def erf(self):
        """Error function."""
        return SymPyExpression(sympy.erf(self._expr))

    def erfc(self):
        """Complementary error function."""
        return SymPyExpression(sympy.erfc(self._expr))

    def subs(self, symbol, value=None):
        """
        Substitute a symbol with a value.

        Parameters
        ----------
        symbol : str or dict
            The symbol to substitute, or a dictionary of substitutions.
        value : SymPyExpression, optional
            The value to substitute with. Not used if symbol is a dictionary.

        Returns
        -------
        SymPyExpression
            The expression with substitutions applied.
        """
        if isinstance(symbol, dict):
            # Convert values to SymPy expressions
            subs_dict = {sympy.Symbol(k): v for k, v in symbol.items()}
            return SymPyExpression(self._expr.subs(subs_dict))
        else:
            return SymPyExpression(self._expr.subs(sympy.Symbol(symbol), value.expr))

    def eval(self):
        """
        Evaluate the expression to a numerical value.

        Returns
        -------
        float
            The numerical value of the expression.

        Raises
        ------
        ValueError
            If the expression cannot be evaluated to a numerical value.
        """
        if not self.is_numeric():
            raise ValueError("Cannot evaluate non-numeric expression")
        return float(self._expr)

    def is_numeric(self):
        """
        Check if the expression is a numerical value.

        Returns
        -------
        bool
            True if the expression is a numerical value, False otherwise.
        """
        return self._expr.is_number

    def series(self, var, point, order):
        """
        Compute the series expansion of the expression around a point.

        Parameters
        ----------
        var : SymPyExpression
            The variable to expand around.
        point : SymPyExpression
            The point to expand around.
        order : int
            The order of the expansion.

        Returns
        -------
        SymPyExpression
            The series expansion.
        """
        return SymPyExpression(self._expr.series(var.expr, point.expr, order))

    def series_to_poly(self, series):
        """
        Convert a series expansion to a polynomial.

        Parameters
        ----------
        series : SymPyExpression
            The series expansion.

        Returns
        -------
        SymPyExpression
            The polynomial.
        """
        # In SymPy, we can use the removeO method to remove the order term
        return SymPyExpression(series.expr.removeO())

    def coeff(self, var, power):
        """
        Extract the coefficient of a variable raised to a power.

        Parameters
        ----------
        var : SymPyExpression
            The variable.
        power : int
            The power.

        Returns
        -------
        SymPyExpression
            The coefficient.
        """
        return SymPyExpression(self._expr.coeff(var.expr, power))

    def evalf(self):
        """
        Evaluate the expression numerically.

        Returns
        -------
        SymPyExpression
            The numerically evaluated expression.
        """
        return SymPyExpression(self._expr.evalf())

    def is_zero(self):
        """
        Check if the expression is zero.

        Returns
        -------
        bool
            True if the expression is zero, False otherwise.
        """
        return self._expr.is_zero

    def to_double(self):
        """
        Convert the expression to a double.

        Returns
        -------
        float
            The expression as a double.

        Raises
        ------
        ValueError
            If the expression cannot be converted to a double.
        """
        if not self.is_numeric():
            raise ValueError("Cannot convert non-numeric expression to double")
        return float(self._expr)

    def to_string(self):
        """
        Convert the expression to a string.

        Returns
        -------
        str
            The string representation of the expression.
        """
        return str(self._expr)

    def to_latex(self):
        """
        Convert the expression to a LaTeX string.

        Returns
        -------
        str
            The LaTeX representation of the expression.
        """
        return sympy.latex(self._expr)

    def to_python(self):
        """
        Convert the expression to a Python string.

        Returns
        -------
        str
            The Python representation of the expression.
        """
        return sympy.python(self._expr)

    def to_cform(self):
        """
        Convert the expression to a C/C++ string.

        Returns
        -------
        str
            The C/C++ representation of the expression.
        """
        return sympy.ccode(self._expr)

    def __str__(self):
        """String representation."""
        return self.to_string()

    def __repr__(self):
        """Representation for debugging."""
        return f"SymPyExpression({self._expr})"


class SymPyFactory:
    """
    Factory for creating SymPy expressions.
    """

    @staticmethod
    def create_symbol(name):
        """
        Create a symbolic variable.

        Parameters
        ----------
        name : str
            The name of the variable.

        Returns
        -------
        SymPyExpression
            The symbolic variable.
        """
        return SymPyExpression(sympy.Symbol(name))

    @staticmethod
    def create_constant(value):
        """
        Create a constant value.

        Parameters
        ----------
        value : float
            The constant value.

        Returns
        -------
        SymPyExpression
            The constant value as a symbolic expression.
        """
        return SymPyExpression(sympy.Float(value))

    @staticmethod
    def create_pi():
        """
        Create the constant pi.

        Returns
        -------
        SymPyExpression
            The constant pi as a symbolic expression.
        """
        return SymPyExpression(sympy.pi)

    @staticmethod
    def create_e():
        """
        Create the constant e.

        Returns
        -------
        SymPyExpression
            The constant e as a symbolic expression.
        """
        return SymPyExpression(sympy.E)

    @staticmethod
    def create_i():
        """
        Create the imaginary unit i.

        Returns
        -------
        SymPyExpression
            The imaginary unit i as a symbolic expression.
        """
        return SymPyExpression(sympy.I)


# Convenience functions
def symbol(name):
    """Create a symbolic variable."""
    return SymPyFactory.create_symbol(name)

def constant(value):
    """Create a constant value."""
    return SymPyFactory.create_constant(value)

def pi():
    """Create the constant pi."""
    return SymPyFactory.create_pi()

def e():
    """Create the constant e."""
    return SymPyFactory.create_e()

def i():
    """Create the imaginary unit i."""
    return SymPyFactory.create_i()

def evaluate_expression(expr, params=None, q_values=None):
    """
    Evaluate a symbolic expression with the given parameters.

    Parameters
    ----------
    expr : SymPyExpression
        The expression to evaluate.
    params : dict, optional
        Dictionary of parameter names and values.
    q_values : array_like, optional
        If provided, evaluate the expression at each q value.

    Returns
    -------
    float or ndarray
        The evaluated expression. If q_values is provided, returns an array
        of the same shape as q_values.
    """
    if params is None:
        params = {}

    if q_values is None:
        # Evaluate at a single point
        subs_expr = expr.subs(params)
        return subs_expr.eval()
    else:
        # Evaluate at multiple q values
        q_array = np.asarray(q_values)
        q_symbol = symbol('q')

        # Create a lambda function for fast evaluation
        subs_expr = expr.subs(params)
        q_func = sympy.lambdify(sympy.Symbol('q'), subs_expr.expr, 'numpy')

        return q_func(q_array)
