"""
SymPy implementation of the symbolic expression interface.
"""

import sympy
import numpy as np
import re
from typing import Dict, Union, Optional

_SYMPY_LOCALS = {
    'besselj': sympy.besselj,
    'dawson': lambda x: sympy.exp(-x**2) * sympy.sqrt(sympy.pi) / 2 * sympy.erfi(x),
    'Integral': sympy.Integral,
    'Si': sympy.Si,
    'hyper': sympy.hyper,
    'gamma': sympy.gamma,
    'struve': sympy.Function('struve'),
    'erf': sympy.erf,
    'erfc': sympy.erfc,
    'erfi': sympy.erfi,
    'Abs': sympy.Abs,
    'sqrt': sympy.sqrt,
    'exp': sympy.exp,
    'log': sympy.log,
    'sin': sympy.sin,
    'cos': sympy.cos,
    'tan': sympy.tan,
    'asin': sympy.asin,
    'acos': sympy.acos,
    'atan': sympy.atan,
    'sinh': sympy.sinh,
    'cosh': sympy.cosh,
    'tanh': sympy.tanh,
    'pi': sympy.pi,
    'E': sympy.E,
    'I': sympy.I,
}

class SymPyExpression:
    """Implementation of the symbolic interface using SymPy"""
    def __init__(self, expr):
        if isinstance(expr, str):
            # Check if this is a GiNaC expression string
            if any(ginac_func in expr for ginac_func in ['besselJ0', 'besselJ1', 'dawsonF']):
                # Convert GiNaC function names to SymPy equivalents
                expr = self._convert_ginac_to_sympy(expr)

            # Parse the expression with SymPy
            try:
                self.expr = sympy.sympify(expr.replace('^', '**'), locals=_SYMPY_LOCALS)
            except Exception as e:
                print(f"Warning: Could not parse expression: {expr}")
                print(f"Error: {e}")
                self.expr = sympy.Symbol(expr)
        elif isinstance(expr, sympy.Expr):
            self.expr = expr
        elif isinstance(expr, SymPyExpression):
            self.expr = expr.expr
        else:
            try:
                self.expr = sympy.sympify(expr)
            except:
                self.expr = expr

    def _convert_ginac_to_sympy(self, expr_str):
        """Convert GiNaC function names to SymPy equivalents"""
        # Function name mappings
        replacements = {
            'besselJ0': 'besselj(0, ',
            'besselJ1': 'besselj(1, ',
            'dawsonF': 'exp(-x**2) * sqrt(pi)/2 * erfi(x)',
            # Add more mappings as needed
        }

        # Apply replacements
        result = expr_str
        for ginac_name, sympy_name in replacements.items():
            if ginac_name == 'dawsonF':
                # Special case for dawson function
                result = re.sub(r'dawsonF\(([^)]+)\)',
                               lambda m: sympy_name.replace('x', m.group(1)),
                               result)
            else:
                # Regular function replacement
                result = result.replace(f"{ginac_name}(", f"{sympy_name}")

        return result

    def add(self, other):
        return SymPyExpression(self.expr + other.expr)

    def subtract(self, other):
        return SymPyExpression(self.expr - other.expr)

    def multiply(self, other):
        return SymPyExpression(self.expr * other.expr)

    def divide(self, other):
        return SymPyExpression(self.expr / other.expr)

    def negate(self):
        return SymPyExpression(-self.expr)

    def pow(self, exponent):
        if isinstance(exponent, (int, float)):
            exp = sympy.Float(exponent)
        else:
            exp = exponent.expr
        return SymPyExpression(self.expr ** exp)

    # Basic functions
    def exp(self):
        return SymPyExpression(sympy.exp(self.expr))

    def log(self):
        return SymPyExpression(sympy.log(self.expr))

    def sin(self):
        return SymPyExpression(sympy.sin(self.expr))

    def cos(self):
        return SymPyExpression(sympy.cos(self.expr))

    def tan(self):
        return SymPyExpression(sympy.tan(self.expr))

    def asin(self):
        return SymPyExpression(sympy.asin(self.expr))

    def acos(self):
        return SymPyExpression(sympy.acos(self.expr))

    def atan(self):
        return SymPyExpression(sympy.atan(self.expr))

    def sinh(self):
        return SymPyExpression(sympy.sinh(self.expr))

    def cosh(self):
        return SymPyExpression(sympy.cosh(self.expr))

    def tanh(self):
        return SymPyExpression(sympy.tanh(self.expr))

    def sqrt(self):
        return SymPyExpression(sympy.sqrt(self.expr))

    def abs(self):
        return SymPyExpression(sympy.Abs(self.expr))

    # Special functions for scattering
    def bessel_j0(self):
        return SymPyExpression(sympy.besselj(0, self.expr))

    def bessel_j1(self):
        return SymPyExpression(sympy.besselj(1, self.expr))

    def dawson(self):
        # Dawson function F(x) = exp(-x²)∫₀ˣexp(t²)dt
        x = self.expr
        return SymPyExpression(sympy.exp(-x**2) * sympy.sqrt(sympy.pi)/2 * sympy.erfi(x))

    def erf(self):
        return SymPyExpression(sympy.erf(self.expr))

    def erfc(self):
        return SymPyExpression(sympy.erfc(self.expr))

    # Substitution and evaluation
    def subs(self, symbol: str, value) -> 'SymPyExpression':
        if isinstance(value, SymPyExpression):
            value = value.expr
        elif isinstance(value, (int, float)):
            value = sympy.Float(value)
        return SymPyExpression(self.expr.subs(symbol, value))

    def eval(self) -> float:
        try:
            result = complex(self.expr.evalf())
            if abs(result.imag) < 1e-10:  # Small imaginary parts are considered numerical error
                return float(result.real)
            return float(abs(result))  # Return magnitude for complex results
        except Exception:
            return float(self.expr)  # Handle special cases

    def is_numeric(self) -> bool:
        return self.expr.is_number

    # Series expansion
    def series(self, var, point, order: int) -> 'SymPyExpression':
        if isinstance(var, SymPyExpression):
            var = var.expr
        if isinstance(point, SymPyExpression):
            point = point.expr
        return SymPyExpression(self.expr.series(var, point, order).removeO())

    def series_to_poly(self, series) -> 'SymPyExpression':
        if isinstance(series, SymPyExpression):
            series = series.expr
        return SymPyExpression(sympy.Poly(series))

    def coeff(self, var, power: int) -> 'SymPyExpression':
        if isinstance(var, SymPyExpression):
            var = var.expr
        return SymPyExpression(self.expr.coeff(var, power))

    def evalf(self) -> 'SymPyExpression':
        return SymPyExpression(self.expr.evalf())

    def is_zero(self) -> bool:
        return self.expr.is_zero

    def to_double(self) -> float:
        return float(self.expr.evalf())

    # String representations
    def __str__(self) -> str:
        return str(self.expr)

    def to_latex(self) -> str:
        return sympy.latex(self.expr)

    def to_python(self) -> str:
        return str(self.expr)

    def to_cform(self) -> str:
        return sympy.ccode(self.expr)

class SymPyFactory:
    """Factory for creating SymPy-based symbolic expressions"""
    @staticmethod
    def createSymbol(name: str) -> SymPyExpression:
        return SymPyExpression(sympy.Symbol(name))

    @staticmethod
    def createConstant(value: float) -> SymPyExpression:
        return SymPyExpression(sympy.Float(value))

    @staticmethod
    def createPi() -> SymPyExpression:
        return SymPyExpression(sympy.pi)

    @staticmethod
    def createE() -> SymPyExpression:
        return SymPyExpression(sympy.E)

    @staticmethod
    def createI() -> SymPyExpression:
        return SymPyExpression(sympy.I)

# Global factory instance
_factory = None

def get_factory():
    global _factory
    if _factory is None:
        _factory = SymPyFactory()
    return _factory

def set_factory(factory):
    global _factory
    _factory = factory

def to_latex(expr) -> str:
    """Convert an expression to LaTeX format."""
    if isinstance(expr, SymPyExpression):
        return expr.to_latex()
    return str(expr)

def to_python(expr) -> str:
    """Convert an expression to Python code."""
    if isinstance(expr, SymPyExpression):
        return expr.to_python()
    return str(expr)

def evaluate_expression(expr, params: Dict[str, float], q_values: np.ndarray) -> np.ndarray:
    """
    Evaluate a symbolic expression for given parameters and q values.

    Args:
        expr: The symbolic expression to evaluate
        params: Dictionary mapping parameter names to values
        q_values: Array of q values to evaluate at

    Returns:
        Array of evaluated values
    """
    if not isinstance(expr, SymPyExpression):
        raise TypeError("Expression must be a SymPyExpression")

    # Create a function that substitutes q and evaluates
    q_symbol = sympy.Symbol('q')
    expr_with_params = expr.expr
    for name, value in params.items():
        expr_with_params = expr_with_params.subs(name, value)

    # Convert to a lambda function for fast evaluation
    func = sympy.lambdify(q_symbol, expr_with_params, modules=['numpy'])

    # Evaluate for all q values
    return func(q_values)
