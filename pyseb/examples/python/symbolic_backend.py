#!/usr/bin/env python3
"""Construct a symbolic expression directly and convert it to SymPy."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import pyseb


def main():
    print("Available backends:", ", ".join(pyseb.available_backends()))

    pyseb.set_backend("portable")
    x = pyseb.symbol("x")
    expr = pyseb.sin(2.0 * x + 3.0) + pyseb.exp(x / 3.0) + x.bessel_j0()

    sympy_expr = pyseb.to_sympy(expr)
    value = float(sympy_expr.subs("x", 1.25).evalf())

    print("Expression:")
    print(sympy_expr)
    print("LaTeX:")
    print(pyseb.to_latex(pyseb.SymPyExpression(sympy_expr)))
    print(f"Value at x=1.25: {value:.12g}")


if __name__ == "__main__":
    main()
