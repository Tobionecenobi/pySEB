#!/usr/bin/env python3
"""Show a SEB expression that contains a symbolic integral."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import pyseb
import sympy
import numpy as np
from scipy import integrate


def main():
    pyseb.set_backend("portable")

    world = pyseb.World()
    world.Add("SolidCylinder", "cylinder")

    form_factor = pyseb.to_sympy(world.FormFactor("cylinder"))
    integrals = sorted(form_factor.atoms(sympy.Integral), key=str)

    print("Solid cylinder symbolic integral")
    print("Symbolic form factor:")
    print(form_factor)
    print("\nIntegral terms:")
    for integral in integrals:
        print(integral)

    print("\nSymPy keeps these as Integral objects.")
    print("Use .doit(), lambdify, or a numeric integration strategy when you want values.")

    # --- Numeric evaluation via lambdify + scipy.integrate.quad ---
    # Physical parameters
    params = {
        "beta_cylinder": 1.0,
        "R_cylinder": 30.0,   # radius in Å
        "L_cylinder": 100.0,  # length in Å
    }

    # There is one integral; extract integrand and integration variable
    (the_integral,) = integrals
    integrand_expr, (t_var, t_lo, t_hi) = the_integral.function, the_integral.limits[0]

    q_sym = sympy.Symbol("q")

    # Substitute physical parameters into the integrand
    integrand_numeric = integrand_expr.subs(params)

    # Build a fast numerical integrand: f(t, q)
    integrand_fn = sympy.lambdify(
        (t_var, q_sym),
        integrand_numeric,
        modules=["scipy", "numpy"],
    )

    def evaluate_ff(q_val):
        result, _ = integrate.quad(integrand_fn, float(t_lo), float(t_hi), args=(q_val,))
        return result

    q_values = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2]
    print(f"\nNumerical form factor  R={params['R_cylinder']} Å  L={params['L_cylinder']} Å")
    print(f"{'q (1/Å)':>10}  {'I(q)':>14}")
    print("-" * 28)
    for q_val in q_values:
        print(f"{q_val:>10.4f}  {evaluate_ff(q_val):>14.6g}")


if __name__ == "__main__":
    main()

