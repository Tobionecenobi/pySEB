#!/usr/bin/env python3
"""Build and evaluate a diblock copolymer form factor."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import numpy as np
import pyseb
import sympy


def main():
    world = pyseb.World()

    graph_id = world.Add("GaussianPolymer", "poly1")
    world.Link("GaussianPolymer", "poly2.end1", "poly1.end2")
    world.Add(graph_id, "diblockcopolymer")

    sympy_form_factor = pyseb.to_sympy(world.FormFactor("diblockcopolymer"))
    params = {
        "beta_poly1": 1.0,
        "beta_poly2": 0.8,
        "Rg_poly1": 5.0,
        "Rg_poly2": 8.0,
    }

    q_values = np.array([0.01, 0.05, 0.1, 0.2])
    q = sympy.Symbol("q")
    numeric_expr = sympy_form_factor.subs(params)
    intensity_fn = sympy.lambdify(q, numeric_expr, modules=["numpy"])
    intensity = np.asarray(intensity_fn(q_values), dtype=float)

    print("Diblock copolymer form factor")
    print("Symbolic form factor:")
    print(sympy_form_factor)
    print("q,I(q)")
    for q_value, value in zip(q_values, intensity):
        print(f"{q_value:.4g},{value:.8g}")


if __name__ == "__main__":
    main()
