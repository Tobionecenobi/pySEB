#!/usr/bin/env python3
"""Build a nested structure from reusable substructures."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import pyseb
import sympy


def main():
    world = pyseb.World()

    diblock = world.Add("GaussianPolymer", "polyA")
    world.Link("GaussianPolymer", "polyB.end1", "polyA.end2")

    star = world.Add(diblock, "diblock1")
    world.Link(diblock, "diblock2:polyA.end1", "diblock1:polyA.end1")

    chain = world.Add(star, "star1")
    world.Link(star, "star2:diblock1:polyB.end2", "star1:diblock2:polyB.end2")
    world.Add(chain, "nestedchain")

    sympy_form_factor = pyseb.to_sympy(world.FormFactor("nestedchain"))
    params = {
        "beta_polyA": 1.0,
        "beta_polyB": 0.75,
        "Rg_polyA": 5.0,
        "Rg_polyB": 8.0,
    }

    q = sympy.Symbol("q")
    numeric_expr = sympy_form_factor.subs(params)

    print("Nested diblock-star chain form factor")
    print("Symbolic form factor:")
    print(sympy_form_factor)
    for q_value in (0.02, 0.05, 0.1):
        value = float(numeric_expr.subs(q, q_value).evalf())
        print(f"q={q_value:.3f} I(q)={value:.8g}")


if __name__ == "__main__":
    main()
