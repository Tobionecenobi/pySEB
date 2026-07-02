#!/usr/bin/env python3
"""Build and evaluate a small sphere-polymer micelle."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import pyseb
import sympy


def main():
    world = pyseb.World()

    graph_id = world.Add("SolidSphere", "core")
    for index in range(6):
        world.Link(
            "GaussianPolymer",
            f"chain{index}.end1",
            f"core.surface#anchor{index}",
            "chain",
        )

    world.Add(graph_id, "micelle")
    sympy_form_factor = pyseb.to_sympy(world.FormFactor("micelle"))

    params = {
        "beta_core": 1.0,
        "beta_chain": 0.25,
        "R_core": 30.0,
        "Rg_chain": 12.0,
    }

    q = sympy.Symbol("q")
    numeric_expr = sympy_form_factor.subs(params)

    print("Micelle form factor")
    print("Symbolic form factor:")
    print(sympy_form_factor)
    for q_value in (0.02, 0.05, 0.1, 0.2):
        value = float(numeric_expr.subs(q, q_value).evalf())
        print(f"q={q_value:.3f} I(q)={value:.8g}")


if __name__ == "__main__":
    main()
