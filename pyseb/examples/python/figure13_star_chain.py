#!/usr/bin/env python3
"""Build the Figure 13 chain of diblock-copolymer stars."""

from pathlib import Path
import os
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "pyseb-matplotlib"))

import matplotlib
import numpy as np
import pyseb
import sympy

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    world = pyseb.World()

    diblock = world.Add("GaussianPolymer", "A")
    world.Link("GaussianPolymer", "B.end1", "A.end2")

    star = world.Add(diblock, "diblock1")
    world.Link(diblock, "diblock2:A.end1", "diblock1:A.end1")
    world.Link(diblock, "diblock3:A.end1", "diblock1:A.end1")
    world.Link(diblock, "diblock4:A.end1", "diblock1:A.end1")

    chain = world.Add(star, "star1")
    world.Link(star, "star2:diblock1:B.end2", "star1:diblock3:B.end2")
    world.Link(star, "star3:diblock1:B.end2", "star2:diblock3:B.end2")
    world.Link(star, "star4:diblock1:B.end2", "star3:diblock3:B.end2")
    world.Link(star, "star5:diblock1:B.end2", "star4:diblock3:B.end2")
    world.Add(chain, "chain")

    sympy_form_factor = pyseb.to_sympy(world.FormFactor("chain"))
    params = {
        "Rg_A": 1.0,
        "Rg_B": 1.0,
        "beta_A": 1.0,
        "beta_B": 1.0,
    }

    q_values = np.array([0.01, 0.05, 0.1, 0.5, 1.0])
    plot_q_values = np.logspace(-2, 1.7, 400)
    q = sympy.Symbol("q")
    numeric_expr = sympy_form_factor.subs(params)
    intensity_fn = sympy.lambdify(q, numeric_expr, modules=["numpy"])
    intensity = np.asarray(intensity_fn(q_values), dtype=float)
    plot_intensity = np.asarray(intensity_fn(plot_q_values), dtype=float)

    print("Figure 13 star-chain form factor")
    print("Structure: five stars, four diblocks per star")
    print("Symbolic form factor:")
    print(sympy_form_factor)
    print("q,I(q)")
    for q_value, value in zip(q_values, intensity):
        print(f"{q_value:.4g},{value:.8g}")

    fig, ax = plt.subplots()
    ax.loglog(plot_q_values, plot_intensity)
    ax.set_xlabel("q")
    ax.set_ylabel("I(q)")
    ax.set_title("Figure 13 star-chain form factor")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    #plt.savefig("figure13_star_chain.png", dpi=300)


if __name__ == "__main__":
    main()
