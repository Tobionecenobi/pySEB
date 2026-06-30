"""
Example of creating a micelle structure and calculating its form factor.
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pyseb
import numpy as np
import matplotlib.pyplot as plt

def main():
    # Create a world
    world = pyseb.World()

    # Spherical core
    g = world.Add("SolidSphere", "sphere")

    # Add polymers to the surface of the sphere
    N = 10
    for i in range(N):
        name = f"poly{i}.end1"
        ref = f"sphere.surface#r{i}"
        world.Link("GaussianPolymer", name, ref, "poly")

    # Define micelle structure
    g1 = world.Add(g, "micelle")

    # Request the symbolic expression for the form factor of this structure
    form_factor = world.FormFactor("micelle")

    # Print out form factor expression
    print("Form Factor:")
    print(form_factor)
    print()

    # Print out LaTeX expression for form factor
    print("Form Factor (LaTeX):")
    print(pyseb.to_latex(form_factor))
    print()

    # Evaluate the form factor for a range of q values
    params = {
        "beta_sphere": 1.0,
        "beta_poly": 0.5,
        "R_sphere": 50.0,
        "Rg_poly": 20.0
    }

    q_values = np.logspace(-3, 0, 100)
    intensity = pyseb.evaluate_expression(world, form_factor, params, q_values)

    # Plot the form factor
    plt.figure(figsize=(10, 6))
    plt.loglog(q_values, intensity)
    plt.xlabel('q')
    plt.ylabel('I(q)')
    plt.title('Form Factor of Micelle')
    plt.grid(True, which="both", ls="-")
    plt.savefig('micelle.png')
    plt.show()

if __name__ == "__main__":
    main()
