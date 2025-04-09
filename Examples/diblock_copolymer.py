"""
Example of creating a diblock copolymer and calculating its form factor.
"""

import pyseb
import numpy as np
import matplotlib.pyplot as plt

def main():
    # Create a world
    world = pyseb.World()

    # Add a single polymer, named 'poly1'
    g = world.Add("GaussianPolymer", "poly1")

    # Add a second polymer, named 'poly2'. Poly2's end1 is linked to poly1's end2 forming a diblock copolymer
    world.Link("GaussianPolymer", "poly2.end1", "poly1.end2")

    # Wrap the g structure naming it "diblockcopolymer"
    g2 = world.Add(g, "diblockcopolymer")

    # Request the symbolic expression for the form factor of this structure
    form_factor = world.FormFactor("diblockcopolymer")

    # Print out form factor expression
    print("Form Factor:")
    print(form_factor)
    print()

    # Print out LaTeX expression for form factor
    print("Form Factor (LaTeX):")
    print(pyseb.to_latex(form_factor))
    print()

    # Print out Python expression for form factor
    print("Form Factor (Python):")
    print(pyseb.to_python(form_factor))
    print()

    # Evaluate the form factor for a range of q values
    params = {
        "beta_poly1": 1.0,
        "beta_poly2": 1.0,
        "Rg_poly1": 10.0,
        "Rg_poly2": 20.0
    }

    q_values = np.logspace(-3, 0, 100)
    intensity = pyseb.evaluate_expression(world, form_factor, params, q_values)

    # Plot the form factor
    plt.figure(figsize=(10, 6))
    plt.loglog(q_values, intensity)
    plt.xlabel('q')
    plt.ylabel('I(q)')
    plt.title('Form Factor of Diblock Copolymer')
    plt.grid(True, which="both", ls="-")
    plt.savefig('diblock_copolymer.png')
    plt.show()

if __name__ == "__main__":
    main()
