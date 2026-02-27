# PySEB: Python Bindings for the Scattering Equation Builder

PySEB provides Python bindings for the [Scattering Equation Builder (SEB)](https://github.com/Tobionecenobi/SEB) library, a C++ library that allows you to build structures from sub-units such as spheres, polymers, rods, etc., and obtain their scattering properties symbolically.

## Installation

### Prerequisites

Before installing PySEB, you need to install the following dependencies:

1. A C++11 compliant compiler
2. [GiNaC](https://www.ginac.de/Download.html) symbolic manipulation library for C++
3. [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)
4. [CMake](https://cmake.org/) (>= 3.18)
5. [Python](https://www.python.org/) (>= 3.7)

#### Ubuntu/Debian

```bash
sudo apt update
sudo apt install libcln-dev libgsl-dev libginac-dev cmake python3-dev
```

#### macOS

```bash
brew install cln gsl ginac cmake python
```

#### Windows

For Windows, we recommend using the Windows Subsystem for Linux (WSL) or installing the dependencies using [MSYS2](https://www.msys2.org/).

### Installing PySEB

Once you have installed the prerequisites, you can install PySEB using pip:

```bash
pip install pyseb
```

Or from source:

```bash
git clone https://github.com/Tobionecenobi/SEB.git
cd SEB
pip install .
```

## Usage

### Evaluation API note

In Python, numerical evaluation is intentionally handled by SymPy utilities.
Use `pyseb.evaluate_expression(...)` for substitutions and numeric evaluation.

- Preferred: `pyseb.evaluate_expression(world, expr, params, q_values)`
- Avoid relying on `World.Evaluate(...)` from Python code

This keeps C++ focused on building symbolic expressions while Python/SymPy handles evaluation.

Here's a simple example of how to use PySEB to create a diblock copolymer and calculate its form factor:

```python
import pyseb
import numpy as np
import matplotlib.pyplot as plt

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
print(form_factor)

# Print out LaTeX expression for form factor
print(pyseb.to_latex(form_factor))

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
plt.show()
```

## Documentation

For more information on how to use PySEB, please refer to the [SEB documentation](https://github.com/Tobionecenobi/SEB/wiki).

## License

PySEB is released under the MIT License. See the LICENSE file for more details.
