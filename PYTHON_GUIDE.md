# PySEB: Comprehensive Guide

This guide provides detailed instructions on how to install, use, and extend the Python bindings for the Scattering Equation Builder (SEB) library.

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

With MSYS2:

```bash
pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-gsl mingw-w64-x86_64-cln mingw-w64-x86_64-ginac mingw-w64-x86_64-cmake mingw-w64-x86_64-python
```

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

## Basic Usage

### Creating a World

The first step in using PySEB is to create a `World` object, which is the main container for all structures and subunits:

```python
import pyseb

# Create a world
world = pyseb.World()
```

### Adding Subunits

You can add subunits to the world using the `Add` method:

```python
# Add a Gaussian polymer
g = world.Add("GaussianPolymer", "poly1")

# Add a solid sphere
s = world.Add("SolidSphere", "sphere1")
```

### Linking Subunits

You can link subunits together using the `Link` method:

```python
# Link poly2's end1 to poly1's end2
world.Link("GaussianPolymer", "poly2.end1", "poly1.end2")

# Link poly3's end1 to a point on the sphere's surface
world.Link("GaussianPolymer", "poly3.end1", "sphere1.surface#r0")
```

### Creating Structures

You can create a structure from a collection of linked subunits using the `Add` method with a GraphID:

```python
# Create a structure from the linked subunits
structure = world.Add(g, "my_structure")
```

### Calculating Scattering Properties

You can calculate various scattering properties of a structure:

```python
# Calculate the form factor
form_factor = world.FormFactor("my_structure")

# Calculate the unnormalized form factor
form_factor_unnorm = world.FormFactor_Unnormalized("my_structure")

# Calculate the form factor amplitude relative to a reference point
form_factor_amp = world.FormFactorAmplitude("my_structure:poly1.end1")

# Calculate the phase factor between two reference points
phase_factor = world.PhaseFactor("my_structure:poly1.end1", "my_structure:poly2.end2")
```

### Converting Expressions

You can convert expressions to different formats:

```python
# Convert to Python
python_expr = pyseb.to_python(form_factor)

# Convert to LaTeX
latex_expr = pyseb.to_latex(form_factor)

# Convert to C/C++
c_expr = pyseb.to_cform(form_factor)
```

### Evaluating Expressions

You can evaluate expressions with specific parameter values:

```python
import numpy as np

# Set parameter values
params = {
    "beta_poly1": 1.0,
    "beta_poly2": 1.0,
    "Rg_poly1": 10.0,
    "Rg_poly2": 20.0
}

# Evaluate at a single q value
intensity = pyseb.evaluate_expression(world, form_factor, params, 0.1)

# Evaluate at multiple q values
q_values = np.logspace(-3, 0, 100)
intensity_array = pyseb.evaluate_expression(world, form_factor, params, q_values)
```

### Plotting Results

You can plot the results using matplotlib:

```python
import matplotlib.pyplot as plt

# Plot the form factor
plt.figure(figsize=(10, 6))
plt.loglog(q_values, intensity_array)
plt.xlabel('q')
plt.ylabel('I(q)')
plt.title('Form Factor')
plt.grid(True, which="both", ls="-")
plt.show()
```

## Advanced Usage

### Available Subunits

PySEB provides the following subunits:

- `GaussianPolymer`: A polymer chain with Gaussian statistics
- `GaussianLoop`: A closed loop polymer with Gaussian statistics
- `ThinCircle`: A thin circular ring
- `ThinRod`: A thin rod with uniform scattering density
- `ThinDisk`: A thin disk with uniform scattering density
- `ThinSphericalShell`: A thin spherical shell with uniform scattering density
- `SolidSphere`: A solid sphere with uniform scattering density
- `SolidSphericalShell`: A solid spherical shell with uniform scattering density
- `SolidCylinder`: A solid cylinder with uniform scattering density
- `Point`: A point scatterer

### Reference Points

Each subunit has specific reference points that can be used for linking. Here are some examples:

- `GaussianPolymer`: `end1`, `end2`
- `SolidSphere`: `center`, `surface#r{i}` (where `{i}` is an index)
- `ThinRod`: `end1`, `end2`, `contour#{i}` (where `{i}` is an index)

### Parameters

Each subunit has specific parameters that can be set. Here are some examples:

- `GaussianPolymer`: `Rg_poly{name}` (radius of gyration), `beta_poly{name}` (scattering length density)
- `SolidSphere`: `R_sphere{name}` (radius), `beta_sphere{name}` (scattering length density)

### Error Handling

PySEB will raise exceptions if there are errors in your code. You can catch these exceptions using a try-except block:

```python
try:
    # Your PySEB code here
except Exception as e:
    print(f"Error: {e}")
```

## Examples

### Diblock Copolymer

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

### Micelle

```python
import pyseb
import numpy as np
import matplotlib.pyplot as plt

# Create a world
world = pyseb.World()

# Spherical core
g = world.Add("SolidSphere", "sphere")

# Add polymers to the surface of the sphere
N = 100
for i in range(N):
    name = f"poly{i}.end1"
    ref = f"sphere.surface#r{i}"
    world.Link("GaussianPolymer", name, ref, "poly")

# Define micelle structure
g1 = world.Add(g, "micelle")

# Request the symbolic expression for the form factor of this structure
form_factor = world.FormFactor("micelle")

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
plt.show()
```

## Troubleshooting

### Installation Issues

If you encounter issues during installation, try the following:

1. Make sure all dependencies are installed correctly
2. Check that your compiler supports C++11
3. Try installing from source instead of using pip
4. If you're using Windows, try using WSL or MSYS2

### Runtime Issues

If you encounter issues during runtime, try the following:

1. Check that you're using the correct reference points for linking
2. Make sure all parameters are set correctly
3. Check for typos in structure and subunit names
4. Try simplifying your structure to isolate the issue

## Contributing

If you'd like to contribute to PySEB, please follow these steps:

1. Fork the repository
2. Create a new branch for your feature
3. Add your feature or fix
4. Write tests for your changes
5. Submit a pull request

## License

PySEB is released under the MIT License. See the LICENSE file for more details.
