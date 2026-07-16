# PySEB Python Guide

PySEB provides Python bindings for SEB. The Python package is a thin layer over
the C++ `seb` core and `seb-symbolic` expression/backend library.

Supported Python versions for packaged wheels are currently 3.9, 3.10, 3.11,
and 3.12. First-release wheels target Linux x86_64, Windows AMD64, and macOS
x86_64 plus Apple Silicon, and CI currently builds and repairs all of these
successfully (macOS x86_64 requires macOS 14+, macOS arm64 requires macOS
15+, both due to Homebrew's bundled GSL). iPhone/iPad iOS wheels are not part
of the first release.

## Install From Source

Create a virtual environment and install the package in editable mode:

```bash
python3 -m venv venv
source venv/bin/activate
python -m pip install -U pip setuptools wheel
python -m pip install -e .
```

The Python packaging path enables the CMake `BUILD_PYTHON` option and builds the
`_pyseb` extension through scikit-build.

Source installs on platforms without a prebuilt wheel require Python 3.9+,
CMake, a C++ compiler, pybind11, GSL, and the platform's Python development
headers/tools.

## Developer Build

To build the Python extension without installing the package:

```bash
make python
```

or directly:

```bash
cmake -S . -B build-python -G Ninja -DBUILD_PYTHON=ON -DSEB_ENABLE_GINAC_BACKEND=OFF
cmake --build build-python
```

When running directly from the repository, make sure Python can find both the
repository package and the freshly built extension:

```bash
PYTHONPATH="$PWD/build-python/pyseb:$PWD" python -c "import pyseb; print(pyseb.available_backends())"
```

## Basic Usage

```python
import numpy as np
import pyseb

world = pyseb.World()
graph_id = world.Add("GaussianPolymer", "poly1")
world.Link("GaussianPolymer", "poly2.end1", "poly1.end2")
world.Add(graph_id, "diblockcopolymer")

form_factor = world.FormFactor("diblockcopolymer")
sympy_expr = pyseb.to_sympy(form_factor)

params = {
    "beta_poly1": 1.0,
    "beta_poly2": 1.0,
    "Rg_poly1": 10.0,
    "Rg_poly2": 20.0,
}

q_values = np.logspace(-3, 0, 100)
intensity = pyseb.evaluate_expression(world, form_factor, params, q_values)
```

## Numerical Subunits

`World` can evaluate structures containing both analytic subunits and numerical
subunits. `DebyeSphereCloud` is the first concrete numerical model:

```python
import pyseb

cloud = pyseb.DebyeSphereCloud([
    pyseb.SphereScatterer(0.0, 0.0, 0.0, 1.0, 2.0),
    pyseb.SphereScatterer(3.0, 0.0, 0.0, 1.0, -0.5),
])
cloud.addReferencePoint("left", 0.0, 0.0, 0.0)

world = pyseb.World()
world.Add(cloud, "cloud")

q_values = [0.0, 0.05, 0.1, 0.2]
form_factor = world.EvaluateFormFactor("cloud", {}, q_values)
rg2 = world.EvaluateRadiusOfGyration2("cloud", {})
```

`NumericalSubunit` also exposes callbacks for form factors, form-factor
amplitudes, phase factors, total excess scattering length, and Guinier data.
This supports geometries whose scattering functions require numerical
integration while preserving the usual symbolic placeholders in symbolic
expressions.

## PDB Atom Clouds

The built-in PDB pipeline is separated into parser, parameter profile, and
cloud builder layers:

```text
PDBParser -> AtomRecord[] -> AtomParameterProfile -> DebyeSphereCloud
```

For example:

```python
import pyseb

profile = pyseb.AtomParameterProfile()
profile.set_element("C", radius=1.70, beta=6.0)
profile.set_element("N", radius=1.55, beta=7.0)
profile.set_element("O", radius=1.52, beta=8.0)
profile.set_element("S", radius=1.80, beta=16.0)

cloud = pyseb.StructureCloudLoader.load_pdb("protein.pdb", profile)
world = pyseb.World()
world.Add(cloud, "protein")

form_factor = world.EvaluateFormFactor("protein", {}, [0.0, 0.05, 0.1])
rg2 = world.EvaluateRadiusOfGyration2("protein", {})
```

Named atom references can be supplied through
`AtomCloudBuildOptions.reference_atom_serials`. They can then be used to attach
analytic subunits, such as connecting `GaussianPolymer.end1` to a terminal
protein atom. The combined structure is evaluated by the same numerical
`World` traversal.

## Examples

Packaged examples live in `pyseb/examples/python/`:

```bash
python pyseb/examples/python/sympy_example.py
python pyseb/examples/python/backend_simplification_comparison.py
python pyseb/examples/python/diblock_evaluation.py
python pyseb/examples/python/figure13_star_chain.py
python pyseb/examples/python/micelle_evaluation.py
python pyseb/examples/python/nested_structures.py
python pyseb/examples/python/symbolic_backend.py
python pyseb/examples/python/debye_sphere_cloud.py
python pyseb/examples/python/crambin_form_factor.py
python pyseb/examples/python/crambin_polymer_form_factor.py
```

These cover SymPy export, numeric evaluation over `q`, a small micelle model,
nested structures, the Figure 13 star-chain model, and direct symbolic
expression construction. The Crambin examples download the official 1CRN PDB
file, verify its checksum, and remove the temporary file after evaluation.

## Symbolic Backends

Python users normally interact with the `sympy` adapter:

```python
pyseb.set_backend("sympy")
```

Internally this selects the portable C++ backend and converts exported
expressions into real SymPy expressions in Python. GiNaC is optional and is not
required for the default Python build.

## Layout

```text
pyseb/
  __init__.py
  symbolic.py
  utils.py
  bindings/          # pybind11 C++ bindings
  examples/python/   # Python examples
  tests/python/      # Python smoke tests
```

## Tests

Run the Python smoke tests against a local extension build:

```bash
make python
PYTHONPATH="$PWD/build-python/pyseb:$PWD" python -m unittest discover -s pyseb/tests/python
```
