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
```

These cover SymPy export, numeric evaluation over `q`, a small micelle model,
nested structures, the Figure 13 star-chain model, and direct symbolic
expression construction.

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
