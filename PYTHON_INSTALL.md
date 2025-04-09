# Building and Installing PySEB

This document provides instructions on how to build and install the Python bindings for the Scattering Equation Builder (SEB) library.

## Prerequisites

Before building PySEB, you need to install the following dependencies:

1. A C++11 compliant compiler
2. [GiNaC](https://www.ginac.de/Download.html) symbolic manipulation library for C++
3. [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)
4. [CMake](https://cmake.org/) (>= 3.18)
5. [Python](https://www.python.org/) (>= 3.7)
6. [pybind11](https://pybind11.readthedocs.io/) (>= 2.6.0)
7. [scikit-build](https://scikit-build.readthedocs.io/) (>= 0.11.0)

## Building from Source

### Step 1: Clone the Repository

```bash
git clone https://github.com/Tobionecenobi/SEB.git
cd SEB
```

### Step 2: Install Python Dependencies

```bash
pip install scikit-build pybind11 numpy matplotlib
```

### Step 3: Build and Install

```bash
pip install .
```

This will build the C++ library and Python bindings, and install the package in your Python environment.

## Development Installation

If you're developing PySEB, you can install it in development mode:

```bash
pip install -e .
```

This will install the package in development mode, so changes to the Python code will be immediately available without reinstalling.

## Building the Wheel

To build a wheel for distribution:

```bash
pip install build
python -m build
```

This will create a wheel file in the `dist` directory, which can be installed with pip.

## Testing the Installation

After installing PySEB, you can test it by running one of the example scripts:

```bash
python examples/diblock_copolymer.py
```

## Troubleshooting

### CMake Can't Find GiNaC or GSL

If CMake can't find GiNaC or GSL, you can specify their locations manually:

```bash
cmake -DCMAKE_PREFIX_PATH=/path/to/ginac;/path/to/gsl ..
```

### Compiler Errors

If you encounter compiler errors, make sure your compiler supports C++11 and that all dependencies are installed correctly.

### Python Import Errors

If you encounter import errors when trying to use PySEB, make sure the package is installed correctly:

```bash
pip list | grep pyseb
```

If it's not listed, try reinstalling:

```bash
pip install .
```

## Package Structure

The PySEB package has the following structure:

```
pyseb/
├── __init__.py       # Package initialization
├── subunits.py       # Subunit classes
├── utils.py          # Utility functions
└── _pyseb.so         # C++ extension module

src/
└── bindings.cpp      # C++ bindings code

examples/
├── diblock_copolymer.py  # Example of a diblock copolymer
└── micelle.py            # Example of a micelle
```

## Building Documentation

To build the documentation:

```bash
pip install sphinx sphinx_rtd_theme
cd docs
make html
```

The documentation will be available in the `docs/build/html` directory.
