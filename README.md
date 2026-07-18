# pySEB

[![Documentation](https://github.com/Tobionecenobi/pySEB/actions/workflows/documentation.yml/badge.svg)](https://github.com/Tobionecenobi/pySEB/actions/workflows/documentation.yml)
[![Packaging](https://github.com/Tobionecenobi/pySEB/actions/workflows/packaging.yml/badge.svg)](https://github.com/Tobionecenobi/pySEB/actions/workflows/packaging.yml)

pySEB is the Python interface to the Scattering Equation Builder (SEB), a C++
library for constructing composite structures from reusable subunits and
deriving their small-angle-scattering expressions.

Build structures from polymers, rods, disks, spheres, cylinders, numerical
coordinate clouds, and PDB-derived atom clouds. Connect them at named or
distributed reference points, then obtain symbolic expressions or evaluate
scattering quantities directly.

[Documentation](https://tobionecenobi.github.io/pySEB/) ·
[Tutorials](https://tobionecenobi.github.io/pySEB/tutorials/) ·
[Reference wiki](https://github.com/Tobionecenobi/pySEB/wiki) ·
[Build guides](https://github.com/Tobionecenobi/pySEB/wiki/Build-Guides)

## Install

Create a virtual environment and install the Python package:

```bash
python3 -m venv .venv-pyseb
source .venv-pyseb/bin/activate
python -m pip install --upgrade pip
python -m pip install pyseb
```

Windows PowerShell activation is documented in the
[installation guides](https://github.com/Tobionecenobi/pySEB/wiki/Build-Guides).

## First model

This example joins two Gaussian polymers and evaluates their normalized form
factor:

```python
import pyseb

world = pyseb.World("FirstModel")
graph = world.Add("GaussianPolymer", "blockA")
world.Link("GaussianPolymer", "blockB.end1", "blockA.end2")
world.Add(graph, "diblock")

parameters = {
    "beta_blockA": 1.0,
    "beta_blockB": 0.8,
    "Rg_blockA": 5.0,
    "Rg_blockB": 8.0,
}
q = [0.001, 0.01, 0.1]
print(world.EvaluateFormFactor("diblock", parameters, q))
```

The complete, executable Python and C++ walkthrough is the
[diblock tutorial](https://tobionecenobi.github.io/pySEB/tutorials/diblock-copolymer.html).

## Documentation map

| Need | Documentation |
| --- | --- |
| Learn by building models | [Executable tutorials](https://tobionecenobi.github.io/pySEB/tutorials/) |
| Understand worlds, graphs, names, and links | [Core concepts](https://github.com/Tobionecenobi/pySEB/wiki/Core-Concepts) |
| Browse supported models | [Subunit catalogue](https://github.com/Tobionecenobi/pySEB/wiki/Subunit-Catalogue) |
| Use Python or C++ APIs | [Python API](https://github.com/Tobionecenobi/pySEB/wiki/Python-API-Reference) · [C++ API](https://github.com/Tobionecenobi/pySEB/wiki/Cpp-API-Reference) |
| Install or build from source | [Linux](https://github.com/Tobionecenobi/pySEB/wiki/Building-on-Linux) · [macOS](https://github.com/Tobionecenobi/pySEB/wiki/Building-on-macOS) · [Windows](https://github.com/Tobionecenobi/pySEB/wiki/Building-on-Windows) |
| Understand the formalism | [Theory](https://github.com/Tobionecenobi/pySEB/wiki/Theory-Overview) · [Publications](https://github.com/Tobionecenobi/pySEB/wiki/Publications-and-Citation) |
| Contribute code | [Repository architecture](https://github.com/Tobionecenobi/pySEB/wiki/Repository-Architecture) · [Adding a subunit](https://github.com/Tobionecenobi/pySEB/wiki/Adding-a-Subunit) |

## Interfaces and backends

- Portable C++14 symbolic construction and numerical evaluation.
- Optional GiNaC symbolic backend for C++.
- Python 3.9–3.12 bindings with native SymPy conversion.
- Export to Python, SymPy, LaTeX, and C-oriented representations.
- Analytic, callback-defined numerical, Debye-cloud, and PDB workflows.

See the [symbolic backend guide](https://github.com/Tobionecenobi/pySEB/wiki/Symbolic-Backends)
for capabilities and selection.

## Example structures

| Linear chains | Dendrimers |
| --- | --- |
| ![Linear chains assembled from different subunits](https://raw.githubusercontent.com/Tobionecenobi/pySEB/main/PaperFigs/fig10.png) | ![Recursively assembled dendrimers](https://raw.githubusercontent.com/Tobionecenobi/pySEB/main/PaperFigs/fig11.png) |
| Surface-linked structures | Diblock-star chain |
| ![Rods and polymers attached to geometric surfaces](https://raw.githubusercontent.com/Tobionecenobi/pySEB/main/PaperFigs/fig12.png) | ![Five linked diblock-copolymer stars](https://raw.githubusercontent.com/Tobionecenobi/pySEB/main/PaperFigs/fig13.png) |

The [Pages overview](https://tobionecenobi.github.io/pySEB/) explains these
figures and links them to current executable tutorials and reference material.

## Development

The canonical build uses CMake 3.18+, C++14, Ninja, and GSL. Development build,
test, validation, release, and contribution procedures live in the wiki:

- [Project status and roadmap](https://github.com/Tobionecenobi/pySEB/wiki/Project-Status-and-Roadmap)
- [Testing and validation](https://github.com/Tobionecenobi/pySEB/wiki/Testing-and-Validation)
- [Repository architecture](https://github.com/Tobionecenobi/pySEB/wiki/Repository-Architecture)
- [Troubleshooting](https://github.com/Tobionecenobi/pySEB/wiki/Troubleshooting-and-FAQ)

## Citation and license

Use the copyable BibTeX and DOI links in
[Publications and Citation](https://github.com/Tobionecenobi/pySEB/wiki/Publications-and-Citation).
pySEB is distributed under the [GPL-3.0 license](LICENSE).
