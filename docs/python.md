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

## Executable Tutorial

The Quarto/Jupyter tutorial
[`docs/tutorials/pyseb_tutorial.qmd`](tutorials/pyseb_tutorial.qmd) covers
symbolic polymers, numerical callbacks, Debye sphere clouds, PDB import, and a
mixed protein--polymer structure.

Render or preview it with:

```bash
PYTHONPATH="$PWD/build-python:$PWD" \
  quarto preview docs/tutorials/pyseb_tutorial.qmd
```

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

## File-defined analytic subunits

Create and use a custom analytic model without rebuilding the extension:

```python
world = pyseb.World("custom")
world.add_subunit_file("MyModel.pyseb.yaml", name="particle")
```

Models can instead be registered explicitly by file or directory and then
constructed through their namespaced ID. See the
[file-defined subunit guide](file-defined-subunits.qmd) for the schema,
expression grammar, validation cases, C++ API, and registry rules.

## World import and export

Complete worlds can be saved as versioned `.pyseb-world.yaml` files. The file
stores subunit instances, reusable structures, graph membership, links, and
root graphs; numerical evaluation parameters remain separate.

```python
world.save("model.pyseb-world.yaml")
loaded = pyseb.load_world("model.pyseb-world.yaml")

yaml_text = pyseb.world_to_yaml(world)
loaded = pyseb.world_from_yaml(yaml_text)
```

Bundled catalogue models are available automatically. External file-defined
models can be supplied while loading:

```python
loaded = pyseb.load_world(
    "model.pyseb-world.yaml",
    model_files=["CustomModel.pyseb.yaml"],
)
```

The command-line helpers validate or print the canonical representation:

```bash
pyseb-world validate model.pyseb-world.yaml
pyseb-world inspect model.pyseb-world.yaml
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

`SolidCylinder` and the integral-backed terms of `ThinDisk` are YAML catalogue
models evaluated through the same `World` API. Their symbolic equations retain
explicit integrals, while numerical evaluation uses the GSL method declared in
the model file. NumPy-friendly helpers accept scalar or array-like `q` values:

~~~python
import numpy as np
import pyseb

world = pyseb.World()
world.Add("SolidCylinder", "cylinder")
parameters = {
    "beta_cylinder": 1.0,
    "R_cylinder": 30.0,
    "L_cylinder": 100.0,
}
q = np.logspace(-3, -1, 100)

form_factor = pyseb.evaluate_form_factor(
    world, "cylinder", parameters, q
)
amplitude = pyseb.evaluate_form_factor_amplitude(
    world, "cylinder.center", parameters, q
)
phase = pyseb.evaluate_phase_factor(
    world,
    "cylinder.hull",
    "cylinder.ends",
    parameters,
    q,
)
~~~

Contributors can express finite one-dimensional integrals with named YAML
`integrals` blocks. Multidimensional, nested, improper, data-driven, or
specialized integration remains available through `IntegratedSubunit` and
normalized numerical callbacks. Missing callbacks automatically fall back to
the analytical `SubUnit` implementation, so `World` requires no subunit-type
checks.

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

## Atomic solution-SAXS clouds

`AtomicSAXSCloud` replaces constant point weights with q-dependent
Cromer--Mann atomic form factors and can subtract a Fraser-style Gaussian for
the solvent displaced by every atomic group:

```python
import pyseb

profile = pyseb.AtomicSAXSParameterProfile.protein_defaults()

options = pyseb.AtomicSAXSOptions()
options.solvent_electron_density = 0.334  # electrons / A^3
options.excluded_volume_scale = 1.0
options.include_excluded_solvent = True

references = pyseb.AtomicSAXSCloudBuildOptions()
references.use_implicit_hydrogens = True
references.calculate_solvent_accessibility = True
references.sasa_probe_radius = 1.4       # A
references.sasa_point_count = 960
references.reference_atom_serials = {
    "nterminus": 1,
    "cterminus": 684,
}

cloud = pyseb.AtomicSAXSStructureLoader.load_pdb(
    "protein.pdb",
    profile,
    pyseb.StructureParseOptions(),
    options,
    references,
)
world = pyseb.World("solution-saxs")
world.Add(cloud, "protein")
intensity = world.EvaluateFormFactor(
    "protein", {}, [0.0, 0.05, 0.1, 0.2]
)
```

The bundled protein profile supplies neutral-atom Cromer--Mann coefficients
for H, C, N, O, P, and S and places standard implicit protein hydrogens on
their bonded heavy-atom groups. Their contribution is spherically averaged
as `n*f_H(q)*sinc(q*d)`, using the C-H, N-H, and O-H distances in Table 1
of [Grudinin et al. (2017)](https://doi.org/10.1107/S2059798317005745),
equation (8). The profile's `implicit_hydrogen_distance` is in Angstroms;
zero recovers coincident hydrogens. That table supplies no S-H distance,
so thiols retain the coincident approximation. Explicit hydrogens retain
their actual coordinates and should not also be counted implicitly.
`protein_defaults()` uses the Fraser 1978 displaced
atomic-volume convention. Select the Pontius 1996 protein-core group means or
the legacy Bondi element spheres explicitly:

```python
fraser = pyseb.AtomicSAXSParameterProfile.protein_defaults()
pontius = pyseb.AtomicSAXSParameterProfile.protein_defaults(
    pyseb.AtomicSAXSVolumeModel.Pontius1996
)
bondi = pyseb.AtomicSAXSParameterProfile.protein_defaults(
    pyseb.AtomicSAXSVolumeModel.BondiElement
)
```

The 20 standard amino acids are resolved by residue and atom name so their
implicit hydrogens use the corresponding chemical group. The Pontius source
does not tabulate every group in the compact comparison table; carbonyl C, P,
S, and isolated-H fallbacks therefore remain explicitly Fraser-valued.
Unsupported atom names resolve only when the profile has an explicit element
rule. Add custom atoms with `set_atom`; their `sasa_radius` is used when the
FreeSASA protein classifier cannot assign a radius.

If the parsed structure contains explicit hydrogen records, set
`build_options.use_implicit_hydrogens = False` to avoid counting hydrogens
twice. This also selects the bundled profile's bare-atom excluded volumes,
so explicit H atoms contribute their own displaced solvent only once.
Custom grouped parameters can declare `atomic_excluded_volume` for this mode;
the default `-1` keeps their supplied `excluded_volume` unchanged.
The bundled group map represents standard internal amino-acid groups;
terminal protonation and nonstandard residues should be supplied explicitly
when they matter for the construct.

The implementation is self-contained: FoXS is neither imported nor executed.
The builder uses the bundled FreeSASA 2.1.3 ProtOr radii and one-thread
Shrake--Rupley calculation by default. Positive-occupancy atoms are full
geometric occluders; zero-occupancy atoms are omitted. Each atomic SASA is
stored as the exposure fraction
$A_i/[4\pi(r_i+r_{probe})^2]$. Explicit hydrogen records participate when
present, while implicit hydrogens remain part of their bonded heavy-atom
group.

Hydration remains opt-in. It uses the stored exposure, 30 A3 of hydration
volume per fully exposed group, the requested shell-density contrast in
electrons/A3, and a normalized water form factor for q dependence:

```python
options.include_hydration = True
options.hydration_shell_contrast = 0.0135  # electrons / A^3
```

The default `AtomCentered` model is an exposure-weighted approximation: it has no
separate spatial shell coordinates. It is not equivalent to Pepsi-SAXS's
surrounding hydration shell, and the same density contrast is not a calibrated
mapping between the two models. A fixed-coordinate comparison on titin and
crambin found roughly 0.6–1.9% RMS differences in normalized solution scattering
without hydration, increasing to 7.5–9.5% at 5% shell contrast. See the
[comparison report](design/atomic-saxs-pepsi.md) and reproducible benchmark
before interpreting hydration parameters as physical shell properties.

### Spatial hydration shell

Select a separate spatial shell when building the cloud:

```python
options.hydration_model = pyseb.AtomicSAXSHydrationModel.SpatialShell
options.include_hydration = True
options.hydration_shell_contrast = 0.0167  # excess electrons / A^3
build = pyseb.AtomicSAXSCloudBuildOptions()
build.hydration_shell_inner_radius = 3.0   # A from any active atomic center
build.hydration_shell_thickness = 3.0      # A
build.hydration_shell_spacing = 3.0        # A; refine to check convergence
cloud = pyseb.AtomicSAXSStructureLoader.load_pdb(
    "protein.pdb", profile, pyseb.StructureParseOptions(), options, build
)
```

Grid points lie outside every inner-radius atomic sphere and within an
outer-radius sphere of at least one atom. Each retained point represents one
voxel of volume `spacing**3`, counted once even where outer spheres overlap.
Positive-occupancy atoms are full geometric occluders; zero-occupancy atoms
are omitted from shell construction. The grid is centered on the active
atomic centroid, so translating the input translates the same grid. Grid
orientation can affect finite-resolution results; refine spacing to assess
discretization error. Enclosed cavities that satisfy the same distance rule
are included; the model does not classify solvent connectivity.

A point contributes `contrast * voxel_volume * normalized_water_form_factor(q)`.
The shell participates in atom-shell and shell-shell interference and in
attachment amplitudes and Guinier moments. Atomic SASA and
`hydration_volume_per_exposed_atom` are not used by this model. Shell points
are included in `get_scatterers()`/`scatterer_count()` and have positive
`hydration_volume`; atomic entries have zero. The reserved `center` reference
remains the atomic centroid, independent of shell sampling.

Choosing `SpatialShell` builds points even when `include_hydration` is false,
allowing contrast and hydration toggles without rebuilding distances. Switching
an atom-only cloud to enabled spatial hydration raises an error: rebuild with
`SpatialShell` first. Changing shell radii or spacing also requires rebuilding.
Custom spatial quadrature can be supplied as zero-vacuum/zero-excluded-volume
scatterers with positive `hydration_volume`. OpenUSD records that volume and
the hydration model alongside the coordinates.

This is a volume-weighted shell, not a calibrated reproduction of Pepsi's
shell normalization. See the [validation and comparison report](design/spatial-hydration.md).

### Reference hydration conventions

The staged [hydration reconciliation](design/hydration-reconciliation.md)
separates shell geometry, water weights, and surface sampling. To reproduce
the Pepsi 3 convention without fitting a molecular scale:

```python
options.hydration_model = pyseb.AtomicSAXSHydrationModel.SpatialShell
options.hydration_water_model = pyseb.AtomicSAXSHydrationWaterModel.BondedWater
options.include_hydration = True
options.hydration_shell_contrast = 0.0167  # Pepsi --dro 5 at bulk density 0.334
build = pyseb.AtomicSAXSCloudBuildOptions()
build.hydration_shell_convention = pyseb.AtomicSAXSHydrationShellConvention.Pepsi3
cloud = pyseb.AtomicSAXSStructureLoader.load_pdb(
    "protein.pdb", profile, pyseb.StructureParseOptions(), options, build
)
```

This convention uses a bounding-box grid with half-cell offset and the
profile's atomic surface radii in the inner exclusion boundary. Point weights
include the current bulk-density factor during evaluation; `hydration_volume`
stores geometric voxel volumes, and `hydration_uses_solvent_density` marks
Pepsi shell points. Changing solvent density therefore agrees with rebuilding.
The `VolumeDensity` default preserves the physical volume-weighted interpretation.
`BondedWater` adds the finite O-H distance to the water form factor and its
Guinier moment. The remaining small discrepancy with Pepsi uses different
four- versus five-Gaussian atomic factors.

FoXS uses an atom-centered hydration contribution with a *solution-water*
form factor. Select `hydration_water_model =
pyseb.AtomicSAXSHydrationWaterModel.FoXSReciprocal` for its full q-dependent
water factor, or `FoXSFast` for its shared Gaussian approximation. The mapping
is `hydration_shell_contrast * hydration_volume_per_exposed_atom = c2 * fw(0)`;
at default density, `fw(0)` is about 3.50968 electrons. This specifies hydration
weights only, not FoXS's protein form factors or displaced-solvent convention.
The validation adapter uses native exposure fractions as one control and
independent FreeSASA fractions as another.

For surface geometry consistent with FoXS, a profile may supply radii
`(3*excluded_volume/(4*pi))**(1/3)` with `build.sasa_use_profile_radii = True`
and `build.sasa_probe_radius = 1.8`. This bypasses the usual ProtOr classifier;
the default remains ProtOr with a 1.4 A probe. Independent surface quadrature
need not match FoXS's coarser dot sampling exactly, even with identical radii.

SASA is calculated once while constructing the cloud. Changing displaced
solvent or hydration parameters reuses both the exposure and pair-distance
caches.

The solution parameters can be changed without rebuilding pair distances:

```python
cloud.set_solvent_electron_density(0.330)
cloud.set_excluded_volume_scale(0.98)
cloud.set_include_excluded_solvent(False)  # vacuum-only control
```

The current backend evaluates the exact Debye pair sum and therefore scales
quadratically with atom count. This is the reference implementation for
correctness. Adaptive multipole methods such as Pepsi-SAXS can have a much
smaller computational prefactor; a future pySEB fast backend can implement
that strategy behind the same `AtomicSAXSCloud` interface without making
Pepsi-SAXS a runtime dependency.

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
